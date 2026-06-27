# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Evaluation harness for natural-language query translation.

Runs a set of example questions (``nl_query_eval_cases.json``) through the live model
and checks the structured interpretation against per-case expectations, reporting a
translation-accuracy score. Optionally (``--search``) it also resolves each query and
runs it against the database to surface zero-result translations.

Requires ``ANTHROPIC_API_KEY``; ``--search`` additionally requires a configured
database (e.g. ``ORD_INTERFACE_POSTGRES``). Run as::

    python -m ord_interface.api.nl_query_eval [--search] [--model MODEL]
"""

from __future__ import annotations

import argparse
import asyncio
import os
import time
from importlib import resources

import anthropic
import yaml
from ord_schema.logging import get_logger
from pydantic import BaseModel

from ord_interface.api.nl_query import (
    NLQuery,
    build_query_params,
    translate,
)
from ord_interface.api.search import run_query

logger = get_logger(__name__)

# Numeric filters the model must populate only when the question asks for them; an
# unexpected value here is over-extraction and counts as a miss.
_NUMERIC_FIELDS = ("min_yield", "max_yield", "min_conversion", "max_conversion")

# Bound each DB search so a pathologically slow query (e.g. a common-scaffold
# SUBSTRUCTURE match) is reported rather than hanging the whole sweep.
SEARCH_TIMEOUT_SECONDS = 60.0


class ComponentExpectation(BaseModel):
    """Expected component constraint.

    The identifier is matched exactly (case- and whitespace-insensitive), so the model
    must reproduce the specific compound -- e.g. "4-aminophenol", not "aminophenol".
    """

    identifier: str
    target: str
    mode: str


class CaseExpectation(BaseModel):
    """Per-case expectations checked against the model's interpretation."""

    components: list[ComponentExpectation] = []
    min_yield: float | None = None
    max_yield: float | None = None
    min_conversion: float | None = None
    max_conversion: float | None = None


class EvalCase(BaseModel):
    """A single evaluation example: a question and its expected interpretation."""

    question: str
    expect: CaseExpectation


def load_cases() -> list[EvalCase]:
    """Loads the evaluation cases bundled alongside this module."""
    raw = (resources.files("ord_interface.api") / "nl_query_eval_cases.yaml").read_text(
        encoding="utf-8"
    )
    return [EvalCase.model_validate(case) for case in yaml.safe_load(raw)]


def check_interpretation(expect: CaseExpectation, interpretation: NLQuery) -> list[str]:
    """Returns a list of mismatch messages between expectation and interpretation.

    An empty list means the interpretation matched all expectations.

    Args:
        expect: The case's expected interpretation.
        interpretation: The structured query the model produced.

    Returns:
        Human-readable mismatch descriptions; empty if the case passed.
    """
    mismatches = []
    remaining = list(interpretation.components)
    for wanted in expect.components:
        for candidate in remaining:
            if (
                candidate.target == wanted.target
                and candidate.mode == wanted.mode
                and candidate.identifier.strip().lower()
                == wanted.identifier.strip().lower()
            ):
                remaining.remove(candidate)
                break
        else:
            mismatches.append(
                f"missing component {wanted.identifier!r} "
                f"({wanted.target}/{wanted.mode})"
            )
    for extra in remaining:
        mismatches.append(
            f"unexpected component {extra.identifier!r} ({extra.target}/{extra.mode})"
        )
    for field in _NUMERIC_FIELDS:
        wanted_value = getattr(expect, field)
        actual_value = getattr(interpretation, field)
        if wanted_value is None and actual_value is not None:
            mismatches.append(f"unexpected {field}={actual_value}")
        elif wanted_value is not None and actual_value != wanted_value:
            mismatches.append(f"{field}: expected {wanted_value}, got {actual_value}")
    return mismatches


class CaseResult(BaseModel):
    """Outcome of evaluating one case, with a per-phase time breakdown (seconds)."""

    question: str
    passed: bool
    mismatches: list[str]
    interpretation: NLQuery
    num_results: int | None = None
    error: str | None = None
    translate_s: float = 0.0
    resolve_s: float = 0.0
    search_s: float = 0.0


async def evaluate_case(
    case: EvalCase, client: anthropic.AsyncAnthropic, search: bool
) -> CaseResult:
    """Evaluates a single case, optionally executing the resolved query.

    Args:
        case: The case to evaluate.
        client: Anthropic async client.
        search: If True, resolve and run the query against the database.

    Returns:
        The case result, including a per-phase time breakdown and any execution error.
    """
    start = time.perf_counter()
    interpretation = await translate(case.question, client)
    translate_s = time.perf_counter() - start
    mismatches = check_interpretation(case.expect, interpretation)
    num_results = None
    error = None
    resolve_s = 0.0
    search_s = 0.0
    if search:
        try:
            mark = time.perf_counter()
            params, _ = await build_query_params(interpretation)
            resolve_s = time.perf_counter() - mark
            mark = time.perf_counter()
            async with asyncio.timeout(SEARCH_TIMEOUT_SECONDS):
                results = await run_query(params, return_ids=True)
            search_s = time.perf_counter() - mark
            num_results = len(results)
        except Exception as exc:  # Resolution/DB failures shouldn't abort the sweep.
            error = f"{type(exc).__name__}: {exc}"
    return CaseResult(
        question=case.question,
        passed=not mismatches,
        mismatches=mismatches,
        interpretation=interpretation,
        num_results=num_results,
        error=error,
        translate_s=translate_s,
        resolve_s=resolve_s,
        search_s=search_s,
    )


async def run_eval(search: bool, model: str | None) -> list[CaseResult]:
    """Runs the full evaluation sweep over all bundled cases.

    Args:
        search: If True, also execute each query against the database.
        model: Optional model override (sets ORD_NL_QUERY_MODEL for this process).

    Returns:
        One CaseResult per case.
    """
    if model:
        os.environ["ORD_NL_QUERY_MODEL"] = model
    client = anthropic.AsyncAnthropic()
    cases = load_cases()
    results = []
    for index, case in enumerate(cases, start=1):
        logger.info(f"[{index}/{len(cases)}] {case.question}")
        result = await evaluate_case(case, client, search)
        # Stream the timing immediately so a slow sweep still yields data per case.
        logger.info(
            f"    translate={result.translate_s:.1f}s "
            f"resolve={result.resolve_s:.1f}s search={result.search_s:.1f}s "
            f"results={result.num_results}"
        )
        results.append(result)
    return results


def _format_report(results: list[CaseResult]) -> str:
    """Formats the eval results as a readable text report with a summary line."""
    lines = []
    passed = 0
    for result in results:
        status = "PASS" if result.passed else "FAIL"
        if result.passed:
            passed += 1
        lines.append(f"[{status}] {result.question}")
        for mismatch in result.mismatches:
            lines.append(f"         - {mismatch}")
        if result.error:
            lines.append(f"         ! search error: {result.error}")
        elif result.num_results is not None:
            note = " (ZERO RESULTS)" if result.num_results == 0 else ""
            lines.append(f"         -> {result.num_results} reactions{note}")
        total_s = result.translate_s + result.resolve_s + result.search_s
        lines.append(
            f"         time {total_s:5.1f}s "
            f"(translate {result.translate_s:.1f} + "
            f"resolve {result.resolve_s:.1f} + search {result.search_s:.1f})"
        )
    total = len(results)
    rate = passed / total if total else 0.0
    lines.append("")
    lines.append(f"Translation accuracy: {passed}/{total} ({rate:.0%})")
    if any(r.num_results is not None for r in results):
        lines.append(
            "Totals: "
            f"translate {sum(r.translate_s for r in results):.1f}s, "
            f"resolve {sum(r.resolve_s for r in results):.1f}s, "
            f"search {sum(r.search_s for r in results):.1f}s"
        )
    return "\n".join(lines)


def main() -> None:
    """CLI entry point for the evaluation harness."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--search",
        action="store_true",
        help="Resolve and execute each query against the database.",
    )
    parser.add_argument(
        "--model",
        default=None,
        help="Override the translation model (default: ORD_NL_QUERY_MODEL or Haiku).",
    )
    args = parser.parse_args()
    results = asyncio.run(run_eval(search=args.search, model=args.model))
    print(_format_report(results))


if __name__ == "__main__":
    main()
