/**
 * Copyright 2026 Open Reaction Database Project Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * Curated content for the landing page. Edit this file to feature datasets or
 * post news — no other code changes needed.
 */

/**
 * Dataset IDs to feature on the landing page, in display order. When empty,
 * the page falls back to the largest datasets from the live /api/datasets
 * listing.
 */
export const FEATURED_DATASET_IDS: string[] = [];

export interface NewsItem {
  /** Display date, e.g. 'November 2021'. */
  date: string;
  title: string;
  text: string;
  href?: string;
}

/** Most recent first. */
export const NEWS_ITEMS: NewsItem[] = [
  {
    date: '2026',
    title: 'A new way to contribute',
    text:
      'Creating and editing datasets moved to the new ORD app, with a ' +
      'redesigned reaction editor and dataset sharing.',
    href: 'https://app.open-reaction-database.org',
  },
  {
    date: 'July 2023',
    title: 'Data sharing in chemistry: lessons learned',
    text:
      'A J. Chem. Inf. Model. perspective makes the case for mandating ' +
      'structured reaction data, drawing on ORD experience.',
    href: 'https://doi.org/10.1021/acs.jcim.3c00607',
  },
  {
    date: 'November 2021',
    title: 'The Open Reaction Database is published in JACS',
    text:
      'The launch paper describes the schema, infrastructure, and vision ' +
      'behind the ORD.',
    href: 'https://doi.org/10.1021/jacs.1c09820',
  },
];

export const CITATION =
  'Kearnes, S. M.; Maser, M. R.; Wleklinski, M.; Kast, A.; Doyle, A. G.; ' +
  'Dreher, S. D.; Hawkins, J. M.; Jensen, K. F.; Coley, C. W. The Open ' +
  'Reaction Database. J. Am. Chem. Soc. 2021, 143 (45), 18820–18826. ' +
  'https://doi.org/10.1021/jacs.1c09820';
