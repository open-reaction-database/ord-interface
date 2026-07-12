/**
 * Copyright 2023 Open Reaction Database Project Authors
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

import React, { useMemo, useState } from 'react';
import { Link, useLocation } from 'wouter';
import {
  Anchor,
  Badge,
  Button,
  Paper,
  SimpleGrid,
  Skeleton,
  Text,
  TextInput,
  Title,
} from '@mantine/core';
import {
  IconBook2,
  IconBrandGithub,
  IconDatabase,
  IconExternalLink,
  IconFileCertificate,
  IconSearch,
  IconUpload,
} from '@tabler/icons-react';
import { useQuery } from '@tanstack/react-query';
import PaperButton from '../components/PaperButton';
import CopyButton from '../components/CopyButton';
import { api } from '../utils/api';
import { CITATION, FEATURED_DATASET_IDS, NEWS_ITEMS } from '../data/homeContent';
import type { Dataset } from '../types/search';
import logoUrl from '../assets/ord-logo.svg';
import classes from './Home.module.scss';

const NUM_FEATURED = 3;

const Stat: React.FC<{ value: string | null; label: string }> = ({ value, label }) => (
  <div className={classes.stat}>
    {value ? (
      <div className={classes.statValue}>{value}</div>
    ) : (
      <Skeleton
        height={34}
        width={110}
      />
    )}
    <div className={classes.statLabel}>{label}</div>
  </div>
);

const DatasetCard: React.FC<{ dataset: Dataset }> = ({ dataset }) => (
  <Paper
    component={Link}
    href={`/dataset/${dataset.dataset_id}`}
    radius="sm"
    p="lg"
    className={classes.datasetCard}
  >
    <div className={classes.datasetCardHeader}>
      <div className={classes.datasetCardName}>
        {dataset.name || dataset.dataset_id}
      </div>
      <Badge
        className={classes.countBadge}
        variant="light"
        color="gray"
        radius="xl"
      >
        {dataset.num_reactions.toLocaleString()} reactions
      </Badge>
    </div>
    <Text
      className={classes.datasetCardDescription}
      lineClamp={3}
    >
      {dataset.description || dataset.dataset_id}
    </Text>
    <div className={classes.datasetCardFooter}>View dataset →</div>
  </Paper>
);

const Home: React.FC = () => {
  const [, navigate] = useLocation();
  const [quickSearch, setQuickSearch] = useState('');

  const { data: datasets } = useQuery<Dataset[]>({
    queryKey: ['datasets'],
    queryFn: () => api.get<Dataset[]>('/datasets').then(res => res.data),
  });

  const totalReactions = useMemo(
    () => datasets?.reduce((sum, dataset) => sum + dataset.num_reactions, 0),
    [datasets],
  );

  const featured = useMemo(() => {
    if (!datasets?.length) return [];
    if (FEATURED_DATASET_IDS.length) {
      const byId = new Map(datasets.map(dataset => [dataset.dataset_id, dataset]));
      return FEATURED_DATASET_IDS.map(id => byId.get(id)).filter(
        (dataset): dataset is Dataset => dataset !== undefined,
      );
    }
    return [...datasets]
      .sort((a, b) => b.num_reactions - a.num_reactions)
      .slice(0, NUM_FEATURED);
  }, [datasets]);

  // Hands the question to the natural-language search (/ask), which resolves
  // compound names to structures with an LLM before querying.
  const submitQuickSearch = () => {
    const query = quickSearch.trim();
    if (!query) {
      navigate('/ask');
      return;
    }
    const params = new URLSearchParams();
    params.set('q', query);
    navigate(`/ask?${params.toString()}`);
  };

  return (
    <div className={classes.home}>
      {/* Hero */}
      <section className={classes.hero}>
        <img
          className={classes.logo}
          src={logoUrl}
          alt="Open Reaction Database logo"
        />
        <Text className={classes.tagline}>
          An open-access schema and infrastructure for structuring and sharing organic
          reaction data — built to power machine learning in synthesis planning,
          reaction prediction, and experiment design.
        </Text>

        <form
          className={classes.searchRow}
          onSubmit={event => {
            event.preventDefault();
            submitQuickSearch();
          }}
        >
          <TextInput
            className={classes.searchInput}
            size="md"
            radius="sm"
            leftSection={<IconSearch size={18} />}
            placeholder="Ask in plain language, e.g. reactions that make an aryl boronic acid"
            value={quickSearch}
            onChange={event => setQuickSearch(event.currentTarget.value)}
            aria-label="Ask about reactions"
          />
          <Button
            type="submit"
            size="md"
            radius="sm"
          >
            Search
          </Button>
        </form>
        <Anchor
          component={Link}
          href="/search"
          className={classes.advancedLink}
        >
          or use the advanced search →
        </Anchor>

        <div className={classes.stats}>
          <Stat
            value={totalReactions?.toLocaleString() ?? null}
            label="reactions"
          />
          <Stat
            value={datasets ? datasets.length.toLocaleString() : null}
            label="datasets"
          />
          <Stat
            value="CC-BY-SA"
            label="open license"
          />
        </div>
      </section>

      {/* Action tiles */}
      <section className={classes.tiles}>
        <PaperButton
          title="Browse datasets"
          description="Explore every dataset in the ORD"
          icon={<IconDatabase />}
          color="var(--color-blue)"
          href="/browse"
        />
        <PaperButton
          title="Search reactions"
          description="By structure, yield, dataset, or DOI"
          icon={<IconSearch />}
          color="var(--color-green)"
          href="/search"
        />
        <PaperButton
          title="Contribute data"
          description="Add your reactions with the ORD app"
          icon={<IconUpload />}
          color="var(--color-orange)"
          href="https://app.open-reaction-database.org"
        />
      </section>

      {/* Featured datasets */}
      <section className={classes.section}>
        <div className={classes.sectionHeader}>
          <Title order={2}>Featured datasets</Title>
          {datasets && (
            <Badge
              className={classes.countBadge}
              variant="light"
              color="gray"
              radius="xl"
            >
              {datasets.length.toLocaleString()} total
            </Badge>
          )}
        </div>
        <SimpleGrid
          cols={{ base: 1, sm: 3 }}
          spacing="md"
        >
          {featured.length > 0
            ? featured.map(dataset => (
                <DatasetCard
                  key={dataset.dataset_id}
                  dataset={dataset}
                />
              ))
            : Array.from({ length: NUM_FEATURED }, (_, index) => (
                <Skeleton
                  key={index}
                  height={168}
                  radius="sm"
                />
              ))}
        </SimpleGrid>
        <Anchor
          component={Link}
          href="/browse"
          className={classes.sectionFooterLink}
        >
          Browse all datasets →
        </Anchor>
      </section>

      {/* News */}
      <section className={classes.section}>
        <div className={classes.sectionHeader}>
          <Title order={2}>News & updates</Title>
        </div>
        <Paper
          radius="sm"
          className={classes.newsPaper}
        >
          {NEWS_ITEMS.map(item => (
            <div
              key={item.title}
              className={classes.newsRow}
            >
              <div className={classes.newsDate}>{item.date}</div>
              <div>
                <div className={classes.newsTitle}>
                  {item.href ? (
                    <a
                      href={item.href}
                      target="_blank"
                      rel="noopener noreferrer"
                    >
                      {item.title}
                      <IconExternalLink size={14} />
                    </a>
                  ) : (
                    item.title
                  )}
                </div>
                <Text className={classes.newsText}>{item.text}</Text>
              </div>
            </div>
          ))}
        </Paper>
      </section>

      {/* Open by design */}
      <section className={classes.section}>
        <div className={classes.sectionHeader}>
          <Title order={2}>Open by design</Title>
        </div>
        <SimpleGrid
          cols={{ base: 1, sm: 3 }}
          spacing="md"
        >
          <div className={classes.feature}>
            <div
              className={classes.featureChip}
              style={{ backgroundColor: 'var(--color-blue)' }}
            >
              <IconBook2 />
            </div>
            <div>
              <Title order={3}>Open schema</Title>
              <Text className={classes.featureText}>
                A structured, machine-readable representation for reactions — from
                benchtop to flow chemistry.
              </Text>
              <Anchor
                href="https://docs.open-reaction-database.org"
                target="_blank"
                rel="noopener noreferrer"
              >
                Read the docs →
              </Anchor>
            </div>
          </div>
          <div className={classes.feature}>
            <div
              className={classes.featureChip}
              style={{ backgroundColor: 'var(--color-green)' }}
            >
              <IconFileCertificate />
            </div>
            <div>
              <Title order={3}>Open data</Title>
              <Text className={classes.featureText}>
                Every record is freely available under CC-BY-SA and downloadable in
                bulk.
              </Text>
              <Anchor
                href="https://github.com/open-reaction-database/ord-data"
                target="_blank"
                rel="noopener noreferrer"
              >
                Get the data →
              </Anchor>
            </div>
          </div>
          <div className={classes.feature}>
            <div
              className={classes.featureChip}
              style={{ backgroundColor: 'var(--color-orange)' }}
            >
              <IconBrandGithub />
            </div>
            <div>
              <Title order={3}>Open source</Title>
              <Text className={classes.featureText}>
                The schema, this interface, and the contribution tools are all developed
                in the open.
              </Text>
              <Anchor
                href="https://github.com/open-reaction-database"
                target="_blank"
                rel="noopener noreferrer"
              >
                View on GitHub →
              </Anchor>
            </div>
          </div>
        </SimpleGrid>
      </section>

      {/* Citation */}
      <section className={classes.section}>
        <Paper
          radius="sm"
          p="lg"
          className={classes.citation}
        >
          <div>
            <Title order={3}>Citing the ORD</Title>
            <Text className={classes.citationText}>{CITATION}</Text>
          </div>
          <CopyButton
            textToCopy={CITATION}
            buttonText="Copy citation"
          />
        </Paper>
      </section>
    </div>
  );
};

export default Home;
