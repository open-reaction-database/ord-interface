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

import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import DownloadResults from './DownloadResults';

interface SentRequest {
  method: string;
  url: string;
  headers: Record<string, string>;
  body: string;
}

const sent: SentRequest[] = [];

// Minimal stand-in for the XMLHttpRequest the component uses to stream the
// .pb.gz back; the tests assert on what it was asked to send.
class FakeXHR {
  static instances: FakeXHR[] = [];
  onload: (() => void) | null = null;
  responseType = '';
  response: unknown = null;
  status = 200;
  private request: SentRequest = { method: '', url: '', headers: {}, body: '' };

  constructor() {
    FakeXHR.instances.push(this);
  }

  open(method: string, url: string) {
    this.request.method = method;
    this.request.url = url;
  }

  setRequestHeader(name: string, value: string) {
    this.request.headers[name] = value;
  }

  send(body: string) {
    this.request.body = body;
    sent.push(this.request);
  }
}

const clickedLinks: HTMLAnchorElement[] = [];

beforeEach(() => {
  sent.length = 0;
  clickedLinks.length = 0;
  FakeXHR.instances.length = 0;
  vi.stubGlobal('XMLHttpRequest', FakeXHR);
  vi.stubGlobal('URL', {
    ...URL,
    createObjectURL: vi.fn(() => 'blob:results'),
    revokeObjectURL: vi.fn(),
  });
  // jsdom does not navigate, and an un-stubbed click on a download anchor logs
  // a "not implemented" error.
  vi.spyOn(HTMLAnchorElement.prototype, 'click').mockImplementation(function (
    this: HTMLAnchorElement,
  ) {
    clickedLinks.push(this);
  });
});

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

const renderDownload = (reactionIds = ['ord-1', 'ord-2'], show = true) => {
  const onHideDownloadResults = vi.fn();
  render(
    <DownloadResults
      reactionIds={reactionIds}
      showDownloadResults={show}
      onHideDownloadResults={onHideDownloadResults}
    />,
  );
  return { onHideDownloadResults };
};

describe('DownloadResults', () => {
  it('renders nothing while hidden', () => {
    const { container } = render(
      <DownloadResults
        reactionIds={['ord-1']}
        showDownloadResults={false}
        onHideDownloadResults={vi.fn()}
      />,
    );
    expect(container).toBeEmptyDOMElement();
  });

  it('defaults to pb.gz', () => {
    renderDownload();
    expect(screen.getByLabelText('File type:')).toHaveValue('pb.gz');
    expect(
      screen.getByRole('button', { name: 'Download pb.gz file' }),
    ).toBeInTheDocument();
  });

  it('restores the previously chosen file type', () => {
    localStorage.setItem('downloadFileType', 'csv');
    renderDownload();
    expect(screen.getByLabelText('File type:')).toHaveValue('csv');
  });

  it('remembers a newly chosen file type', async () => {
    const user = userEvent.setup();
    renderDownload();

    // The other formats are not implemented yet, so drive the select directly
    // rather than through its disabled options.
    await user.selectOptions(screen.getByLabelText('File type:'), [
      screen.getByRole('option', { name: 'pb.gz' }),
    ]);

    expect(localStorage.getItem('downloadFileType')).toBe('pb.gz');
  });

  it('offers only pb.gz for now', () => {
    renderDownload();
    expect(screen.getByRole('option', { name: 'pb.gz' })).toBeEnabled();
    expect(screen.getByRole('option', { name: /csv/ })).toBeDisabled();
    expect(screen.getByRole('option', { name: /pbtxt/ })).toBeDisabled();
  });

  it('posts the reaction ids to the download endpoint', async () => {
    const user = userEvent.setup();
    renderDownload(['ord-1', 'ord-2']);

    await user.click(screen.getByRole('button', { name: /Download/ }));

    expect(sent).toHaveLength(1);
    expect(sent[0].method).toBe('POST');
    expect(sent[0].url).toBe('/api/download_search_results');
    expect(sent[0].headers['Content-Type']).toBe('application/json');
    expect(JSON.parse(sent[0].body)).toEqual({ reaction_ids: ['ord-1', 'ord-2'] });
  });

  it('saves the response under a .pb.gz filename', async () => {
    const user = userEvent.setup();
    renderDownload();

    await user.click(screen.getByRole('button', { name: /Download/ }));
    const xhr = FakeXHR.instances[0];
    xhr.response = new Blob(['data']);
    xhr.onload?.();

    expect(clickedLinks).toHaveLength(1);
    expect(clickedLinks[0].download).toBe('ord_search_results.pb.gz');
    expect(clickedLinks[0].href).toContain('blob:results');
  });

  it('does not save anything when the download fails', async () => {
    const user = userEvent.setup();
    renderDownload();

    await user.click(screen.getByRole('button', { name: /Download/ }));
    const xhr = FakeXHR.instances[0];
    xhr.status = 500;
    xhr.onload?.();

    expect(clickedLinks).toHaveLength(0);
  });

  it('closes on the modal close control', async () => {
    const user = userEvent.setup();
    const { onHideDownloadResults } = renderDownload();

    await user.click(screen.getByText('✕'));

    expect(onHideDownloadResults).toHaveBeenCalledTimes(1);
  });
});
