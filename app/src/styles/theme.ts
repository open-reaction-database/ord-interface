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

// Mantine theme mirroring ord-app (ui/src/common/styling/theme.ts) so the two
// ORD frontends share one look. Keep the two files in sync until the theme is
// extracted into a shared package.
import { colorsTuple, createTheme } from '@mantine/core';

type MantineColorsTuple = ReturnType<typeof colorsTuple>;
type MutableMantineColorsTuple = {
  -readonly [K in keyof MantineColorsTuple]: MantineColorsTuple[K];
};

const primaryColors = colorsTuple('#3C78D8') as MutableMantineColorsTuple;
primaryColors[9] = '#001926';
const secondaryColors = colorsTuple(['#aab9c5', '#637d92', '#4a5e6d']);

export const theme = createTheme({
  fontFamily: '"Roboto", "Helvetica", "Arial", sans-serif',
  fontSizes: {
    sm: '0.875rem',
    md: '0.875rem',
  },
  lineHeights: {
    md: '1.5',
  },
  headings: {
    sizes: {
      h1: {
        fontSize: '2rem',
        lineHeight: '1.25',
      },
      h2: {
        fontSize: '1.5rem',
        lineHeight: '1.34',
      },
      h3: {
        fontSize: '1rem',
        lineHeight: '1.5',
      },
      h4: {
        fontSize: '0.875rem',
        lineHeight: '1.5',
      },
      h5: {
        fontSize: '0.75rem',
        lineHeight: '1.34',
      },
    },
  },
  colors: {
    primary: primaryColors as MantineColorsTuple,
    secondary: secondaryColors,
  },
  primaryColor: 'primary',
  radius: {
    xs: '2px',
    sm: '6px',
    md: '10px',
    lg: '24px',
    xl: '36px',
  },
  spacing: {
    xs: '0.25rem',
    sm: '0.75rem',
    md: '1rem',
    lg: '1.5rem',
    xl: '2rem',
  },
  components: {
    Paper: {
      defaultProps: {
        withBorder: true,
      },
    },
    Tooltip: {
      styles: {
        tooltip: {
          wordBreak: 'break-all',
          maxWidth: '60dvw',
          overflow: 'hidden',
          textWrap: 'auto',
        },
      },
    },
    Modal: {
      styles: {
        title: {
          fontSize: '24px',
          lineHeight: '32px',
          fontWeight: 600,
        },
      },
    },
    InputWrapper: {
      defaultProps: {
        inputWrapperOrder: ['label', 'input', 'description', 'error'],
      },
      styles: {
        label: {
          color: 'var(--color-text-secondary-2)',
          paddingBottom: '8px',
        },
        error: {
          marginTop: '4px',
        },
      },
    },
    Textarea: {
      styles: {
        input: {
          minHeight: '92px',
        },
      },
    },
    Accordion: {
      styles: {
        control: {
          backgroundColor: '#F2F2F2',
        },
        content: {
          padding: 'var(--mantine-spacing-xs) var(--mantine-spacing-md)',
        },
      },
    },
    Badge: {
      styles: {
        label: {
          textTransform: 'none',
        },
      },
    },
    Drawer: {
      styles: {
        header: {
          borderBottom: '1px solid var(--color-border-1)',
        },
        body: {
          padding: 'var(--mantine-spacing-md) var(--mantine-spacing-lg)',
        },
        title: {
          fontSize: 'var(--mantine-h2-font-size)',
          lineHeight: 'var(--mantine-h2-line-height)',
          fontWeight: 600,
        },
      },
    },
  },
});
