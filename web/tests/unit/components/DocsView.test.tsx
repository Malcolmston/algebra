import { describe, it, expect, vi, beforeEach } from 'vitest';
import { render, screen } from '@testing-library/react';
import { DocsView } from '../../../src/components/DocsView';
import type { DocIndex } from 'go-ui';

// A minimal DocIndex the stubbed fetch returns for DocsApp's doc.json request.
const DOC_INDEX: DocIndex = {
  module: 'github.com/malcolmston/algebra',
  packages: [
    {
      importPath: 'github.com/malcolmston/algebra',
      name: 'algebra',
      synopsis: 'Package algebra is a standard-library-only Go computer-algebra system.',
      doc: 'Package algebra is a standard-library-only Go computer-algebra system.',
      consts: [],
      vars: [],
      types: [
        {
          name: 'Expr',
          signature: 'type Expr interface{}',
          doc: 'Expr is the interface implemented by every node in an expression tree.',
          consts: [],
          vars: [],
          funcs: [],
          methods: [],
        },
      ],
      funcs: [{ name: 'Sym', signature: 'func Sym(name string) Expr', doc: 'Sym returns the symbol with the given name.' }],
    },
  ],
};

describe('DocsView', () => {
  beforeEach(() => {
    // DocsApp fetches doc.json; return the small index.
    global.fetch = vi.fn((input: RequestInfo | URL) => {
      if (String(input).includes('doc.json')) {
        return Promise.resolve({ ok: true, json: () => Promise.resolve(DOC_INDEX) } as Response);
      }
      return new Promise<Response>(() => {});
    }) as unknown as typeof fetch;
  });

  it('renders the inline React API reference from the fetched doc.json', async () => {
    const { container } = render(<DocsView />);
    expect(container.querySelector('#view-docs')).not.toBeNull();
    expect(
      screen.getByRole('heading', { level: 2, name: /API documentation/ }),
    ).toBeInTheDocument();

    // DocsApp fetches asynchronously, then renders the package view + symbols.
    expect(await screen.findByRole('heading', { name: /package algebra/ })).toBeInTheDocument();
    expect(container.querySelector('#sym-Sym'), 'func Sym symbol card').not.toBeNull();
    expect(container.querySelector('#sym-Expr'), 'type Expr symbol card').not.toBeNull();

    // The secondary link to the raw generated static HTML remains.
    expect(screen.getByRole('link', { name: /Open the raw generated HTML/ })).toHaveAttribute('href', './api/');
  });
});
