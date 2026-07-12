import { describe, it, expect, vi, beforeEach } from 'vitest';
import { render, screen } from '@testing-library/react';
import { Hero } from '../../../src/components/Hero';
import { ALGEBRA } from '../../../src/data';

describe('Hero', () => {
  beforeEach(() => {
    // VersionBadge fetches on mount; keep it pending so the hero renders cleanly.
    global.fetch = vi.fn().mockReturnValue(new Promise(() => {}));
  });

  it('renders the name, package path and tagline', () => {
    render(<Hero lib={ALGEBRA} />);
    expect(screen.getByRole('heading', { level: 2, name: /algebra/ })).toBeInTheDocument();
    expect(screen.getByText(ALGEBRA.pkg)).toBeInTheDocument();
    expect(screen.getByText(ALGEBRA.tagline)).toBeInTheDocument();
  });

  it('renders the GitHub link opening in a new tab', () => {
    render(<Hero lib={ALGEBRA} />);
    const github = screen.getByRole('link', { name: /GitHub/ });
    expect(github).toHaveAttribute('href', ALGEBRA.repo);
    expect(github).toHaveAttribute('target', '_blank');
    expect(github).toHaveAttribute('rel', expect.stringContaining('noopener'));
  });
});
