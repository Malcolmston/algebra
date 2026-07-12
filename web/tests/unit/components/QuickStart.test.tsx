import { describe, it, expect } from 'vitest';
import { render, screen } from '@testing-library/react';
import { QuickStart } from '../../../src/components/QuickStart';
import { ALGEBRA } from '../../../src/data';

describe('QuickStart', () => {
  it('renders the Quick start heading and highlighted Go snippet', () => {
    const { container } = render(<QuickStart lib={ALGEBRA} />);
    expect(container.querySelector(`#${ALGEBRA.id}-quick`)).not.toBeNull();
    expect(screen.getByRole('heading', { name: 'Quick start' })).toBeInTheDocument();
    // The snippet mentions algebra.MustParse.
    expect(container.textContent).toContain('algebra.MustParse');
  });
});
