# Benchmark Specification

This document is the contract all language implementations must follow. Deviating from it makes results incomparable.

## Objective

Compare wall-clock time and numerical accuracy of the Franck-Condon I-V simulation across languages and compilers. The goal is an apples-to-apples comparison of the core numerical work, not I/O or startup overhead.

## Parameter Specs

Two parameter specs are provided. All implementations must support both.

| Spec | File | N | Purpose | Typical MATLAB time |
|------|------|---|---------|---------------------|
| **quick** | `benchmark/spec/quick_params.json` | 6 | Fast CI/sanity check, cross-language numerical validation | Minutes |
| **default** | `benchmark/spec/default_params.json` | 15 | Publication-quality I-V curves, real performance benchmark | Hours (MATLAB), minutes (compiled/JIT) |

Use `quick` for development iteration and correctness checks. Use `default` for the official benchmark numbers and publication-quality figures.

All parameters except N are identical between specs. Do not hardcode parameters — always read from the JSON file.

## What Is Timed

The entire I-V sweep loop, from the first bias point to the last. This includes:

- Rate matrix construction (FC matrix elements, sequential + cotunneling rates)
- Steady-state probability solve (linear system 0 = WP)
- Current calculation (I_seq + I_cot)

Start the timer immediately before the first bias point. Stop it immediately after the last current value is computed.

## What Is NOT Timed

- File I/O, JSON parsing, result saving
- Path setup, imports, module loading
- Compilation (for JIT languages, a warm-up run must be completed before timing starts)
- Plotting or visualization

## Required Outputs

Each implementation must produce **four output files** per run:

| File | Format | Purpose |
|------|--------|---------|
| `{language}_{spec}_results.json` | JSON | Full results with metadata |
| `{language}_{spec}_IV.csv` | CSV | I-V data for reuse in other tools |
| `{language}_{spec}_IV.pdf` | PDF | Publication-quality plot (vector) |
| `{language}_{spec}_IV.png` | PNG | Publication-quality plot (300 dpi) |

All outputs go to `benchmark/results/`.

### JSON format

```json
{
  "spec": "default",
  "language": "string",
  "version": "string",
  "wall_time_seconds": 0.0,
  "parameters": { ... },
  "Vsd": [...],
  "I_tol": [...],
  "I_seq": [...],
  "I_cot": [...]
}
```

### CSV format

```
Vsd_V,I_tol_A,I_seq_A,I_cot_A
0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00
3.000000e-03,...,...,...
```

### Plot requirements

The I-V plot must contain three curves on a single axes:

| Curve | Style | Label |
|-------|-------|-------|
| I_total | Filled markers (small) | `$I_{\mathrm{total}}$` |
| I_seq | Solid line | `$I_{\mathrm{seq}}$` |
| I_cot | Dashed line | `$I_{\mathrm{cot}}$` |

Required plot elements:
- X-axis: `$V_{\mathrm{sd}}$ (V)` (LaTeX)
- Y-axis: `$I$ (A)` (LaTeX)
- Title: key parameters (N, lambda, T, vmode)
- Legend with LaTeX labels
- Annotation: tunnel coupling, voltage division, wall time, language
- Figure size: ~12 cm x 9 cm (single-column paper width)
- Font: LaTeX-rendered, 10-12 pt
- Resolution: 300 dpi for PNG

The goal is publication quality — the plot should be directly usable in a journal paper.

## Correctness Validation

All implementations must produce I-V curves matching the MATLAB reference within a relative tolerance of 1e-10. The comparison formula is:

```
max(abs(I_impl - I_ref) ./ max(abs(I_ref), 1e-30)) < 1e-10
```

Run this check before reporting timing results. A fast but wrong implementation is not a valid entry.

## Benchmark Protocol

1. **Cold start**: Fresh process or session, no cached results from a previous run.
2. **Warm-up** (JIT languages only): Complete one full run and discard the result before starting the timer. This accounts for JIT compilation overhead.
3. **Timed runs**: Minimum 3 runs. Report the median wall time.
4. **Hardware**: Report CPU model, RAM, and OS version alongside results.

## Fair Comparison Notes

- **Single-threaded** comparison is the primary metric. This is the apples-to-apples baseline.
- Parallel versions may be reported separately, with thread count noted explicitly.
- No GPU acceleration. The algorithm has sequential dependencies across the bias sweep that make GPU offloading impractical for this problem size.
- Use the same parameter spec for all languages in a comparison. Do not tune N or other parameters to make your language look faster.
