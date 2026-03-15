# Franck-Condon Blockade: Multi-Language I-V Simulation Benchmark

A cross-language benchmark comparing implementations of the Franck-Condon blockade I-V characteristics simulation for single-molecule junctions.

## Physics Background

A single-molecule transistor consists of a molecule bridging two metallic leads, with a gate electrode tuning the molecular energy level. When an electron tunnels onto the molecule, it couples to a vibrational mode (phonon) of frequency ω. This electron-phonon coupling, parameterized by the dimensionless constant λ, shifts the equilibrium nuclear coordinate upon charging. The overlap between vibrational wavefunctions before and after charging is captured by the Franck-Condon matrix elements.

For strong coupling (λ >> 1), the Franck-Condon matrix elements for low-phonon transitions are exponentially suppressed. This creates the "Franck-Condon blockade": even when the molecular level is within the bias window, current is blocked because the zero-phonon transition is forbidden. Current only flows once the bias is large enough to excite high-phonon states with appreciable FC overlap.

The simulation computes I-V curves by solving quantum master equations that include two tunneling mechanisms. Sequential tunneling (first-order in the tunnel coupling Γ) describes electrons hopping through one lead at a time, with the molecule spending real time in the charged state. Cotunneling (second-order in Γ) involves a virtual intermediate state where an electron transfers between leads without the molecule being charged for any measurable time. The steady-state occupation probabilities P are found by solving the rate equation 0 = WP subject to ΣP = 1, then the total current is I = I_seq + I_cot. The implementation follows Koch, von Oppen, and Glazman, PRB 74, 205438 (2006).

## Key Parameters

| Parameter | Symbol | Description | Default |
|-----------|--------|-------------|---------|
| N | N | Number of phonon Fock states | 15 (default) / 6 (quick) |
| vmode | ℏω | Phonon energy | 73 meV |
| alphaL/R | Γ_L,R / ℏω | Tunnel coupling ratios | 0.02 |
| lambda | λ | Electron-phonon coupling | 5 |
| T | T | Temperature | 4.2 K |
| eta | η | Voltage division factor | 0.5 |
| Vg | V_g | Gate voltage | 0 V |
| tau | τ | Phonon relaxation time | ∞ (unequilibrated) |

## Project Structure

```
├── matlab/          # MATLAB implementation (reference)
│   ├── src/         # Core simulation code
│   ├── utils/       # Plotting utilities
│   ├── test/        # Test functions
│   ├── setup_path.m # Path configuration
│   └── run_benchmark.m
├── julia/           # Julia implementation ✅
│   ├── src/         # Core simulation module (FranckCondon.jl)
│   ├── run_benchmark.jl
│   ├── Project.toml
│   └── Manifest.toml
├── fortran/         # Fortran implementation (planned)
├── python/          # Python+Numba implementation (planned)
├── c/               # C+GSL implementation (planned)
├── rust/            # Rust implementation (planned)
├── benchmark/
│   ├── spec/        # Benchmark parameters (quick + default)
│   └── results/     # Outputs: JSON, CSV, PDF, PNG per language
├── README.md
├── BENCHMARK.md     # Benchmark specification
└── AGENTS.md        # AI agent knowledge base
```

## Quick Start (MATLAB)

```matlab
cd matlab
results = run_benchmark('quick');    % Fast validation (N=6, ~30 min)
results = run_benchmark('default');  % Full benchmark (N=15, hours)
results = run_benchmark();           % Same as 'default'
```

## Quick Start (Julia)

```bash
cd julia
julia run_benchmark.jl quick        # Fast validation (N=6, ~5 seconds)
julia run_benchmark.jl default      # Full benchmark (N=15)
julia run_benchmark.jl              # Same as 'default'
```

First run will install dependencies automatically via `Project.toml`. Subsequent runs use the cached environment. The benchmark includes a JIT warm-up pass before timing.

## Outputs

Each run produces four files in `benchmark/results/`:
- `{lang}_{spec}_results.json` — full results with metadata and timing
- `{lang}_{spec}_IV.csv` — I-V data for reuse
- `{lang}_{spec}_IV.pdf` — publication-quality vector plot
- `{lang}_{spec}_IV.png` — publication-quality raster plot (300 dpi)

## Benchmark Results

| Language | Quick (N=6) | Default (N=15) | Speedup vs MATLAB |
|----------|-------------|----------------|---------------------|
| MATLAB | 1844 s | — | 1x (reference) |
| Julia | 5.2 s | — | **358x** |

*Quick spec on Apple Silicon. Julia uses `SpecialFunctions.jl` for native complex digamma/trigamma, replacing MATLAB's symbolic math bottleneck.*

## Language Status

| Language | Status | Notes |
|----------|--------|-------|
| MATLAB | ✅ Reference | Symbolic Math Toolbox required |
| Julia | ✅ Complete | 358x faster (quick spec) |
| Fortran | 🔲 Planned | |
| Python (Numba/JAX) | 🔲 Planned | |
| C (GSL) | 🔲 Planned | |
| Rust | 🔲 Planned | |

## Reference

Koch, J., von Oppen, F., and Glazman, L. I., *Franck-Condon blockade and giant Fano factors in transport through single molecules*, Phys. Rev. B **74**, 205438 (2006).

## Citation

If you use this software in your research, please cite it. Click **"Cite this repository"** on GitHub or use:

```bibtex
@software{ning2026franck-condon-benchmark,
  author       = {Ning, Shanglong},
  title        = {Franck-Condon Benchmark: Multi-Language I-V Simulation for Single-Molecule Junctions},
  year         = {2026},
  url          = {https://github.com/QuantumWhisper/franck-condon-benchmark}
}
```

## License

[MIT License](LICENSE). Copyright (c) 2022-2026 Shanglong Ning.
