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
├── python/          # Python+Numba implementation ✅
│   ├── src/         # Core simulation modules
│   ├── run_benchmark.py
│   └── requirements.txt
├── c/               # C+GSL implementation ✅
│   ├── src/         # Core simulation modules (17 .h/.c files)
│   ├── main.c       # Benchmark runner
│   ├── Makefile     # Build for macOS + Linux
│   ├── cJSON.h/.c   # Vendored JSON library
│   └── c_benchmark  # Compiled binary
├── fortran/         # Fortran implementation ✅
│   ├── src/         # Core simulation modules (14 .f90 files)
│   ├── main.f90     # Benchmark runner
│   └── Makefile     # Build for macOS + Linux
├── rust/            # Rust implementation ✅ + interactive explorer
├── cpp/             # C++ implementation ✅
│   ├── src/         # Core simulation modules (14 .hpp/.cpp files)
│   ├── CMakeLists.txt # CMake build system
│   └── build/       # Build directory
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

## Quick Start (Python)

```bash
cd python
pip install -r requirements.txt
python run_benchmark.py quick       # Fast validation (N=6, ~11 seconds)
python run_benchmark.py default     # Full benchmark (N=15)
python run_benchmark.py             # Same as 'default'
```

Uses `scipy.special.digamma` for complex digamma and a custom Numba-JIT'd trigamma via asymptotic series. No warm-up needed (not JIT at the simulation level).

## Quick Start (C)

```bash
cd c
brew install gsl gnuplot   # macOS; Linux: apt install libgsl-dev gnuplot
make
./c_benchmark quick         # Fast validation (N=6, ~2.3 seconds)
./c_benchmark default       # Full benchmark (N=15)
./c_benchmark               # Same as 'default'
```

Requires GSL (GNU Scientific Library) for QR decomposition and complex digamma. Uses cJSON (vendored) for JSON I/O and gnuplot for PDF/PNG plots. No warm-up needed (compiled, not JIT).

## Quick Start (Rust)

```bash
cd rust
cargo build --release
./target/release/rust_benchmark quick     # Fast validation (N=6, ~0.28 seconds)
./target/release/rust_benchmark default   # Full benchmark (N=15)
./target/release/rust_benchmark           # Same as 'default'
```

Pure Rust implementation with no C/FFI dependencies. Uses `num-complex` for complex arithmetic, `nalgebra` for QR decomposition, and custom asymptotic series for complex digamma/trigamma. Includes the factored digamma precomputation optimization from the C port. Uses gnuplot for PDF/PNG plots.

### Interactive Explorer (Rust TUI)

```bash
cd rust
cargo run --release --bin rust_explorer                      # Fresh session
cargo run --release --bin rust_explorer -- --load session.json  # Resume saved session
```

A terminal UI for real-time parameter exploration with three modes:

- **I-V Curve**: adjust parameters with arrow keys and see the I-V curve update automatically (300ms debounce). Three curves: I_tol (green), I_seq (yellow), I_cot (cyan) rendered with Braille markers.
- **Stability Diagram**: 2D color map of |I(Vg, Vsd)| with Viridis colormap, log-scale colorbar, axis labels, and a progress bar. FCCache and DigammaTable are shared across all Vg points for efficiency.
- **Temperature Diagram**: 2D color map of |I(T, Vsd)| — sweeps temperature (1–50 K default) on the x-axis and bias voltage on the y-axis. Viridis colormap with log-scale colorbar. Useful for visualizing how the Franck-Condon blockade evolves with temperature.

Each mode supports four display types, cycled with `d`: **Current** (I), **Conductance** (dI/dV), **IETS** (d²I/dV²), and **Normalized IETS** ((d²I/dV²)/(dI/dV)). IETS uses 3-point central finite differences; normalized IETS divides by conductance with safe zero-guard.

Sessions can be saved (`s` key) and resumed later (`--load`), preserving all parameters, computed results, and display settings. This avoids re-running expensive stability/temperature diagram computations.

Keybindings: `↑↓` select parameter, `←→` adjust (Shift=fine), `d` cycle display (I → G → IETS → nIETS), `Enter` run, `Tab` switch mode, `e` export CSV, `s` save session, `p` plot, `Esc` cancel, `q` quit.

## Quick Start (C++)

```bash
cd cpp
brew install eigen nlohmann-json libomp gnuplot   # macOS
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/cpp_benchmark quick     # Fast validation (N=6, ~1.2 seconds)
./build/cpp_benchmark default   # Full benchmark (N=15)
./build/cpp_benchmark           # Same as 'default'
```

Modern C++17 implementation with no C library dependencies. Uses Eigen (header-only) for QR decomposition, nlohmann/json (header-only) for JSON I/O, custom asymptotic series for complex digamma/trigamma, and OpenMP for parallel bias-point computation. Uses gnuplot for PDF/PNG plots.

## Quick Start (Fortran)

```bash
cd fortran
make
./fortran_benchmark quick     # Fast validation (N=6, ~1.6 seconds)
./fortran_benchmark default   # Full benchmark (N=15)
./fortran_benchmark           # Same as 'default'
```

Pure Fortran 2008 implementation with no C/FFI dependencies. Uses LAPACK (via macOS Accelerate or `-llapack -lblas`) for linear system solve. Custom complex digamma/trigamma via Rust's optimized 10-term Bernoulli asymptotic series. Includes factored digamma precomputation. Uses gnuplot for PDF/PNG plots.

## Outputs

Each run produces four files in `benchmark/results/`:
- `{lang}_{spec}_results.json` — full results with metadata and timing
- `{lang}_{spec}_IV.csv` — I-V data for reuse
- `{lang}_{spec}_IV.pdf` — publication-quality vector plot
- `{lang}_{spec}_IV.png` — publication-quality raster plot (300 dpi)

## Benchmark Results

| Language | Quick (N=6) | Default (N=15) | Speedup vs MATLAB |
|----------|-------------|----------------|---------------------|
| MATLAB | 1844 s | 11725 s | 1× (reference) |
| Python | 11 s | — ᵇ | **169×** |
| Julia | 5.2 s | 102 s | **115×** |
| C (GSL) | 2.3 s | 17.1 s | **685×** |
| Fortran | 1.6 s | 12.3 s | **953×** |
| C++ | 1.2 s | 8.7 s | **1348×** |
| Rust | 0.28 s | 1.9 s | **6202×** |

Speedup column refers to the default spec (N=15), the primary benchmark. All ports replace MATLAB's symbolic digamma bottleneck with native complex implementations.

**Hardware:**
- Quick spec ᵃ: Apple M4 Max (14 cores: 4E+10P, 32 GPU cores)
- Default spec: MATLAB on Apple M4 Max; all other languages on Apple M5 (10 cores: 4S+6E, 24 GB)

ᵃ Quick-spec speedups vs MATLAB: Python 169×, Julia 358×, C 809×, Fortran 1153×, C++ 1598×, Rust 6638×.
ᵇ Python's default spec (N=15) is too slow for practical benchmarking on consumer hardware (~30 min+ per run estimated). The quick spec (N=6) completes in 11 seconds.

**Numerical accuracy (default spec):** All ports match MATLAB to **6.1×10⁻⁶** max relative error (excluding 2 solver-artifact bias points at Vsd ≈ 0.219 V and 0.438 V where the rate matrix W is ill-conditioned and different linear solvers give different results). Non-MATLAB ports agree with each other to **<5×10⁻¹³**.

## Language Status

| Language | Status | Speedup (Default) | Notes |
|----------|--------|-------------------|-------|
| MATLAB | ✅ Reference | 1× | Symbolic Math Toolbox required |
| Rust | ✅ Complete | **6202×** | Fastest. Rayon parallel + optimized digamma |
| C++ | ✅ Complete | **1348×** | OpenMP parallel + Eigen QR |
| Fortran | ✅ Complete | **953×** | LAPACK (Accelerate) + Rust-derived digamma |
| C (GSL) | ✅ Complete | **685×** | GSL digamma + factored precomputation |
| Julia | ✅ Complete | **115×** | SpecialFunctions.jl digamma |
| Python (Numba) | ✅ Complete | — | 169× on quick spec; default spec too slow |

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
