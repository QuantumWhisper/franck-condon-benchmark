"""Publication-quality I-V characteristic plot.

Generates plots matching the BENCHMARK.md specification:
- 3 curves: I_total (filled markers), I_seq (solid line), I_cot (dashed line)
- LaTeX labels and title
- 12 x 9 cm figure, 300 dpi PNG, vector PDF
- Annotation with tunnel coupling, voltage division, wall time, language
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for headless rendering
import matplotlib.pyplot as plt


def plot_iv(Vsd, I_tol, I_seq, I_cot, N, lambda_, T, vmode, alphaL, alphaR, eta, Vg,
            wall_time, pdf_path, png_path):
    """Generate publication-quality I-V characteristic plot.

    Matches MATLAB plot_IV in run_benchmark.m and Julia plotting.jl.

    Parameters
    ----------
    Vsd : ndarray
        Bias voltages.
    I_tol, I_seq, I_cot : ndarray
        Total, sequential, and cotunneling currents.
    N, lambda_, T, vmode, alphaL, alphaR, eta, Vg : float
        Physical parameters for title/annotation.
    wall_time : float
        Wall time in seconds.
    pdf_path, png_path : str
        Output file paths.
    """
    # Try to use LaTeX rendering, fall back gracefully
    try:
        plt.rcParams.update({
            'text.usetex': True,
            'font.family': 'serif',
            'font.serif': ['Computer Modern Roman'],
        })
        use_latex = True
    except Exception:
        use_latex = False

    # Colors matching MATLAB/Julia reference
    c_tol = (0.1, 0.1, 0.1)       # near-black for total
    c_seq = (0.0, 0.35, 0.75)     # blue for sequential
    c_cot = (0.85, 0.15, 0.15)    # red for cotunneling

    # Figure: 12 x 9 cm (single-column paper width)
    cm = 1.0 / 2.54  # cm to inches
    fig, ax = plt.subplots(1, 1, figsize=(12 * cm, 9 * cm))
    fig.patch.set_facecolor('white')

    # I_total: filled markers
    ax.plot(Vsd, I_tol, 'o',
            color=c_tol, markersize=2.5, markerfacecolor=c_tol,
            label=r'$I_{\mathrm{total}}$' if use_latex else 'I_total')

    # I_seq: solid line
    ax.plot(Vsd, I_seq, '-',
            color=c_seq, linewidth=1.4,
            label=r'$I_{\mathrm{seq}}$' if use_latex else 'I_seq')

    # I_cot: dashed line
    ax.plot(Vsd, I_cot, '--',
            color=c_cot, linewidth=1.4,
            label=r'$I_{\mathrm{cot}}$' if use_latex else 'I_cot')

    # Labels
    if use_latex:
        ax.set_xlabel(r'$V_{\mathrm{sd}}$ (V)', fontsize=12)
        ax.set_ylabel(r'$I$ (A)', fontsize=12)
        title_str = (r'$N{=}%d,\;\lambda{=}%.1f,\;T{=}%.1f$ K,'
                     r'\;$\hbar\omega{=}%d$ meV' % (N, lambda_, T, vmode * 1e3))
        ax.set_title(title_str, fontsize=11)
    else:
        ax.set_xlabel('Vsd (V)', fontsize=12)
        ax.set_ylabel('I (A)', fontsize=12)
        ax.set_title(f'N={N}, λ={lambda_:.1f}, T={T:.1f} K, ħω={vmode*1e3:.0f} meV',
                     fontsize=11)

    # Legend
    lgd = ax.legend(fontsize=10, loc='upper left',
                    framealpha=0.85, edgecolor='black', fancybox=False)

    # Axes styling
    ax.tick_params(labelsize=10)
    ax.spines['top'].set_linewidth(0.6)
    ax.spines['bottom'].set_linewidth(0.6)
    ax.spines['left'].set_linewidth(0.6)
    ax.spines['right'].set_linewidth(0.6)
    ax.grid(True, linewidth=0.5, alpha=0.12)

    # Annotation: tunnel coupling, voltage division, wall time, language
    if use_latex:
        anno_str = (r'$\Gamma_{L,R}/\hbar\omega = %.2f$, $\eta = %.1f$, '
                    r'$V_g = %.1f$ V' '\n'
                    r'Python, $t = %.1f$ s' % (alphaL, eta, Vg, wall_time))
    else:
        anno_str = (f'Γ_L,R/ħω = {alphaL:.2f}, η = {eta:.1f}, '
                    f'Vg = {Vg:.1f} V\n'
                    f'Python, t = {wall_time:.1f} s')

    ax.annotate(anno_str, xy=(0.98, 0.02), xycoords='axes fraction',
                fontsize=8, ha='right', va='bottom',
                bbox=dict(boxstyle='square,pad=0.3', facecolor='white', alpha=0.7,
                          edgecolor='none'))

    fig.tight_layout()

    # Save
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight', facecolor='white')
    fig.savefig(png_path, format='png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print("Plot saved.")
