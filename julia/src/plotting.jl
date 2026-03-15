using CairoMakie
using LaTeXStrings

# Unit conversion for CairoMakie (internal unit = CSS px)
const _inch = 96.0
const _pt = 4.0 / 3.0
const _cm = _inch / 2.54

"""
    plot_iv(Vsd, I_tol, I_seq, I_cot, N, lambda, T, vmode, alphaL, alphaR, eta, Vg,
            wall_time, pdf_path, png_path)

Generate publication-quality I-V characteristic plot.
Matches MATLAB `plot_IV` in run_benchmark.m: 3 curves, LaTeX labels, 12×9 cm.
"""
function plot_iv(Vsd, I_tol, I_seq, I_cot, N, lambda, T, vmode, alphaL, alphaR, eta, Vg,
                 wall_time, pdf_path, png_path)

    # Colors matching MATLAB reference
    c_tol = RGBf(0.1, 0.1, 0.1)       # near-black for total
    c_seq = RGBf(0.0, 0.35, 0.75)     # blue for sequential
    c_cot = RGBf(0.85, 0.15, 0.15)    # red for cotunneling

    # Publication theme
    pub_theme = Theme(
        fontsize = 10 * _pt,
        Axis = (
            spinewidth = 0.6,
            xtickwidth = 0.6,
            ytickwidth = 0.6,
            xgridvisible = true,
            ygridvisible = true,
            xgridwidth = 0.5,
            ygridwidth = 0.5,
            xgridcolor = RGBAf(0, 0, 0, 0.12),
            ygridcolor = RGBAf(0, 0, 0, 0.12),
            xgridstyle = :solid,
            ygridstyle = :solid,
            xticklabelsize = 10 * _pt,
            yticklabelsize = 10 * _pt,
            xlabelsize = 12 * _pt,
            ylabelsize = 12 * _pt,
            titlesize = 11 * _pt,
        ),
        Legend = (
            framewidth = 0.5,
            labelsize = 10 * _pt,
            patchsize = (20 * _pt, 8 * _pt),
            padding = (6, 6, 4, 4),
        ),
    )

    with_theme(pub_theme) do
        fig = Figure(size = (12 * _cm, 9 * _cm))
        ax = Axis(fig[1, 1],
            xlabel = L"$V_{\mathrm{sd}}$ (V)",
            ylabel = L"$I$ (A)",
            title = latexstring(
                "\$N{=}$(N),\\;\\lambda{=}$(lambda),\\;T{=}$(T)\$ K,\\;\$\\hbar\\omega{=}$(round(Int, vmode*1e3))\$ meV"
            ),
        )

        # I_total: filled markers
        scatter!(ax, collect(Vsd), collect(I_tol),
            color = c_tol,
            marker = :circle,
            markersize = 2.5 * _pt,
            strokewidth = 0,
            label = L"$I_{\mathrm{total}}$",
        )

        # I_seq: solid line
        lines!(ax, collect(Vsd), collect(I_seq),
            color = c_seq,
            linewidth = 1.4,
            label = L"$I_{\mathrm{seq}}$",
        )

        # I_cot: dashed line
        lines!(ax, collect(Vsd), collect(I_cot),
            color = c_cot,
            linewidth = 1.4,
            linestyle = :dash,
            label = L"$I_{\mathrm{cot}}$",
        )

        # Legend
        axislegend(ax,
            position = :lt,
            framecolor = :black,
            backgroundcolor = RGBAf(1, 1, 1, 0.85),
        )

        # Annotation: tunnel coupling, voltage division, wall time, language
        anno_str = latexstring(
            "\$\\Gamma_{L,R}/\\hbar\\omega = $(alphaL)\$, " *
            "\$\\eta = $(eta)\$, " *
            "\$V_g = $(Vg)\$ V\n" *
            "Julia, \$t = $(round(wall_time, digits=1))\$ s"
        )
        text!(ax, 0.98, 0.02,
            text = anno_str,
            fontsize = 8 * _pt,
            color = :black,
            align = (:right, :bottom),
            space = :relative,
        )

        # Save
        save(pdf_path, fig)
        save(png_path, fig; px_per_unit = 300 / _inch)

        println("Plot saved.")
    end
end
