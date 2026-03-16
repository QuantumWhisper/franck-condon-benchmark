module fc_plotting
    implicit none
contains

    subroutine plot_iv(csv_path, pdf_path, png_path, N, lambda, T, vmode, &
                       alphaL, alphaR, eta, Vg, wall_time, spec)
        character(len=*), intent(in) :: csv_path, pdf_path, png_path, spec
        integer, intent(in) :: N
        real(8), intent(in) :: lambda, T, vmode, alphaL, alphaR, eta, Vg, wall_time
        integer :: unit_num, ios
        character(len=1024) :: label_str

        ! Write gnuplot script to a temporary file
        open(newunit=unit_num, file='_fc_plot.gp', status='replace', iostat=ios)
        if (ios /= 0) then
            write(0, '(A)') 'Warning: could not create gnuplot script.'
            return
        end if

        write(unit_num, '(A)') "set datafile separator ','"
        write(unit_num, '(A)') "set key left top box opaque"
        write(unit_num, '(A)') "set grid lw 0.5 lc rgb '#cccccc'"
        write(unit_num, '(A)') "set xlabel 'V_{sd} (V)'"
        write(unit_num, '(A)') "set ylabel '{/Italic I} (A)'"
        write(unit_num, '(A,I0,A,F4.1,A,F5.1,A,F5.0,A)') &
            "set title 'N=", N, ", lambda=", lambda, &
            ", T=", T, " K, hbar*omega=", vmode * 1.0d3, " meV'"
        write(label_str, '(A,F5.2,A,F5.2,A,F3.1,A,F4.1,A,F7.2,A)') &
            "set label 1 'Fortran: aL=", alphaL, &
            ", aR=", alphaR, ", eta=", eta, &
            ", Vg=", Vg, " V  t=", wall_time, &
            " s' at graph 0.98,0.03 right front"
        write(unit_num, '(A)') trim(label_str)
        write(unit_num, '(A)') "set style line 1 lc rgb '#111111' pt 7 ps 0.35 lw 1.0"
        write(unit_num, '(A)') "set style line 2 lc rgb '#0059b3' lw 1.5"
        write(unit_num, '(A)') "set style line 3 lc rgb '#cc1f1f' lw 1.5 dt 2"

        write(unit_num, '(A)') "set terminal pdfcairo enhanced font 'Times,10' size 12cm,9cm"
        write(unit_num, '(A)') "set output '"//trim(pdf_path)//"'"
        write(unit_num, '(A)') "plot '"//trim(csv_path)// &
            "' every ::1 using 1:2 with points ls 1 title 'I_{total}'," &
            //" '' every ::1 using 1:3 with lines ls 2 title 'I_{seq}'," &
            //" '' every ::1 using 1:4 with lines ls 3 title 'I_{cot}'"

        write(unit_num, '(A)') "set terminal pngcairo enhanced font 'Times,10' size 12cm,9cm"
        write(unit_num, '(A)') "set output '"//trim(png_path)//"'"
        write(unit_num, '(A)') "replot"
        write(unit_num, '(A)') "unset output"

        close(unit_num)

        ! Execute gnuplot
        call execute_command_line('gnuplot _fc_plot.gp', exitstat=ios)
        if (ios /= 0) then
            write(0, '(A)') 'Warning: gnuplot not available, skipping plot generation.'
        end if

        ! Clean up
        call execute_command_line('rm -f _fc_plot.gp', exitstat=ios)
    end subroutine plot_iv

end module fc_plotting
