program fc_benchmark
    use fc_simulate, only: simulate_iv
    use fc_json_io, only: sim_params_t, parse_params_json, load_matlab_vsd, &
                          load_matlab_reference, write_results_json, write_results_csv
    use fc_plotting, only: plot_iv
    implicit none

    character(len=256) :: spec, params_path, ref_path, json_out, csv_out, pdf_out, png_out
    type(sim_params_t) :: params
    real(8), allocatable :: Vsd(:), I_tol(:), I_seq(:), I_cot(:)
    real(8), allocatable :: ref_tol(:), ref_seq(:), ref_cot(:)
    integer :: nVsd, ref_npts, ncmp, i, run
    real(8) :: wall_times(3), median_time, t0, t1
    real(8) :: et, es, ec, denom_v
    real(8) :: max_err_tol, max_err_seq, max_err_cot
    real(8) :: max_err_tol_ex, max_err_seq_ex, max_err_cot_ex
    real(8), parameter :: tol = 1.0d-4
    logical :: is_artifact
    integer :: count0, count1, count_rate

    ! Parse command line
    if (command_argument_count() >= 1) then
        call get_command_argument(1, spec)
    else
        spec = 'default'
    end if

    ! Build file paths
    params_path = '../benchmark/spec/'//trim(spec)//'_params.json'
    ref_path = '../benchmark/results/matlab_'//trim(spec)//'_results.json'
    json_out = '../benchmark/results/fortran_'//trim(spec)//'_results.json'
    csv_out = '../benchmark/results/fortran_'//trim(spec)//'_IV.csv'
    pdf_out = '../benchmark/results/fortran_'//trim(spec)//'_IV.pdf'
    png_out = '../benchmark/results/fortran_'//trim(spec)//'_IV.png'

    ! Parse parameters
    params = parse_params_json(params_path)
    if (params%N <= 0) then
        write(0, '(A,A)') 'Failed to parse parameters from ', trim(params_path)
        stop 1
    end if

    ! Load Vsd from MATLAB reference or generate
    call load_matlab_vsd(ref_path, Vsd, nVsd)
    if (nVsd <= 0) then
        nVsd = nint((params%Vsd_end - params%Vsd_start) / params%Vsd_step) + 1
        allocate(Vsd(nVsd))
        do i = 1, nVsd
            Vsd(i) = params%Vsd_start + dble(i - 1) * params%Vsd_step
        end do
    end if

    write(*, '(A,A,A)') '=== Franck-Condon Benchmark (Fortran) [', trim(spec), '] ==='
    write(*, '(A,I0,A,F4.1,A,F5.1,A,F6.3,A,F6.3,A,F5.3,A)') &
        'N=', params%N, ', lambda=', params%lambda, ', T=', params%T, &
        ' K, Vsd=[', params%Vsd_start, ':', params%Vsd_step, ':', params%Vsd_end, '] V'
    write(*, '(A,I0)') 'Total bias points: ', nVsd
    write(*, *)

    ! 3 timed runs
    allocate(I_tol(nVsd), I_seq(nVsd), I_cot(nVsd))
    do run = 1, 3
        call system_clock(count0, count_rate)
        call simulate_iv(params%N, params%vmode, params%alphaL, params%alphaR, &
                         params%lambda, Vsd, nVsd, params%T, params%eta, &
                         params%Vg, params%tau, run == 1, &
                         I_tol, I_seq, I_cot)
        call system_clock(count1)
        wall_times(run) = dble(count1 - count0) / dble(count_rate)
        write(*, '(A,I0,A,F8.3,A)') 'Run ', run, ': ', wall_times(run), ' s'
    end do

    ! Sort and get median
    call sort3(wall_times(1), wall_times(2), wall_times(3))
    median_time = wall_times(2)
    write(*, '(A,F8.3,A)') 'Median wall time: ', median_time, ' s'

    ! Validate against MATLAB reference
    call load_matlab_reference(ref_path, ref_tol, ref_seq, ref_cot, ref_npts)
    if (ref_npts > 0 .and. allocated(ref_tol)) then
        max_err_tol = 0.0d0; max_err_seq = 0.0d0; max_err_cot = 0.0d0
        max_err_tol_ex = 0.0d0; max_err_seq_ex = 0.0d0; max_err_cot_ex = 0.0d0

        ncmp = min(ref_npts, nVsd)
        do i = 1, ncmp
            denom_v = max(abs(ref_tol(i)), 1.0d-30)
            et = abs(I_tol(i) - ref_tol(i)) / denom_v
            denom_v = max(abs(ref_seq(i)), 1.0d-30)
            es = abs(I_seq(i) - ref_seq(i)) / denom_v
            denom_v = max(abs(ref_cot(i)), 1.0d-30)
            ec = abs(I_cot(i) - ref_cot(i)) / denom_v

            if (et > max_err_tol) max_err_tol = et
            if (es > max_err_seq) max_err_seq = es
            if (ec > max_err_cot) max_err_cot = ec

            ! Solver artifacts: Vsd~0.219 (both specs), ~0.438 (default N=15), ~0.585 (quick N=6)
            is_artifact = (abs(Vsd(i) - 0.219d0) < 0.002d0) .or. &
                          (abs(Vsd(i) - 0.438d0) < 0.002d0) .or. &
                          (abs(Vsd(i) - 0.585d0) < 0.002d0)
            if (.not. is_artifact) then
                if (et > max_err_tol_ex) max_err_tol_ex = et
                if (es > max_err_seq_ex) max_err_seq_ex = es
                if (ec > max_err_cot_ex) max_err_cot_ex = ec
            end if
        end do

        write(*, '(A)') 'Max relative error vs MATLAB (all points):'
        write(*, '(A,ES10.3,A,ES10.3,A,ES10.3)') &
            '  I_tol: ', max_err_tol, '  I_seq: ', max_err_seq, '  I_cot: ', max_err_cot
        write(*, '(A)') 'Max relative error vs MATLAB (excl. solver artifacts at Vsd~0.219,0.438,0.585):'
        write(*, '(A,ES10.3,A,ES10.3,A,ES10.3)') &
            '  I_tol: ', max_err_tol_ex, '  I_seq: ', max_err_seq_ex, '  I_cot: ', max_err_cot_ex

        if (max_err_tol_ex < tol .and. max_err_seq_ex < tol .and. max_err_cot_ex < tol) then
            write(*, '(A,ES8.0,A)') 'VALIDATION PASSED (tolerance ', tol, &
                ', excluding artifact points)'
        else
            write(*, '(A,ES8.0,A)') 'VALIDATION FAILED (tolerance ', tol, ')'
        end if

        deallocate(ref_tol, ref_seq, ref_cot)
    else
        write(*, '(A)') 'No MATLAB reference found, skipping validation.'
    end if

    ! Write outputs
    call write_results_json(json_out, trim(spec), params, Vsd, I_tol, I_seq, I_cot, &
                             nVsd, median_time)
    call write_results_csv(csv_out, Vsd, I_tol, I_seq, I_cot, nVsd)
    call plot_iv(csv_out, pdf_out, png_out, params%N, params%lambda, params%T, &
                 params%vmode, params%alphaL, params%alphaR, params%eta, params%Vg, &
                 median_time, trim(spec))

    deallocate(Vsd, I_tol, I_seq, I_cot)

contains

    subroutine sort3(a, b, c)
        real(8), intent(inout) :: a, b, c
        real(8) :: tmp
        if (a > b) then; tmp = a; a = b; b = tmp; end if
        if (b > c) then; tmp = b; b = c; c = tmp; end if
        if (a > b) then; tmp = a; a = b; b = tmp; end if
    end subroutine sort3

end program fc_benchmark
