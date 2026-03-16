module fc_simulate
    use fc_constants, only: ELEMENTARY_CHARGE
    use fc_fc_matrix, only: fc_cache_t, fc_cache_init
    use fc_rate, only: rate_store_t, rate_store_init, rate_store_free, calculate_all_rateW
    use fc_current, only: current_from_rate_equations
    implicit none
contains

    subroutine simulate_iv(N, vmode, alphaL, alphaR, lambda, Vsd_vec, nVsd, &
                           T, eta, Vg, tau, verbose, &
                           I_tol_out, I_seq_out, I_cot_out)
        integer, intent(in) :: N, nVsd
        real(8), intent(in) :: vmode, alphaL, alphaR, lambda, T, eta, Vg, tau
        real(8), intent(in) :: Vsd_vec(nVsd)
        logical, intent(in) :: verbose
        real(8), intent(out) :: I_tol_out(nVsd), I_seq_out(nVsd), I_cot_out(nVsd)
        type(fc_cache_t) :: fc
        type(rate_store_t) :: store
        real(8) :: v, s_sign, cr_tol, cr_seq, cr_cot
        integer :: vv
        character(len=80) :: msg

        call fc_cache_init(fc, lambda)

        do vv = 1, nVsd
            v = Vsd_vec(vv)
            if (verbose) then
                write(msg, '(A,I4,A,I4,A,F7.4,A)') &
                    char(13)//'Bias point ', vv, '/', nVsd, ' (Vsd = ', v, ' V)'
                write(0, '(A)', advance='no') trim(msg)
            end if

            call rate_store_init(store, N)
            call calculate_all_rateW(store, N, vmode, alphaL, alphaR, lambda, &
                                      v, T, eta, 1, Vg, fc)
            call calculate_all_rateW(store, N, vmode, alphaL, alphaR, lambda, &
                                      v, T, eta, -1, Vg, fc)

            call current_from_rate_equations(v, N, vmode, alphaL, alphaR, lambda, &
                                             T, eta, Vg, tau, store, &
                                             cr_tol, cr_seq, cr_cot)

            if (v /= 0.0d0) then
                if (v > 0.0d0) then
                    s_sign = -1.0d0
                else
                    s_sign = 1.0d0
                end if
                cr_tol = cr_tol * s_sign
                cr_seq = cr_seq * s_sign
                cr_cot = cr_cot * s_sign
            end if

            I_tol_out(vv) = cr_tol * ELEMENTARY_CHARGE
            I_seq_out(vv) = cr_seq * ELEMENTARY_CHARGE
            I_cot_out(vv) = cr_cot * ELEMENTARY_CHARGE

            call rate_store_free(store)
        end do

        if (verbose) then
            write(0, *)
        end if
    end subroutine simulate_iv

end module fc_simulate
