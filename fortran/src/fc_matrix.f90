module fc_fc_matrix
    use fc_laguerre, only: laguerre_L
    implicit none

    integer, parameter :: FC_MAX_N = 256

    type :: fc_cache_t
        real(8) :: cache(0:FC_MAX_N-1, 0:FC_MAX_N-1)
        logical :: valid(0:FC_MAX_N-1, 0:FC_MAX_N-1)
        real(8) :: lambda
    end type fc_cache_t

contains

    subroutine fc_cache_init(fc, lambda)
        type(fc_cache_t), intent(out) :: fc
        real(8), intent(in) :: lambda
        fc%lambda = lambda
        fc%cache = 0.0d0
        fc%valid = .false.
    end subroutine fc_cache_init

    pure function fc_matrix_single(q1, q2, lambda) result(M)
        integer, intent(in) :: q1, q2
        real(8), intent(in) :: lambda
        real(8) :: M
        integer :: q, Q_val
        real(8) :: sign_factor, L_val, log_coeff

        q = min(q1, q2)
        Q_val = max(q1, q2)

        if (q1 == 0 .and. q2 == 0) then
            M = exp(-lambda * lambda / 2.0d0)
            return
        end if

        sign_factor = 1.0d0
        if (q2 < q1 .and. mod(q1 - q2, 2) /= 0) then
            sign_factor = -1.0d0
        end if

        L_val = laguerre_L(q, Q_val - q, lambda * lambda)

        if (lambda == 0.0d0) then
            if (q1 == q2) then
                M = L_val
            else
                M = 0.0d0
            end if
            return
        end if

        log_coeff = dble(Q_val - q) * log(lambda) &
                  - lambda * lambda / 2.0d0 &
                  + 0.5d0 * (log_gamma(dble(q) + 1.0d0) - log_gamma(dble(Q_val) + 1.0d0))

        M = sign_factor * exp(log_coeff) * L_val
        if (isnan(M) .or. abs(M) > huge(1.0d0)) then
            M = 0.0d0
        end if
    end function fc_matrix_single

    ! Pre-populate all FC cache entries up to max_q (thread-safe after this call)
    subroutine fc_cache_populate(fc, max_q)
        type(fc_cache_t), intent(inout) :: fc
        integer, intent(in) :: max_q
        integer :: q1, q2, n

        n = min(max_q, FC_MAX_N)
        do q2 = 0, n - 1
            do q1 = 0, n - 1
                if (.not. fc%valid(q1, q2)) then
                    fc%cache(q1, q2) = fc_matrix_single(q1, q2, fc%lambda)
                    fc%valid(q1, q2) = .true.
                end if
            end do
        end do
    end subroutine fc_cache_populate

    function fc_cache_get(fc, q1, q2) result(val)
        type(fc_cache_t), intent(inout) :: fc
        integer, intent(in) :: q1, q2
        real(8) :: val

        if (q1 < 0 .or. q2 < 0 .or. q1 >= FC_MAX_N .or. q2 >= FC_MAX_N) then
            val = 0.0d0
            return
        end if

        if (.not. fc%valid(q1, q2)) then
            fc%cache(q1, q2) = fc_matrix_single(q1, q2, fc%lambda)
            fc%valid(q1, q2) = .true.
        end if
        val = fc%cache(q1, q2)
    end function fc_cache_get

end module fc_fc_matrix
