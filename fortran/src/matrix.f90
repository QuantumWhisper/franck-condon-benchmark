module fc_matrix_w
    use fc_constants, only: KB_EV
    use fc_rate, only: rate_store_t, rateW_lead_from_store
    implicit none
contains

    pure function peq_(q, vmode, T) result(p)
        integer, intent(in) :: q
        real(8), intent(in) :: vmode, T
        real(8) :: p
        real(8) :: beta
        beta = 1.0d0 / (KB_EV * T)
        p = exp(-dble(q) * vmode * beta) * (1.0d0 - exp(-vmode * beta))
    end function peq_

    function sigma_W(store, N, n1, n2, q1) result(s)
        type(rate_store_t), intent(in) :: store
        integer, intent(in) :: N, n1, n2, q1
        real(8) :: s
        integer :: q2
        s = 0.0d0
        do q2 = 0, N - 1
            s = s + rateW_lead_from_store(store, n1, n2, q1, q2)
        end do
    end function sigma_W

    subroutine generate_matrix_W(store, N, vmode, T, tau, W)
        type(rate_store_t), intent(in) :: store
        integer, intent(in) :: N
        real(8), intent(in) :: vmode, T, tau
        real(8), intent(out) :: W(2*N, 2*N)
        integer :: ii, jj, modjj, dim, row
        real(8) :: inv_tau, val

        dim = 2 * N
        inv_tau = 1.0d0 / tau  ! IEEE: 1/Inf = 0

        ! Upper block: empty state (n=0), rows 1..N
        do ii = 1, N
            do jj = 1, dim
                modjj = mod(jj - 1, N)  ! 0-based phonon index
                if (jj <= N) then
                    ! Empty-to-empty transitions
                    if ((ii - 1) == modjj) then
                        val = rateW_lead_from_store(store, 0, 0, modjj, modjj) &
                            - sigma_W(store, N, 0, 0, modjj) &
                            - sigma_W(store, N, 0, 1, modjj) &
                            - inv_tau + peq_(modjj, vmode, T) * inv_tau
                    else
                        val = rateW_lead_from_store(store, 0, 0, modjj, ii - 1) &
                            + peq_(ii - 1, vmode, T) * inv_tau
                    end if
                else
                    ! Charged-to-empty transitions
                    val = rateW_lead_from_store(store, 1, 0, modjj, ii - 1)
                end if
                W(ii, jj) = val
            end do
        end do

        ! Lower block: charged state (n=1), rows N+1..2N
        do ii = 1, N
            row = N + ii
            do jj = 1, dim
                modjj = mod(jj - 1, N)  ! 0-based phonon index
                if (jj <= N) then
                    ! Empty-to-charged transitions
                    val = rateW_lead_from_store(store, 0, 1, modjj, ii - 1)
                else
                    ! Charged-to-charged transitions
                    if (modjj == (ii - 1)) then
                        val = rateW_lead_from_store(store, 1, 1, modjj, ii - 1) &
                            - sigma_W(store, N, 1, 0, ii - 1) &
                            - sigma_W(store, N, 1, 1, ii - 1) &
                            - inv_tau + peq_(ii - 1, vmode, T) * inv_tau
                    else
                        val = rateW_lead_from_store(store, 1, 1, modjj, ii - 1) &
                            + peq_(ii - 1, vmode, T) * inv_tau
                    end if
                end if
                W(row, jj) = val
            end do
        end do
    end subroutine generate_matrix_W

end module fc_matrix_w
