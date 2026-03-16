module fc_current
    use fc_rate, only: rate_store_t, rateW_from_store
    use fc_matrix_w, only: generate_matrix_W
    use fc_solver, only: solve_steady_state
    implicit none
contains

    subroutine current_from_rate_equations(Vsd, N, vmode, alphaL, alphaR, lambda, &
                                           T, eta, Vg, tau, store, &
                                           I_tol, I_seq, I_cot)
        real(8), intent(in) :: Vsd, vmode, alphaL, alphaR, lambda, T, eta, Vg, tau
        integer, intent(in) :: N
        type(rate_store_t), intent(in) :: store
        real(8), intent(out) :: I_tol, I_seq, I_cot
        real(8), allocatable :: W(:,:), P(:)
        real(8) :: I_seq0, I_seq1, sum_diff, w_R, w_L, w_RL, w_LR
        integer :: dim, q1, q2, nn, leadR, leadL

        dim = 2 * N
        leadR = -1
        leadL = 1

        allocate(W(dim, dim), P(dim))

        call generate_matrix_W(store, N, vmode, T, tau, W)
        call solve_steady_state(W, P, N)

        ! I_seq0: n=0->1 (empty to charged), diffW = w_R - w_L
        I_seq0 = 0.0d0
        do q1 = 0, N - 1
            sum_diff = 0.0d0
            do q2 = 0, N - 1
                w_R = rateW_from_store(store, 0, 1, q1, q2, leadR)
                w_L = rateW_from_store(store, 0, 1, q1, q2, leadL)
                sum_diff = sum_diff + (w_R - w_L)
            end do
            I_seq0 = I_seq0 + P(q1 + 1) * sum_diff
        end do

        ! I_seq1: n=1->0 (charged to empty), diffW = w_L - w_R
        I_seq1 = 0.0d0
        do q1 = 0, N - 1
            sum_diff = 0.0d0
            do q2 = 0, N - 1
                w_R = rateW_from_store(store, 1, 0, q1, q2, leadR)
                w_L = rateW_from_store(store, 1, 0, q1, q2, leadL)
                sum_diff = sum_diff + (w_L - w_R)
            end do
            I_seq1 = I_seq1 + P(N + q1 + 1) * sum_diff
        end do

        I_seq = I_seq0 + I_seq1

        ! I_cot: cotunneling (n->n), diffW = w_RL - w_LR
        I_cot = 0.0d0
        do nn = 0, 1
            do q1 = 0, N - 1
                sum_diff = 0.0d0
                do q2 = 0, N - 1
                    w_RL = rateW_from_store(store, nn, nn, q1, q2, leadR)
                    w_LR = rateW_from_store(store, nn, nn, q1, q2, leadL)
                    sum_diff = sum_diff + (w_RL - w_LR)
                end do
                I_cot = I_cot + P(nn * N + q1 + 1) * sum_diff
            end do
        end do

        I_tol = I_seq + I_cot

        deallocate(W, P)
    end subroutine current_from_rate_equations

end module fc_current
