module fc_regularized
    use fc_constants, only: KB_EV, PI
    use fc_digamma, only: digamma_c, trigamma_c
    use fc_fermi_bose, only: bose_fcn
    implicit none
contains

    ! Regularized I integral with factored digamma precomputation
    ! Key optimization: 4*N calls instead of 4*N*N
    subroutine regularized_I(E1, E2, epsilon1, n_eps1, epsilon2, n_eps2, T, out)
        real(8), intent(in) :: E1, E2, T
        integer, intent(in) :: n_eps1, n_eps2
        real(8), intent(in) :: epsilon1(n_eps1), epsilon2(n_eps2)
        real(8), intent(out) :: out(n_eps1, n_eps2)
        real(8) :: beta, bose_val, coeff, denom, val
        complex(8) :: a1, a3, a2, a4
        complex(8), allocatable :: row_diff(:), col_diff(:)
        integer :: i, j

        if (n_eps1 <= 0 .or. n_eps2 <= 0) return

        beta = 1.0d0 / (KB_EV * T)
        bose_val = bose_fcn(E2 - E1, T)
        coeff = beta / (2.0d0 * PI)

        allocate(row_diff(n_eps1), col_diff(n_eps2))

        ! Precompute row-only digamma differences
        do i = 1, n_eps1
            a1 = cmplx(0.5d0, coeff * (E2 - epsilon1(i)), kind=8)
            a3 = cmplx(0.5d0, coeff * (E1 - epsilon1(i)), kind=8)
            row_diff(i) = digamma_c(a1) - digamma_c(a3)
        end do

        ! Precompute column-only digamma differences
        do j = 1, n_eps2
            a2 = cmplx(0.5d0, -coeff * (E2 - epsilon2(j)), kind=8)
            a4 = cmplx(0.5d0, -coeff * (E1 - epsilon2(j)), kind=8)
            col_diff(j) = digamma_c(a2) - digamma_c(a4)
        end do

        ! Compute matrix entries
        do j = 1, n_eps2
            do i = 1, n_eps1
                denom = epsilon1(i) - epsilon2(j)
                val = bose_val / denom * real(row_diff(i) - col_diff(j))
                if (val /= val .or. abs(val) > huge(1.0d0)) then
                    out(i, j) = 0.0d0
                else
                    out(i, j) = val
                end if
            end do
        end do

        deallocate(row_diff, col_diff)
    end subroutine regularized_I

    ! Regularized J integral (vector E2)
    subroutine regularized_J(E1, E2, n_E2, epsilon, n_eps, T, out)
        real(8), intent(in) :: E1, T
        integer, intent(in) :: n_E2, n_eps
        real(8), intent(in) :: E2(n_E2), epsilon(n_eps)
        real(8), intent(out) :: out(n_eps, n_E2)
        real(8) :: beta, coeff, val
        real(8), allocatable :: bose_vals(:)
        complex(8), allocatable :: trig_a2(:)
        complex(8) :: a1, a2
        integer :: i, j

        if (n_E2 <= 0 .or. n_eps <= 0) return

        beta = 1.0d0 / (KB_EV * T)
        coeff = beta / (2.0d0 * PI)

        allocate(bose_vals(n_E2), trig_a2(n_eps))

        ! Precompute Bose factors
        do j = 1, n_E2
            bose_vals(j) = bose_fcn(E2(j) - E1, T)
        end do

        ! Precompute trigamma(a2) per row
        do i = 1, n_eps
            a2 = cmplx(0.5d0, coeff * (E1 - epsilon(i)), kind=8)
            trig_a2(i) = trigamma_c(a2)
        end do

        ! Compute matrix entries
        do j = 1, n_E2
            do i = 1, n_eps
                a1 = cmplx(0.5d0, coeff * (E2(j) - epsilon(i)), kind=8)
                val = coeff * bose_vals(j) * aimag(trigamma_c(a1) - trig_a2(i))
                if (val /= val .or. abs(val) > huge(1.0d0)) then
                    out(i, j) = 0.0d0
                else
                    out(i, j) = val
                end if
            end do
        end do

        deallocate(bose_vals, trig_a2)
    end subroutine regularized_J

    ! Regularized J integral (matrix epsilon) for n=1->1 cotunneling
    subroutine regularized_J_matrix(E1, E2, n_E2, epsilon_matrix, eps_rows, eps_cols, T, out)
        real(8), intent(in) :: E1, T
        integer, intent(in) :: n_E2, eps_rows, eps_cols
        real(8), intent(in) :: E2(n_E2)
        real(8), intent(in) :: epsilon_matrix(eps_rows, eps_cols)
        real(8), intent(out) :: out(eps_cols, n_E2)
        real(8) :: beta, coeff, eps, val
        real(8), allocatable :: bose_vals(:)
        complex(8) :: a1, a2
        integer :: i, j

        if (n_E2 <= 0 .or. eps_rows <= 0 .or. eps_cols <= 0) return

        beta = 1.0d0 / (KB_EV * T)
        coeff = beta / (2.0d0 * PI)

        allocate(bose_vals(n_E2))
        do j = 1, n_E2
            bose_vals(j) = bose_fcn(E2(j) - E1, T)
        end do

        do j = 1, n_E2
            do i = 1, eps_cols
                if (j > eps_rows) then
                    out(i, j) = 0.0d0
                else
                    eps = epsilon_matrix(j, i)
                    a1 = cmplx(0.5d0, coeff * (E2(j) - eps), kind=8)
                    a2 = cmplx(0.5d0, coeff * (E1 - eps), kind=8)
                    val = coeff * bose_vals(j) * aimag(trigamma_c(a1) - trigamma_c(a2))
                    if (val /= val .or. abs(val) > huge(1.0d0)) then
                        out(i, j) = 0.0d0
                    else
                        out(i, j) = val
                    end if
                end if
            end do
        end do

        deallocate(bose_vals)
    end subroutine regularized_J_matrix

end module fc_regularized
