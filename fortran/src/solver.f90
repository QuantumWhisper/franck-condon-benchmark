module fc_solver
    implicit none
contains

    subroutine solve_steady_state(W, P, N)
        integer, intent(in) :: N
        real(8), intent(in) :: W(2*N, 2*N)
        real(8), intent(out) :: P(2*N)
        integer :: dim, i, j, info
        real(8), allocatable :: A(:,:), d(:)
        integer, allocatable :: ipiv(:)
        real(8) :: total, v

        dim = 2 * N
        if (dim <= 0) return

        allocate(A(dim, dim), d(dim), ipiv(dim))

        ! Copy W to A (column-major), replace last row with normalization
        do j = 1, dim
            do i = 1, dim
                if (i == dim) then
                    A(i, j) = 1.0d0
                else
                    ! W is stored row-major in the caller: W(row, col)
                    A(i, j) = W(i, j)
                end if
            end do
        end do

        ! RHS: all zeros except last element = 1
        d = 0.0d0
        d(dim) = 1.0d0

        ! Solve via LAPACK dgesv (LU factorization)
        call dgesv(dim, 1, A, dim, ipiv, d, dim, info)

        ! Post-process: clamp negatives, renormalize
        total = 0.0d0
        do i = 1, dim
            v = d(i)
            if (v /= v .or. v < 0.0d0) then  ! NaN check: v /= v
                v = 0.0d0
            end if
            P(i) = v
            total = total + v
        end do

        if (total > 0.0d0) then
            P = P / total
        else
            P = 0.0d0
            P(dim) = 1.0d0
        end if

        deallocate(A, d, ipiv)
    end subroutine solve_steady_state

end module fc_solver
