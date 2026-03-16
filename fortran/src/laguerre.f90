module fc_laguerre
    implicit none
contains

    pure function laguerre_L(n, alpha, x) result(L)
        integer, intent(in) :: n, alpha
        real(8), intent(in) :: x
        real(8) :: L
        real(8) :: L_prev, L_next
        integer :: k

        if (n <= 0) then
            L = 1.0d0
            return
        end if

        L_prev = 1.0d0
        L = 1.0d0 + dble(alpha) - x
        if (n == 1) return

        do k = 1, n - 1
            L_next = ((2.0d0 * dble(k) + 1.0d0 + dble(alpha) - x) * L &
                     - (dble(k) + dble(alpha)) * L_prev) / (dble(k) + 1.0d0)
            L_prev = L
            L = L_next
        end do
    end function laguerre_L

end module fc_laguerre
