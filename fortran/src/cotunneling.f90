module fc_cotunneling
    use fc_fc_matrix, only: fc_cache_t, fc_cache_get
    use fc_regularized, only: regularized_I, regularized_J, regularized_J_matrix
    implicit none
contains

    pure function sanitize(x) result(s)
        real(8), intent(in) :: x
        real(8) :: s
        if (x /= x .or. abs(x) > huge(1.0d0)) then
            s = 0.0d0
        else
            s = x
        end if
    end function sanitize

    function rel_diff_log10(a, b) result(r)
        real(8), intent(in) :: a, b
        real(8) :: r
        complex(8) :: z
        z = cmplx(b / a, 0.0d0, kind=8)
        r = abs(log(z)) / log(10.0d0)
    end function rel_diff_log10

    function rel_diff_log(a, b) result(r)
        real(8), intent(in) :: a, b
        real(8) :: r
        complex(8) :: z
        z = cmplx(b / a, 0.0d0, kind=8)
        r = abs(log(z))
    end function rel_diff_log

    ! ================================================================
    ! m_sumMMr: inner function for n=0->0 cotunneling single-sum
    ! ================================================================
    subroutine m_sumMMr(NN, q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                        epsilond, T, fc, out_arr)
        integer, intent(in) :: NN, q1, nq2
        integer, intent(in) :: q2_vec(nq2)
        real(8), intent(in) :: lambda, muL, muR, vmode, epsilond, T
        type(fc_cache_t), intent(inout) :: fc
        real(8), intent(out) :: out_arr(nq2)
        real(8) :: MM_sq(NN, nq2), Jr(NN, nq2), E2_vec(nq2), eps_vec(NN)
        real(8) :: fc_q1r, fc_q2r, prod, sum_val, term
        integer :: r, j

        if (NN <= 0 .or. nq2 <= 0) then
            out_arr = 0.0d0
            return
        end if

        ! Precompute energy arrays
        do j = 1, nq2
            E2_vec(j) = muR - dble(q1 - q2_vec(j)) * vmode
        end do
        do r = 1, NN
            eps_vec(r) = epsilond - dble(q1 - (r-1)) * vmode
        end do

        ! Compute FC matrix products (squared)
        do r = 1, NN
            fc_q1r = fc_cache_get(fc, q1, r-1)
            do j = 1, nq2
                fc_q2r = fc_cache_get(fc, q2_vec(j), r-1)
                prod = fc_q2r * fc_q1r
                MM_sq(r, j) = prod * prod
            end do
        end do

        ! Call regularized_J
        call regularized_J(muL, E2_vec, nq2, eps_vec, NN, T, Jr)

        ! Sum over r
        do j = 1, nq2
            sum_val = 0.0d0
            do r = 1, NN
                term = MM_sq(r, j) * Jr(r, j)
                sum_val = sum_val + sanitize(term)
            end do
            out_arr(j) = sanitize(sum_val)
        end do
    end subroutine m_sumMMr

    ! ================================================================
    ! m_sumMMr11: inner function for n=1->1 cotunneling single-sum
    ! ================================================================
    subroutine m_sumMMr11(NN, q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                          epsilond, T, fc, out_arr)
        integer, intent(in) :: NN, q1, nq2
        integer, intent(in) :: q2_vec(nq2)
        real(8), intent(in) :: lambda, muL, muR, vmode, epsilond, T
        type(fc_cache_t), intent(inout) :: fc
        real(8), intent(out) :: out_arr(nq2)
        real(8) :: MM_sq(NN, nq2), Jr(NN, nq2), E2_vec(nq2), eps_matrix(nq2, NN)
        real(8) :: fc_q1r, fc_q2r, prod, sum_val, term
        integer :: r, j

        if (NN <= 0 .or. nq2 <= 0) then
            out_arr = 0.0d0
            return
        end if

        ! Precompute E2 vector
        do j = 1, nq2
            E2_vec(j) = muR - dble(q1 - q2_vec(j)) * vmode
        end do

        ! Compute FC products and epsilon matrix
        do r = 1, NN
            fc_q1r = fc_cache_get(fc, q1, r-1)
            do j = 1, nq2
                fc_q2r = fc_cache_get(fc, q2_vec(j), r-1)
                prod = fc_q2r * fc_q1r
                MM_sq(r, j) = prod * prod
                eps_matrix(j, r) = epsilond + dble(q2_vec(j) - (r-1)) * vmode
            end do
        end do

        ! Call regularized_J_matrix
        call regularized_J_matrix(muL, E2_vec, nq2, eps_matrix, nq2, NN, T, Jr)

        ! Sum over r
        do j = 1, nq2
            sum_val = 0.0d0
            do r = 1, NN
                term = MM_sq(r, j) * Jr(r, j)
                sum_val = sum_val + sanitize(term)
            end do
            out_arr(j) = sanitize(sum_val)
        end do
    end subroutine m_sumMMr11

    ! ================================================================
    ! m_sumMMMMrs: inner function for n=0->0 cotunneling double-sum
    ! ================================================================
    function m_sumMMMMrs(NN, q1, q2, lambda, muL, muR, vmode, &
                         epsilond, T, fc) result(total)
        integer, intent(in) :: NN, q1, q2
        real(8), intent(in) :: lambda, muL, muR, vmode, epsilond, T
        type(fc_cache_t), intent(inout) :: fc
        real(8) :: total
        real(8) :: Irs(NN, NN), eps1(NN), eps2(NN)
        real(8) :: E2, fc_q2r, fc_q1r, fc_q2s, fc_q1s, mmmm, term
        integer :: r, s

        total = 0.0d0
        if (NN <= 0) return

        do r = 1, NN
            eps1(r) = epsilond - dble(q1 - (r-1)) * vmode
            eps2(r) = epsilond - dble(q1 - (r-1)) * vmode
        end do

        E2 = muR - dble(q1 - q2) * vmode
        call regularized_I(muL, E2, eps1, NN, eps2, NN, T, Irs)

        do r = 1, NN
            fc_q2r = fc_cache_get(fc, q2, r-1)
            fc_q1r = fc_cache_get(fc, q1, r-1)
            do s = 1, NN
                if (r == s) cycle
                fc_q2s = fc_cache_get(fc, q2, s-1)
                fc_q1s = fc_cache_get(fc, q1, s-1)
                mmmm = fc_q2r * fc_q1r * fc_q2s * fc_q1s
                term = mmmm * Irs(r, s)
                total = total + sanitize(term)
            end do
        end do

        total = sanitize(total)
    end function m_sumMMMMrs

    ! ================================================================
    ! m_sumMMMMrs11: inner function for n=1->1 cotunneling double-sum
    ! ================================================================
    function m_sumMMMMrs11(NN, q1, q2, lambda, muL, muR, vmode, &
                           epsilond, T, fc) result(total)
        integer, intent(in) :: NN, q1, q2
        real(8), intent(in) :: lambda, muL, muR, vmode, epsilond, T
        type(fc_cache_t), intent(inout) :: fc
        real(8) :: total
        real(8) :: Irs(NN, NN), eps1(NN), eps2(NN)
        real(8) :: E2, fc_q2r, fc_q1r, fc_q2s, fc_q1s, mmmm, term
        integer :: r, s

        total = 0.0d0
        if (NN <= 0) return

        do r = 1, NN
            eps1(r) = epsilond + dble(q2 - (r-1)) * vmode
            eps2(r) = epsilond + dble(q2 - (r-1)) * vmode
        end do

        E2 = muR - dble(q1 - q2) * vmode
        call regularized_I(muL, E2, eps1, NN, eps2, NN, T, Irs)

        do r = 1, NN
            fc_q2r = fc_cache_get(fc, q2, r-1)
            fc_q1r = fc_cache_get(fc, q1, r-1)
            do s = 1, NN
                if (r == s) cycle
                fc_q2s = fc_cache_get(fc, q2, s-1)
                fc_q1s = fc_cache_get(fc, q1, s-1)
                mmmm = fc_q2r * fc_q1r * fc_q2s * fc_q1s
                term = mmmm * Irs(r, s)
                total = total + sanitize(term)
            end do
        end do

        total = sanitize(total)
    end function m_sumMMMMrs11

    ! ================================================================
    ! sumMMr: convergence wrapper for m_sumMMr (n=0->0 single-sum)
    ! ================================================================
    subroutine sumMMr_conv(q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                           epsilond, T, fc, out_arr)
        integer, intent(in) :: q1, nq2
        integer, intent(in) :: q2_vec(nq2)
        real(8), intent(in) :: lambda, muL, muR, vmode, epsilond, T
        type(fc_cache_t), intent(inout) :: fc
        real(8), intent(out) :: out_arr(nq2)
        integer :: tempN, step, j, k, nloc
        real(8), parameter :: epsilon_conv = 1.0d-14
        real(8) :: temp_tol(nq2), temp_tol2(nq2), relative_diff(nq2)
        integer :: loc_indices(nq2), q2_subset(nq2)
        real(8) :: partial(nq2)

        if (nq2 <= 0) return

        tempN = nint(lambda**2.2d0 * 3.0d0)

        call m_sumMMr(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                      epsilond, T, fc, temp_tol)
        tempN = tempN + 5
        call m_sumMMr(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                      epsilond, T, fc, temp_tol2)

        do j = 1, nq2
            relative_diff(j) = rel_diff_log10(temp_tol(j), temp_tol2(j))
        end do

        do while (.true.)
            nloc = 0
            do j = 1, nq2
                if (relative_diff(j) > epsilon_conv) then
                    nloc = nloc + 1
                    loc_indices(nloc) = j
                end if
            end do
            if (nloc == 0) exit

            step = nint(dble(tempN) * 0.5d0)
            if (step < 10) step = 10
            if (step > 20) step = 20
            tempN = tempN + step

            temp_tol = temp_tol2

            do k = 1, nloc
                q2_subset(k) = q2_vec(loc_indices(k))
            end do

            call m_sumMMr(tempN, q1, q2_subset(1:nloc), nloc, lambda, muL, muR, &
                          vmode, epsilond, T, fc, partial(1:nloc))

            do k = 1, nloc
                temp_tol2(loc_indices(k)) = partial(k)
            end do

            do j = 1, nq2
                relative_diff(j) = rel_diff_log10(temp_tol(j), temp_tol2(j))
            end do
        end do

        do j = 1, nq2
            out_arr(j) = sanitize(temp_tol2(j))
        end do
    end subroutine sumMMr_conv

    ! ================================================================
    ! sumMMr11: convergence wrapper for m_sumMMr11 (n=1->1 single-sum)
    ! ================================================================
    subroutine sumMMr11_conv(q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                             epsilond, T, fc, out_arr)
        integer, intent(in) :: q1, nq2
        integer, intent(in) :: q2_vec(nq2)
        real(8), intent(in) :: lambda, muL, muR, vmode, epsilond, T
        type(fc_cache_t), intent(inout) :: fc
        real(8), intent(out) :: out_arr(nq2)
        integer :: tempN, step, j, k, nloc
        real(8), parameter :: epsilon_conv = 1.0d-14
        real(8) :: temp_tol(nq2), temp_tol2(nq2), relative_diff(nq2)
        integer :: loc_indices(nq2), q2_subset(nq2)
        real(8) :: partial(nq2)

        if (nq2 <= 0) return

        tempN = nint(lambda**2.2d0 * 3.0d0)

        call m_sumMMr11(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                        epsilond, T, fc, temp_tol)
        tempN = tempN + 5
        call m_sumMMr11(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                        epsilond, T, fc, temp_tol2)

        do j = 1, nq2
            relative_diff(j) = rel_diff_log10(temp_tol(j), temp_tol2(j))
        end do

        do while (.true.)
            nloc = 0
            do j = 1, nq2
                if (relative_diff(j) > epsilon_conv) then
                    nloc = nloc + 1
                    loc_indices(nloc) = j
                end if
            end do
            if (nloc == 0) exit

            step = nint(dble(tempN) * 0.5d0)
            if (step < 10) step = 10
            if (step > 20) step = 20
            tempN = tempN + step

            temp_tol = temp_tol2

            do k = 1, nloc
                q2_subset(k) = q2_vec(loc_indices(k))
            end do

            call m_sumMMr11(tempN, q1, q2_subset(1:nloc), nloc, lambda, muL, muR, &
                            vmode, epsilond, T, fc, partial(1:nloc))

            do k = 1, nloc
                temp_tol2(loc_indices(k)) = partial(k)
            end do

            do j = 1, nq2
                relative_diff(j) = rel_diff_log10(temp_tol(j), temp_tol2(j))
            end do
        end do

        do j = 1, nq2
            out_arr(j) = sanitize(temp_tol2(j))
        end do
    end subroutine sumMMr11_conv

    ! ================================================================
    ! sumMMMMrs: convergence wrapper for m_sumMMMMrs (n=0->0 double-sum)
    ! ================================================================
    subroutine sumMMMMrs_conv(q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                              epsilond, T, fc, out_arr)
        integer, intent(in) :: q1, nq2
        integer, intent(in) :: q2_vec(nq2)
        real(8), intent(in) :: lambda, muL, muR, vmode, epsilond, T
        type(fc_cache_t), intent(inout) :: fc
        real(8), intent(out) :: out_arr(nq2)
        integer :: idx, q2, tempN, step
        real(8), parameter :: epsilon_conv = 1.0d-14
        real(8) :: temp_tol_s, temp_tol2_s, rd

        do idx = 1, nq2
            q2 = q2_vec(idx)
            tempN = nint(lambda * lambda * 4.0d0)

            temp_tol_s = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, &
                                     epsilond, T, fc)
            tempN = tempN + 5
            temp_tol2_s = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, &
                                      epsilond, T, fc)

            rd = rel_diff_log(temp_tol_s, temp_tol2_s)

            do while (rd == rd .and. rd > epsilon_conv)
                temp_tol_s = temp_tol2_s
                step = nint(dble(tempN) * 0.5d0)
                if (step < 20) step = 20
                if (step > 40) step = 40
                tempN = tempN + step
                temp_tol2_s = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, &
                                          epsilond, T, fc)
                rd = rel_diff_log(temp_tol_s, temp_tol2_s)
            end do

            out_arr(idx) = sanitize(temp_tol2_s)
        end do
    end subroutine sumMMMMrs_conv

    ! ================================================================
    ! sumMMMMrs11: convergence wrapper for m_sumMMMMrs11 (n=1->1 double-sum)
    ! ================================================================
    subroutine sumMMMMrs11_conv(q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                                epsilond, T, fc, out_arr)
        integer, intent(in) :: q1, nq2
        integer, intent(in) :: q2_vec(nq2)
        real(8), intent(in) :: lambda, muL, muR, vmode, epsilond, T
        type(fc_cache_t), intent(inout) :: fc
        real(8), intent(out) :: out_arr(nq2)
        integer :: idx, q2, tempN, secondN, step
        real(8), parameter :: epsilon_conv = 1.0d-14
        real(8) :: temp_tol_s, temp_tol2_s, rd

        do idx = 1, nq2
            q2 = q2_vec(idx)
            tempN = nint(lambda * lambda * 4.0d0)

            temp_tol_s = m_sumMMMMrs11(tempN, q1, q2, lambda, muL, muR, vmode, &
                                        epsilond, T, fc)
            tempN = tempN + 5
            secondN = tempN - 1
            if (secondN < 2) secondN = 2
            temp_tol2_s = m_sumMMMMrs11(secondN, q1, q2, lambda, muL, muR, vmode, &
                                         epsilond, T, fc)

            rd = rel_diff_log(temp_tol_s, temp_tol2_s)

            do while (rd == rd .and. rd > epsilon_conv)
                step = nint(dble(tempN) * 0.5d0)
                if (step < 20) step = 20
                if (step > 40) step = 40
                tempN = tempN + step
                temp_tol_s = temp_tol2_s
                temp_tol2_s = m_sumMMMMrs11(tempN, q1, q2, lambda, muL, muR, vmode, &
                                             epsilond, T, fc)
                rd = rel_diff_log(temp_tol_s, temp_tol2_s)
            end do

            out_arr(idx) = sanitize(temp_tol2_s)
        end do
    end subroutine sumMMMMrs11_conv

end module fc_cotunneling
