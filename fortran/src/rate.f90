module fc_rate
    use fc_constants, only: HBAR_EV, PI
    use fc_fermi_bose, only: fermi
    use fc_fc_matrix, only: fc_cache_t, fc_cache_get
    use fc_cotunneling, only: sumMMr_conv, sumMMr11_conv, sumMMMMrs_conv, sumMMMMrs11_conv
    implicit none

    type :: rate_store_t
        integer :: N
        real(8), allocatable :: data(:)  ! flat array: 2*2*N*2*N
    end type rate_store_t

contains

    subroutine rate_store_init(store, N)
        type(rate_store_t), intent(out) :: store
        integer, intent(in) :: N
        store%N = N
        allocate(store%data(4 * N * 2 * N))
        store%data = 0.0d0
    end subroutine rate_store_init

    subroutine rate_store_free(store)
        type(rate_store_t), intent(inout) :: store
        if (allocated(store%data)) deallocate(store%data)
        store%N = 0
    end subroutine rate_store_free

    pure function lead_to_idx(lead) result(idx)
        integer, intent(in) :: lead
        integer :: idx
        if (lead == 1) then
            idx = 0
        else
            idx = 1
        end if
    end function lead_to_idx

    ! Flat index into rate store (0-based inputs, 1-based output for Fortran array)
    pure function rate_idx(N, n1, n2, q1, lead_idx, q2) result(idx)
        integer, intent(in) :: N, n1, n2, q1, lead_idx, q2
        integer :: idx
        idx = ((((n1 * 2 + n2) * N + q1) * 2 + lead_idx) * N + q2) + 1
    end function rate_idx

    pure function sanitize_rate(x) result(s)
        real(8), intent(in) :: x
        real(8) :: s
        if (x /= x .or. abs(x) > huge(1.0d0)) then
            s = 0.0d0
        else
            s = x
        end if
    end function sanitize_rate

    pure function spin_degeneracy(n1, n2) result(s)
        integer, intent(in) :: n1, n2
        real(8) :: s
        if (n1 == 0) then
            s = 2.0d0
        else if (n2 == 0) then
            s = 1.0d0
        else
            s = 2.0d0
        end if
    end function spin_degeneracy

    subroutine m_rateW(n1, n2, q1, q2_vec, nq2, vmode, alphaL, alphaR, lambda, &
                       Vsd, T, eta, lead, Vg, fc, out_arr)
        integer, intent(in) :: n1, n2, q1, nq2, lead
        integer, intent(in) :: q2_vec(nq2)
        real(8), intent(in) :: vmode, alphaL, alphaR, lambda, Vsd, T, eta, Vg
        type(fc_cache_t), intent(inout) :: fc
        real(8), intent(out) :: out_arr(nq2)
        real(8) :: gammaL, gammaR, epsilond, muL, muR
        real(8) :: s, gamma_val, mu_val, fc_val, f_val, val, prefac
        real(8), allocatable :: single_sum(:), double_sum(:)
        integer :: i, q2

        gammaL = alphaL * vmode
        gammaR = alphaR * vmode
        epsilond = 0.0d0 + Vg
        muL = eta * Vsd
        muR = -(1.0d0 - eta) * Vsd

        s = spin_degeneracy(n1, n2)
        if (lead == 1) then
            gamma_val = gammaL
            mu_val = muL
        else
            gamma_val = gammaR
            mu_val = muR
        end if

        if (n1 == 1 .and. n2 == 0) then
            ! n=1->0: electron leaves
            do i = 1, nq2
                q2 = q2_vec(i)
                fc_val = fc_cache_get(fc, q1, q2)
                f_val = fermi(epsilond - dble(q2 - q1) * vmode, mu_val, T)
                val = s * gamma_val / HBAR_EV * fc_val * fc_val * (1.0d0 - f_val)
                out_arr(i) = sanitize_rate(val)
            end do

        else if (n1 == 0 .and. n2 == 1) then
            ! n=0->1: electron enters
            do i = 1, nq2
                q2 = q2_vec(i)
                fc_val = fc_cache_get(fc, q1, q2)
                f_val = fermi(epsilond + dble(q2 - q1) * vmode, mu_val, T)
                val = s * gamma_val / HBAR_EV * fc_val * fc_val * f_val
                out_arr(i) = sanitize_rate(val)
            end do

        else if (n1 == 0 .and. n2 == 0) then
            ! Cotunneling n=0->0
            allocate(single_sum(nq2), double_sum(nq2))
            if (lead == 1) then
                call sumMMr_conv(q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                                 epsilond, T, fc, single_sum)
                call sumMMMMrs_conv(q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                                    epsilond, T, fc, double_sum)
            else
                call sumMMr_conv(q1, q2_vec, nq2, lambda, muR, muL, vmode, &
                                 epsilond, T, fc, single_sum)
                call sumMMMMrs_conv(q1, q2_vec, nq2, lambda, muR, muL, vmode, &
                                    epsilond, T, fc, double_sum)
            end if
            prefac = s / (2.0d0 * PI * HBAR_EV) * gammaL * gammaR
            do i = 1, nq2
                out_arr(i) = sanitize_rate(prefac * (single_sum(i) + double_sum(i)))
            end do
            deallocate(single_sum, double_sum)

        else
            ! Cotunneling n=1->1
            allocate(single_sum(nq2), double_sum(nq2))
            if (lead == 1) then
                call sumMMr11_conv(q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                                   epsilond, T, fc, single_sum)
                call sumMMMMrs11_conv(q1, q2_vec, nq2, lambda, muL, muR, vmode, &
                                      epsilond, T, fc, double_sum)
            else
                call sumMMr11_conv(q1, q2_vec, nq2, lambda, muR, muL, vmode, &
                                   epsilond, T, fc, single_sum)
                call sumMMMMrs11_conv(q1, q2_vec, nq2, lambda, muR, muL, vmode, &
                                      epsilond, T, fc, double_sum)
            end if
            prefac = s / (2.0d0 * PI * HBAR_EV) * gammaL * gammaR
            do i = 1, nq2
                out_arr(i) = sanitize_rate(prefac * (single_sum(i) + double_sum(i)))
            end do
            deallocate(single_sum, double_sum)
        end if
    end subroutine m_rateW

    subroutine calculate_all_rateW(store, N, vmode, alphaL, alphaR, lambda, &
                                    Vsd, T, eta, lead, Vg, fc)
        type(rate_store_t), intent(inout) :: store
        integer, intent(in) :: N, lead
        real(8), intent(in) :: vmode, alphaL, alphaR, lambda, Vsd, T, eta, Vg
        type(fc_cache_t), intent(inout) :: fc
        integer :: n1, n2, q1, q2, lidx
        integer, allocatable :: q2_vec(:)
        real(8), allocatable :: out_arr(:)

        allocate(q2_vec(N), out_arr(N))
        do q2 = 1, N
            q2_vec(q2) = q2 - 1  ! 0-based q2 values
        end do

        lidx = lead_to_idx(lead)
        do n1 = 0, 1
            do n2 = 0, 1
                do q1 = 0, N - 1
                    call m_rateW(n1, n2, q1, q2_vec, N, vmode, alphaL, alphaR, &
                                 lambda, Vsd, T, eta, lead, Vg, fc, out_arr)
                    do q2 = 0, N - 1
                        store%data(rate_idx(N, n1, n2, q1, lidx, q2)) = out_arr(q2 + 1)
                    end do
                end do
            end do
        end do

        deallocate(q2_vec, out_arr)
    end subroutine calculate_all_rateW

    function rateW_from_store(store, n1, n2, q1, q2, lead) result(w)
        type(rate_store_t), intent(in) :: store
        integer, intent(in) :: n1, n2, q1, q2, lead
        real(8) :: w
        integer :: lidx
        lidx = lead_to_idx(lead)
        w = store%data(rate_idx(store%N, n1, n2, q1, lidx, q2))
    end function rateW_from_store

    function rateW_lead_from_store(store, n1, n2, q1, q2) result(w)
        type(rate_store_t), intent(in) :: store
        integer, intent(in) :: n1, n2, q1, q2
        real(8) :: w
        w = rateW_from_store(store, n1, n2, q1, q2, 1) &
          + rateW_from_store(store, n1, n2, q1, q2, -1)
    end function rateW_lead_from_store

end module fc_rate
