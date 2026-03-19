module fc_digamma
    use fc_constants, only: PI
    implicit none

    ! Even Bernoulli numbers B_2, B_4, ..., B_20 (10 terms)
    ! Used directly by trigamma asymptotic series
    real(8), parameter :: BERNOULLI_EVEN(10) = (/ &
        1.0d0 / 6.0d0, &
        -1.0d0 / 30.0d0, &
        1.0d0 / 42.0d0, &
        -1.0d0 / 30.0d0, &
        5.0d0 / 66.0d0, &
        -691.0d0 / 2730.0d0, &
        7.0d0 / 6.0d0, &
        -3617.0d0 / 510.0d0, &
        43867.0d0 / 798.0d0, &
        -174611.0d0 / 330.0d0 /)

    ! Pre-computed digamma coefficients: B_{2(k+1)} / (2*(k+1)) for k=0..9
    real(8), parameter :: DIGAMMA_COEFF(10) = (/ &
        1.0d0 / 6.0d0 / 2.0d0, &
        -1.0d0 / 30.0d0 / 4.0d0, &
        1.0d0 / 42.0d0 / 6.0d0, &
        -1.0d0 / 30.0d0 / 8.0d0, &
        5.0d0 / 66.0d0 / 10.0d0, &
        -691.0d0 / 2730.0d0 / 12.0d0, &
        7.0d0 / 6.0d0 / 14.0d0, &
        -3617.0d0 / 510.0d0 / 16.0d0, &
        43867.0d0 / 798.0d0 / 18.0d0, &
        -174611.0d0 / 330.0d0 / 20.0d0 /)

contains

    ! Fast-path digamma: 5-term asymptotic, no reflection/recurrence
    ! Valid for Re(z) > 0 and |z|^2 > 900 (>99% of calls at T=4.2K)
    function digamma_asymptotic5(z) result(res)
        complex(8), intent(in) :: z
        complex(8) :: res
        complex(8) :: inv_z, inv_z_sq, inv_power
        integer :: k

        inv_z = (1.0d0, 0.0d0) / z
        inv_z_sq = inv_z * inv_z
        res = log(z) - inv_z * 0.5d0
        inv_power = inv_z_sq
        do k = 1, 5
            res = res - inv_power * DIGAMMA_COEFF(k)
            inv_power = inv_power * inv_z_sq
        end do
    end function digamma_asymptotic5

    ! Fast-path trigamma: 5-term asymptotic, no reflection/recurrence
    ! Valid for Re(z) > 0 and |z|^2 > 900
    function trigamma_asymptotic5(z) result(res)
        complex(8), intent(in) :: z
        complex(8) :: res
        complex(8) :: iz, iz2, power
        integer :: k

        iz = (1.0d0, 0.0d0) / z
        iz2 = iz * iz
        res = iz + iz2 * 0.5d0
        power = iz2 * iz
        do k = 1, 5
            res = res + power * BERNOULLI_EVEN(k)
            power = power * iz2
        end do
    end function trigamma_asymptotic5

    ! Dispatching digamma: fast-path if |z|^2 > 900 and Re(z) > 0, else full
    function fast_digamma(z) result(res)
        complex(8), intent(in) :: z
        complex(8) :: res
        real(8) :: norm_sq

        norm_sq = real(z)**2 + aimag(z)**2
        if (real(z) > 0.0d0 .and. norm_sq > 900.0d0) then
            res = digamma_asymptotic5(z)
        else
            res = digamma_c(z)
        end if
    end function fast_digamma

    ! Dispatching trigamma: fast-path if |z|^2 > 900 and Re(z) > 0, else full
    function fast_trigamma(z) result(res)
        complex(8), intent(in) :: z
        complex(8) :: res
        real(8) :: norm_sq

        norm_sq = real(z)**2 + aimag(z)**2
        if (real(z) > 0.0d0 .and. norm_sq > 900.0d0) then
            res = trigamma_asymptotic5(z)
        else
            res = trigamma_c(z)
        end if
    end function fast_trigamma

    ! Complex digamma function psi(z) via asymptotic series
    ! Optimized: 10-term Bernoulli, threshold |z|>=10, multiply-instead-of-divide
    function digamma_c(z) result(res)
        complex(8), intent(in) :: z
        complex(8) :: res
        complex(8) :: z_work, inv_z, inv_z_sq, inv_power, pi_z
        logical :: reflection
        integer :: k
        real(8) :: norm_sq

        z_work = z
        res = cmplx(0.0d0, 0.0d0, kind=8)
        reflection = .false.

        ! Step 1: Reflection formula for Re(z) <= 0
        if (real(z_work) <= 0.0d0) then
            reflection = .true.
            z_work = cmplx(1.0d0 - real(z_work), -aimag(z_work), kind=8)
        end if

        ! Step 2: Recurrence shift until |z|^2 >= 100
        norm_sq = real(z_work)**2 + aimag(z_work)**2
        do while (norm_sq < 100.0d0)
            res = res - (1.0d0, 0.0d0) / z_work
            z_work = z_work + (1.0d0, 0.0d0)
            norm_sq = real(z_work)**2 + aimag(z_work)**2
        end do

        ! Step 3: Asymptotic expansion
        inv_z = (1.0d0, 0.0d0) / z_work
        res = res + log(z_work) - inv_z * 0.5d0

        inv_z_sq = inv_z * inv_z
        inv_power = inv_z_sq
        do k = 1, 10
            res = res - inv_power * DIGAMMA_COEFF(k)
            inv_power = inv_power * inv_z_sq
        end do

        ! Step 4: Reflection correction
        if (reflection) then
            pi_z = z * PI
            res = res - cos(pi_z) / sin(pi_z) * PI
        end if
    end function digamma_c

    ! Complex trigamma function psi'(z) via asymptotic series
    ! 10-term Bernoulli, threshold |z|>=10
    function trigamma_c(z) result(res)
        complex(8), intent(in) :: z
        complex(8) :: res
        complex(8) :: z_work, iz, iz2, power, sinval, ratio
        logical :: reflection
        integer :: k
        real(8) :: norm_sq

        z_work = z
        res = cmplx(0.0d0, 0.0d0, kind=8)
        reflection = .false.

        ! Step 1: Reflection formula for Re(z) <= 0
        if (real(z_work) <= 0.0d0) then
            reflection = .true.
            z_work = cmplx(1.0d0 - real(z_work), -aimag(z_work), kind=8)
        end if

        ! Step 2: Recurrence shift until |z|^2 >= 100
        norm_sq = real(z_work)**2 + aimag(z_work)**2
        do while (norm_sq < 100.0d0)
            res = res + (1.0d0, 0.0d0) / (z_work * z_work)
            z_work = z_work + (1.0d0, 0.0d0)
            norm_sq = real(z_work)**2 + aimag(z_work)**2
        end do

        ! Step 3: Asymptotic expansion
        iz = (1.0d0, 0.0d0) / z_work
        iz2 = iz * iz
        res = res + iz + iz2 * 0.5d0

        power = iz2 * iz
        do k = 1, 10
            res = res + power * BERNOULLI_EVEN(k)
            power = power * iz2
        end do

        ! Step 4: Reflection correction
        if (reflection) then
            sinval = sin(z * PI)
            ratio = (1.0d0, 0.0d0) / sinval * PI
            res = ratio * ratio - res
        end if
    end function trigamma_c

end module fc_digamma
