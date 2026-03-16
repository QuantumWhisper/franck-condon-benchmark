module fc_fermi_bose
    use fc_constants, only: KB_EV
    implicit none
contains

    pure function fermi(x, mu, T) result(f)
        real(8), intent(in) :: x, mu, T
        real(8) :: f
        f = 1.0d0 / (exp((x - mu) / (KB_EV * T)) + 1.0d0)
    end function fermi

    pure function bose_fcn(x, T) result(b)
        real(8), intent(in) :: x, T
        real(8) :: b
        real(8) :: beta
        beta = 1.0d0 / (KB_EV * T)
        b = 1.0d0 / (exp(x * beta) - 1.0d0)
    end function bose_fcn

end module fc_fermi_bose
