module fc_json_io
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    implicit none

    type :: sim_params_t
        integer :: N
        real(8) :: vmode, alphaL, alphaR, lambda, T, eta, Vg, tau
        real(8) :: Vsd_start, Vsd_end, Vsd_step
    end type sim_params_t

contains

    ! Extract a numeric value after a key in JSON text
    function json_get_number(text, key, default_val) result(val)
        character(len=*), intent(in) :: text, key
        real(8), intent(in) :: default_val
        real(8) :: val
        integer :: pos, colon_pos, end_pos, ios
        character(len=256) :: num_str

        val = default_val
        pos = index(text, '"'//trim(key)//'"')
        if (pos == 0) return

        colon_pos = index(text(pos:), ':')
        if (colon_pos == 0) return
        colon_pos = pos + colon_pos  ! absolute position after ':'

        ! Find start of number (skip whitespace)
        do while (colon_pos <= len(text) .and. &
                  (text(colon_pos:colon_pos) == ' ' .or. text(colon_pos:colon_pos) == char(9)))
            colon_pos = colon_pos + 1
        end do

        ! Check for string value "Inf"
        if (text(colon_pos:colon_pos) == '"') then
            if (index(text(colon_pos:), 'Inf') > 0) then
                val = ieee_value(1.0d0, ieee_positive_inf)
            end if
            return
        end if

        ! Find end of number
        end_pos = colon_pos
        do while (end_pos <= len(text))
            if (text(end_pos:end_pos) == ',' .or. text(end_pos:end_pos) == '}' .or. &
                text(end_pos:end_pos) == char(10) .or. text(end_pos:end_pos) == char(13)) exit
            end_pos = end_pos + 1
        end do

        num_str = adjustl(text(colon_pos:end_pos-1))
        read(num_str, *, iostat=ios) val
        if (ios /= 0) val = default_val
    end function json_get_number

    ! Extract an integer value after a key
    function json_get_int(text, key, default_val) result(val)
        character(len=*), intent(in) :: text, key
        integer, intent(in) :: default_val
        integer :: val
        val = nint(json_get_number(text, key, dble(default_val)))
    end function json_get_int

    ! Parse benchmark spec JSON file
    function parse_params_json(filepath) result(params)
        character(len=*), intent(in) :: filepath
        type(sim_params_t) :: params
        character(len=:), allocatable :: text
        integer :: unit_num, fsize, ios

        params%N = 0
        params%tau = huge(1.0d0)

        ! Read file into string
        open(newunit=unit_num, file=filepath, status='old', access='stream', &
             form='unformatted', iostat=ios)
        if (ios /= 0) return

        inquire(unit=unit_num, size=fsize)
        allocate(character(len=fsize) :: text)
        read(unit_num, iostat=ios) text
        close(unit_num)
        if (ios /= 0) then
            deallocate(text)
            return
        end if

        params%N = json_get_int(text, 'N', 0)
        params%vmode = json_get_number(text, 'vmode', 0.0d0)
        params%alphaL = json_get_number(text, 'alphaL', 0.0d0)
        params%alphaR = json_get_number(text, 'alphaR', 0.0d0)
        params%lambda = json_get_number(text, 'lambda', 0.0d0)
        params%T = json_get_number(text, 'T', 0.0d0)
        params%eta = json_get_number(text, 'eta', 0.0d0)
        params%Vg = json_get_number(text, 'Vg', 0.0d0)
        params%tau = json_get_number(text, 'tau', huge(1.0d0))
        params%Vsd_start = json_get_number(text, 'Vsd_start', 0.0d0)
        params%Vsd_end = json_get_number(text, 'Vsd_end', 0.0d0)
        params%Vsd_step = json_get_number(text, 'Vsd_step', 0.0d0)

        deallocate(text)
    end function parse_params_json

    ! Load Vsd array from MATLAB reference JSON
    subroutine load_matlab_vsd(filepath, Vsd, nVsd)
        character(len=*), intent(in) :: filepath
        real(8), allocatable, intent(out) :: Vsd(:)
        integer, intent(out) :: nVsd
        character(len=:), allocatable :: text
        integer :: unit_num, fsize, ios, pos, start_pos, end_pos, i, cnt
        real(8) :: val
        character(len=64) :: num_str

        nVsd = 0
        open(newunit=unit_num, file=filepath, status='old', access='stream', &
             form='unformatted', iostat=ios)
        if (ios /= 0) return

        inquire(unit=unit_num, size=fsize)
        allocate(character(len=fsize) :: text)
        read(unit_num, iostat=ios) text
        close(unit_num)
        if (ios /= 0) then
            deallocate(text)
            return
        end if

        ! Find "Vsd": [
        pos = index(text, '"Vsd"')
        if (pos == 0) then
            deallocate(text)
            return
        end if
        pos = index(text(pos:), '[')
        if (pos == 0) then
            deallocate(text)
            return
        end if
        ! Adjust to absolute position
        start_pos = index(text, '"Vsd"')
        pos = start_pos + index(text(start_pos:), '[') - 1

        ! Count elements
        cnt = 0
        i = pos + 1
        do while (i <= len(text))
            if (text(i:i) == ']') exit
            if (text(i:i) == ',' .or. text(i:i) == '[') then
                i = i + 1
                cycle
            end if
            if (text(i:i) == ' ' .or. text(i:i) == char(10) .or. &
                text(i:i) == char(13) .or. text(i:i) == char(9)) then
                i = i + 1
                cycle
            end if
            ! Found start of a number
            cnt = cnt + 1
            ! Skip to end of number
            do while (i <= len(text))
                if (text(i:i) == ',' .or. text(i:i) == ']' .or. &
                    text(i:i) == ' ' .or. text(i:i) == char(10) .or. &
                    text(i:i) == char(13)) exit
                i = i + 1
            end do
        end do

        if (cnt == 0) then
            deallocate(text)
            return
        end if

        nVsd = cnt
        allocate(Vsd(nVsd))

        ! Parse elements
        cnt = 0
        i = pos + 1
        do while (i <= len(text) .and. cnt < nVsd)
            if (text(i:i) == ']') exit
            if (text(i:i) == ',' .or. text(i:i) == ' ' .or. &
                text(i:i) == char(10) .or. text(i:i) == char(13) .or. &
                text(i:i) == char(9)) then
                i = i + 1
                cycle
            end if
            ! Found start of number
            start_pos = i
            do while (i <= len(text))
                if (text(i:i) == ',' .or. text(i:i) == ']' .or. &
                    text(i:i) == ' ' .or. text(i:i) == char(10) .or. &
                    text(i:i) == char(13)) exit
                i = i + 1
            end do
            end_pos = i - 1
            num_str = text(start_pos:end_pos)
            cnt = cnt + 1
            read(num_str, *, iostat=ios) val
            if (ios == 0) then
                Vsd(cnt) = val
            else
                Vsd(cnt) = 0.0d0
            end if
        end do

        deallocate(text)
    end subroutine load_matlab_vsd

    ! Load MATLAB reference I-V arrays
    subroutine load_matlab_reference(filepath, ref_tol, ref_seq, ref_cot, npts)
        character(len=*), intent(in) :: filepath
        real(8), allocatable, intent(out) :: ref_tol(:), ref_seq(:), ref_cot(:)
        integer, intent(out) :: npts
        ! Use the same array parser for each field
        real(8), allocatable :: tmp(:)
        integer :: n1, n2, n3

        npts = 0
        call load_json_array(filepath, 'I_tol', ref_tol, n1)
        call load_json_array(filepath, 'I_seq', ref_seq, n2)
        call load_json_array(filepath, 'I_cot', ref_cot, n3)

        npts = min(n1, n2, n3)
    end subroutine load_matlab_reference

    ! Generic JSON array parser for a named field
    subroutine load_json_array(filepath, field_name, arr, npts)
        character(len=*), intent(in) :: filepath, field_name
        real(8), allocatable, intent(out) :: arr(:)
        integer, intent(out) :: npts
        character(len=:), allocatable :: text
        integer :: unit_num, fsize, ios, pos, start_pos, i, cnt, end_pos
        character(len=64) :: num_str
        real(8) :: val

        npts = 0
        open(newunit=unit_num, file=filepath, status='old', access='stream', &
             form='unformatted', iostat=ios)
        if (ios /= 0) return

        inquire(unit=unit_num, size=fsize)
        allocate(character(len=fsize) :: text)
        read(unit_num, iostat=ios) text
        close(unit_num)
        if (ios /= 0) then
            deallocate(text)
            return
        end if

        pos = index(text, '"'//trim(field_name)//'"')
        if (pos == 0) then
            deallocate(text)
            return
        end if
        pos = pos + index(text(pos:), '[') - 1

        ! Count elements
        cnt = 0
        i = pos + 1
        do while (i <= len(text))
            if (text(i:i) == ']') exit
            if (text(i:i) == ',' .or. text(i:i) == ' ' .or. &
                text(i:i) == char(10) .or. text(i:i) == char(13) .or. &
                text(i:i) == char(9)) then
                i = i + 1
                cycle
            end if
            cnt = cnt + 1
            do while (i <= len(text))
                if (text(i:i) == ',' .or. text(i:i) == ']' .or. &
                    text(i:i) == ' ' .or. text(i:i) == char(10) .or. &
                    text(i:i) == char(13)) exit
                i = i + 1
            end do
        end do

        if (cnt == 0) then
            deallocate(text)
            return
        end if

        npts = cnt
        allocate(arr(npts))

        cnt = 0
        i = pos + 1
        do while (i <= len(text) .and. cnt < npts)
            if (text(i:i) == ']') exit
            if (text(i:i) == ',' .or. text(i:i) == ' ' .or. &
                text(i:i) == char(10) .or. text(i:i) == char(13) .or. &
                text(i:i) == char(9)) then
                i = i + 1
                cycle
            end if
            start_pos = i
            do while (i <= len(text))
                if (text(i:i) == ',' .or. text(i:i) == ']' .or. &
                    text(i:i) == ' ' .or. text(i:i) == char(10) .or. &
                    text(i:i) == char(13)) exit
                i = i + 1
            end do
            end_pos = i - 1
            cnt = cnt + 1
            num_str = text(start_pos:end_pos)
            read(num_str, *, iostat=ios) val
            if (ios == 0) then
                arr(cnt) = val
            else
                arr(cnt) = 0.0d0
            end if
        end do

        deallocate(text)
    end subroutine load_json_array

    ! Write JSON results file
    subroutine write_results_json(filepath, spec, params, Vsd, I_tol, I_seq, I_cot, &
                                   nVsd, wall_time)
        character(len=*), intent(in) :: filepath, spec
        type(sim_params_t), intent(in) :: params
        integer, intent(in) :: nVsd
        real(8), intent(in) :: Vsd(nVsd), I_tol(nVsd), I_seq(nVsd), I_cot(nVsd)
        real(8), intent(in) :: wall_time
        integer :: unit_num, ios, i
        character(len=20) :: timestamp
        integer :: dt(8)

        call date_and_time(values=dt)
        write(timestamp, '(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') &
            dt(1), '-', dt(2), '-', dt(3), ' ', dt(5), ':', dt(6), ':', dt(7)

        open(newunit=unit_num, file=filepath, status='replace', iostat=ios)
        if (ios /= 0) return

        write(unit_num, '(A)') '{'
        write(unit_num, '(A)') '  "spec": "'//trim(spec)//'",'
        write(unit_num, '(A)') '  "language": "Fortran",'
        write(unit_num, '(A)') '  "version": "gfortran",'
        write(unit_num, '(A,ES23.17,A)') '  "wall_time_seconds": ', wall_time, ','
        write(unit_num, '(A)') '  "parameters": {'
        write(unit_num, '(A,I0,A)') '    "N": ', params%N, ','
        write(unit_num, '(A,ES23.17,A)') '    "vmode": ', params%vmode, ','
        write(unit_num, '(A,ES23.17,A)') '    "alphaL": ', params%alphaL, ','
        write(unit_num, '(A,ES23.17,A)') '    "alphaR": ', params%alphaR, ','
        write(unit_num, '(A,ES23.17,A)') '    "lambda": ', params%lambda, ','
        write(unit_num, '(A,ES23.17,A)') '    "T": ', params%T, ','
        write(unit_num, '(A,ES23.17,A)') '    "eta": ', params%eta, ','
        write(unit_num, '(A,ES23.17,A)') '    "Vg": ', params%Vg, ','
        if (params%tau > 1.0d30) then
            write(unit_num, '(A)') '    "tau": "Inf"'
        else
            write(unit_num, '(A,ES23.17)') '    "tau": ', params%tau
        end if
        write(unit_num, '(A)') '  },'
        write(unit_num, '(A)') '  "bias_sweep": {'
        write(unit_num, '(A,ES23.17,A)') '    "Vsd_start": ', params%Vsd_start, ','
        write(unit_num, '(A,ES23.17,A)') '    "Vsd_end": ', params%Vsd_end, ','
        write(unit_num, '(A,ES23.17)') '    "Vsd_step": ', params%Vsd_step
        write(unit_num, '(A)') '  },'

        call write_json_array_f(unit_num, 'Vsd', Vsd, nVsd, .true.)
        call write_json_array_f(unit_num, 'I_tol', I_tol, nVsd, .true.)
        call write_json_array_f(unit_num, 'I_seq', I_seq, nVsd, .true.)
        call write_json_array_f(unit_num, 'I_cot', I_cot, nVsd, .true.)
        write(unit_num, '(A)') '  "timestamp": "'//trim(timestamp)//'"'
        write(unit_num, '(A)') '}'

        close(unit_num)
    end subroutine write_results_json

    subroutine write_json_array_f(unit_num, name, arr, n, comma)
        integer, intent(in) :: unit_num, n
        character(len=*), intent(in) :: name
        real(8), intent(in) :: arr(n)
        logical, intent(in) :: comma
        integer :: i

        write(unit_num, '(A)') '  "'//trim(name)//'": ['
        do i = 1, n
            if (i < n) then
                write(unit_num, '(A,ES23.17,A)') '    ', arr(i), ','
            else
                write(unit_num, '(A,ES23.17)') '    ', arr(i)
            end if
        end do
        if (comma) then
            write(unit_num, '(A)') '  ],'
        else
            write(unit_num, '(A)') '  ]'
        end if
    end subroutine write_json_array_f

    ! Write CSV results file
    subroutine write_results_csv(filepath, Vsd, I_tol, I_seq, I_cot, nVsd)
        character(len=*), intent(in) :: filepath
        integer, intent(in) :: nVsd
        real(8), intent(in) :: Vsd(nVsd), I_tol(nVsd), I_seq(nVsd), I_cot(nVsd)
        integer :: unit_num, ios, i

        open(newunit=unit_num, file=filepath, status='replace', iostat=ios)
        if (ios /= 0) return

        write(unit_num, '(A)') 'Vsd_V,I_tol_A,I_seq_A,I_cot_A'
        do i = 1, nVsd
            write(unit_num, '(ES13.6,A,ES13.6,A,ES13.6,A,ES13.6)') &
                Vsd(i), ',', I_tol(i), ',', I_seq(i), ',', I_cot(i)
        end do

        close(unit_num)
    end subroutine write_results_csv

end module fc_json_io
