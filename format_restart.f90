program format_restart

implicit none

character(len=1) :: in_name, option
!character(len=32) :: out_name

real*8,     dimension(:,:,:), allocatable   :: KS_eigenvalue1
real*8,     dimension(:,:,:,:), allocatable :: KS_eigenvector1

real*8,     dimension(:,:,:), allocatable   :: occ_numbers1
integer :: n_spin1 = 2

integer :: i_basis1, i_states1, i_spin1, n_basis1, n_states1, n_test, line_count

call getarg(1, in_name)
call getarg(2, option)

if (option.eq."i") then

    open(file = "info.f"//in_name, unit = 25, status = 'old', form = 'formatted')

    do line_count = 1, 3
        read(25,*) ! move line pointer to actual line with information 
    end do
        read(25,*) n_basis1, n_states1, n_spin1
    close(unit=25)

    allocate(KS_eigenvector1(n_basis1,n_states1,n_spin1,1))
    allocate(KS_eigenvalue1(n_states1,n_spin1,1))
    allocate(occ_numbers1(n_states1,n_spin1,1))

    open(file = "restart.f"//in_name, unit = 7, status = 'old', form = 'unformatted')

    read(7) n_test

    do i_basis1 = 1, n_basis1
        do i_states1 = 1, n_states1
            read(7) (KS_eigenvector1(i_basis1,i_states1,i_spin1,1),i_spin1 = 1, n_spin1)
        end do
    end do

    do i_states1 = 1, n_states1
        do i_spin1 = 1, n_spin1
            read(7) KS_eigenvalue1(i_states1,i_spin1,1), occ_numbers1(i_states1,i_spin1,1)
        end do
    end do

    close(unit=7)


    ! now write formatted output
    open(file = "restart.form.f"//in_name, unit = 8, status = 'new', form = 'formatted')

    write(8,*) n_test

    do i_basis1 = 1, n_basis1
        do i_states1 = 1, n_states1
            write(8,*) (KS_eigenvector1(i_basis1,i_states1,i_spin1,1),i_spin1 = 1, n_spin1)
        end do
    end do

    do i_states1 = 1, n_states1
        do i_spin1 = 1, n_spin1
            write(8,*) KS_eigenvalue1(i_states1,i_spin1,1), occ_numbers1(i_states1,i_spin1,1)
        end do
    end do

elseif (option.eq."o") then
    open(file = "info.f"//in_name, unit = 25, status = 'old', form = 'formatted')

    do line_count = 1, 3
        read(25,*) ! move line pointer to actual line with information 
    end do
        read(25,*) n_basis1, n_states1, n_spin1
    close(unit=25)

    allocate(KS_eigenvector1(n_basis1,n_states1,n_spin1,1))
    allocate(KS_eigenvalue1(n_states1,n_spin1,1))
    allocate(occ_numbers1(n_states1,n_spin1,1))

    open(file = "restart.form.f"//in_name, unit = 7, status = 'old', form = 'formatted')

    read(7,*) n_test

    do i_basis1 = 1, n_basis1
        do i_states1 = 1, n_states1
            read(7,*) (KS_eigenvector1(i_basis1,i_states1,i_spin1,1),i_spin1 = 1, n_spin1)
        end do
    end do

    do i_states1 = 1, n_states1
        do i_spin1 = 1, n_spin1
            read(7,*) KS_eigenvalue1(i_states1,i_spin1,1), occ_numbers1(i_states1,i_spin1,1)
        end do
    end do

    close(unit=7)


    ! now write unformatted output
    open(file = "restart.f"//in_name, unit = 8, status = 'new', form = 'unformatted')

    write(8) n_test

    do i_basis1 = 1, n_basis1
        do i_states1 = 1, n_states1
            write(8) (KS_eigenvector1(i_basis1,i_states1,i_spin1,1),i_spin1 = 1, n_spin1)
        end do
    end do

    do i_states1 = 1, n_states1
        do i_spin1 = 1, n_spin1
            write(8) KS_eigenvalue1(i_states1,i_spin1,1), occ_numbers1(i_states1,i_spin1,1)
        end do
    end do
end if

end program format_restart
