!  TEST_NAIVE - test the naive forward substitution algorithm for a 3x3 problem

    program test_naive
    use naive
    use, intrinsic:: iso_fortran_env
    implicit none
    integer:: j, ierr
    namelist/matrixdata/ a, b

    real(DP):: a(3,3), b(3)

    read(INPUT_UNIT,matrixdata)
    call forward_naive(a, 3, b)
    if(ierr.eq.0) then
       write(OUTPUT_UNIT,*) b
    else
       write(OUTPUT_UNIT,*) 'error: ', ierr
    endif  ! or END IF
    stop
    end program test_naive
