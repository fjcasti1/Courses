program testForward
    use substitution
    use, intrinsic:: iso_fortran_env
    implicit none
    integer:: j, ierr
    namelist/matrixdataL/ L, b

    real(DP):: L(3,3), b(3)
    ierr=0
    read(INPUT_UNIT,matrixdataL)
    call forward_subs(L, 3, b, ierr)
    if(ierr.eq.0) then
       write(OUTPUT_UNIT,*) b
    else
       write(OUTPUT_UNIT,*) 'error: ', ierr
    endif  ! or END IF
    stop
end program testForward
