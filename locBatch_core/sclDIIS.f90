subroutine ctrlDIIS(mode)
    
    integer :: mode
end subroutine ctrlDIIS

subroutine initDIIS(mcf, dimn, mode)
    implicit none

    double precision, dimension(dimn, dimn, dimn, dimn) :: AO2eIND, MO2eIND
    double precision, dimension(dimn, dimn) :: mcf, orth, trans, seed, store, curt, orthPrior
    double precision, dimension(dimn) :: eivl
    double precision, parameter :: absthrsh = 1.0e-4
    integer :: i, j, k, l, m, n, dimn, iter
    logical :: isSymm

    call AOIntglVal(dimn, AO2eIND)
    call MOIntglVal(dimn, mcf, AO2eIND, MO2eIND)
    
    do i = 1, dimn
        do j = 1, dimn
            trans(j, i) = MO2eIND(j, i, i, i)
        end do
    end do

    call dgemm('T', 'N', dimn, dimn, dimn, 1.d0, trans, dimn, trans, dimn, 0.d0, seed, dimn)
    rdTrans(seed, dimn, eivl, store)
    call dgemm('N', 'N', dimn, dimn, dimn, 1.d0, trans, dimn, store, dimn, 0.d0, orthPrior, dimn)

    do while(.not.(iter .and. isSymm(trans, dimn, absthrsh)) .and. iter .lt. iterlimit)
        do i = 1, dimn
            do j = 1, dimn
                trans(j, i) = loLcl2eInt(j, i, dimn, MO2eIND, orth)
            end do
        end do
        call dgemm('T', 'N', dimn, dimn, dimn, 1.d0, trans, dimn, trans, dimn, 0.d0, seed, dimn)
        rdTrans(seed, dimn, eivl, store)
        call dgemm('N', 'N', dimn, dimn, dimn, 1.d0, trans, dimn, store, dimn, 0.d0, curt, dimn)
        call dgemm('N', 'N', dimn, dimn, dimn, 1.d0, curt, dimn, orthPrior, dimn, 0.d0, orth, dimn)
        orthPrior = orth
    end do

end subroutine initDIIS

subroutine rdTrans(seed, dimn, eivl, store)
    implicit none

    double precision, dimension(dimn, dimn) :: seed, eivx, store, intmd
    double precision, dimension(dimn) :: eivl
    integer :: i, j, dimn, status

    call construct(seed, dimn, eivl, eivx, status)

    do i = 1, dimn
        do j = i + 1, dimn
            store(i, i) = 1/sqrt(eivl(i))
            store(i, j) = 0.d0
            store(j, i) = 0.d0
        end do
    end do
    
    call dgemm('N', 'N', dimn, dimn, dimn, 1.d0, eivx, dimn, store, dimn, 0.d0, intmd, dimn)
    call dgemm('N', 'T', dimn, dimn, dimn, 1.d0, intmd, dimn, eivx, dimn, 0.d0, store, dimn)
end subroutine rdTrans

logical function isSymm(scMat, dimn, absthrsh)
    implicit none

    double precision, dimension(dimn, dimn) :: scMat
    double precision :: absthrsh
    integer :: i, j, dimn
    
    isSymm = .true.
    do i = 1, dimn
        do j = i + 1, dimn
            if(scMat(i, j) .gt. scMat(j, i) + absthrsh  &
               .or. scMat(i, j) .lt. scMat(j, i) - absthrsh)
               isSymm = .false.
        end do
    end do
end function isSymm

double precision function loLcl2eInt(j, i, dimn, MO2eIND, orth)  result(elval)
    implicit none
    
    double precision, dimension(dimn, dimn, dimn, dimn) :: MO2eIND
    double precision, dimension(dimn, dimn) :: orth
    double precision :: elval
    integer :: i, j, k, l, m, n, dimn

    elval = 0.d0

    do k = 1, dimn
        do l = 1, dimn
            do m = 1, dimn
                do n = 1, dimn
                    elval = elval + orth(k, i)*orth(l, j)*orth*(m, j)*orth(n, j)*MO2eIND(k, l, m, n)
                end do
            end do
        end do
    end do
end function loLcl2eInt