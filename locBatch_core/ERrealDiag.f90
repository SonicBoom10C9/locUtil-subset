module matOut
    implicit none
        integer, parameter, private :: errunit = 1000
        integer, parameter, private :: stdunit = 1001
        integer, parameter, private :: auxunit = 1002
    contains
    subroutine n_nOutput(stdmat, rDim, cDim, ldR, ldC, mode, lCaption)
        implicit none

        double precision, dimension(rdim, cdim) :: stdmat
        integer :: rDim, cDim, ldR, ldC, mode, i, j
        character (len = 64) :: lCaption

        print "(A)", trim(lcaption)//":"
        do i = ldR, rDim
            write (mode + 1000, '(*(ES12.4E4))') (stdmat(i, j), j = ldC, cDim)
        end do
    end subroutine n_nOutput
end module matOut

subroutine erDiagCall(debug)
    implicit none

    double precision, dimension(:, :), allocatable :: mcf
    double precision, dimension(:), allocatable :: redorb
    integer :: dimn, mode
    logical :: debug

    call copyInp(mcf, dimn)
    
    if(.not. allocated(mcf))        allocate(mcf(dimn, dimn))
    if(.not. allocated(redorb))     allocate(redorb(dimn))
            
    if(debug)    then
        !   call iterDiagUniv(mcf, dimn, redorb, debug)
    else
        call iterDiagReal(mcf, dimn, redorb, debug)
    endif

    contains
    subroutine copyInp(mcf, dimn)
        use gaussprog, only: mocoef, ntot
        
        implicit none

        double precision, dimension(:, :), allocatable	:: mcf
        integer :: dimn

        if(ntot .gt. 0) then
            dimn = ntot
	    !   print*,'Copied ntot to dimn',dimn
        else
        !	Do something, write to an error file
        endif

        if(.not. allocated(mcf))	allocate(mcf(dimn, dimn))

        mcf = mocoef

        !	#ifdef _RDBG
        !	#endif
    end subroutine copyInp
end subroutine erDiagCall

subroutine iterDiagReal(mcf, dimn, redorb, debug)
    use omp_lib
    
    implicit none
    
    double precision, dimension(dimn, dimn, dimn, dimn) :: AO2eIND, MO2eIND, lcl2eIND
    double precision, dimension(dimn, dimn) :: mcf, orth, seed, asmat, genr
    double precision, dimension(dimn) :: eivl, eivlPrior, redorb, diagmon
    double precision, parameter :: absThrsh = 1.0e-4, eps = 1.0e-00
    double precision :: infas, intglEval, diagonal, evalTrans
    integer :: dimn, i, j, k, l, m, n, index, status, iter
    logical :: debug, rangeChk, converged
    integer, parameter :: iterLimit = 1e6

    !   Cache values for AO 2e integrals
    call AOIntglVal(dimn, AO2eIND)
    !   Cache values for MO 2e integrals
    call MOIntglVal(dimn, mcf, AO2eIND, MO2eIND)
    return
    ! Initial seed generator
    do m = 1, dimn
        do n = 1, dimn
            ! Load element into AS matrix for seeding
            asmat(m, n) = MO2eIND(n, n, n, m) - MO2eIND(m, n, m, m)
        end do
    end do
    
    ! Defaults call mmul for dimn indices, seed = asmat*asmat
    call dgemm('N', 'N', dimn, dimn, dimn, 1.d0, asmat, dimn, asmat, dimn, 0.d0, seed, dimn)
    ! Sanity check, δD = ΣεF² > 0
    if(debug .and. .not. rangeChk(seed, eps, dimn, absThrsh))   then
        ! Do something
        return
    endif

    call construct(seed, dimn, eivl, genr, status)

    debug = .true.
    if(debug)   then
        if(.not. status) then
            print *, "Eigendecomposition routine reports no errors"
        else
            print *, "Something went wrong. Run debug backtrace for memory allocation"
        endif
    
        print *, "Eigenvectors (column ordered): "
        do i = 1, dimn
            write (*, 201) (genr(i, j), j = 1, dimn)
        end do
    
        print *, "Eigenvalues: "
        write(*, 201) (eivl(i), i = 1, dimn)
        
        201    format(*(ES12.4E4))
    endif
    debug = .false.

    iter = 0

    do while(.not.(iter .and. converged(eivl, eivlPrior, dimn, absThrsh)) .and. iter .lt. iterlimit)
        ! Copy for convergence check before caluclating new eigenvalues
        eivlPrior = eivl
        ! Solve the eigenvalue problem for eigenvalues and generator matrix
        call construct(seed, dimn, eivl, genr, status)
        ! Build the orthogonal matrix from these values and the genr matrix
        call generate(eps, eivl, genr, dimn, orth)   

        207    format(*(ES12.4E4))

        !   Deprecated call to form localised orbitals
        call form(mcf, dimn, orth, eps, redorb)
        !   Form localised 2e Integral values
        call lcl2eInt(dimn, orth, MO2eIND, lcl2eIND)

        do j = 1, dimn
            do k = 1, dimn
                ! Recalculate AS matrix elements from new orbitals for seed
                asmat(j, k) = lcl2eIND(k, k, k, j) - lcl2eIND(j, k, j, j)
            end do
        end do
        ! Update seed
        call dgemm('N', 'N', dimn, dimn, dimn, 1.d0, asmat, dimn, asmat, dimn, 0.d0, seed, dimn)
        iter = iter + 1
        !   print *, iter
    end do
    
    print *, "Symmetric matrix reduction: "
    do i = 1, dimn
        write (*, 204) (seed(i, j), j = 1, dimn)
    end do

    print *, "Antisymmetric matrix reduction: "
    do i = 1, dimn
        write (*, 204) (asmat(i, j), j = 1, dimn)
    end do

    diagonal = 0.d0
    print *, "Final ER diagonalisation Functional: "
    do i = 1, dimn
        diagonal = diagonal + lcl2eIND(i, i, i, i)
    end do
    write (*, 204) diagonal
    print *, "Iterations to Complete: ", iter
    
    204    format(*(ES12.4E4))
end subroutine iterDiagReal

logical function converged(vx, vxPrior, dimn, absThrsh)
    implicit none

    double precision, dimension(dimn) :: vx, vxPrior
    double precision :: absThrsh
    integer :: dimn, i, j

    converged = .true.
    do i = 1, dimn
        if((abs(vx(i)) - abs(vxPrior(i)) .gt. absThrsh) .or. (abs(vx(i)) - abs(vxPrior(i)) .lt. -absThrsh))   then
            converged = .false.
            exit
        endif
    end do
end function converged

logical function rangeChk(seed, eps, dimn, absthresh)
    implicit none

    double precision, dimension(dimn, dimn) :: seed
    double precision :: eps, absthresh
    integer :: dimn, m, n

    do m = 1, dimn-1
        do n = m + 1, dimn
            ! Quite impossible. This is just here to sheild against invalid ε
            if(eps*(seed(m, n)**2) .lt. 0) then
                rangeChk = .false.
                return
            endif
        end do
    end do

    rangeChk = .true.
    return
end function rangeChk

subroutine construct(seed, dimn, eivl, genr, status)
    implicit none

    double precision, dimension(dimn, dimn) :: seed, genr
    double precision, dimension(dimn) :: eivl
    double precision, dimension(:), allocatable :: work
    integer, dimension(:), allocatable :: support, iwork
    integer :: dimn, lwork, liwork, status, ilaenv, eict
    double precision :: dlamch
    
    allocate(support(dimn*dimn), work(1), iwork(1))
    
    liwork = -1; lwork = -1

    ! Get optimal working parameters for work, iwork
    call    dsyevr('V', 'A', 'U', dimn, seed, dimn, 0.d0, 0.d0, 0, 0, dlamch('S'),   &
            eict, eivl, genr, dimn, support, work, lwork, iwork, liwork, status)

    lwork = int(work(1));   liwork = iwork(1)
    deallocate(work, iwork)
    allocate(work(lwork), iwork(liwork))
    !   print *, size(work), size(iwork)

    ! Solve the eigenvalue problem. This is achieved through LAPACK
    ! ?syevr, present in the man pages. Briefly, it reduces input
    ! through iterative approximation of the tridiagonal form,
    ! and its reduction through scope detection of eigenvectors.
    ! Here is chosen the upper t.d. matrix, but doesn't matter,
    ! it is the calling convention in ilaenv for block detection.
    call    dsyevr('V', 'A', 'U', dimn, seed, dimn, 0.d0, 0.d0, 0, 0, dlamch('S'),   &
            eict, eivl, genr, dimn, support, work, lwork, iwork, liwork, status)

    !   deallocate(support, work , iwork)
end subroutine construct

subroutine generate(eps, eivl, genr, dimn, orth)
    implicit none

    double precision, dimension(dimn, dimn) :: genr, orth
    double precision, dimension(dimn) :: eivl
    double precision :: eps
    integer :: dimn, m, n, k

    do m = 1, dimn
        do n = 1, dimn
            orth(m, n) = 0.d0
            do k = 1, dimn
                ! Compute orthogonal matrix element from generator
                orth(m, n) = orth(m, n) + exp(eps*eivl(k))*genr(m, k)*genr(k, n)
            end do
        end do
    end do
end subroutine generate

subroutine form(mcf, dimn, orth, eps, redorb)
    implicit none

    double precision, dimension(dimn, dimn) :: mcf, orth
    double precision, dimension(dimn) :: redorb
    double precision :: eps
    integer :: dimn, k, n, i

    do k = 1, dimn
        redorb(k) = 0.d0
        do n = 1, dimn
            ! If debugging, orth(n, k) varies with ε, change parameter
            ! Calculate the reduced orbital of calling super-iteration
            do i = 1, dimn
                redorb(k) = redorb(k) + mcf(n, i)*orth(n, k)
            end do
        end do
    end do
end subroutine form

!   Obsolete
double precision function evalTrans(ldA, ldB, dimn, mcf, orth)  result(element)
    implicit none

    integer :: ldA, ldB, dimn, a, b, c, d, k, l, m, n
    double precision, dimension(dimn, dimn)  :: mcf, orth
    double precision, dimension(dimn, dimn, dimn, dimn) :: AO2eIND
    double precision :: element, expr, intglEval
    
    element = 0.d0
    do a = 1, dimn
        do b = 1, dimn
            do c = 1, dimn
                do d = 1, dimn
                    
                    expr = 0.d0
                    do k = 1, dimn
                        do l = 1, dimn
                            do m = 1, dimn
                                do n = 1, dimn
                                    expr = expr + mcf(a, k)*mcf(b, l)*mcf(c, m)*mcf(d, n)*AO2eIND(k, l, m, n)
                                end do
                            end do	
                        end do
                    end do
                    element = element + orth(b, ldA)*orth(d, ldB)*expr*(orth(a, ldB)*orth(c, ldB) - orth(a, ldA)*orth(c, ldA))
                
                end do
            end do
        end do
    end do
end function evalTrans

subroutine AOIntglVal(dimn, AO2eIND)
    implicit none

    double precision, dimension(dimn, dimn, dimn, dimn) :: AO2eIND
    double precision :: intglEval
    integer :: dimn, i, j, k, l

    do i = 1, dimn
        do j = 1, dimn
            do k = 1, dimn
                do l = 1, dimn
                    AO2eIND(i, j, k, l) = intglEval(i, j, k, l)
                end do
            end do
        end do
    end do
end subroutine AOIntglVal

subroutine MOIntglVal(dimn, mcf, AO2eIND, MO2eIND)
    implicit none

    double precision, dimension(dimn, dimn, dimn, dimn) :: MO2eIND, AO2eIND
    double precision, dimension(dimn, dimn) :: mcf
    double precision :: elval
    integer :: dimn, i, j, k, l, p, q, r, s

    do i = 1, dimn
        do j = 1, dimn
            do k = 1, dimn
                do l = 1, dimn
                    
                    elval = 0.d0
                    do p = 1, dimn
                        do q = 1, dimn
                            do r = 1, dimn
                                do s = 1, dimn
                                    elval = elval + mcf(i, p)*mcf(j, q)*mcf(k, r)*mcf(l, s)*AO2eIND(p, q, r, s)
                                end do
                            end do
                        end do
                    end do
                    MO2eInd(i, j, k, l) = elval

                end do
            end do
        end do
    end do
end subroutine MOIntglVal

subroutine lcl2eInt(dimn, orth, MO2eIND, lcl2eIND)
    double precision, dimension(dimn, dimn, dimn, dimn) :: MO2eIND, lcl2eIND
    double precision, dimension(dimn, dimn) :: orth
    double precision :: elval
    integer :: dimn, i, j, k, l, p, q, r, s

    do p = 1, dimn
        do q = 1, dimn
            do r = 1, dimn
                do s = 1, dimn
                    
                    elval = 0.d0
                    do i = 1, dimn
                        do j = 1, dimn
                            do k = 1, dimn
                                do l = 1, dimn
                                    elval = elval + MO2eIND(i, j, k, l)*orth(i, p)*orth(j, q)*orth(k, r)*orth(l, s)
                                end do
                            end do
                        end do
                    end do
                    lcl2eIND(p, q, r, s) = elval
                    
                end do
            end do
        end do
    end do
end subroutine lcl2eInt

double precision function intglEval(aoL1, aoL2, aoR1, aoR2)    result(int2e)
    use genvar
    use formmd
    use gaussprog
    use fock_test
    use hf_intgls
    
    implicit none
    
    double precision(8) :: Fopt, Flarge, os_2e, fval
    integer :: com(8), flag, tang, aoL1, aoL2, aoR1, aoR2

    tcut = dble(nx4); mvalnmax4 = maxval(nmax4)
    int2e = 0.d0
    
    do ii = 1, 3
        e1t(ii) = parami(aoL1, aoR1, 3 + ii - 1);   e2t(ii) = parami(aoL1, aoR1, 6 + ii - 1)
        e3t(ii) = parami(aoL2, aoR2, 3 + ii - 1);   e4t(ii) = parami(aoL2, aoR2, 6 + ii - 1)
        a(ii) = nucl(parami(aoL1, aoR1, 1), ii);    b(ii) = nucl(parami(aoL1, aoR1, 2), ii)
        c(ii) = nucl(parami(aoL2, aoR2, 1), ii);    d(ii) = nucl(parami(aoL2, aoR2, 2), ii)
    end do

    tang = sum(e1t) + sum(e2t) + sum(e3t) + sum(e4t)
    call comp4(e1t, e2t, e3t, e4t, flag, com)
    ll = 0 

    do ii = 1, parami(aoL1, aoR1, 9)
        do jj = 1, parami(aoL2, aoR2, 9)
            a1 = paramr(aoL1, aoR1, ii, 1); a2 = paramr(aoL1, aoR1, ii, 2)
            a3 = paramr(aoL2, aoR2, jj, 1); a4 = paramr(aoL2, aoR2, jj, 2)
            aL = a1 + a2; aR = a3 + a4; rho = (aL*aR)/(aL + aR); 
            ll = ll + 1
            
            ovlpprod = (ovl2e(aoL1, aoR1, ii)*ovl2e(aoL2, aoR2, jj))
            p(:) = ((a1*a(:)) + (a2*b(:)))/aL
            q(:) = ((a3*c(:)) + (a4*d(:)))/aR
            vec = p - q
            T = rho*dot_product(vec, vec)
            ! s, p, d basis
            
            if (T <= tcut) then
                nm = nmax4(floor(T)) + 5
                if (nm > mvalnmax4) nm = mvalnmax4
                fval = Fopt(tmax, T, nmax, del, gammaf, nginv4, gaminv4, nm)
            end if
            
            if (T > tcut) fval = Flarge(tmax, T, nmax, del, gammaf)
            int2e = int2e + os_2e(tang, flag, com, ovlpprod, a1, a2, a3, a4, a, b, c, d, fval)
        end do
    end do
end function intglEval