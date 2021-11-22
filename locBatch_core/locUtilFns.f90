subroutine randInit(dimn, mcf)
	implicit none

	double precision, dimension(dimn, dimn) :: mcf
	double precision, dimension(dimn*3) :: work
	double precision, dimension(4) :: dseed
	integer, dimension(4) :: iseed
	integer :: dimn, i, info

	write(*, '(A31)', advance = 'no')	"Generating random multiplier..."
	
	do i = 1, dimn
		call cpu_time(dseed(i))
		call sleep(1)
    end do
    
    !   Use generated orth matrix and normalise for initial guess

	iseed = int(dseed)
	iseed(4) = 2*iseed(4) - 1
	!	write(*, '(I4)')	(iseed(i), i = 1, 4)

	call dlaror('L', 'N', dimn, dimn, mcf, dimn, iseed, work, info)
	write(*, '(A4)')	"done"
end subroutine randInit

subroutine gaussVxOutNS(stdMat, rDim, cDim, ldR, ldC)
	implicit none

	double precision, dimension(rDim, cDim) :: stdMat
	double precision, dimension(rDim*cDim) :: sprint
	integer :: rDim, cDim, ldR, ldC, i, j, k

	k = 1
    do i = ldR, rDim
        do j = ldC, cDim
            sprint(k) = stdMat(i, j)
            k = k + 1
        end do
	end do

	!	Gaussian Standard output format
    do i = 0, rDim**2/5 - 1
		write(*, 274)   (sprint(j), j = 5*i + 1, 5*(i + 1))
	end do
	write(*, 274) sprint(rDim*cDim)
	print *, ""

	274    format(*(ES16.8E4))
end subroutine gaussVxOutNS

subroutine gaussVxOut(stdMat, dimn)
	implicit none

	double precision, dimension(dimn, dimn) :: stdMat
	double precision, dimension(dimn*dimn) :: sprint
	integer :: dimn, i, j, k

	k = 1
    do i = 1, dimn
        do j = 1, dimn
            sprint(k) = stdMat(i, j)
            k = k + 1
        end do
    end do
    
    !   print *, size(sprint)

    !	Gaussian Standard output format
    k = mod(dimn**2, 5)

    do i = 1, dimn**2 - k, 5
        write(*, 274)   (sprint(j), j = i, i + 4)
    end do
    write(*, 274)   (sprint (j), j = i, dimn**2)
    
    
	print *, ""

	274    format(*(ES16.8E4))
end subroutine gaussVxOut

subroutine AOMmtIntglVal(dimn, AO1eIND, swid)
	use externAccess

	implicit none

	double precision, dimension(dimn, dimn, 3) :: AO1eIND
    integer :: dimn, i, j
    logical :: swid
	
    if(.not. swid)  then
    !	Copy pre-evaluated integral values
	    AO1eIND = moments
	else
	!	For completeness, no need to use this if just copying
		do i = 1, dimn
			do j = 1, dimn
				AO1eIND(i, j, 1) = moments(j, i, 1)
				AO1eIND(i, j, 2) = moments(j, i, 2)
				AO1eIND(i, j, 3) = moments(j, i, 3)
			end do
        end do
    endif
end subroutine AOMmtIntglVal

subroutine TransMmtIntglVal(dimn, start, lim, trans, trans1eIND, Cur1eIND, Cur1eINDrd, swid)
    use omp_lib
    use externAccess

	implicit none

	double precision, dimension(dimn, dimn, 3) :: trans1eIND, Cur1eIND
	double precision, dimension(dimn, dimn), optional :: Cur1eINDrd
	double precision, dimension(dimn, dimn) :: trans
    double precision, dimension(3) :: lval
    double precision :: xval, yval, zval
	integer :: dimn, i, j, k, l, execlim, start, lim
	logical :: swid

	!	Simple basis transformation for moment integrals
	!	"trans" MUST be an orthogonal coefficient matrix, preferably normalised
    !	TBM: Generalise for n-D instead of 3-D, if ever necessary 
    
    execlim = omp_get_num_procs()

	if(swid)	then
		do i = start, lim
			do j = start, lim
                xval = 0.d0;	yval = 0.d0;	zval = 0.d0
                do k = 1, dimn
					do l = 1, dimn
						xval = xval + trans(k, i)*trans(l, j)*trans1eIND(k, l, 1)
						yval = yval + trans(k, i)*trans(l, j)*trans1eIND(k, l, 2)
                        zval = zval + trans(k, i)*trans(l, j)*trans1eIND(k, l, 3)
					end do
                end do
                Cur1eIND(i, j, 1) = xval
                Cur1eIND(i, j, 2) = yval
                Cur1eIND(i, j, 3) = zval
				!	Store scalval for quick calculation, optional
				if(present(Cur1eINDrd))	Cur1eINDrd(i, j) = dsqrt(Cur1eIND(i, j, 1)**2 + Cur1eIND(i, j, 2)**2 + Cur1eIND(i, j, 3)**2)
			end do
		end do
	else
		do i = start, lim
			do j = start, lim
                xval = 0.d0;	yval = 0.d0;	zval = 0.d0
                do k = 1, dimn
					do l = 1, dimn
						xval = xval + trans(i, k)*trans(j, l)*trans1eIND(k, l, 1)
						yval = yval + trans(i, k)*trans(j, l)*trans1eIND(k, l, 2)
						zval = zval + trans(i, k)*trans(j, l)*trans1eIND(k, l, 3)
					end do
                end do
                Cur1eIND(i, j, 1) = xval
                Cur1eIND(i, j, 2) = yval
                Cur1eIND(i, j, 3) = zval
				!	Store scalval for quick calculation, optional
				if(present(Cur1eINDrd))	Cur1eINDrd(i, j) = dsqrt(xval**2 + yval**2 + zval**2)
			end do
		end do
	end if
end subroutine TransMmtIntglVal

subroutine TransMmtIntglValPl(dimn, start, lim, trans, trans1eIND, Cur1eIND, Cur1eINDrd, swid)
    use omp_lib    
    use externAccess

	implicit none

	double precision, dimension(dimn, dimn, 3) :: trans1eIND, Cur1eIND
	double precision, dimension(dimn, dimn), optional :: Cur1eINDrd
	double precision, dimension(dimn, dimn) :: trans
	double precision :: xval, yval, zval, xvi, yvi, zvi
	integer :: dimn, i, j, k, l, execlim, lim, start
	logical :: swid

	!	Simple basis transformation for moment integrals
	!	"trans" MUST be an orthogonal coefficient matrix, preferably normalised
    !	TBM: Generalise for n-D instead of 3-D, if ever necessary 
    
    execlim = omp_get_num_procs()

	if(swid)	then
		do i = start, lim
			do j = start, lim
                xval = 0.d0;	yval = 0.d0;	zval = 0.d0
                xvi = 0.d0;     yvi = 0.d0;     zvi = 0.d0
                !$OMP parallel shared(xval, yval, zval) private(xvi, yvi, zvi) num_threads(execlim)
                !$OMP do schedule(dynamic)
                do k = 1, dimn
					do l = 1, dimn
						xvi = xvi + trans(k, i)*trans(l, j)*trans1eIND(k, l, 1)
						yvi = yvi + trans(k, i)*trans(l, j)*trans1eIND(k, l, 2)
						zvi = zvi + trans(k, i)*trans(l, j)*trans1eIND(k, l, 3)
					end do
                end do
                !$OMP end do
                !$OMP critical
				xval = xval + xvi
                yval = yval + yvi
                zval = zval + zvi
                !$OMP end critical
                !$OMP end parallel
                Cur1eIND(i, j, :) = (/xval, yval, zval/)
				!	Store scalval for quick calculation, optional
				if(present(Cur1eINDrd))	Cur1eINDrd(i, j) = dsqrt(xval**2 + yval**2 + zval**2)
			end do
		end do
	else
		do i = start, lim
			do j = start, lim
                xval = 0.d0;	yval = 0.d0;	zval = 0.d0
                xvi = 0.d0;     yvi = 0.d0;     zvi = 0.d0
                !$OMP parallel shared(xval, yval, zval) private(xvi, yvi, zvi) num_threads(execlim)
                !$OMP do schedule(dynamic)
				do k = 1, dimn
					do l = 1, dimn
						xvi = xvi + trans(i, k)*trans(j, l)*trans1eIND(k, l, 1)
						yvi = yvi + trans(i, k)*trans(j, l)*trans1eIND(k, l, 2)
						zvi = zvi + trans(i, k)*trans(j, l)*trans1eIND(k, l, 3)
					end do
                end do
                !$OMP end do
                !$OMP critical
				xval = xval + xvi
                yval = yval + yvi
                zval = zval + zvi
                !$OMP end critical
                !$OMP end parallel
                Cur1eIND(i, j, :) = (/xval, yval, zval/)
				if(present(Cur1eINDrd))	Cur1eINDrd(i, j) = sqrt(xval**2 + yval**2 + zval**2)
			end do
		end do
	end if
end subroutine TransMmtIntglValPl

logical function isConverged(dimn, Lcl1eIND, absThrsh)
	implicit none

	double precision, dimension(dimn, dimn, 3) :: Lcl1eIND
	double precision ::	rKl, absThrsh
	integer :: dimn, i, j

	isConverged = .false.

	rKl = 0.d0
	do i = 1, dimn
		do j = i + 1, dimn
			!	Convergence to D(\phi), linear gradient only
			rKl = rKl + dot_product((Lcl1eIND(i, i, :) - Lcl1eIND(j, j, :)), Lcl1eIND(i, j, :))
		end do
	end do
	rKl = 4*dsqrt(rKl)
	print *, rKl
	
	! Determine if gradient for line min. is within margin of error
	if(rKl .lt. absThrsh .and. rKl .gt. -absThrsh) isConverged = .true.
end function isConverged

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

    ! Solve the eigenvalue problem. This is achieved through BLAS
    ! ?syevr, present in the man pages. Briefly, it reduces input
    ! through iterative approximation of the tridiagonal form,
    ! and its reduction through scope detection of eigenvectors.
    ! Here is chosen the upper t.d. matrix, but doesn't matter,
    ! it is the calling convention in ilaenv for block detection.
    call    dsyevr('V', 'A', 'U', dimn, seed, dimn, 0.d0, 0.d0, 0, 0, dlamch('S'),   &
            eict, eivl, genr, dimn, support, work, lwork, iwork, liwork, status)

    !   deallocate(support, work, iwork)
end subroutine construct

logical function isOrth(inpMat, dimn, absThrsh)
	implicit none
    double precision, dimension(dimn, dimn) :: inpMat, resMat
    double precision :: absThrsh
    integer :: dimn, i, j

    call dgemm('N', 'T', dimn, dimn, dimn, 1.d0, inpMat, dimn, inpMat, dimn, 0.d0, resMat, dimn)

	isOrth = .true.
	do i = 1, dimn
		do j = i + 1, dimn
			if(resMat(i, i) .gt. (1 + absThrsh) .or. resMat(i, i) .lt. (1 - absThrsh)	&
			  .or. resMat(i, j) .gt. absThrsh .or. resMat(i, j) .lt. -absThrsh)	then
				isOrth = .false.
			endif
		end do
	end do
end function isOrth

subroutine dbOut(stdMat, rDim, cDim, ldR, ldC)
	implicit none

    double precision, dimension(rdim, cdim) :: stdmat
    integer :: rDim, cDim, ldR, ldC, i, j
    character (len = 64) :: lCaption

    do i = ldR, rDim
        write (*, '(*(ES16.8E4))') (stdmat(i, j), j = ldC, cDim)
    end do
    100 format(*(A))
end subroutine dbOut

logical function scalConverged(sc, scPrior, absThrsh)
    implicit none

    double precision :: sc, scPrior
    double precision :: absThrsh

    scalConverged = .true.
	if((abs(sc) - abs(scPrior) .gt. absThrsh) .or. (abs(sc) - abs(scPrior) .lt. -absThrsh))   then
		scalConverged = .false.
	endif
end function scalConverged

subroutine rstTrMat(uplo, stdMat, dimn)
    implicit none

    double precision, dimension(dimn, dimn) :: stdMat
    integer :: i, j, dimn
    character :: uplo

    select case(uplo)
    !   Null upper triangular elements
    case('L', 'l')
        do i = 1, dimn
            do j = i + 1, dimn
                stdMat(i, j) = 0.d0
            end do
        end do
    !   Null lower triangular elements
    case('U', 'u')
        do i = 1, dimn
            do j = 1, i - 1
                stdMat(i, j) = 0.d0
            end do
        end do
    !   Null off-diagonal  elements
    case('A', 'a')
        do i = 1, dimn
            do j = i + 1, dimn
                stdMat(i, j) = 0.d0
                stdMat(j, i) = 0.d0
            end do
        end do
    end select
end subroutine rstTrMat

logical function isSymm(scMat, dimn, absthrsh)
    implicit none

    double precision, dimension(dimn, dimn) :: scMat
    double precision :: absthrsh
    integer :: i, j, dimn
    
    isSymm = .true.
    do i = 1, dimn
        do j = i + 1, dimn
            !   Probably a good idea to tweak this for rough symmetry
            if(scMat(i, j) .gt. scMat(j, i) + absthrsh  &
               .or. scMat(i, j) .lt. scMat(j, i) - absthrsh)    then
               isSymm = .false.
            endif
        end do
    end do
end function isSymm

subroutine AOIntglVal(dimn, AO2eIND)
    use omp_lib
    
    implicit none

    double precision, dimension(dimn, dimn, dimn, dimn) :: AO2eIND
    double precision, dimension(:, :, :, :), allocatable :: str2eIND
    double precision :: intglEval
    integer :: dimn, i, j, k, l, execlim

    allocate(str2eIND(dimn, dimn, dimn, dimn))

    execlim = omp_get_num_procs()

    !$omp parallel num_threads(execlim)
    !$omp do schedule(static)
    do i = 1, dimn
        do j = 1, dimn
            do k = 1, dimn
                do l = 1, dimn
                    str2eIND(i, j, k, l) = intglEval(i, j, k, l)
                end do
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

    AO2eIND = str2eIND
end subroutine AOIntglVal

subroutine MOIntglVal(dimn, mcf, AO2eIND, MO2eIND)
    implicit none

    double precision, dimension(dimn, dimn, dimn, dimn) :: MO2eIND, AO2eIND
    double precision, dimension(dimn, dimn) :: mcf
    double precision :: elval, elpriv
    integer :: dimn, i, j, k, l, p, q, r, s, execlim, omp_get_num_procs

    execlim = omp_get_num_procs()

    do i = 1, dimn
        do j = 1, dimn
            write(*, '(2I4)') i, j
            do k = 1, dimn
                do l = 1, dimn
                    
                    elval = 0.d0
                    elpriv = 0.d0
                    
                    !$OMP parallel shared(elval) private(elpriv) num_threads(execlim)
                    !$OMP do schedule(dynamic)
                    do p = 1, dimn
                        do q = 1, dimn
                            do r = 1, dimn
                                do s = 1, dimn
                                    elpriv = elpriv + mcf(i, p)*mcf(j, q)*mcf(k, r)*mcf(l, s)*AO2eIND(p, q, r, s)
                                end do
                            end do
                        end do
                    end do
                    !$OMP end do
                    !$OMP critical
                    elval = elval + elpriv
                    !$OMP end critical
                    !$OMP end parallel
                    MO2eInd(i, j, k, l) = elval

                end do
            end do
        end do
    end do
end subroutine MOIntglVal

subroutine Intgl2eTransPS(dimn, mcf, AO2eIND, MO2eIND)
    implicit none

    double precision, dimension(dimn, dimn, dimn, dimn) :: MO2eIND, AO2eIND, str1_2e, str2_2e, str3_2e
    double precision, dimension(dimn, dimn) :: mcf
    double precision :: elval, elpriv
    integer :: dimn, i, j, k, l, p, q, r, s, execlim, omp_get_num_procs

    execlim = omp_get_num_procs()

    MO2eIND(:, :, :, :) = 0.d0

    do p = 1, dimn
        do q = 1, dimn
            do r = 1, dimn
                do s = 1, dimn
                    do l = 1, dimn
                        str1_2e(p, q, r, l) = str1_2e(p, q, r, l) + mcf(l, s)*AO2eIND(p, q, r, s)
                    end do
                end do
                !   str2eIND = MO2eIND
                !   MO2eIND(:, :, :, :) = 0.d0
                do k = 1, dimn
                    do l = 1, dimn
                        str2_2e(p, q, k, l) = str2_2e(p, q, k, l) + mcf(k, r)*str1_2e(p, q, r, l)
                    end do
                end do
            end do
            !   str2eIND = MO2eIND
            !   MO2eIND(:, :, :, :) = 0.d0
            do j = 1, dimn
                do k = 1, dimn
                    do l = 1, dimn
                        str3_2e(p, j, k, l) = str3_2e(p, j, k, l) + mcf(j, q)*str2_2e(p, q, k, l)
                    end do
                end do
            end do
        end do
        !   str2eIND = MO2eIND
        !   MO2eIND(:, :, :, :) = 0.d0
        do i = 1, dimn
            do j = 1, dimn
                do k = 1, dimn
                    do l = 1, dimn
                        MO2eIND(i, j, k, l) = MO2eIND(i, j, k, l) + mcf(i, p)*str3_2e(p, j, k, l)
                    end do
                end do
            end do
        end do
    end do
end subroutine Intgl2eTransPS

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

subroutine updateIntglsPl(i, j, dimn, verifMat, Lcl1eIND, Lcl1eINDrd, Str1eIND, Str1eINDrd)
	use omp_lib

	implicit none

	double precision, dimension(dimn, dimn, 3) :: Lcl1eIND, Str1eIND
	double precision, dimension(dimn, dimn) :: verifMat, Lcl1eINDrd, Str1eINDrd
	double precision :: xpriv, ypriv, zpriv, xv, yv, zv
	integer :: i, j, k, l, dimn, execlim

	Str1eIND(:, :, :) = 0.d0;	Str1eINDrd(:, :) = 0.d0
	xpriv = 0.d0;	ypriv = 0.d0;	zpriv = 0.d0
	xv = 0.d0;	yv = 0.d0;	zv = 0.d0

	execlim = omp_get_num_procs()

	!$omp parallel shared(xv, yv, zv) private(xpriv, ypriv, zpriv)
	!$omp do schedule(dynamic)
	do k = 1, dimn
		do l = 1, dimn
			!	Compare with ovlp_mo output from trans1
			xpriv = xpriv + verifMat(k, i)*verifMat(l, j)*Lcl1eIND(k, l, 1)
			ypriv = ypriv + verifMat(k, i)*verifMat(l, j)*Lcl1eIND(k, l, 2)
			zpriv = zpriv + verifMat(k, i)*verifMat(l, j)*Lcl1eIND(k, l, 3)
		end do
	end do
	!$omp end do
	!$omp critical
	xv = xv + xpriv
	yv = yv + ypriv
	zv = zv + zpriv
	!$omp end critical
	!$omp end parallel
	Str1eIND(i, j, :) = (/xv, yv, zv/)

	Str1eINDrd(i, j) = dsqrt(Str1eIND(i, j, 1)**2 + Str1eIND(i, j, 2)**2 + Str1eIND(i, j, 3)**2)
end subroutine updateIntglsPl

subroutine updateIntgls(i, j, dimn, verifMat, Lcl1eIND, Lcl1eINDrd, Str1eIND, Str1eINDrd)
	implicit none

	double precision, dimension(dimn, dimn, 3) :: Lcl1eIND, Str1eIND
	double precision, dimension(dimn, dimn) :: verifMat, Lcl1eINDrd, Str1eINDrd
	double precision :: xpriv, ypriv, zpriv, xv, yv, zv
	integer :: i, j, k, l, dimn, execlim

	Str1eIND(:, :, :) = 0.d0;	Str1eINDrd(:, :) = 0.d0
	xv = 0.d0;	yv = 0.d0;	zv = 0.d0

	do k = 1, dimn
		do l = 1, dimn
			!	Compare with ovlp_mo output from trans1
			xv = xv + verifMat(i, k)*verifMat(j, l)*Lcl1eIND(k, l, 1)
			yv = yv + verifMat(i, k)*verifMat(j, l)*Lcl1eIND(k, l, 2)
			zv = zv + verifMat(i, k)*verifMat(j, l)*Lcl1eIND(k, l, 3)
		end do
	end do

	Str1eIND(i, j, :) = (/xv, yv, zv/)

	Str1eINDrd(i, j) = dsqrt(Str1eIND(i, j, 1)**2 + Str1eIND(i, j, 2)**2 + Str1eIND(i, j, 3)**2)
end subroutine updateIntgls