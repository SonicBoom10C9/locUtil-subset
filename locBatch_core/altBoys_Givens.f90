subroutine byDiag()
	use omp_lib
	use boysParams
	use externAccess

	implicit none

	!  Calls gaussin() from input_setup.f90 - ensure in $PWD

	double precision, dimension(:, :), allocatable :: mcf
	double precision :: BoysFnl, dlamch, tcoll
	double precision, parameter :: absthrsh = 1e-8
	integer :: dimn, iter, boys, mode, occDimn, verify
	logical :: dbg, rndS

	print *, "======== Localisation calculations begin here ========"
	!  Read MO coeffs using gaussin()
	call copyInp(mcf, dimn)
	!	if(.not. allocated(mcf))	allocate(mcf(dimn, dimn))

	!	0 = Occupied only, 1 = Virtual Only
	!	2 = Occupied + Virtual (output of occupied -> input of virtual)
	!	Anything else = complete basis
	write (*, *, advance='no')	"Reading mode: "
	read(*, *)	mode

	dbg = .false.	
	print *, "Print iterative coefficients? (T/F) [F]: "
	read (*, *, iostat=verify)	dbg

	rndS = .false.
	!	Warning: Do not use a very high precision when doing this,
	!	It causes a runaway iterative condition
	!	Max recommended threshold is 1e-8
	print *, "Subject coeffs to a random orthogonal transformation? Not recommended. (T/F) [F]: "
	read(*, *, iostat=verify)	rndS
	
	iter = boys(mcf, occDimn, dimn, absthrsh, BoysFnl, mode, dbg, rndS)
	
	print *, ""
	write(*, 100) " Boys functional and iterations: ", BoysFnl, ", ", iter

	100	format(A32, ES12.4E4, A2, I8)

	contains
	subroutine copyInp(mcf, dimn)
		use gaussprog, only: mocoef, ntot
		
		implicit none

		interface
			subroutine gaussin()
			end subroutine gaussin
		end interface

		double precision, dimension(:, :), allocatable, intent(out)	:: mcf
        integer, intent(out) :: dimn
        integer :: i, j

		if(ntot .gt. 0) then
			dimn = ntot
		else
		!	Do something, write to an error file. Up to the user.
		endif

		if(.not. allocated(mcf))	allocate(mcf(dimn, dimn))

        mcf = mocoef
	end subroutine copyInp
end subroutine byDiag

integer function boys(mcf, occDimn, dimn, absThrsh, BoysFnl, mode, dbg, rndS)	result(iter)
	use boysParams
	use externAccess

	implicit none
	
	double precision, dimension(dimn, dimn, 3) :: AO1eIND, MO1eIND, Lcl1eIND
	double precision, dimension(dimn, dimn) :: AO1eINDrd, MO1eINDrd, Lcl1eINDrd
    double precision, dimension(dimn, dimn) :: mcf, orth, orthPrior, gemmResult, verifMat
    double precision, dimension(dimn*dimn) :: sprint
	double precision, dimension(dimn) :: eival, long, lat
	double precision :: BoysFnl, maxFnl, iterGain, iterGainPrev, absThrsh, gamma, A_ij, B_ij, cgamma, sgamma
	integer :: i, j, k, l, occDimn, iter, dimn, mode, start, end, spliter
	!	Clamp at static convergence point for now, iterlimit ~ 100 for 1e-6 (remember to check signflip)
	integer, parameter :: iterlimit = 1e4
	logical :: scalConverged, isOrth, orthTest, isConverged, swid, dbg, rndS

	!	Leave true unless orbital coefficients are transposed
	swid = .true.		
	!	Get this from input, use as is for now in LiH STO-3G
	occDimn = 2

	!	Set localisation bounds
	select case(mode)
	case(0, 2)
		start = 1
		end = occDimn
	case(1)
		start = occDimn + 1
		end = dimn
	case default
		start = 1
		end = dimn
	end select

	call AOMmtIntglVal(dimn, AO1eIND)
	call TransMmtIntglVal(dimn, mcf, AO1eIND, MO1eIND, MO1eINDrd, swid)
	if(rndS)	call randInit(dimn, mcf)

	BoysFnl = 0.d0
	iter = 0

	orth(:, :) = 0.d0
	
	do i = 1, dimn
		BoysFnl = BoysFnl + MO1eINDrd(i, i)**2
		orth(i, i) = 1.d0
	end do
	
	if(mode - 2)	write(*, '(A34, ES12.4E4)') " Initial Boys SOS Function value: ", BoysFnl
	print *, "------------------------------------------------------"

	orthPrior = orth
	gemmResult = orthPrior
	verifMat = mcf
	Lcl1eIND = MO1eIND
	Lcl1eINDrd = MO1eINDrd

	do while(.not.(iter .and. scalConverged(iterGain, iterGainPrev, absThrsh)) .and. iter .lt. iterlimit)
		do i = start, end - 1
			!	Strict upper half, invertible
			do j = i + 1, end
                !	orthPrior = gemmResult
                
                A_ij = dot_product(Lcl1eIND(i, j, :), Lcl1eIND(i, j, :)) - 0.25*dot_product(	&
					(Lcl1eIND(i, i, :) - Lcl1eIND(j, j, :)), (Lcl1eIND(i, i, :) - Lcl1eIND(j, j, :))	&
				)
                B_ij = dot_product((Lcl1eIND(i, i, :) - Lcl1eIND(j, j, :)), Lcl1eIND(i, j, :))

                if(A_ij**2 + B_ij**2 .lt. 1D-12)	cycle

                gamma = 0.25*sign(1d0, B_ij)*acos(-A_ij/dsqrt(A_ij**2 + B_ij**2))

                long = cos(gamma)*mcf(i, :) + sin(gamma)*mcf(j, :) 
                lat = -sin(gamma)*mcf(i, :) + cos(gamma)*mcf(j, :)
				verifMat(i, :) = long
				verifMat(j, :) = lat
			end do
		end do
		
		if(dbg)	then
			print *, "Coeffcients after iteration: ", iter + 1
			call gaussVxOut(verifMat, dimn, dimn, 1, 1)
		endif

		call TransMmtIntglVal(dimn, verifMat, MO1eIND, Lcl1eIND, Lcl1eINDrd, swid)

		print *, isConverged(dimn, Lcl1eIND, absThrsh)
		
		if(scalConverged(iterGain, iterGainPrev, 1e-4) .and. mode .ne. 2)	then
			print *, "Note: Bouncing off threshold. Convergence was achieved a while back, consider lowering precision."
			print *, "If it is sufficiently low, this was an already mostly localised basis"
		endif

		maxFnl = BoysFnl
		BoysFnl = 0.d0
		do i = 1, dimn
			BoysFnl = BoysFnl + Lcl1eINDrd(i, i)**2
		end do
		
		iterGainPrev = iterGain
		iterGain = BoysFnl - maxFnl
        iter = iter + 1
        
        !   print *, iterGain
	end do
	!	Orthogonality test, use for verification if generation method changes
	!	print *, "Orthogonality test: ", isOrth(gemmResult, dimn, absThrsh)

	!	call dgemm('N', 'N', dimn, dimn, dimn, 1.d0, gemmResult, dimn, mcf, dimn, 0.d0, verifMat, dimn)

	if(mode .eq. 2)	then
		spliter = boys(verifMat, occDimn, dimn, absThrsh, BoysFnl, 1, dbg, .false.)
		iter = iter + spliter
		write(*, '(A35, I4)'), " Secondary localisation iterations: ", spliter
		return
	endif

	print *, ""
	if(mode .eq. 0)	then
		print *, "Occupied-only localised coefficients: "
	else	
		print *, "Final coefficients: "
	endif

	call gaussVxOut(verifMat, dimn, dimn, 1, 1)
	204    format(*(ES16.8E4))
end function boys

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

	iseed = int(dseed)
	iseed(4) = 2*iseed(4) - 1
	!	write(*, '(I4)')	(iseed(i), i = 1, 4)

	call dlaror('L', 'N', dimn, dimn, mcf, dimn, iseed, work, info)
	write(*, '(A4)')	"done"
end subroutine randInit

subroutine gaussVxOut(stdMat, rDim, cDim, ldR, ldC)
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
end subroutine gaussVxOut

subroutine AOMmtIntglVal(dimn, AO1eIND)
	use externAccess

	implicit none

	double precision, dimension(dimn, dimn, 3) :: AO1eIND
	integer :: dimn, i, j
	
	!	Copy pre-evaluated integral values
	AO1eIND = moments
	
	!	For completeness, no need to use this if just copying
	!	do i = 1, dimn
	!		do j = 1, dimn
	!			AO1eIND(i, j, 1) = moments(i, j, 1)
	!			AO1eIND(i, j, 2) = moments(i, j, 2)
	!			AO1eIND(i, j, 3) = moments(i, j, 3)
	!		end do
	!	end do
end subroutine AOMmtIntglVal

subroutine TransMmtIntglVal(dimn, trans, trans1eIND, Cur1eIND, Cur1eINDrd, swid)
	use externAccess

	implicit none

	double precision, dimension(dimn, dimn, 3) :: trans1eIND, Cur1eIND
	double precision, dimension(dimn, dimn), optional :: Cur1eINDrd
	double precision, dimension(dimn, dimn) :: trans
	double precision :: xval, yval, zval
	integer :: dimn, i, j, k, l
	logical :: swid

	!	Simple basis transformation for moment integrals
	!	"trans" MUST be an orthogonal coefficient matrix, preferably normalised
	!	TBM: Generalise for n-D instead of 3-D, if ever necessary 

	if(swid)	then
		do i = 1, dimn
			do j = 1, dimn
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
				if(present(Cur1eINDrd))	Cur1eINDrd(i, j) = sqrt(xval**2 + yval**2 + zval**2)
			end do
		end do
	else
		do i = 1, dimn
			do j = 1, dimn
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
				if(present(Cur1eINDrd))	Cur1eINDrd(i, j) = sqrt(xval**2 + yval**2 + zval**2)
			end do
		end do
	end if
end subroutine TransMmtIntglVal

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