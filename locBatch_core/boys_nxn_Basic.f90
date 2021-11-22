module boysParams
	type locSet
		double precision :: A_ij, B_ij, cfgamma, sfgamma, cgamma, sgamma, gamma, gain, count
		logical :: thrshmet
	end type locset
	double precision, parameter :: pi = 3.14159265358979323846264338327950D+00
end module boysParams

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
    double precision, dimension(rdim, cdim) :: stdmat
    integer :: rDim, cDim, ldR, ldC, i, j
    character (len = 64) :: lCaption

    do i = ldR, rDim
        write (*, '(*(ES12.4E4))') (stdmat(i, j), j = ldC, cDim)
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

subroutine byDiag()
	use omp_lib
	use boysParams
	use externAccess

	implicit none

	!  Calls gaussin() from input_setup.f90 - ensure in $PWD

	type(locSet), dimension(:, :), allocatable :: localiser
	double precision, dimension(:, :), allocatable :: mcf
	double precision, parameter :: absThrsh = 1.0e-2
	double precision :: BoysFnl
	integer :: dimn, iter, boys
	
	!  Read MO coeffs using gaussin()
	call copyInp(mcf, dimn)
	
	if(.not. allocated(mcf))	allocate(mcf(dimn, dimn))
	allocate(localiser(dimn, dimn))

	iter = Boys(mcf, dimn, localiser, absThrsh, BoysFnl)

	write(*, 100) "Boys functional and iterations: ", BoysFnl, ", ", iter

	100	format(A31, ES12.4E4, A2, I8)

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

		if(ntot .gt. 0) then
			dimn = ntot
		else
		!	Do something, write to an error file. Up to the user.
		endif

		if(.not. allocated(mcf))	allocate(mcf(dimn, dimn))

		mcf = mocoef
	end subroutine copyInp
end subroutine byDiag

integer function boys(mcf, dimn, localiser, absThrsh, BoysFnl)	result(iter)
	use boysParams
	use externAccess

	implicit none
	
	double precision, dimension(dimn, dimn, 3) :: AO1eIND, MO1eIND, Lcl1eIND
	double precision, dimension(dimn, dimn) :: mcf, orth, orthPrior, gemmResult, verifMat, MO1eINDrd, Lcl1eINDrd
	double precision, dimension(dimn) :: eival
	double precision :: BoysFnl, maxFnl, iterGain, absThrsh
	type(locSet), dimension(dimn, dimn) :: localiser
	integer :: i, j, k, l, dimn, iter
	!	Clamp at known convergence point for now
	integer, parameter :: iterlimit = 1e1
	logical :: scalConverged, isOrth, orthTest, isConverged
	
	call AOMmtIntglVal(dimn, AO1eIND)
	call TransMmtIntglVal(dimn, mcf, AO1eIND, MO1eIND, MO1eINDrd)
	call initBoysParams(dimn, MO1eIND, MO1eINDrd, localiser, orth, orthPrior)

	BoysFnl = 0.d0
	iter = 0
	
	do i = 1, dimn
		BoysFnl = BoysFnl + MO1eINDrd(i, i)**2
	end do
	
	write(*, '(A24, ES12.4E4)') "Initial Boys Functional: ", BoysFnl

	gemmResult = orthPrior

	do while(.not.(iter .and. scalConverged(BoysFnl, BoysFnl - iterGain, absThrsh)) .and. iter .lt. iterlimit)
		do i = 1, dimn
			!	Strict upper half, invertible
			do j = i + 1, dimn
				orthPrior = gemmResult
				!	Check if localisation threshold is met
				if(localiser(i, j)%thrshmet)	iterGain = iterGain + localiser(i, j)%gain
				!	Update integral values for this iteration
				call TransMmtIntglVal(dimn, orth, MO1eIND, Lcl1eIND, Lcl1eINDrd)
				!	Update values of Aij, Bij from integral values in current basis
				localiser(i, j)%A_ij = Lcl1eINDrd(i, j)**2 - 0.25*dot_product(	&
					(Lcl1eIND(i, i, :) - Lcl1eIND(j, j, :)),	&
					(Lcl1eIND(i, i, :) - Lcl1eIND(j, j, :))		&
				)**2
				
				localiser(i, j)%B_ij = dot_product(	&
					(Lcl1eIND(i, i, :) - Lcl1eIND(j, j, :)),	&
					Lcl1eIND(i, j, :)	&
				)
				!	Determine rotation angles
				localiser(i, j)%cfgamma	= -localiser(i, j)%A_ij/(sqrt(localiser(i, j)%A_ij**2 + localiser(i, j)%B_ij**2))
            	localiser(i, j)%sfgamma	=  localiser(i, j)%B_ij/(sqrt(localiser(i, j)%A_ij**2 + localiser(i, j)%B_ij**2))
				
				localiser(i, j)%count = localiser(i, j)%count + 1
				!	Leave filter alone, unless necessary
                !	localiser(i, j)%thrshmet = .true.
				
				!	localiser(i, j)%cgamma = sqrt((1 + sqrt((1 + localiser(i, j)%cfgamma**2)/2.0))/2.0)
				!	localiser(i, j)%sgamma = sqrt(1.d0 - localiser(i, j)%cgamma**2)
				
				localiser(i, j)%gamma = 0.25*atan(localiser(i, j)%sfgamma/localiser(i, j)%cfgamma)
				!	print *, "Rotation: ", localiser(i, j)%gamma
			
				localiser(i, j)%cgamma = cos(localiser(i, j)%gamma)
				localiser(i, j)%sgamma = sin(localiser(i, j)%gamma)
                
                orth(i, i) =  localiser(i, j)%cgamma
				orth(i, j) =  localiser(i, j)%sgamma
				orth(j, i) = -localiser(i, j)%sgamma
                orth(j, j) =  localiser(i, j)%cgamma
                
				call dgemm('N', 'N', dimn, dimn, dimn, 1.d0, orth, dimn, orthPrior, dimn, 0.d0, gemmResult, dimn)
				!	Reset values for next iteration
				orth(i, i) = 1.d0
				orth(i, j) = 0.d0
				orth(j, i) = 0.d0
				orth(j, j) = 1.d0
				!	print *, iterGain
			end do
        end do
        !	print *, "Orthogonality test: ", isOrth(gemmResult, dimn, absThrsh)
		maxFnl = BoysFnl
		!	print *, "Orth matrix for iter: ", iter + 1
		!	call dbOut(orthPrior, dimn, dimn, 1, 1)
		!	print *, "Converged: ", isConverged(dimn, Lcl1eIND, absThrsh)
		!	print *, "Orth test: ", isOrth(orthPrior, dimn, absThrsh)

		!	Equation 6, (1974) Boys functional
		do i = 1, dimn
			BoysFnl = BoysFnl + Lcl1eINDrd(i, i)**2
		end do
		iterGain = BoysFnl - maxFnl
		iter = iter + 1
	end do
	!	Orthogonality test, use for verification if generation method changes
	!	print *, "Orthogonality test: ", isOrth(gemmResult, dimn, absThrsh)

	call dgemm('N', 'T', dimn, dimn, dimn, 1.d0, gemmResult, dimn, mcf, dimn, 0.d0, verifMat, dimn)

	print *, "Localised Boys Coeffs"
	call dbOut(verifMat, dimn, dimn, 1, 1)

	204    format(*(ES12.4E4))
end function boys

subroutine initBoysParams(dimn, MO1eIND, MO1eINDrd, localiser, orth, orthPrior)
	use boysParams
	
	implicit none

	double precision, dimension(dimn, dimn, 3)	:: MO1eIND
	double precision, dimension(dimn, dimn) :: orth, orthPrior, gemmResult, MO1eINDrd
	type(locSet), dimension(dimn, dimn) :: localiser
	integer :: dimn, i, j

	do i = 1, dimn
		orth(i, i) = 1.d0
	end do

	gemmResult = orth

	!	Sweep F:{i, j} -> (0:n(n-1)/2) => (i, j>i) A_ij, B_ij; n = dimn
	do i = 1, dimn
		!	ΔL:(i .eq. j) = 0, for obvious reasons, consider (i, j + i) -> n(n - 1)/2
		do j = i + 1, dimn
			orthPrior = gemmResult
			
			localiser(i, j)%A_ij = MO1eINDrd(i, j)**2 - 0.25*dot_product(	&
				(MO1eIND(i, i, :) - MO1eIND(j, j, :)),	&
				(MO1eIND(i, i, :) - MO1eIND(j, j, :))	&
			)**2
			localiser(i, j)%B_ij = dot_product(	&
				(MO1eIND(i, i, :) - MO1eIND(j, j, :)),	& 
				MO1eIND(i, j, :)	&
			)
			
			localiser(i, j)%cfgamma	= -localiser(i, j)%A_ij/(sqrt(localiser(i, j)%A_ij**2 + localiser(i, j)%B_ij**2))
            localiser(i, j)%sfgamma	=  localiser(i, j)%B_ij/(sqrt(localiser(i, j)%A_ij**2 + localiser(i, j)%B_ij**2))
            
			localiser(i, j)%count = 0
			!	Set all to localise, for initial iteration
			localiser(i, j)%thrshmet = .true.
			
			!	localiser(i, j)%cgamma = sqrt((1 + sqrt((1 + localiser(i, j)%cfgamma**2)/2.0))/2.0)
			!	localiser(i, j)%sgamma = sqrt(1.d0 - localiser(i, j)%cgamma**2)
			localiser(i, j)%gamma = 0.25*atan(localiser(i, j)%sfgamma/localiser(i, j)%cfgamma)
			
			localiser(i, j)%cgamma = cos(localiser(i, j)%gamma)
			localiser(i, j)%sgamma = sin(localiser(i, j)%gamma)

            orth(i, i) =  localiser(i, j)%cgamma
            orth(i, j) =  localiser(i, j)%sgamma
            orth(j, i) = -localiser(i, j)%sgamma
			orth(j, j) =  localiser(i, j)%cgamma
			
			call dgemm('N', 'N', dimn, dimn, dimn, 1.d0, orth, dimn, orthPrior, dimn, 0.d0, gemmResult, dimn)
			
			orth(i, i) = 1.d0
			orth(i, j) = 0.d0
			orth(j, i) = 0.d0
			orth(j, j) = 1.d0
		end do
	end do
end subroutine initBoysParams

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

subroutine TransMmtIntglVal(dimn, trans, trans1eIND, Cur1eIND, Cur1eINDrd)
	use externAccess

	implicit none

	double precision, dimension(dimn, dimn, 3) :: trans1eIND, Cur1eIND
	double precision, dimension(dimn, dimn), optional :: Cur1eINDrd
	double precision, dimension(dimn, dimn) :: trans
	double precision :: xval, yval, zval
	integer :: dimn, i, j, k, l

	!	Simple basis transformation for moment integrals
	!	"trans" MUST be an orthogonal coefficient matrix, preferably normalised
	!	TBM: Generalise for n-D instead of 3-D, if ever necessary 

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
end subroutine TransMmtIntglVal

logical function isConverged(dimn, Lcl1eIND, absThrsh)
	double precision, dimension(dimn, dimn, 3) :: Lcl1eIND
	double precision ::	rKl, absThrsh
	integer :: dimn, i, j

	isConverged = .false.

	rKl = 0.d0
	do i = 1, dimn
		do j = i + 1, dimn
			!	Convergence to D(\phi), linear gradient only
			rKl = rKl + dot_product(	&
				(Lcl1eIND(i, i, :) - Lcl1eIND(j, j, :)),	&
				Lcl1eIND(i, j, :)	&
			)**2
		end do
	end do
	rKl = 4*sqrt(rKl)
	
	! Determine if gradient for line min. is within margin of error
	if(rKl .lt. absThrsh .and. rKl .gt. -absThrsh) isConverged = .true.
end function isConverged