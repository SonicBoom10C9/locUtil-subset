subroutine eiDiag()
	use omp_lib
	use params

	implicit none

	!  Calls gaussin() from input_setup.f90 - ensure in $PWD

	double precision, dimension(:, :), allocatable :: mcf
	double precision, parameter :: absThrsh = 1e-12
	double precision :: ERfunl
	integer :: dimn, iter, jacobi, occDimn, mode, verify
	logical :: dbg
	
	!  Read MO coeffs using gaussin()
	call copyInp(mcf, dimn, occDimn)
	!	0 = Occupied only, 1 = Virtual Only
	!	2 = Occupied + Virtual (output of occupied -> input of virtual)
	!	Anything else = complete basis
	write (*, *, advance='no')	"Reading mode: "
	read(*, *)	mode

	dbg = .false.	
	print *, "Print iterative coefficients? (T/F) [F]: "
	read (*, *, iostat=verify)	dbg
	
	!	if(.not. allocated(mcf))	allocate(mcf(dimn, dimn))
    
	iter = jacobi(mcf, occDimn, dimn, absThrsh, ERfunl, mode, dbg)

	write(*, 100) "ER functional and iterations: ", ERfunl, ", ", iter

	100	format(A31, ES12.4E4, A2, I8)

	contains
	subroutine copyInp(mcf, dimn, occDimn)
		use gaussprog, only: mocoef, ntot, nocc
		
		implicit none

		interface
			subroutine gaussin()
			end subroutine gaussin
		end interface

		double precision, dimension(:, :), allocatable, intent(out)	:: mcf
		integer, intent(out) :: dimn, occDimn
		integer :: i, j

		if(ntot .gt. 0) then
			dimn = ntot
			occDimn = nocc
		else
		!	Do something, write to an error file
		endif

		if(.not. allocated(mcf))	allocate(mcf(dimn, dimn))

		mcf = mocoef
	end subroutine copyInp
end subroutine eiDiag

integer function jacobi(mcf, occDimn, dimn, absThrsh, ERFunl, mode, dbg)	result(iter)
	use params

	implicit none
	
	double precision, dimension(dimn, dimn, dimn, dimn) :: AO2eIND, MO2eIND, lcl2eIND, str2eIND
	double precision, dimension(dimn, dimn) :: mcf, verifMat, R_1, R_2
	double precision, dimension(dimn) :: eival, long, lat
	double precision :: ERFunl, iterGain, absThrsh, A_ij, B_ij, gamma, gain, count, R1prev, R2prev, cfg, sfg, fracGain
	integer :: i, j, k, l, occDimn, dimn, iter, mode, start, end, spliter
	integer, parameter :: iterlimit = 1e0
	logical :: scalConverged, isOrth, dbg

	!	Get this from input, use as is for now in LiH STO-3G
	!	occDimn = 2

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
	
	call AOIntglVal(dimn, AO2eIND)
	print *, "AO Evaluation complete"
	call Intgl2eTransPS(dimn, mcf, AO2eIND, MO2eIND)
	print *, "MO integral evaluation complete"
	!	call initParams(dimn, MO2eIND, R_1, R_2)
	
	ERFunl = 0.d0
	do i = 1, dimn
		ERFunl = ERFunl + MO2eIND(i, i, i, i)
	end do
	
	write(*, 1000) "Initial ER Functional: ", ERFunl

	1000	format(A24, ES12.4E4)

	iter = 0
	verifMat = mcf
	lcl2eIND = MO2eIND

	!	Convergence is determined from if the value of the functional sits inside the threshold.
	do while(iter .lt. iterlimit .and. .not.(iter .and. scalConverged(ERFunl, ERFunl - iterGain, absThrsh)))
		iterGain = 0.d0
		!	ERFunl = 0.d0
		do i = start, end - 1
			!	Strict upper half, invertible
			do j = i + 1, end
				!	Simplified expressions, refer Theor Chim Acta (1992) 86:149-165, Section 2.12, 2.13, 4.1
				!	Now that they are initialised, R_1 and R_2 are calculated as transformations, (2.12e) and (2.12f)
				!	This is understood from Section 2.3, Case 1, points (i) and (ii)
				!	R1prev = R_1(i, j)
				!	R2prev = R_2(i, j)

				!	cfg = (-R_1(i, j)/(dsqrt(R_1(i, j)**2 + R_2(i, j)**2)))
				!	sfg = (R_2(i, j)/(dsqrt(R_1(i, j)**2 + R_2(i, j)**2)))
				
				!	R_1(i, j) = R1prev*cfg + R2prev*sfg
				!	R_2(i, j) = R2prev*cfg - R1prev*sfg
				
				!	Expanded as A_ij and B_ij for readability. Refer ER (1963)
				!	A_ij = -2*R_1i, j);	B_ij = 2*R_2(i, j)
				
				A_ij = lcl2eIND(i, i, j, j) - 0.25d0*(lcl2eIND(i, i, i, i) + lcl2eIND(j, j, j, j) - 2.d0*lcl2eIND(i, j, i, j))
				B_ij = lcl2eIND(i, i, j, i) - lcl2eIND(i, j, j, j)

				if(A_ij**2 + B_ij**2 .lt. 1D-12)	cycle

				gamma = 0.25*atan(-B_ij/A_ij)
				
				iterGain = iterGain + A_ij + dsqrt(A_ij**2 + B_ij**2)
                
				!	call dgemm('N', 'N', dimn, dimn, dimn, 1.d0, orth, dimn, orthPrior, dimn, 0.d0, gemmResult, dimn)

				long = cos(gamma)*verifMat(i, :) + sin(gamma)*verifMat(j, :) 
                lat = -sin(gamma)*verifMat(i, :) + cos(gamma)*verifMat(j, :)
				verifMat(i, :) = long
				verifMat(j, :) = lat
				!	Reset values for next iteration
			end do
		end do

		ERFunl = ERFunl + iterGain

		call Intgl2eTransPS(dimn, verifMat, lcl2eIND, str2eIND)
		!	do i = 1, dimn
			!	ERFunl = ERFunl + str2eIND(i, i, i, i)
		!	end do
		lcl2eIND = str2eIND

		!	print *, iterGain
	
		if(dbg)	then
			print *, "Coeffcients after iteration: ", iter + 1
			call gaussVxOut(verifMat, dimn, dimn, 1, 1)
		endif

		iter = iter + 1

		!	call lcl2eInt(dimn, verifMat, MO2eIND, lcl2eIND)
	end do
	!	Orthogonality test, use for verification if generation method changes
	!	print *, "Orthogonality test: ", isOrth(gemmResult, dimn, absThrsh)

	!	call dgemm('N', 'T', dimn, dimn, dimn, 1.d0, gemmResult, dimn, mcf, dimn, 0.d0, verifMat, dimn)

	if(mode .eq. 2)	then
		spliter = jacobi(verifMat, occDimn, dimn, absThrsh, ERFunl, 1, dbg)
		iter = iter + spliter
		write(*, '(A35, I4)'), " Secondary localisation iterations: ", spliter
		return
	endif

	print *, "Localised ER Coeffs"
    call gaussVxOut(verifMat, dimn, dimn, 1, 1)

	204    format(*(ES12.4E4))
end function jacobi

subroutine initParams(dimn, MO2eIND, R_1, R_2)
	use params

	implicit none

	double precision, dimension(dimn, dimn, dimn, dimn)	:: MO2eIND
	double precision, dimension(dimn, dimn) :: R_1, R_2
	integer :: dimn, i, j, k, l

	!	Sweep F:{i, j} -> (0:n(n-1)/2) => (i, j>i) A_ij, B_ij; n = dimn
	do i = 1, dimn - 1
		!	ΔL:(i .eq. j) = 0, for obvious reasons, consider (i, j + i) -> n(n - 1)/2
		do j = i + 1, dimn
			!	Expressions for R_1 (2.13e), R_2 (2.13f), first from the MO integral values
			R_1(i, j) = 0.125d0*(MO2eIND(i, i, i, i) + MO2eIND(j, j, j, j)	&
								- 2*MO2eIND(i, j, i, j) - 4*MO2eIND(i, i, j, j))
			R_2(i, j) = 0.5d0*(MO2eIND(i, i, i, j) - MO2eIND(j, i, j, j))
		end do
	end do
end subroutine initParams