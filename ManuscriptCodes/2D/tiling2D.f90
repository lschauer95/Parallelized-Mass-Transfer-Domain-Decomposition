program Tiling
implicit none

real, allocatable       :: divisors(:), ratiocheck(:), refratio(:) ! divisors vector
integer                 :: num_cores, allocated_cores
integer                 :: nnx, nny, factor1, factor2, i, num_divisors, idxRatio, count
real                    :: width, height, domratio, dx, dy, paddist
logical                 :: wide



! ENTER YOUR INFORMATION HERE !
! ==================================================== !
real,   parameter       :: xleftlim     = 0, xrightlim  = 2000 ! global x domain limits
real,   parameter       :: ybottomlim   = 0, ytoplim  = 1000 ! global y domain limits
real,   parameter       :: dt           = 0.1, D0 = 1, pad_const = 3 ! simulation parameters
real,   parameter       :: eff = 0.75 ! user-defined efficiency threshold
allocated_cores = 1006
! ==================================================== !
! NOW SIT BACK AND RELAX -- KEEP ALL HANDS AND FEET INSIDE THE VEHICLE AT ALL TIMES !



height = ytoplim - ybottomlim
width  = xrightlim - xleftlim

if (height >= width) then
    domratio = height/width
    wide = .false.
else
    domratio = width/height
    wide = .true.
endif

! Find divisors of the number of cores and pack them into an array
! This will be used to evenly distribute the rows and columns of the decomposition

Write(*,*) "Calculating your best tiling method..."



do num_cores = allocated_cores, 1, -1

    allocate(divisors(1:num_cores), ratiocheck(1:num_cores))

    divisors = 0

    do i = 1, num_cores
        if(mod(num_cores,i) == 0) then
            divisors(i) = i
        else
            divisors(i) = 0
        endif
    enddo

    divisors = pack(divisors,divisors/=0)

    num_divisors = size(divisors,1)
    ratiocheck = 0


    if (mod(num_divisors,2) == 0) then
	    do i = 1, (num_divisors/2)
            ratiocheck(i) = divisors(num_divisors - i + 1)/divisors(i)
	    enddo
	else
        do i = 1, floor(real(num_divisors/2))
            ratiocheck(i) = divisors(num_divisors - i + 1)/divisors(i)
        enddo
        ratiocheck(floor(real(num_divisors/2))+1) = 1
	endif

    ratiocheck = pack(ratiocheck,ratiocheck/=0)

    idxRatio = minloc(abs(ratiocheck-domratio),1)

    factor1 = divisors(idxRatio)
    factor2 = num_cores/factor1

    if (wide .eqv. .true.) then
        nny = factor1
    else
        nny = factor2
    endif
    nnx = num_cores/nny


    dx = width/real(nnx,8)
    dy = height/real(nny,8)

    paddist = pad_const * sqrt(4 * D0 * dt)


    count = count + 1

    deallocate(divisors, ratiocheck)

    if ((paddist/min(dx,dy)) <= (0.5*((1.0/eff)**(1.0/2.0)-1.0))) then
        exit
    endif

enddo 


    Write(*,*) "Given the",allocated_cores,"cores you have allocated, we think you should use",&
    num_cores,"cores, tiling the domain with nnx =",nnx,"and nny =",nny



end program Tiling
