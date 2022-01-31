program Tiling
implicit none

real, allocatable       :: divisors(:), newdivisors(:), ratiocheck(:), refratio(:) ! divisors vector
integer                 :: num_cores, new_num_cores, allocated_cores, step
integer                 :: nnx, nny, nnz, factor1, factor2, factor3, i, num_divisors, idxRatio, idxMin, count
real                    :: width, height, depth, domratio, smallratio, smallratnum, dx, dy,dz, paddist, dimlengths(3)
logical                 :: wide, tall, deep, xy, xz, yz


                           ! ENTER YOUR SIMULATION INFORMATION HERE !
! ============================================================================================== !

real,   parameter       :: xleftlim     = 0, xrightlim  = 200 ! global x domain limits
real,   parameter       :: ybottomlim   = 0, ytoplim    = 200 ! global y domain limits
real,   parameter       :: zbacklim    = 0, zfirstlim   = 200 ! global z domain limits
real,   parameter       :: dt           = 0.1 ! time step
real,   parameter       :: D0 = 1             ! constant coefficient for total diffusion in the system
real,   parameter       :: pad_const = 3      ! number of standard deviations for colocation prob. search
real,   parameter       :: eff = 0.5         ! user-defined desired efficiency

allocated_cores = 1006 ! total number of cores the user has access to

! ============================================================================================== !
        ! NOW SIT BACK AND RELAX -- KEEP HANDS AND FEET INSIDE THE VEHICLE AT ALL TIMES !

dimlengths(1)  = xrightlim - xleftlim
dimlengths(2)  = ytoplim - ybottomlim
dimlengths(3)  = zfirstlim - zbacklim

if (dimlengths(2) >= dimlengths(1)) then

    if (dimlengths(3) >= dimlengths(1)) then
        yz = .true.
    else
        xy = .true.
    endif

else if (dimlengths(1) >= dimlengths(2)) then

    if(dimlengths(3) >= dimlengths(2)) then
        xz = .true.
    else
        xy = .true.
    endif

endif


if (yz .eqv. .true.) then

    if (dimlengths(3) >= dimlengths(2)) then
        domratio = dimlengths(3)/dimlengths(2)
        smallratio = dimlengths(1)/dimlengths(2)
    else if (dimlengths(2) >= dimlengths(3)) then
        domratio = dimlengths(2)/dimlengths(3)
        smallratio = dimlengths(1)/dimlengths(3)
    endif

else if (xy .eqv. .true.) then

    if (dimlengths(2) >= dimlengths(1)) then
        domratio = dimlengths(2)/dimlengths(1)
        smallratio = dimlengths(3)/dimlengths(1)
    else if (dimlengths(1) >= dimlengths(2)) then
        domratio = dimlengths(1)/dimlengths(2)
        smallratio = dimlengths(3)/dimlengths(2)
    endif

else if (xz .eqv. .true.) then

    if (dimlengths(1) >= dimlengths(3)) then
        domratio = dimlengths(1)/dimlengths(3)
        smallratio = dimlengths(2)/dimlengths(3)
    else if (dimlengths(3) >= dimlengths(1)) then
        domratio = dimlengths(3)/dimlengths(1)
        smallratio = dimlengths(2)/dimlengths(1)
    endif

endif

smallratnum = (1.0/domratio)*floor((real(allocated_cores)**(1.0/3.0)))

!Write(*,*) "smallratnum is",smallratnum

Write(*,*) "Calculating your best tiling method..."


step = 0

do num_cores = allocated_cores, 1, -1

    
    step = step + 1
    allocate(divisors(1:num_cores), newdivisors(1:num_cores), ratiocheck(1:num_cores))

    newdivisors = 0
    divisors = 0
    factor1 = 0
    factor2 = 0
    factor3 = 0

    do i = 1, num_cores
        if(mod(num_cores,i) == 0) then
            divisors(i) = i
        else
            divisors(i) = 0
        endif
    enddo

    newdivisors = pack(divisors,divisors/=0)

    idxMin = minloc(abs(newdivisors-smallratnum),1)

    factor3 = newdivisors(idxMin)

    new_num_cores = num_cores/factor3

    divisors = 0

! Find divisors of the number of cores and pack them into an array

    do i = 1, new_num_cores
        if(mod(new_num_cores,i) == 0) then
            divisors(i) = i
        else
            divisors(i) = 0
        endif
    enddo

    divisors = pack(divisors,divisors/=0)

    num_divisors = size(divisors,1)
    ratiocheck = 0

! Create an array of ratios between factor pairs.
! "else" case is used for perfect squares.

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

! Now find the pair of divisors whose ratio most closely resembles the domain's dimension's ratio 

    idxRatio = minloc(abs(ratiocheck-domratio),1)

! By construction, factor1 will always be smaller.

    factor1 = divisors(idxRatio)
    factor2 = new_num_cores/factor1

! Now that the factors have been determined, set the smaller factor in the shorter direction of the domain.

    if (xy .eqv. .true.) then

        if (factor1 <= factor3) then
            nnz = factor1

            if (dimlengths(1) <= dimlengths(2)) then
                nnx = factor3
                nny = factor2
            else if (dimlengths(2) <= dimlengths(1)) then
                nny = factor3
                nnx = factor2
            endif

        else if (factor3 <= factor1) then
            nnz = factor3

            if (dimlengths(1) <= dimlengths(2)) then
                nnx = factor1
                nny = factor2
            else if (dimlengths(2) <= dimlengths(1)) then
                nny = factor1
                nnx = factor2
            endif

        endif


    else if (xz .eqv. .true.) then

        if (factor1 <= factor3) then
            nny = factor1

            if (dimlengths(1) <= dimlengths(3)) then
                nnx = factor3
                nnz = factor2
            else if (dimlengths(3) <= dimlengths(1)) then
                nnz = factor3
                nnx = factor2
            endif

        else if (factor3 <= factor1) then
            nny = factor3

            if (dimlengths(1) <= dimlengths(3)) then
                nnx = factor1
                nnz = factor2
            else if (dimlengths(3) <= dimlengths(1)) then
                nnz = factor1
                nnx = factor2
            endif

        endif

    else if (yz .eqv. .true.) then

        if (factor1 <= factor3) then
            nnx = factor1

            if (dimlengths(2) <= dimlengths(3)) then
                nny = factor3
                nnz = factor2
            else if (dimlengths(3) <= dimlengths(2)) then
                nnz = factor3
                nny = factor2
            endif

        else if (factor3 <= factor1) then
            nnx = factor3

            if (dimlengths(2) <= dimlengths(3)) then
                nny = factor1
                nnz = factor2
            else if (dimlengths(3) <= dimlengths(2)) then
                nnz = factor1
                nny = factor2
            endif

        endif

    endif
    
    dx = dimlengths(1)/real(nnx,8)
    dy = dimlengths(2)/real(nny,8)
    dz = dimlengths(3)/real(nnz,8)

    paddist = pad_const * sqrt(4 * D0 * dt)

    deallocate(divisors, newdivisors, ratiocheck)


    ! Check to ensure that the paddist is not more than 20% of any grid box (check grid smallest dimension)
    ! If condition is violated, program will take away a core and try a factorization at the next biggest core.
    ! This will repeat until a suitable tiling is found. 

    if ((paddist/min(dx,dy,dz)) <= (0.5*((1.0/eff)**(1.0/3.0)-1.0))) then
        exit
    endif


    ! Can uncomment and print at top of loop to keep track of failed core numbers
    ! count = count + 1

enddo 


    Write(*,*) "Given the",allocated_cores,"cores you have allocated, we think you should use",&
    num_cores,"cores, tiling the domain with nnx =",nnx,", nny =",nny,", and nnz =",nnz



end program Tiling
