program DDC_2D_mpi_paper
use mod_DDC_2D_mpi_paper
implicit none

! ==============================================================================
!                           SIMULATION PARAMETERS
! ==============================================================================

real(pDki), parameter :: xleftlim  = 0.0_pDki, xrightlim  = 1000.0_pDki ! global x domain limits
real(pDki), parameter :: ybottomlim= 0.0_pDki, ytoplim  = 1000.0_pDki ! global y domain limits
real(pDki), parameter :: xmidpt    = (xrightlim + xleftlim) / 2.0_pDki ! x midpoint of the domain (used for the initial condition)
real(pDki), parameter :: ymidpt    = (ytoplim + ybottomlim) / 2.0_pDki ! y midpoint of the domain (used for the initial condition)
real(pDki), parameter :: D0        = 1e0_pDki ! coefficient of total diffusion in the system
real(pDki), parameter :: kappa_RW  = 0.5_pDki ! amount of diffusion to be simulated by random walks
real(pDki), parameter :: DRW       = kappa_RW * D0 ! diffusion coefficient for the random walks
real(pDki), parameter :: DMT       = (1.0_pDki - kappa_RW) * D0 ! amount of diffusion to be simulated by mass transfers
real(pDki), parameter :: beta      = 1 ! beta parameter, encoding the bandwidth of the mass transfer kernel
real(pDki), parameter :: maxtime   = 1e1_pDki ! simulation end time (used for RMSE calculation)
real(pDki), parameter :: cut_const = 6.0_pDki ! multiplier for cutdist
real(pDki), parameter :: pad_const = 6.0_pDki ! multiplier for paddist
integer,    parameter :: Num_DDCMT = 1 ! corresponds to the 2 mass transfer schemes (i.e., SPaM = 1, PLSS = 2)
integer,    parameter :: DDC_WhichWay = 2 ! corresponds to the 2 decomposition schemes (i.e., Vertical Slices = 1, Checkerboard = 2)
integer,    parameter :: Np                = 10e6 ! number of particles in the simulation
real(pDki), parameter :: dt                = 1e-1_pDki ! time step size
! real(pDki)            :: dt ! time step size 


! ==============================================================================
!                             GENERAL VARIABLES
! ==============================================================================
! integer(kind=8)                 :: random
integer                         :: i, Np_idx, ens, tstep, DDC_MTflag ! loop iterators
integer                         :: errcode, ierr ! error codes
real(pDki)                      :: time ! time step
integer                         :: Nsteps ! number of time steps to take

type(ParticleType), allocatable :: p(:) ! ParticleType array (these are the particles)

integer                         :: nnx, nny ! number of horizontal/vertical grid nodes
integer                         :: NDom ! number of subdomains (set below to be the number of cores in the MPI pool)
integer                         :: rem ! modulo remainder (used for load balancing, below)
integer                         :: coreNp, pVecSize ! number of particles initially on a core and the (larger) size the particle vector will be allocated to
real(pDki)                      :: Dlims(4) ! subdomain boundary limits

real(pDki)                      :: tic, toc ! timer variables for tracking run time
real(pDki)                      :: error, totError ! RMSE at final time (core local, and total)
real(pDki),         allocatable :: analytical(:) !analytical solution array
real(pDki),         allocatable :: errvec(:) ! vector for error calculation
real(pDki)                      :: dx, dy, width, height, minDTcheck ! horizontal/vertical grid spacing
real(pDki)                      :: cutdist ! number of standard deviations for kD search
real(pDki)                      :: paddist ! number of standard deviations for subdomain pad
real(pDki),           parameter :: desiredeff = 0.75_pDki ! how efficient the user wants to be
real(pDki)                      :: eff ! how efficient the user can expect to be
real(pDki)                      :: proc ! suggested number of processors

real(pDki)                      :: ds, factor, denom, spillX, spillY ! these are used for co-location probability calculations
integer,            allocatable :: idx(:), idxActive(:) ! these are used to get the indices of active particles and idxActive is used as an argument to the particle array
integer                         :: Nactive ! number of active particles
integer                         :: xmultiplier, ymultiplier ! for the local domain limits
logical                         :: mask, broadcast ! for finding the particle to give the initial mass 
integer                         :: jumpedTotSum, ghostTotSum, totSendSum ! various numbers of interesting sends in swapDomains

! ==============================================================================
!                             READ/WRITE VARIABLES
! ==============================================================================

! unit number and filename for writing the run time arrays
integer                 :: sizewrite
integer,      parameter :: uTimeDDCp1       = 11
character(*), parameter :: timeDDCNamep1    = 'PDDC_timesp1.dat'

! unit number and filename for writing the dt arrays
integer,      parameter :: udtDDC       = 12
character(*), parameter :: dtDDCName    = 'PDDC_dt.dat'

! unit numbers and filenames for writing the MSE error numbers
integer,      parameter :: uErrDDCpSize1        = 21
character(*), parameter :: ErrDDCNamepSize1     = 'PDDC_Error1.dat'

! unit number and filename for writing the error arrays
integer,      parameter :: uErrVecDDC        = 24
character(*), parameter :: ErrVecDDCName     = 'PDDC_Error_Vec.dat'

! ==============================================================================
! Plotting Variables (uncomment this block if you want to plot positions/masses)
! ==============================================================================

! particle array to hold all particles on master processor
type(ParticleType), allocatable :: masterP(:)
integer                         :: num, begin ! used for indexing in the masterP array

! unit numbers and filenames for writing the mass/location arrays
integer,      parameter :: uLocX     = 34
integer,      parameter :: uLocY     = 35
integer,      parameter :: uMass    = 36
character(*), parameter :: locXName  = 'locsX.dat', locYName = 'locsY.dat'
character(*), parameter :: massName = 'mass.dat'

! ==============================================================================

! initialize the random seed right off the bat
call init_random_seed()

! initialize the mpi session
call mpi_init(ierror)
call mpi_comm_size(mpi_comm_world, num_cores, ierror)
call mpi_comm_rank(mpi_comm_world, my_rank, ierror)
call mpi_get_processor_name(procname, namelen, ierror)

! create an mpi real type to correspond to pDki
call create_MPIRealType(D0, realType)
! create and commit a derived mpi type to correspond to the ParticleType
call build_derived_pType(mpi_pType)


! ! print summary information to screen
! if (my_rank == master) then
!     write (*, *)
!     print *, '==================== Run summary information ===================='
!     write (*, "((i2, 1x), 'particle count(s),')") Np
!     write (*, "((i2, 1x), 'DDC MT mode(s)')") Num_DDCMT
!     print *, '================================================================='
!     write (*, *)
! endif

! ============================== Wake up the CPU ===============================
if (my_rank == master) then
    print *, 'Wait a moment, this will make the CPU work out a bit first...'
    ! Write(*,*) "Size of random is",huge(random)
    ! Write(*,*) "Size of int(Np,8) is",huge(int(Np,8))
    write (*, *)
endif

! Start timer for setup time
tic = mpi_wtime()


select case (DDC_WhichWay)
    case (1) ! Vertical Slices

        nny = 1
        nnx = num_cores

    case (2) ! Checkerboard

        ! INPUT SIMULATION PARAMETERS INTO TILING2D.F90 TO GET TILING SUGGESTIONS TO PUT HERE
        ! Remember the nny*nnx must equal the number of cores with which you submit the slurm job
        nny = 1
        nnx = 1

    case default
        print *, '**** ERROR: Invalid transfer mode ****'

end select


! grid spacing
dx = (xrightlim-xleftlim)/real(nnx,pDki)
dy = (ytoplim-ybottomlim)/real(nny,pDki)

height = ytoplim - ybottomlim
width  = xrightlim - xleftlim
minDTcheck = beta*(width*height)/(2.0_pDki*D0*real(Np,pDki))
! dt = beta*(width*height)/(2.0_pDki*D0*real(Np,pDki))


if (my_rank == master) then
    print *, '--- Starting runs now ---'
    write (*, *)
    if (DDC_WhichWay == 1) then
        Write(*,*) "This is the vertical slices DDC method."
    else
        Write(*,*) "This is the checkerboard DDC method."
        Write(*,*) "minDT cutoff is",minDTcheck
        Write(*,*) "dt for this simulation is",dt
        Write(*,*) "Total simulation time is",maxtime
    endif

    if (dt >= minDTcheck) then
        Write(*,*) "Chosen time step is sufficiently large for the simulation parameters."
    else
        Write(*,*) "A larger time step is suggested for the given simulation parameter."
        Write(*,*) "We suggest using a time step larger than",minDTcheck
    endif

endif

    ! try and load balance the initial particles per core
    rem = mod(Np, num_cores)
    coreNp = Np / num_cores
    if (my_rank < rem) then
        coreNp = coreNp + 1
    endif

    ! find multipliers for determining local subdomain limits
    xmultiplier = floor(real(my_rank,pDki)/real(nny,pDki))
    ymultiplier = mod(my_rank,nny)


    ! each core's horizontal domain limits
    Dlims(1) = real(xmultiplier, pDki) * dx
    Dlims(2) = Dlims(1) + dx


    ! now calculate each core's vertical domain limits
    select case(DDC_WhichWay)

    !vertical slices
    case(1)

        Dlims(3) = ybottomlim
        Dlims(4) = ytoplim


    !checkerboard    
    case(2)

        Dlims(3) = real(ymultiplier, pDki) * dy
        Dlims(4) = real(ymultiplier + 1, pDki) * dy

    end select


    ! allocate the local p vectors some factor larger, to account for ghost
    ! particles and more particles in one domain than the other after random walks start

    ! NOTE: this is currently an ad hoc solution and makes them unnecessarily large
        ! there is likely a better way to do this. however, everything in the
        ! Engdahl, et al., "Accelerating and Parallelizing... "
        ! paper will definitely run using this
        if (num_cores == 1) then

            pVecSize = Np

        else

            ! pVecSize = ceiling(coreNp + (2.0_pDki * paddist) / (Dlims(2) - Dlims(1)) * 1.2_pDki * coreNp)
            pVecSize = 4*coreNp

        endif


    ! number of subdomains is equal to the number of cores
    NDom = num_cores

    allocate(p(pVecSize), idx(pVecSize), idxActive(pVecSize))
   

    ! initialize all variables for good measure each run
    ! masterP%mass = 0.0_pDki
    p%mass = 0.0_pDki
    p%loc(1) = 0.0_pDki
    p%loc(2) = 0.0_pDki
    error = 0.0_pDki
    totError = 0.0_pDki
    denom = 0.0_pDki
    factor = 0.0_pDki
    cutdist = 0.0_pDki
    paddist = 0.0_pDki
    Nsteps = 0
    time = 0
    Nactive = 0
    spillX = 0
    spillY = 0

    p%active = .false.
    p%ghost(1) = .false.
    p%ghost(2) = .false.
    p%ghost(3) = .false.
    p%ghost(4) = .false.
    p%ghost(5) = .false.
    p%ghost(6) = .false.
    p%ghost(7) = .false.
    p%ghost(8) = .false.
    p%jumped(1) = .false.
    p%jumped(2) = .false.
    p%jumped(3) = .false.
    p%jumped(4) = .false.
    p%jumped(5) = .false.
    p%jumped(6) = .false.
    p%jumped(7) = .false.
    p%jumped(8) = .false.

    ! master array of indices, to be used for logical indexing
    idx = (/(i, i = 1, pVecSize)/)


    ! number of time steps
    Nsteps    = maxtime/real(dt,pDki)

    ! search distance for MT algorithm (\psi in manuscript)
    cutdist   = cut_const * sqrt(2.0_pDki * DMT * dt)
    paddist   = pad_const * sqrt(2.0_pDki * DMT * dt)
    
    eff = num_cores**(-1.0_pDki)*((sqrt(real(num_cores,pDki))**(-1.0_pDki)+(2.0_pDki*cutdist)/width))**(-2.0_pDki)
    proc = (desiredeff)**(-1.0_pDki)*(((1.0_pDki-sqrt(desiredeff))*width)/(2.0_pDki*cutdist))**(2.0_pDki)
    

    ! these are for the co-location probability calculations
        ! while ds and factor are not strictly necessary, they are included to
        ! match the "Accelerating and Parallelizing... " manuscript
    ds     = (xrightlim - xleftlim) / Np
    factor = ds / sqrt(4.0_pDki * pi * DMT * dt)
    denom  = (1 / beta) * (-4.0_pDki) * DMT * dt

    do DDC_MTflag = 1, Num_DDCMT ! mass-transfer mode loop (SPaM = 1, PLSS = 2)

                ! Each core randomly spaces all of its particles in its local subdomain limits
                ! and then imposes the given initial condition 

                call InitialSpacingAndIC(p, Dlims, coreNp, xmidpt, ymidpt, nnx, nny, spillX, spillY)


                ! if you are using delta initial condition in InitialSpacingAndIC, uncomment the two lines below
                ! to share the initial spill location with all cores for error calculation later

                ! call MPI_Bcast(spillX, 1, mpi_double_precision, (floor(nnx/real(2,pDki))*nny + floor(nny/real(2,pDki))), &
                !     mpi_comm_world, ierror)
                ! call MPI_Bcast(spillY, 1, mpi_double_precision, (floor(nnx/real(2,pDki))*nny + floor(nny/real(2,pDki))), &
                !     mpi_comm_world, ierror)    


             ! end the clock for setup time
             toc = mpi_wtime()
             if (my_rank == master) then

                Write(*,*) "Setup time is",toc-tic
                Write(*,*) "The user can expect to be",eff,"percent efficient with these simulation parameters."
                Write(*,*) "To remain",desiredeff,"efficient, we suggest not using more than",proc,"processors."

             !    ! If you want to see where the spill is located
                !   Write(*,*) "spillX is",spillX
                !   Write(*,*) "spillY is",spillY

             endif


            ! ============================================================================
            !  ***Uncomment this block if you want to plot the initial positions/masses***
            ! ============================================================================

            !! write out the plot information with initial == true
            !! This takes a lot of time with large particle numbers


            !! get indices of the currently active particles

            !  if (my_rank == master) then
            !      allocate(masterP(1 : Np))
            !      masterP(1 : coreNp) = p(1 : coreNp)
            !  endif

            ! idxActive = 0
            ! Nactive = count(p%active)
            ! idxActive(1 : Nactive) = pack(idx, p%active)

            !         begin = Nactive + 1
            !         do i = 1, num_cores - 1

            !             if (my_rank == i) then
            !                 call mpi_send(Nactive, 1, mpi_integer, master, tag, mpi_comm_world, ierror)
            !             endif
            !             if (my_rank == master) then
            !                 call mpi_recv(num, 1, mpi_integer, i, tag, mpi_comm_world, status, ierror)
            !             endif
            !             if (my_rank == i) then
            !                 call mpi_send(p(idxActive(1 : Nactive)), Nactive, mpi_pType, master, tag, mpi_comm_world, ierror)
            !             endif
            !             if (my_rank == master) then
            !                 call mpi_recv(masterP(begin : begin + num - 1), num, mpi_pType,&
            !                               i, tag, mpi_comm_world, status, ierror)
            !             endif
            !           begin = begin + num

            !         enddo

            !         if (my_rank == master) then

            !                  call write_plot( uLocX,  locXName,  masterP%loc(1), .true.)
            !                  call write_plot( uLocY,  locYName,  masterP%loc(2), .true.)
            !                  call write_plot( uMass, massName, masterP%mass, .true.)
            !                 Write(*,*) "Total mass at initial write is",sum(masterP%mass)

            !                 ! Write(*,*) "After initial write"

            !         endif

            !  if (my_rank == master) then
            !      deallocate(masterP)
            !  endif

             ! ==================================================================
             ! ==================================================================


            ! start the clock for recording run time
            tic = mpi_wtime()
            
            if(my_rank == master) then
                ! Write (*,*) "cutdist is",cutdist
                ! Write (*,*) "Nactive is",Nactive
                Write (*,*) "Nsteps is",Nsteps
            endif
            

            do tstep = 1, Nsteps ! time stepping loop

                ! print time step if you want to track real-time progress
                if (my_rank == master) then
                    ! if (mod(tstep,10) == 0) then
                        Write(*,*) "tstep is",tstep
                    ! endif
                endif

                ! get indices of the currently active particles
                idxActive = 0
                Nactive = count(p%active)
                idxActive(1 : Nactive) = pack(idx, p%active)

                ! random walk diffusion
                   call diffuse(Nactive, DRW, dt, idxActive(1 : Nactive), p)

                ! reflecting boundary conditions

                   call reflectLow(xleftlim, idxActive(1 : Nactive), p, 1)
                   call reflectLow(ybottomlim, idxActive(1 : Nactive), p, 2)

                   call reflectHigh(xrightlim, idxActive(1 : Nactive), p, 1)
                   call reflectHigh(ytoplim, idxActive(1 : Nactive), p, 2)


                ! swap any particles that have crossed subdomain boundaries or
                ! are ghost particles to the neighboring subdomain/processor
                if (num_cores > 1) then
                    call swapDomains(pVecSize, Nactive, idxActive, Dlims, paddist, p, nnx, &
                                nny, DDC_WhichWay, tstep, jumpedTotSum, ghostTotSum, totSendSum)
                endif

                ! now that particles have been swapped between subdomains,
                ! recalculate the indices of currently active particles for MT algorithm
                idxActive = 0
                Nactive = count(p%active)
                idxActive(1 : Nactive) = pack(idx, p%active)
                
               

                ! do the mass transfers
                select case (DDC_MTflag)
                    case (1) ! SPaM
                        call massTrans_kDMat(Nactive, idxActive, cutdist, denom, p, beta)
                    case (2) ! PLSS
                        call massTrans_kDLoop(Nactive, idxActive, cutdist, factor, denom, p)
                    case default
                        print *, '**** ERROR: Invalid transfer mode ****'
                end select

                ! deactivate the ghost particles, now that mass transfers have occurred
                if (num_cores > 1) then
                    where (p%ghost(1)) p%active = .false.
                    where (p%ghost(2)) p%active = .false.
                    where (p%ghost(3)) p%active = .false.
                    where (p%ghost(4)) p%active = .false.
                    where (p%ghost(5)) p%active = .false.
                    where (p%ghost(6)) p%active = .false.
                    where (p%ghost(7)) p%active = .false.
                    where (p%ghost(8)) p%active = .false.
                endif


                ! ================================================================================
                !***Uncomment this block if you want to plot positions/masses at each time step***
                ! ================================================================================

                ! send all the particles to master for plotting write-out

            !     if (mod(tstep,20) == 0) then
            !   !  if (tstep == Nsteps) then

            !         if (my_rank == master) then
            !          allocate(masterP(1 : Np))
            !         endif

            !         idxActive = 0
            !         Nactive = count(p%active)
            !         idxActive(1 : Nactive) = pack(idx, p%active)


            !          if (my_rank == master) then
            !              masterP(1 : Nactive) = p(idxActive(1 : Nactive))
            !          endif
            !          !===============================================================

            !         begin = Nactive + 1
            !         do i = 1, num_cores - 1

            !             if (my_rank == i) then
            !                 call mpi_send(Nactive, 1, mpi_integer, master, tag, mpi_comm_world, ierror)
            !             endif
            !             if (my_rank == master) then
            !                 call mpi_recv(num, 1, mpi_integer, i, tag, mpi_comm_world, status, ierror)
            !             endif
            !             if (my_rank == i) then
            !                 call mpi_send(p(idxActive(1 : Nactive)), Nactive, mpi_pType, master, tag, mpi_comm_world, ierror)
            !             endif
            !             if (my_rank == master) then
            !                 call mpi_recv(masterP(begin : begin + num - 1), num, mpi_pType,&
            !                               i, tag, mpi_comm_world, status, ierror)
            !             endif
            !           begin = begin + num

            !         enddo

            !         if (my_rank == master) then

            !                  call write_plot( uLocX,  locXName,  masterP%loc(1), .false.)
            !                  call write_plot( uLocY,  locYName,  masterP%loc(2), .false.)
            !                  call write_plot( uMass, massName, masterP%mass, .false.)
            !                 !  Write(*,*) "Done with write at time step",tstep
            !                  Write(*,*) "Total mass at time step",tstep,"is",sum(masterP%mass)

            !         endif

            !         if (my_rank == master .and. tstep /= Nsteps) then
            !              deallocate(masterP)
            !         endif

            !     endif
                 
                ! ==============================================================
                ! ==============================================================

            enddo ! time stepping loop

            ! if (my_rank == master) then
            !     Write(*,*) "After time loop"
            ! endif

            ! stock the clock for recording run time
            toc = mpi_wtime()

            ! write to screen some summary info and the run time
            if (my_rank == master) then
                select case (DDC_MTflag)
                    case (1) ! SPaM
                        write (*, "('SPaM', (i4, 1x), '(MT mode ', (i1, 1x), '): N = ', (i10, 1x))") NDom, DDC_MTflag, Np
                        write (*, *) 'run time = ', toc - tic
                    case (2) ! PLSS
                        write (*, "('PLSS', (i4, 1x), '(MT mode ', (i1, 1x), '): N = ', (i7, 1x),&
                               &' Ensemble number = ', (i2, 1x))") NDom, DDC_MTflag, Np, ens
                        write (*, *) 'run time = ', toc - tic
                    case default
                        print *, '**** ERROR: Invalid transfer mode ****'
                end select
            endif

            ! write the run time information to file
            ! if (my_rank == master) then
            !             call write_time(uTimeDDCp1, timeDDCNamep1, toc - tic,&
            !                         .true.)
            ! endif  



            ! compute error from analytic solution in parallel 

            ! first compute individual "MSEs" on each core
            ! NOTE: these are divided by the TOTAL particle number so are not
                ! actually local mean-squared-errors


            ! should be equal to maxtime, unless mod(maxtime,dt) is not an integer
            time = real(Nsteps,pDki) * dt

            if (my_rank == master) then
                ! Write(*,*) "Total simulation time is",time
                ! Write(*,*) "dt for this simulation is",dt
                allocate(errvec(1:num_cores))
            endif

            ! determine active particles for error calculation
            idxActive = 0
            Nactive = count(p%active)
            idxActive(1 : Nactive) = pack(idx, p%active)


            ! allocate local analytical solution vector
            allocate(analytical(pVecSize))
            analytical = 0.0_pDki


            ! Heaviside analytic solution
            ! analytical = 0.5_pDki * erfc(-(p(idxActive(1 : Nactive))%loc(1) - xmidpt)/sqrt(4.0_pDki * D0 * time))


            ! Delta function analytic solution for a single spill 
            analytical = exp(-1.0_pDki*((p(idxActive(1 : Nactive))%loc(1)-spillX)**2+&
                (p(idxActive(1 : Nactive))%loc(2)-spillY)**2)/(4.0_pDki*D0*time))/(4.0_pDki*pi*D0*time)
                
            analytical = analytical/sum(analytical)


            ! local "MSEs"
            error = sum((analytical-p(idxActive(1 : Nactive))%mass)**2)
            
            ! Write(*,*) "The error on rank",my_rank,"is",error

            if(my_rank == master) then
                errvec(1) = error
            endif
            
            do i = 1, num_cores-1
            
                if (my_rank == i) then
                        call mpi_send(error, 1, mpi_double_precision, master, tag, mpi_comm_world, ierror)
                endif
                
                if (my_rank == master) then
                        call mpi_recv(errvec(i+1), 1, mpi_double_precision, i, tag, mpi_comm_world, status, ierror)
                endif
            
            enddo


            ! sum all the individual core's "MSEs" on lead core
            ! call mpi_reduce(error, totError, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierror)



            if (my_rank == master) then

                ! Write(*,*) "global error sum before square root is",totError

                ! totError = sqrt(totError)
                
                totError = sqrt(sum(errvec)/real(Np,pDki))

                ! print the error to screen
                ! print *,'totError = ',totError
                ! print *,'================================'

                ! write the error and dt information to file if desired

                !         call write_error(uErrDDCpSize1, ErrDDCNamepSize1, totError, .true.)
                !         call write_error(udtDDC, dtDDCName, dt, .true.)


            endif



    enddo ! mass-transfer mode loop

    deallocate(p, idx, idxActive, analytical)


call mpi_finalize(ierror)

end program DDC_2D_mpi_paper