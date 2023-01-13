module mod_DDC_2D_mpi_paper
use mpi
use kdtree2_module
implicit none

! ==============================================================================
!                              GLOBAL VARIABLES
! ==============================================================================

integer,    parameter :: spkind = kind(1.0), dpkind = kind(1.0d0)
! choose single or double precision, below
integer,    parameter :: pDki = dpkind
real(pDki), parameter :: pi = 4.0_pDki * atan(1.0_pDki)
integer               :: dim = 2 ! NOTE: this is the hard-coded 2 spatial dimension
real(pDki), parameter :: eps = 10**(-15)

! ==============================================================================
!                                 MPI VARIABLES
! ==============================================================================

integer                             :: ierror, num_cores, my_rank, namelen
integer, parameter                  :: master = 0
integer, parameter                  :: tag = 8273
integer, dimension(mpi_status_size) :: status
character(mpi_max_processor_name)   :: procname

! these are custom MPI types
integer                             :: realType  ! this will correspond to pDki, above
integer                             :: mpi_pType ! this will correspond to the derived ParticleType, below

! ==============================================================================

! derived type for particles
type ParticleType
    real(pDki) :: loc(2) ! real-valued spatial location
    real(pDki) :: mass ! real-valued mass
    integer    :: box ! which box the particle is in 
    logical    :: active ! whether a given particle is currently active
    logical    :: ghost(8) ! indicates whether a particle is a ghost particle in a corresponding direction.
    logical    :: jumped(8) ! indicates whether the particle jumped subdomains and which way it went.

                            ! The logicals for both ghost(8) and jumped(8) corresopnd to a direction based on the 
                            ! value of the 1 in the vector. The positions of the 1 correspond to the following send
                            ! directions: L, R, D, U, DL, DR, UL, UR
end type ParticleType

! sparse matrix type in coordinate format
type SparseMatType
    integer(kind=8)    :: row
    integer(kind=8)    :: col
    real(pDki)         :: val
end type SparseMatType

! holds the results of the kD tree fixed radius search
type kDRS_ResultType
    integer                 :: num    ! number of nearest neighbors (size of idx and rad)
    integer,    allocatable :: idx(:) ! indices of nearest neighbors
    real(pDki), allocatable :: rad(:) ! distances to nearest neighbors
end type kDRS_ResultType

! all the subroutines are below
contains


! subroutine to initialize the random number generator seed from clock time
! source: https://gcc.gnu.org/onlinedocs/gcc-4.2.1/gfortran/RANDOM_005fSEED.html
subroutine init_random_seed()
    integer              :: i, n, clock
    integer, allocatable :: seed(:)

    call random_seed(size = n)
    allocate (seed(n))
    call system_clock(count = clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
end subroutine init_random_seed


! function to generate a linearly-spaced array with n points ranging [low, high]
function linspace(low, high, n) result(vec)
    real(pDki), intent(in   ) :: low, high
    integer,    intent(in   ) :: n
    real(pDki)                :: vec(n)
    real(pDki)                :: stride
    integer                   :: i

    if (n < 2) then
        print *, 'ERROR: linspace requires n > 1'
        call exit
    else
        vec(1) = low
        vec(n) = high
        stride = (high - low) / real(n - 1, pDki)
        do i = 2, n - 1
            vec(i) = low + real(i - 1, pDki) * stride
        enddo
    endif
end function linspace


! create an MPI real type to correspond to the chosen real kind
subroutine create_MPIRealType(x, realType)
    real(pDki), intent(in   ) :: x
    integer,    intent(  out) :: realType

    ! ========================== LOCAL VARIABLES ===============================
    integer :: p, r ! precision and range
    integer :: ierror

    p = precision(x)
    r = range(x)

    CALL MPI_Type_create_f90_real(p, r, realType, ierror)
end subroutine create_MPIRealType


! create and commit the ParticleType as an MPI type
subroutine build_derived_pType(mesg_mpi_pType)
    integer, intent(  out) :: mesg_mpi_pType

    ! ========================== LOCAL VARIABLES ===============================
    integer, parameter        :: number = 6 ! number of fields in the ParticleType
    integer                   :: block_lengths(number)
    integer(mpi_address_kind) :: displacements(number)
    integer                   :: typelist(number)
    real(pDki)                :: r
    logical                   :: log
    integer                   :: errcode, ierr, disp1, disp2, disp3, disp4, disp5, disp6

    typelist = (/ realType, realType, mpi_integer, mpi_logical, mpi_logical, mpi_logical /)
    block_lengths = (/ 2, 1, 1, 1, 8, 8 /)

    disp1 = 0
    disp2 = int(2,mpi_address_kind)*sizeof(r)
    disp3 = disp2 + int(sizeof(r),mpi_address_kind)
    disp4 = disp3 + int(sizeof(ierr),mpi_address_kind)
    disp5 = disp4 + sizeof(log)
    disp6 = disp5 + int(8,mpi_address_kind)*sizeof(log)

    displacements = (/ disp1, disp2, disp3, disp4, disp5, disp6 /)

    call mpi_type_create_struct(number, block_lengths, displacements, &
                                typelist, mesg_mpi_pType, ierr)
    if (ierr /= 0 ) then
        print *, 'Error in type create: ', ierr
        call mpi_abort(mpi_comm_world, errcode, ierr)
    endif

    call mpi_type_commit(mesg_mpi_pType, ierr)
    if (ierr /= 0 ) then
        print *, 'Error in type commit: ', ierr
        call mpi_abort(mpi_comm_world, errcode, ierr)
    endif

end subroutine build_derived_pType


! moves particles via random walk diffusion
subroutine diffuse(np, D, dt, idxActive, p)
    integer,            intent(in   ) :: np ! number of particles
    real(pDki),         intent(in   ) :: D, dt ! diffusion coefficient and time step
    integer,            intent(in   ) :: idxActive(np) ! array containing the indices of the active particles
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! ========================== LOCAL VARIABLES ===============================
    real(pDki)                        :: normvec1(np), normvec2(np) ! vector which will hold Normal(0, 1) values

    ! call N(0, 1) generator
    ! call init_random_seed()

    call box_mullerp(np, normvec1)
    call box_mullerp(np, normvec2)

    p(idxActive)%loc(1) = p(idxActive)%loc(1) + sqrt(2.0_pDki * D * dt) * normvec1
    p(idxActive)%loc(2) = p(idxActive)%loc(2) + sqrt(2.0_pDki * D * dt) * normvec2
end subroutine diffuse


! this subroutine uses the Box-Muller transform to generate N(0,1)
! random numbers from U(0,1) numbers
! https://goo.gl/DQgmMu
! Note: this polar formulation seems to be consistently ~20% faster than the
! version that uses trig functions
! source for polar version (and standard version):
! https://www.taygeta.com/random/gaussian.html
subroutine box_mullerp(n, z)
    integer,    intent(in   ) :: n ! size of random vector to be generated
    real(pDki), intent(  out) :: z(n)

    ! ========================== LOCAL VARIABLES ===============================
    integer                   :: j
    real(pDki)                :: w, x1, x2
    real(pDki)                :: rand(2)

    ! initialize the random seed, just in case
    ! call init_random_seed()

    do j = 1, n/2
        w = 1.0_pDki
        do while (w >= 1.0_pDki)
            call random_number(rand)
            x1 = 2.0_pDki * rand(1) - 1.0_pDki
            x2 = 2.0_pDki * rand(2) - 1.0_pDki
            w = x1**2 + x2**2
        enddo
        w = sqrt((-2.0_pDki * log(w)) / w)
        z(2 * j - 1 : 2 * j) = (/x1 * w, x2 * w/)
    enddo

    if (mod(n, 2) /= 0) then
        w = 1.0_pDki
        do while (w >= 1.0_pDki)
            call random_number(rand)
            x1 = 2.0_pDki * rand(1) - 1.0_pDki
            x2 = 2.0_pDki * rand(2) - 1.0_pDki
            w = x1**2 + x2**2
        enddo
        w = sqrt((-2.0_pDki * log(w)) / w)
        z(n) = x1 * w
    endif
end subroutine box_mullerp


! reflective lower boundary
subroutine reflectLow(low, idxActive, p, case)
    real(pDki),         intent(in   ) :: low ! lower spatial boundary
    integer,            intent(in   ) :: idxActive(:)! indices of the active particles
    integer,            intent(in   ) :: case ! whether you are reflecting x or y values
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! if a particle has exited the boundary, bounce it back the distance
    ! by which it overstepped
    where (p(idxActive)%loc(case) < low) 
        p(idxActive)%loc(case) = 2.0_pDki * low - p(idxActive)%loc(case)
    endwhere
end subroutine reflectLow


! reflective upper boundary
subroutine reflectHigh(high, idxActive, p, case)
    real(pDki),         intent(in   ) :: high ! upper spatial boundary
    integer,            intent(in   ) :: idxActive(:)! indices of the active particles
    integer,            intent(in   ) :: case ! whether you are reflecting x or y values
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! if a particle has exited the boundary, bounce it back the distance
    ! by which it overstepped
    where (p(idxActive)%loc(case) > high) 
        p(idxActive)%loc(case) = 2.0_pDki * high - p(idxActive)%loc(case)
    endwhere
end subroutine reflectHigh


! this is the big communications step between cores/subdomains
! any particle that random-walked across subdomain boundaries is sent to the
! neighboring subdomain, and any particle that is within paddist of the boundary
! is also sent to the neighboring subdomain as a ghost particle
subroutine swapDomains(pVecSize, Nactive, idxActive, Dlims, paddist, p, nnx, nny, DDC_WhichWay, tstep,&
                        jumpedTotSum, ghostTotSum, totSendSum)
    integer,            intent(in   ) :: pVecSize ! total size of the local particle array
    integer,            intent(in   ) :: DDC_WhichWay ! which DDC method is being used
    integer,            intent(in   ) :: nny, nnx ! bring in the vert/horiz nodes
    integer,            intent(in   ) :: tstep ! current time step
    integer,            intent(inout) :: jumpedTotSum, ghostTotSum, totSendSum
    integer,            intent(inout) :: Nactive ! number of active particles
    integer,            intent(inout) :: idxActive(:) ! indices of the active particles
    real(pDki),         intent(in   ) :: Dlims(4) ! subdomain boundary limits
    real(pDki),         intent(in   ) :: paddist ! ghost particle padding distance
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! ========================== LOCAL VARIABLES ===============================
    integer, allocatable              :: idx(:) ! these are all used for bookkeeping in sending and receiving (hopefully what they are should be apparent by usage)
    integer                           :: active(Nactive)
    integer                           :: idxBig(pVecSize)
    integer                           :: numSendLeft, numSendRight, numSendUp, numSendDown, numSendDR, &
                                        numSendDL, numSendUR, numSendUL, recvSizeRight, recvSizeLeft, &
                                        recvSizeDown, recvSizeUp, recvSizeDL, recvSizeDR, recvSizeUL, recvSizeUR
    integer                           :: nJumpedRight, nJumpedLeft, nJumpedDown, nJumpedUp, nJumpedUL, &
                                        nJumpedUR, nJumpedDL, nJumpedDR, totRecv, nGhostRight, nGhostLeft, &
                                        nGhostDown, nGhostUp, nGhostDR, nGhostDL, nGhostUR, nGhostUL, totRecvLR, &
                                        totRecvLRU, totRecvCard, totRecvDiag1, totRecvDiag2, totRecvDiag3, &
                                        jumpedTot, ghostTot, totSend
    integer                           :: i, notActive
    real(pDki)                        :: ghostTotAv, jumpedTotAv, totSendAv
    type(ParticleType), allocatable   :: tmp_pLeft(:), tmp_pRight(:), tmp_pDown(:), tmp_pUp(:),&
                                        tmp_pDR(:), tmp_pDL(:), tmp_pUR(:), tmp_pUL(:) ! particle arrays for sending
    type(ParticleType), allocatable   :: recv_pLeft(:), recv_pRight(:), recv_pDown(:),&
                                        recv_pUp(:), recv_pUR(:), recv_pUL(:), recv_pDR(:),&
                                        recv_pDL(:) ! particle arrays for receiving

    ! initialize all of the relevant working info
    p%jumped(1) = .false.
    p%jumped(2) = .false.
    p%jumped(3) = .false.
    p%jumped(4) = .false.
    p%jumped(5) = .false.
    p%jumped(6) = .false.
    p%jumped(7) = .false.
    p%jumped(8) = .false.

    p%ghost(1) = .false.
    p%ghost(2) = .false.
    p%ghost(3) = .false.
    p%ghost(4) = .false.
    p%ghost(5) = .false.
    p%ghost(6) = .false.
    p%ghost(7) = .false.
    p%ghost(8) = .false.

    numSendLeft = 0
    numSendRight = 0
    numSendUp = 0
    numSendDown = 0
    numSendUL = 0
    numSendUR = 0
    numSendDL = 0
    numSendDR = 0
    recvSizeRight = 0
    recvSizeLeft = 0
    recvSizeDown = 0
    recvSizeUp = 0
    recvSizeUR = 0
    recvSizeUL = 0
    recvSizeDL = 0
    recvSizeDR = 0
    nJumpedRight = 0
    nJumpedLeft = 0
    nJumpedDown = 0
    nJumpedUp = 0
    nGhostRight = 0
    nGhostLeft = 0
    nGhostDown = 0
    nGhostUp = 0
    nGhostUL = 0
    nGhostUR = 0
    nGhostDL = 0
    nGhostDR = 0
    idxBig = (/(i, i = 1, pVecSize)/)

    ! this array will be used as an argument in the particle array
    allocate(idx(Nactive),recv_pLeft(5*Nactive))
    active = idxActive(1 : Nactive)

    ! NOTE: subdomain boundaries are [lower, upper)

    select case(DDC_WhichWay)


    !vertical slices
    case(1)

    ! first send particles to the left
    if (my_rank > 0) then

        allocate(tmp_pLeft(5*Nactive))
        ! if it is outside of the lower boundary, tag it to be sent to the left
        where (p(active)%loc(1) < Dlims(1))
            p(active)%jumped(1) = .true.
        endwhere

        nJumpedLeft = count(p(active)%jumped(1))
        idx(1 : nJumpedLeft) = pack(active, p(active)%jumped(1))
        ! put it in a temporary array for sending
        if (nJumpedLeft > 0) then
            tmp_pLeft(1 : nJumpedLeft) = p(idx(1 : nJumpedLeft))
        endif
        ! tag the particles that will be sent as ghost particles
        where (p(active)%loc(1) >= Dlims(1) .and. &
                p(active)%loc(1) <= (Dlims(1) + paddist))
            p(active)%ghost(1) = .true.
        endwhere

        ! where (p(active)%ghost(1) .eqv. .true. .and. p(active)%mass < eps)
        !     p(active)%ghost(1) = .false.
        ! endwhere

        idx = 0
        nGhostLeft = count(p(active)%ghost(1))
        idx(1 : nGhostLeft) = pack(active, p(active)%ghost(1))

        numSendLeft = nJumpedLeft + nGhostLeft
        ! add the ghost particles to the temporary particles array
        tmp_pLeft(nJumpedLeft + 1 : numSendLeft) = p(idx(1 : nGhostLeft))

        ! turn off the ghost particle indicator for the particles that are staying
        p%ghost(1) = .false.
        ! deactivate the particles that left the domain
        where (p(active)%jumped(1)) 
            p(active)%active = .false.
        endwhere

        ! send the number of particles to be sent to the left
        call mpi_send(numSendLeft, 1, mpi_integer, my_rank - 1, tag + my_rank - 1,&
                      mpi_comm_world, ierror)
    endif

    if (my_rank < num_cores - 1) then
        ! receive the number of particles to be received from the right
        call mpi_recv(recvSizeRight, 1, mpi_integer, my_rank + 1, tag + my_rank, &
                      mpi_comm_world, status, ierror)
    endif

    ! if there are > 0 particles to be sent, send them
    if (my_rank > 0 .and. numSendLeft > 0) then
        call mpi_send(tmp_pLeft(1 : numSendLeft), numSendLeft, mpi_pType, my_rank - 1,&
                      tag + my_rank - 1, mpi_comm_world, ierror)
        deallocate(tmp_pLeft)
    endif
    
    
    ! if there are > 0 particles to be received, receive them
    if (my_rank < num_cores - 1 .and. recvSizeRight > 0) then
        allocate(recv_pRight(5*Nactive))
        call mpi_recv(recv_pRight(1 : recvSizeRight), recvSizeRight, mpi_pType,&
                      my_rank + 1, tag + my_rank, mpi_comm_world, status, ierror)

        ! Write(*,*) "Receiving on the right on core",my_rank,"is",recv_pRight(1 : recvSizeRight)
    endif

    ! now send particles to the right (this works the same as above, and, thus, is
    ! not commented)
    if (my_rank < num_cores - 1) then
        allocate(tmp_pRight(5*Nactive))
        where (p(active)%loc(1) >= Dlims(2))
            p(active)%jumped(2) = .true.
        endwhere
        nJumpedRight = count(p(active)%jumped(2))
        idx(1 : nJumpedRight) = pack(active, p(active)%jumped(2))
        if (nJumpedRight > 0) then
            tmp_pRight(1 : nJumpedRight) = p(idx(1 : nJumpedRight))
        endif

        where (p(active)%loc(1) < Dlims(2) .and. &
                p(active)%loc(1) >= (Dlims(2) - paddist))
            p(active)%ghost(2) = .true.
        endwhere

        ! where (p(active)%ghost(2) .eqv. .true. .and. p(active)%mass < eps)
        !     p(active)%ghost(2) = .false.
        ! endwhere

        nGhostRight = count(p(active)%ghost(2))
        idx(1 : nGhostRight) = pack(active, p(active)%ghost(2))
        numSendRight = nJumpedRight + nGhostRight
        tmp_pRight(nJumpedRight + 1 : numSendRight) = p(idx(1 : nGhostRight))

        p%ghost(2) = .false.
        where (p(active)%jumped(2)) 
            p(active)%active = .false.
        endwhere

        call mpi_send(numSendRight, 1, mpi_integer, my_rank + 1, tag + my_rank + 1,&
                      mpi_comm_world, ierror)
    endif

    if (my_rank > 0) then
        call mpi_recv(recvSizeLeft, 1, mpi_integer, my_rank - 1, tag + my_rank, &
                      mpi_comm_world, status, ierror)
    endif

    if (my_rank < num_cores - 1 .and. numSendRight > 0) then
        call mpi_send(tmp_pRight(1 : numSendRight), numSendRight, mpi_pType, my_rank + 1,&
                      tag + my_rank + 1, mpi_comm_world, ierror)
        deallocate(tmp_pRight)
    endif
    if (my_rank > 0 .and. recvSizeLeft > 0) then
        call mpi_recv(recv_pLeft(1 : recvSizeLeft), recvSizeLeft, mpi_pType,&
                      my_rank - 1, tag + my_rank, mpi_comm_world, status, ierror)
    endif

    !checkerboard method
    case(2)
    ! first send particles to the left (up left and down left too)
    if (my_rank > (nny - 1)) then

        allocate(tmp_pLeft(Nactive),tmp_pDL(Nactive/2),tmp_pUL(Nactive/2))

        ! if it is outside of the lower boundary, tag it to be sent to the left
        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) > Dlims(3) .and. p(active)%loc(2) < Dlims(4))
            p(active)%jumped(1) = .true.
        endwhere

        idx = 0
        nJumpedLeft = count(p(active)%jumped(1))
        idx(1 : nJumpedLeft) = pack(active, p(active)%jumped(1))
        ! put it in a temporary array for sending
        tmp_pLeft(1 : nJumpedLeft) = p(idx(1 : nJumpedLeft))


        ! below left boundary and below bottom boundary
        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) < Dlims(3))
            p(active)%jumped(5) = .true.
        endwhere

        idx = 0
        nJumpedDL = count(p(active)%jumped(5))
        idx(1 : nJumpedDL) = pack(active, p(active)%jumped(5))
        ! put it in a temporary array for sending
        tmp_pDL(1 : nJumpedDL) = p(idx(1 : nJumpedDL))

        !below left boundary and above top boundary
        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) > Dlims(4))
            p(active)%jumped(7) = .true.
        endwhere

        ! determine which particles were just tagged for UL
        idx = 0
        nJumpedUL = count(p(active)%jumped(7))
        idx(1 : nJumpedUL) = pack(active, p(active)%jumped(7))
        ! put them in a temporary array for sending
        tmp_pUL(1 : nJumpedUL) = p(idx(1 : nJumpedUL))




        ! tag the particles that will be sent as ghost particles to the left
        where (p(active)%loc(1) >= Dlims(1) .and. p(active)%loc(1) < (Dlims(1) + paddist))

            p(active)%ghost(1) = .true.

        endwhere
        nGhostLeft = count(p(active)%ghost(1))
        idx(1 : nGhostLeft) = pack(active, p(active)%ghost(1))
        numSendLeft = nJumpedLeft + nGhostLeft
        ! add the ghost particles to the temporary particles array
        tmp_pLeft(nJumpedLeft + 1 : numSendLeft) = p(idx(1 : nGhostLeft))


        ! tagging DL ghosts
        where (p(idx(1 : nGhostLeft))%loc(2) <= Dlims(3) + paddist)

            p(idx(1 : nGhostLeft))%ghost(5) = .true.

        endwhere

        ! tagging UL ghosts
        where (p(idx(1 : nGhostLeft))%loc(2) >= Dlims(4) - paddist)

            p(idx(1 : nGhostLeft))%ghost(7) = .true.

        endwhere

        ! determine where DL ghosts are
        idx = 0
        nGhostDL = count(p(active)%ghost(5))
        idx(1 : nGhostDL) = pack(active, p(active)%ghost(5))
        numSendDL = nJumpedDL + nGhostDL
        ! add DL ghosts to temp particle array 
        tmp_pDL(nJumpedDL + 1 : numSendDL) = p(idx(1 : nGhostDL))

        ! determine where UL ghosts are
        idx = 0
        nGhostUL = count(p(active)%ghost(7))
        idx(1 : nGhostUL) = pack(active, p(active)%ghost(7))
        numSendUL = nJumpedUL + nGhostUL
        ! add UL ghosts to temp array 
        tmp_pUL(nJumpedUL + 1 : numSendUL) = p(idx(1 : nGhostUL))


        ! turn off the ghost particle indicator for the particles that are staying
        p%ghost(1) = .false.
        p%ghost(5) = .false.
        p%ghost(7) = .false.

        ! deactivate the particles that left the domain
        where (p(active)%jumped(1)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(5)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(7)) 
            p(active)%active = .false.
        endwhere

        ! send the number of particles to be sent to the left
        call mpi_send(numSendLeft, 1, mpi_integer, my_rank - nny, tag + my_rank - nny,&
                      mpi_comm_world, ierror)
    endif

    if (my_rank < num_cores - nny) then
        ! receive the number of particles to be received from the right
        call mpi_recv(recvSizeRight, 1, mpi_integer, my_rank + nny, tag + my_rank, &
                      mpi_comm_world, status, ierror)
    endif

    ! if there are > 0 particles to be sent, send them
    if (my_rank > nny - 1 .and. numSendLeft > 0) then
        call mpi_send(tmp_pLeft(1 : numSendLeft), numSendLeft, mpi_pType, my_rank - nny,&
                      tag + my_rank - nny, mpi_comm_world, ierror)
        deallocate(tmp_pLeft)
    endif


    ! if there are > 0 particles to be received, receive them
    if (my_rank < num_cores - nny .and. recvSizeRight > 0) then
        allocate(recv_pRight(Nactive))
        call mpi_recv(recv_pRight(1 : recvSizeRight), recvSizeRight, mpi_pType,&
                      my_rank + nny, tag + my_rank, mpi_comm_world, status, ierror)
    endif


    ! now send particles to the right (this works the same as above, and thus is
    ! not commented)
    if (my_rank < num_cores - nny) then

        allocate(tmp_pRight(Nactive),tmp_pUR(Nactive/2),tmp_pDR(Nactive/2))

        where (p(active)%loc(1) > Dlims(2) .and. p(active)%loc(2) < Dlims(4) .and. p(active)%loc(2) > Dlims(3))
            p(active)%jumped(2) = .true.
        endwhere
        nJumpedRight = count(p(active)%jumped(2))
        idx(1 : nJumpedRight) = pack(active, p(active)%jumped(2))
        tmp_pRight(1 : nJumpedRight) = p(idx(1 : nJumpedRight))


        ! above right boundary and below bottom boundary
        where (p(active)%loc(1) > Dlims(2) .and. p(active)%loc(2) < Dlims(3))
            p(active)%jumped(6) = .true.
        endwhere

        idx = 0
        nJumpedDR = count(p(active)%jumped(6))
        idx(1 : nJumpedDR) = pack(active, p(active)%jumped(6))
        ! put it in a temporary array for sending
        tmp_pDR(1 : nJumpedDR) = p(idx(1 : nJumpedDR))

        !below lower boundary and above top boundary
        where (p(active)%loc(1) > Dlims(2) .and. p(active)%loc(2) > Dlims(4))
            p(active)%jumped(8) = .true.
        endwhere

        idx = 0
        nJumpedUR = count(p(active)%jumped(8))
        idx(1 : nJumpedUR) = pack(active, p(active)%jumped(8))
        ! put it in a temporary array for sending
        tmp_pUR(1 : nJumpedUR) = p(idx(1 : nJumpedUR))


        where (p(active)%loc(1) < Dlims(2) .and. p(active)%loc(1) >= (Dlims(2) - paddist))

            p(active)%ghost(2) = .true.
        
        endwhere
        nGhostRight = count(p(active)%ghost(2))
        idx(1 : nGhostRight) = pack(active, p(active)%ghost(2))
        numSendRight = nJumpedRight + nGhostRight
        tmp_pRight(nJumpedRight + 1 : numSendRight) = p(idx(1 : nGhostRight))

        where (p(idx(1 : nGhostRight))%loc(2) <= Dlims(3) + paddist)

            p(idx(1 : nGhostRight))%ghost(6) = .true.

        endwhere

        where (p(idx(1 : nGhostRight))%loc(2) >= Dlims(4) - paddist)

            p(idx(1 : nGhostRight))%ghost(8) = .true.

        endwhere

        idx = 0
        nGhostDR = count(p(active)%ghost(6))
        idx(1 : nGhostDR) = pack(active, p(active)%ghost(6))
        numSendDR = nJumpedDR + nGhostDR
        tmp_pDR(nJumpedDR + 1 : numSendDR) = p(idx(1 : nGhostDR))

        idx = 0
        nGhostUR = count(p(active)%ghost(8))
        idx(1 : nGhostUR) = pack(active, p(active)%ghost(8))
        numSendUR = nJumpedUR + nGhostUR
        tmp_pUR(nJumpedUR + 1 : numSendUR) = p(idx(1 : nGhostUR))

        p%ghost(2) = .false.
        p%ghost(6) = .false.
        p%ghost(8) = .false.
        where (p(active)%jumped(2)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(6)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(8)) 
            p(active)%active = .false.
        endwhere

        call mpi_send(numSendRight, 1, mpi_integer, my_rank + nny, tag + my_rank + nny,&
                      mpi_comm_world, ierror)
    endif

    if (my_rank > nny - 1) then
        call mpi_recv(recvSizeLeft, 1, mpi_integer, my_rank - nny, tag + my_rank, &
                      mpi_comm_world, status, ierror)
    endif

    if (my_rank < num_cores - nny .and. numSendRight > 0) then
        call mpi_send(tmp_pRight(1 : numSendRight), numSendRight, mpi_pType, my_rank + nny,&
                      tag + my_rank + nny, mpi_comm_world, ierror)
        deallocate(tmp_pRight)
    endif

    if (my_rank > nny - 1 .and. recvSizeLeft > 0) then
        call mpi_recv(recv_pLeft(1 : recvSizeLeft), recvSizeLeft, mpi_pType,&
                      my_rank - nny, tag + my_rank, mpi_comm_world, status, ierror)
    endif


    ! ! now send particles to up/down

        ! ! first we send the particles down

        if (mod(my_rank,nny) > 0) then

            allocate(tmp_pDown(Nactive))

            where (p(active)%loc(2) < Dlims(3) .and. p(active)%loc(1) > Dlims(1) .and. p(active)%loc(1) < Dlims(2))
                p(active)%jumped(3) = .true.
            endwhere

            idx = 0
            nJumpedDown = count(p(active)%jumped(3))
            idx(1 : nJumpedDown) = pack(active, p(active)%jumped(3))
            tmp_pDown(1 : nJumpedDown) = p(idx(1 : nJumpedDown))

            where (p(active)%loc(2) >= Dlims(3) .and. p(active)%loc(2) <= (Dlims(3) + paddist))

                p(active)%ghost(3) = .true.

            endwhere

            idx = 0
            nGhostDown = count(p(active)%ghost(3))
            idx(1 : nGhostDown) = pack(active, p(active)%ghost(3))
            numSendDown = nJumpedDown + nGhostDown
            tmp_pDown(nJumpedDown + 1 : numSendDown) = p(idx(1 : nGhostDown))

            p%ghost(3) = .false.
            where (p(active)%jumped(3)) 
                p(active)%active = .false.
            endwhere

            call mpi_send(numSendDown, 1, mpi_integer, my_rank - 1, tag + my_rank - 1,&
                          mpi_comm_world, ierror)

            ! send diagonal particles if necessary

            ! down and to the right
            if (my_rank < num_cores - nny) then
                call mpi_send(numSendDR, 1, mpi_integer, my_rank + nny - 1, tag + my_rank + nny - 1,&
                          mpi_comm_world, ierror)
            endif

            ! down and to the left
            if (my_rank > (nny - 1)) then
                call mpi_send(numSendDL, 1, mpi_integer, my_rank - nny - 1, tag + my_rank - nny - 1,&
                          mpi_comm_world, ierror)
            endif


        endif

        ! recv number from above
        if (mod(my_rank,nny) /= (nny - 1)) then
            call mpi_recv(recvSizeUp, 1, mpi_integer, my_rank + 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif


        !recv number from UR
        if (mod(my_rank,nny) /= (nny - 1) .and. my_rank < (num_cores - nny)) then
            call mpi_recv(recvSizeUR, 1, mpi_integer, my_rank + nny + 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif

        !recv number from UL
        if (mod(my_rank,nny) /= (nny - 1) .and. my_rank > (nny - 1)) then
            call mpi_recv(recvSizeUL, 1, mpi_integer, my_rank - nny + 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif

        !send array down
        if (mod(my_rank,nny) > 0 .and. numSendDown > 0) then
            call mpi_send(tmp_pDown(1 : numSendDown), numSendDown, mpi_pType, my_rank - 1,&
                          tag + my_rank - 1, mpi_comm_world, ierror)
            deallocate(tmp_pDown)
        endif

        !recv array from above
        if (mod(my_rank,nny) /= (nny - 1) .and. recvSizeUp > 0) then
            allocate(recv_pUp(Nactive))
            call mpi_recv(recv_pUp(1 : recvSizeUp), recvSizeUp, mpi_pType,&
                          my_rank + 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        !send array DR
        if (mod(my_rank,nny) > 0 .and. my_rank < (num_cores - nny) .and. numSendDR > 0) then
            call mpi_send(tmp_pDR(1 : numSendDR), numSendDR, mpi_pType, my_rank + nny - 1, tag + my_rank + nny - 1,&
                          mpi_comm_world, ierror)
            deallocate(tmp_pDR)
        endif 

        !recv array from UL
        if (mod(my_rank,nny) /= (nny-1) .and. my_rank > (nny - 1) .and. recvSizeUL > 0) then
            allocate(recv_pUL(Nactive/2))
            call mpi_recv(recv_pUL(1 : recvSizeUL), recvSizeUL, mpi_pType,&
                          my_rank - nny + 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        !send array DL
        if (mod(my_rank,nny) > 0 .and. my_rank > (nny - 1) .and. numSendDL > 0) then
            call mpi_send(tmp_pDL(1 : numSendDL), numSendDL, mpi_pType, my_rank - nny - 1, tag + my_rank - nny - 1,&
                          mpi_comm_world, ierror)
            deallocate(tmp_pDL)
        endif

        !recv array UR
        if (mod(my_rank,nny) /= (nny-1) .and. my_rank < (num_cores - nny) .and. recvSizeUR > 0) then
            allocate(recv_pUR(Nactive/2))
            call mpi_recv(recv_pUR(1 : recvSizeUR), recvSizeUR, mpi_pType,&
                          my_rank + nny + 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        ! ! now send particles up


        if (mod(my_rank,nny) /= (nny - 1)) then

            allocate(tmp_pUp(Nactive))

            where (p(active)%loc(2) > Dlims(4) .and. p(active)%loc(1) > Dlims(1) .and. p(active)%loc(1) < Dlims(2))
                p(active)%jumped(4) = .true.
            endwhere
            idx = 0
            nJumpedUp = count(p(active)%jumped(4))
            idx(1 : nJumpedUp) = pack(active, p(active)%jumped(4))
            tmp_pUp(1 : nJumpedUp) = p(idx(1 : nJumpedUp))

            where (p(active)%loc(2) < Dlims(4) .and. p(active)%loc(2) >= (Dlims(4) - paddist))

                p(active)%ghost(4) = .true.

            endwhere
            idx = 0
            nGhostUp = count(p(active)%ghost(4))
            idx(1 : nGhostUp) = pack(active, p(active)%ghost(4))
            numSendUp = nJumpedUp + nGhostUp
            tmp_pUp(nJumpedUp + 1 : numSendUp) = p(idx(1 : nGhostUp))

            p%ghost(4) = .false.
            where (p(active)%jumped(4)) 
                p(active)%active = .false.
            endwhere

            call mpi_send(numSendUp, 1, mpi_integer, my_rank + 1, tag + my_rank + 1,&
                          mpi_comm_world, ierror)

            !send size UR
            if (my_rank < num_cores - nny) then
                call mpi_send(numSendUR, 1, mpi_integer, my_rank + nny + 1, tag + my_rank + nny + 1,&
                          mpi_comm_world, ierror)
            endif

            !send size UL
            if (my_rank > (nny - 1)) then
                call mpi_send(numSendUL, 1, mpi_integer, my_rank - nny + 1, tag + my_rank - nny + 1,&
                          mpi_comm_world, ierror)
            endif
        endif

        !recv number from below
        if (mod(my_rank,nny) /= 0) then
            call mpi_recv(recvSizeDown, 1, mpi_integer, my_rank - 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif

        !send array up
        if (mod(my_rank,nny) /= nny-1 .and. numSendUp > 0) then
            call mpi_send(tmp_pUp(1 : numSendUp), numSendUp, mpi_pType, my_rank + 1,&
                          tag + my_rank + 1, mpi_comm_world, ierror)
            deallocate(tmp_pUp)
        endif    

        !recv array from below
        if (mod(my_rank,nny) /= 0 .and. recvSizeDown > 0) then
            allocate(recv_pDown(Nactive))
            call mpi_recv(recv_pDown(1 : recvSizeDown), recvSizeDown, mpi_pType,&
                          my_rank - 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        !recv size DL
        if (mod(my_rank,nny) /= 0 .and. my_rank > nny - 1) then
            call mpi_recv(recvSizeDL, 1, mpi_integer, my_rank - nny - 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif

        !recv size DR
        if (mod(my_rank,nny) /= 0 .and. my_rank < (num_cores - nny)) then
            call mpi_recv(recvSizeDR, 1, mpi_integer, my_rank + nny - 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif

        !send array UR
        if (mod(my_rank,nny) /= (nny-1) .and. my_rank < (num_cores - nny) .and. numSendUR > 0) then
            call mpi_send(tmp_pUR(1 : numSendUR), numSendUR, mpi_pType, my_rank + nny + 1,&
                          tag + my_rank + nny + 1, mpi_comm_world, ierror)
            deallocate(tmp_pUR)
        endif    

        !recv array DL
        if (mod(my_rank,nny) /= 0 .and. my_rank > (nny - 1) .and. recvSizeDL > 0) then
            allocate(recv_pDL(Nactive/2))
            call mpi_recv(recv_pDL(1 : recvSizeDL), recvSizeDL, mpi_pType,&
                          my_rank - nny - 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        !send array UL
        if (mod(my_rank,nny) /= (nny-1) .and. my_rank > (nny - 1) .and. numSendUL > 0) then
            call mpi_send(tmp_pUL(1 : numSendUL), numSendUL, mpi_pType, my_rank - nny + 1,&
                          tag + my_rank - nny + 1, mpi_comm_world, ierror)
            deallocate(tmp_pUL)
        endif

        !recv array DR
        if (mod(my_rank,nny) /= 0 .and. my_rank < (num_cores - nny) .and. recvSizeDR > 0) then
            allocate(recv_pDR(Nactive/2))
            call mpi_recv(recv_pDR(1 : recvSizeDR), recvSizeDR, mpi_pType,&
                          my_rank + nny - 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif




    end select


    deallocate(idx)
    allocate(idx(pVecSize))

    ! idx will show where to fill the newly received particles into the non-active spaces in "p"
    idx = 0
    idx = pack(idxBig, .not. p%active)

    select case(DDC_WhichWay)

    !vertical slices
    case(1)

        totRecv = recvSizeLeft + recvSizeRight

    !checkerboard
    case(2)
        totRecvLR = recvSizeLeft + recvSizeRight
        totRecvLRU = totRecvLR + recvSizeUp
        totRecvCard = totRecvLRU + recvSizeDown
        totRecvDiag1 = totRecvCard + recvSizeDL
        totRecvDiag2 = totRecvDiag1 + recvSizeDR
        totRecvDiag3 = totRecvDiag2 + recvSizeUL
        totRecv = totRecvDiag3 + recvSizeUR


        ! can comment some of this off and on to get relevant send information for MPI performance

        ! if (my_rank == (floor(nnx/real(2,pDki))*nny + floor(nny/real(2,pDki)))) then
        !     jumpedTot = nJumpedLeft+nJumpedRight+nJumpedUp+nJumpedDown+nJumpedDL+&
        !             nJumpedDR+nJumpedUL+nJumpedUR

        !     ghostTot = nGhostLeft+nGhostRight+nGhostDown+nGhostUp+nGhostDL+&
        !             nGhostDR+nGhostUL+nGhostUR

        !     totSend = ghostTot+jumpedTot

        !     jumpedTotSum = jumpedTotSum + jumpedTot

        !     jumpedTotAv = real(jumpedTotSum,pDki)/real(tstep,pDki)

        !     ghostTotSum = ghostTotSum + ghostTot

        !     ghostTotAv = real(ghostTotSum,pDki)/real(tstep,pDki)

        !     totSendSum = totSendSum + totSend

        !     totSendAv = real(totSendSum,pDki)/real(tstep,pDki)

        !     ! Write(*,*) "jumpedTot is",jumpedTot
        !     ! Write(*,*) "ghostTot is",ghostTot
        !     ! Write(*,*) "totSend is",totSend
        !     ! Write(*,*) "totRecv is",totRecv

        !     ! if(tstep == 100) then

        !     !     Write(*,*) "Average particles jumped is",jumpedTotAv
        !     !     Write(*,*) "Average ghost particles is",ghostTotAv
        !     !     ! Write(*,*) "Average particles sent is",totSendAv

        !     ! endif

        ! endif

    end select


    ! if > 0 particles were received, put them in the main particles array,
    ! first filling them into the empty spaces of the non-active particles

    select case(DDC_WhichWay)

    case(1)
        if (totRecv > 0) then
            recv_pLeft(recvSizeLeft + 1 : totRecv) = recv_pRight(1 : recvSizeRight)
            do i = 1, totRecv
                p(idx(i)) = recv_pLeft(i)
                p(idx(i))%active = .true.
            enddo            
        endif

        if (my_rank < num_cores - 1 .and. recvSizeRight > 0) then
            deallocate(recv_pRight)
        endif

        deallocate(recv_pLeft)

    case(2)

        ! put all received particles into current particle array for the next MT step
        if (totRecv > 0) then
            recv_pLeft(recvSizeLeft + 1 : totRecvLR) = recv_pRight(1 : recvSizeRight)
            if (nny > 1) then
                recv_pLeft(totRecvLR + 1 : totRecvLRU) = recv_pUp(1 : recvSizeUp)
                recv_pLeft(totRecvLRU + 1 : totRecvCard) = recv_pDown(1 : recvSizeDown)
                recv_pLeft(totRecvCard + 1 : totRecvDiag1) = recv_pDL(1 : recvSizeDL)
                recv_pLeft(totRecvDiag1 + 1 : totRecvDiag2) = recv_pDR(1 : recvSizeDR)
                recv_pLeft(totRecvDiag2 + 1 : totRecvDiag3) = recv_pUL(1 : recvSizeUL)
                recv_pLeft(totRecvDiag3 + 1 : totRecv) = recv_pUR(1 : recvSizeUR)
            endif


            ! dynamically deallocate arrays if they were allocated above
            if (my_rank < num_cores - nny .and. recvSizeRight > 0) then
                    deallocate(recv_pRight)
            endif
            if (mod(my_rank,nny) /= (nny - 1) .and. recvSizeUp > 0) then
                    deallocate(recv_pUp)
            endif
            if (mod(my_rank,nny) /= (nny-1) .and. my_rank > (nny - 1) &
                .and. recvSizeUL > 0) then
                    deallocate(recv_pUL)
            endif
            if (mod(my_rank,nny) /= (nny-1) .and. my_rank < (num_cores - nny) &
                .and. recvSizeUR > 0) then
                    deallocate(recv_pUR)
            endif
            if (mod(my_rank,nny) /= 0 .and. recvSizeDown > 0) then
                    deallocate(recv_pDown)
            endif
            if (mod(my_rank,nny) /= 0 .and. my_rank > (nny - 1) &
                .and. recvSizeDL > 0) then
                    deallocate(recv_pDL)
            endif
            if (mod(my_rank,nny) /= 0 .and. my_rank < (num_cores - nny) &
                .and. recvSizeDR > 0) then
                    deallocate(recv_pDR)
            endif

            ! activate all received particles for MT step 
            ! All received particles will be tagged as a ghost, so they will be easily deactivated after MT step. 
            do i = 1, totRecv
                p(idx(i)) = recv_pLeft(i)
                p(idx(i))%active = .true.
            enddo
            deallocate(recv_pLeft)
        endif

    end select

    deallocate(idx)
end subroutine swapDomains

! write plotting information to file
subroutine write_plot(uname, filename, vec, initial, N, Nsteps)
    integer,           intent(in   ) :: uname ! unit number to write to
    character(*),      intent(in   ) :: filename ! filename to write to
    real(pDki),        intent(in   ) :: vec(:) ! the vector to be written (either mass or location)
    logical,           intent(in   ) :: initial ! true if writing the initial time step
    integer, optional, intent(in   ) :: N, Nsteps ! number of particles and number of time steps, for the header

    ! if it is the initial time step, write the header indicating the shape of the the array
    if (initial) then
        open (unit = uname, file = filename, action = 'write', status='replace')
        write (uname, *) vec
        close (unit = uname, status='keep')
    else
        open (unit = uname, file = filename, action = 'write', status='old', access='append')
        write (uname, *) vec
        close (unit = uname, status='keep')
    endif

end subroutine write_plot

! write run time to file for the current loop realization
subroutine write_time(uname, filename, time, initial, Num_ens, Np_ens, Np_list)
    integer,           intent(in   ) :: uname ! unit number to write to
    character(*),      intent(in   ) :: filename ! filename to write to
    real(pDki),        intent(in   ) :: time ! run time for the current realization
    logical,           intent(in   ) :: initial ! true if writing the for the first time => write the header
    integer, optional, intent(in   ) :: Num_ens, Np_ens ! number of runs in the ensemble (for averaging) and number of members in the particle number ensemble
    integer, optional, intent(in   ) :: Np_list(:) ! the number of particles in each member of the particle number ensemble

    ! if it is the initial time step, write the header indicating the shape of the the array
    if (initial) then
        open (unit = uname, file = filename, action = 'write',status = 'replace')
        write (uname, *) time
        close (unit = uname, status='keep')
    else
         open (unit = uname, file = filename, action = 'write', status='old', access='append')
         write (uname, *) time
         close (unit = uname, status='keep')
    endif
end subroutine write_time

! write error to file for the current loop realization
subroutine write_error(uname, filename, error, initial, Num_ens, Np_ens, Np_list)
    integer,           intent(in   ) :: uname ! unit number to write to
    character(*),      intent(in   ) :: filename ! filename to write to
    real(pDki),        intent(in   ) :: error ! error for the current realization
    logical,           intent(in   ) :: initial ! true if writing the for the first time => write the header
    integer, optional, intent(in   ) :: Num_ens, Np_ens ! number of runs in the ensemble (for averaging) and number of members in the particle number ensemble
    integer, optional, intent(in   ) :: Np_list(:) ! the number of particles in each member of the particle number ensemble

    ! if it is the initial time step, write the header indicating the shape of the the array
    if (initial) then
        open (unit = uname, file = filename, action = 'write', status ='replace')
        write(uname, *) error
        ! write (uname, *) Num_ens, Np_ens
        ! write (uname, *) Np_list
        close (unit = uname, status='keep')
    else
        open (unit = uname, file = filename, action = 'write', status='old', access='append')
        write (uname, *) error
        close (unit = uname, status='keep')
    endif
end subroutine write_error

! solves y = A * x, when A is in sparse coordinate format
! source: https://www.it.uu.se/education/phd_studies/phd_courses/pasc/lecture-1
subroutine SP_matVecMult(A, x, n, y)
    type(SparseMatType), intent(in   ) :: A(:)
    real(pDki),          intent(in   ) :: x(:)
    integer(kind=8),     intent(in   ) :: n ! number of entries in A
    real(pDki),          intent(  out) :: y(:)
    integer                            :: i

    y = 0.0_pDki

    do i = 1, n
        y(A(i)%row) = y(A(i)%row) + A(i)%val * x(A(i)%col)
    enddo
end subroutine SP_matVecMult

! this is the sparse matrix-based mass transfer algorithm (SPaM)
subroutine massTrans_kDMat(n, idxActive, cutdist, denom, p, beta)
    integer,            intent(in   ) :: n ! number of particles
    integer,            intent(in   ) :: idxActive(:) ! indices of the active particles
    real(pDki),         intent(in   ) :: cutdist ! cutoff distance for the kD tree fixed-radius search
    real(pDki),         intent(in   ) :: denom ! denominator of the exponential in the co-location probability density
    type(ParticleType), intent(inout) :: p(:) ! particle array
    real(pDki),         intent(in   ) :: beta ! beta parameter, encoding the
   ! bandwidth of the mass-transfer kernel
    ! ========================== LOCAL VARIABLES ===============================
    type(SparseMatType), allocatable :: Emat(:) ! mass transfer matrix
    integer(kind=8)                  :: start(n), finish(n)
    integer(kind=8)                  :: Nclose ! used for building the mass transfer matrix
    real(pDki)                       :: tmpmass(n) ! temporary array for holding particle masses
    type(ParticleType)               :: tmp_p(n) ! temporary particle array for dealing with ghost particles
    integer                          :: idx(n) ! indexing array
    integer                          :: i
    integer                          :: nNotGhost, idxNotGhost(n), idxNotGhostTmp(n) ! indexing array for non-ghost particles
    logical                          :: logNotGhost(n) ! logical array for non-ghost particles

    tmp_p = p(idxActive(1 : n))

    ! build the pairwise distance matrix
    call build_DistmatSparse(n, cutdist, tmp_p, Emat, start, finish, Nclose)

    ! build the matrix of co-location probabilities
    call build_PmatSparse(n, denom, start, finish, Emat)

    ! build the mass transfer matrix
    call build_EmatSparse(n, Nclose, Emat, beta)

    tmpmass = tmp_p%mass

    ! conduct the mass transfers with sparse matrix-vector multiplication
    call SP_matVecMult(Emat, tmpmass, Nclose, tmp_p%mass)

    ! only change the masses of non-ghost particles
    idx = (/(i, i = 1, n)/)
    logNotGhost = tmp_p%ghost(1) .or. tmp_p%ghost(2) .or. tmp_p%ghost(3) .or. tmp_p%ghost(4) .or.&
                    tmp_p%ghost(5) .or. tmp_p%ghost(6) .or. tmp_p%ghost(7) .or. tmp_p%ghost(8)
    logNotGhost = .not. logNotGhost
    nNotGhost = count(logNotGhost)
    idxNotGhostTmp(1 : nNotGhost) = pack(idx, logNotGhost)
    idxNotGhost(1 : nNotGhost) = pack(idxActive(1 : n), logNotGhost)

    p(idxNotGhost(1 : nNotGhost))%mass = tmp_p(idxNotGhostTmp(1 : nNotGhost))%mass

    deallocate(Emat)
end subroutine massTrans_kDMat


! this is the looping sparse mass transfer algorithm (PLSS)
! portions of this work the same as in massTrans_kDMat(), and thus are not commented
subroutine massTrans_kDLoop(n, idxActive, cutdist, factor, denom, p)
    integer,            intent(in   ) :: n ! number of particles
    integer,            intent(in   ) :: idxActive(:) ! indices of the active particles
    real(pDki),         intent(in   ) :: cutdist ! cutoff distance for the kD tree fixed-radius search
    real(pDki),         intent(in   ) :: factor ! constant multiplying the co-location probability density
    real(pDki),         intent(in   ) :: denom ! denominator of the exponential in the co-location probability density
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! ========================== LOCAL VARIABLES ===============================
    type(kDRS_ResultType), allocatable :: neighbors(:) ! holds the results of the kD tree fixed radius search
    integer                            :: i, j ! loop iterators
    integer                            :: length
    integer, allocatable               :: idx(:) ! temporary array for holding indices of nearby particles (according to kD search)
    real(pDki), allocatable            :: rad(:) ! temporary array for holding distances to nearby particles (according to kD search)
    real(pDki), allocatable            :: Ptot(:) ! used to normalize the probabilities to sum to 1
    real(pDki)                         :: rescale ! used to normalize the probabilities to sum to 1
    real(pDki)                         :: vs_ds ! normalized co-location probability for the i-j particle pair
    real(pDki)                         :: dm ! mass differential for the i-j particle pair
    integer                            :: numTrans ! these next few are used for indexing arguments
    integer                            :: jpart
    logical, allocatable               :: logTrans(:)
    integer, allocatable               :: idxBig(:), idxTrans(:) ! used as indexing arguments to the particle array
    type(ParticleType)                 :: tmp_p(n) ! see above for what everything below this point is
    integer                            :: idxtmp(n)
    integer                            :: nNotGhost, idxNotGhost(n), idxNotGhostTmp(n)
    logical                            :: logNotGhost(n)

    tmp_p = p(idxActive(1 : n))

    ! conduct the kD tree fixed-radius search
    call fixedRadSearch(n, cutdist, tmp_p, neighbors)

    ! the following vectors are made as large as they possibly need to be to avoid
    ! having to allocate/deallocate repeatedly
    length = maxval(neighbors%num)
    allocate(idx(length), rad(length), Ptot(length), logTrans(length),&
             idxBig(length), idxTrans(length))

    idxBig = (/(i, i = 1, length)/)
    idxTrans = 0

    ! nested loop for conducting the mass transfers
    do i = 1, n
        ! we are only doing mass transfers for particles with an index greater
        ! than i, in order to avoid doing them twice
        logTrans(1 : neighbors(i)%num) = neighbors(i)%idx > i
        numTrans = count(logTrans(1 : neighbors(i)%num))
        if (numTrans > 0) then
            idxTrans(1 : numTrans) = pack(idxBig(1 : neighbors(i)%num), logTrans(1 : neighbors(i)%num))

            ! get the relevant indices and radii
            idx(1 : numTrans) = neighbors(i)%idx(idxTrans(1 : numTrans))
            rad(1 : numTrans) = neighbors(i)%rad(idxTrans(1 : numTrans))

            ! normalize the probabilities to sum to 1
            Ptot(1 : neighbors(i)%num) = factor * exp(neighbors(i)%rad / denom)
            rescale = sum(Ptot(1 : neighbors(i)%num))
            Ptot(1 : numTrans) = Ptot(idxTrans(1 : numTrans))

            ! make the indexing array
            idxTrans(1 : numTrans) = idxBig(1 : numTrans)

            ! do the transfers in a randomized order to avoid numerical artifacts
            call shuffleInt(idxTrans(1 : numTrans))

            ! do the pairwise mass transfers
            do j = 1, numTrans
                jpart = idx(idxTrans(j))
                vs_ds = Ptot(idxTrans(j)) / rescale
                dm = 0.5_pDki * (tmp_p(i)%mass - tmp_p(jpart)%mass) * vs_ds
                ! don't allow negative masses
                tmp_p(i)%mass     = maxval((/     tmp_p(i)%mass - dm, 0.0_pDki /))
                tmp_p(jpart)%mass = maxval((/ tmp_p(jpart)%mass + dm, 0.0_pDki /))
            enddo
            idxTrans = 0
        endif
    enddo

    ! only change the masses of non-ghost particles
    idxtmp = (/(i, i = 1, n)/)
    logNotGhost = tmp_p%ghost(1) .or. tmp_p%ghost(2) .or. tmp_p%ghost(3) .or. tmp_p%ghost(4) .or.&
        tmp_p%ghost(5) .or. tmp_p%ghost(6) .or. tmp_p%ghost(7) .or. tmp_p%ghost(8)
    logNotGhost = .not. logNotGhost
    nNotGhost = count(logNotGhost)
    idxNotGhostTmp(1 : nNotGhost) = pack(idxtmp(1 : n), logNotGhost)
    idxNotGhost(1 : nNotGhost) = pack(idxActive(1 : n), logNotGhost)

    p(idxNotGhost(1 : nNotGhost))%mass = tmp_p(idxNotGhostTmp(1 : nNotGhost))%mass

    deallocate(idx, rad, Ptot, logTrans, idxTrans)
end subroutine massTrans_kDLoop

! this subroutine builds the (squared) distance matrix
! NOTE: hard-coded for 1D
subroutine build_DistmatSparse(n, cutdist, p, Distmat, start, finish, Nclose)
    integer,                          intent(in   ) :: n ! number of particles
    real(pDki),                       intent(in   ) :: cutdist ! cutoff distance for the kD tree fixed-radius search
    type(ParticleType),               intent(in   ) :: p(:) ! particle array
    type(SparseMatType), allocatable, intent(  out) :: Distmat(:) ! sparse distance matrix
    integer(kind=8),                  intent(  out) :: start(n), finish(n) ! indices (in the Distmat vectors) for the start and finish of each column of Distmat
    integer(kind=8),                  intent(  out) :: Nclose ! total number of neighbor particles found by the kD search (i.e, the length of the vectors in Distmat)

    ! ========================== LOCAL VARIABLES ===============================
    type(kDRS_ResultType), allocatable :: neighbors(:) ! holds the results of the kD tree fixed radius search
    integer                            :: i ! loop iterator
    integer(kind=8)                    :: tmpstart ! temporary variable to handle the n+1^th calculation of start in the loop

    ! conduct the kD tree fixed-radius search
    ! NOTE: this returns squared distances between particles

    call fixedRadSearch(n, cutdist, p, neighbors)

    ! allocate Distmat to have length = total number of neighbors found by the kD search
    Nclose = sum(int(neighbors%num,8))
    
    allocate(Distmat(Nclose))

    ! fill in Distmat
    tmpstart = 1
    do i = 1, n
        start(i) = tmpstart
        finish(i) = start(i) - 1 + int(neighbors(i)%num,8)
        Distmat(start(i) : finish(i))%col = i
        Distmat(start(i) : finish(i))%row = neighbors(i)%idx
        Distmat(start(i) : finish(i))%val = real(neighbors(i)%rad, pDki)
        tmpstart = finish(i) + 1
    enddo
    
    deallocate (neighbors)
end subroutine build_DistmatSparse

! perform the kD tree fixed-radius search
subroutine fixedRadSearch(n, cutdist, p, neighbors)
    integer,                            intent(in   ) :: n ! number of particles
    real(pDki),                         intent(in   ) :: cutdist ! cutoff distance for the kD tree fixed-radius search
    type(ParticleType),                 intent(in   ) :: p(:) ! particle array
    type(kDRS_ResultType), allocatable, intent(  out) :: neighbors(:) ! results of the fixed-radius search

    ! ========================== LOCAL VARIABLES ===============================
    type(kdtree2), pointer          :: tree ! this is the KD tree
    real(kdkind)                    :: locs(dim,n), locs1(n), locs2(n) ! locations array to be passed to kD tree
    real(kdkind)                    :: r2 ! value of squared search radius for kD tree

    ! convert particle locations to kdkind
    ! locs = real(p%loc(1), kdkind)
    locs1 = real(p%loc(1), kdkind)
    locs2 = real(p%loc(2), kdkind)
    locs(1,:) = locs1
    locs(2,:) = locs2

    ! squared search cutoff distance in kdkind
    r2   = real(cutdist, kdkind)**2

    ! build the KD tree and search it
    call maketree(tree, 2, n, locs)

    allocate (neighbors(n))

    ! this finds the closest mobile particles to each immobile particle
    ! NOTE: this search returns the SQUARED distance between two particles
    call FRsearch(n, tree, r2, neighbors)
    call kdtree2_destroy(tree)
end subroutine fixedRadSearch

! build the sparse matrix of co-location probabilities
subroutine build_PmatSparse(n, denom, start, finish, Pmat)
    integer,             intent(in   ) :: n ! number of entries in the matrix
    real(pDki),          intent(in   ) :: denom ! denominator of the exponential in the co-location probability density
    integer(kind=8),     intent(in   ) :: start(n), finish(n) ! indices (in the Distmat vectors) for the start and finish of each column of Distmat
    type(SparseMatType), intent(inout) :: Pmat(:) ! sparse matrix of co-location probabilities

    ! ========================== LOCAL VARIABLES ===============================
    integer :: i

    ! calculate the encounter probabilities
    ! NOTE: the leading constant term is neglected here, due to the normalization
    ! NOTE: Pmat%val is already squared separation distance because of the kD search
    Pmat%val = exp(Pmat%val / denom) / (-denom * pi)

    ! normalize the columns of Pmat to sum to 1
    ! do i = 1, n
    !    Pmat(start(i) : finish(i))%val = Pmat(start(i) : finish(i))%val /&
    !                                     sum(Pmat(start(i) : finish(i))%val)
    ! enddo
end subroutine build_PmatSparse

! build the sparse mass-transfer matrix
! \vec Emat = \vec I - beta * [diag(\vec Pmat \vec 1)  - \vec Pmat]
!           = \vec I - beta * [diag(rowsum(\vec Pmat)) - \vec Pmat]
!           = \vec I - beta * diag(rowsum(\vec Pmat)) + beta * \vec Pmat
subroutine build_EmatSparse(n, Nclose, Emat, beta)
    integer,             intent(in   ) :: n ! total number of particles
    integer(kind=8),     intent(in   ) :: Nclose ! total number of neighbor particles found by the kD search (i.e, the length of the vectors in Distmat)
    type(SparseMatType), intent(inout) :: Emat(:) ! sparse mass-transfer matrix
    real(pDki),          intent(in   ) :: beta ! beta parameter
    ! ========================== LOCAL VARIABLES ===============================
    integer(kind=8)    :: i, j
    integer(kind=8)    :: diag(n) ! linear indices of the diagonal elements of Pmat
    real(pDki)         :: rowsum(n), colsum(n) ! arrays holding the row/col sums of Pmat

    ! compute the rowsums of Pmat
    rowsum = 0.0_pDki
    colsum = 0.0_pDki
    do i = 1, Nclose
        rowsum(Emat(i)%row) = rowsum(Emat(i)%row) + Emat(i)%val
        colsum(Emat(i)%col) = colsum(Emat(i)%col) + Emat(i)%val
        if (Emat(i)%row == Emat(i)%col) then
            diag(Emat(i)%row) = i
        endif
    enddo

    do i = 1, Nclose
        Emat(i)%val = Emat(i)%val / ((rowsum(Emat(i)%row) + colsum(Emat(i)%col)) / 2.0_pDki)
    enddo

    rowsum = 0.0_pDki
    do i = 1, Nclose
        rowsum(Emat(i)%row) = rowsum(Emat(i)%row) + Emat(i)%val
    enddo

    ! this is the I - beta * diag(P * 1) step
    rowsum = 1.0_pDki - beta * rowsum

    ! this is the beta * \vec Pmat step
    Emat%val = beta * Emat%val

    ! finally, add them together
    Emat(diag)%val =  Emat(diag)%val + rowsum
end subroutine build_EmatSparse

! randomly shuffles an array of integers
! source: https://www.rosettacode.org/wiki/Knuth_shuffle#Fortran
subroutine shuffleInt(a)
    integer, intent(inout) :: a(:)
    integer                :: i, randpos, temp
    real                   :: r

    do i = size(a), 2, -1
        call random_number(r)
        randpos = int(r * i) + 1
        temp = a(randpos)
        a(randpos) = a(i)
        a(i) = temp
    enddo
end subroutine shuffleInt

! this builds a KD tree
subroutine maketree(tree2, d, n, locs)
    type(kdtree2), pointer, intent(  out) :: tree2 ! this is the KD tree
    integer,                intent(in   ) :: d, n ! number of spatial dimensions, number of particles
    real(kdkind),           intent(in   ) :: locs(d, n) ! location array for particles, with dimension d x n (number of spatial dimensions x number of particles)

    ! build the tree
    tree2 => kdtree2_create(locs, dim = d, sort = .false., rearrange = .true.)
        ! currently don't see a need to sort, as false appears to be quicker
        ! rearrange = true also appears to be quicker
end subroutine maketree

! this searches an already built KD tree
subroutine FRsearch(n, tree, r2, neighbors)
    integer,                intent(in   ) :: n ! number of particles
    type(kdtree2), pointer, intent(in   ) :: tree ! the KD tree
    real(kdkind),           intent(in   ) :: r2 ! squared search radius
    type(kDRS_ResultType),  intent(  out) :: neighbors(:) ! holds the results of the kD tree fixed radius search

    ! ========================== LOCAL VARIABLES ===============================
    integer                           :: i
    type(kdtree2_result), allocatable :: results(:) ! results array from KD tree module
    integer                           :: num_alloc

    ! allocate results as big as it could possibly be
    ! there's probably a more memory-efficient way to do this
    num_alloc = n
    allocate (results(num_alloc))

    ! loop over all particles
    do i = 1, n
        ! the type of search used here finds all the particles within
        ! squared distance r2 from the i^th particle in the list
        ! the hard-coded 0 is the 'correlation time' of the search
        ! as far as i can tell, this means a particle, itself, is included the
        ! idx list, while 1 would leave the particle, itself, out of the list
        call kdtree2_r_nearest_around_point(tree, i, 0, r2, neighbors(i)%num, num_alloc, results)

        ! allocate these based on how many nearby particles were found
        allocate (neighbors(i)%idx(neighbors(i)%num), neighbors(i)%rad(neighbors(i)%num))

        ! put the results in the neighbors array
        neighbors(i)%idx = results(1 : neighbors(i)%num)%idx
        neighbors(i)%rad = results(1 : neighbors(i)%num)%dis
    enddo

    deallocate (results)
end subroutine FRsearch


! this subroutine initially places the particles for Heaviside IC and random initial positions
subroutine InitialSpacingAndIC(p, Dlims, coreNp, xmidpt, ymidpt, nnx, nny, spillX, spillY)

    real(pDki),         intent(in   ) :: Dlims(4) ! subdomain boundary limits
    integer,            intent(in   ) :: coreNp ! number of particles in local subdomain
    integer,            intent(in   ) :: nnx, nny ! grid specifications
    real(pDki),         intent(in   ) :: xmidpt, ymidpt ! midpoint of domain for IC
    real(pDki),         intent(  out) :: spillX, spillY ! location of spill
    type(ParticleType), intent(inout) :: p(:) ! particle array


            ! call random_number(p%loc(1))
            call random_number(p%loc(2))
            ! p(1 : coreNp)%loc(1) = Dlims(1) + (Dlims(2) - Dlims(1)) * p(1 : coreNp)%loc(1)
            p(1 : coreNp)%loc(1) = linspace(Dlims(1),Dlims(2),coreNp)
            p(1 : coreNp)%loc(2) = Dlims(3) + (Dlims(4) - Dlims(3)) * p(1 : coreNp)%loc(2)  

            ! initialize the particle type variables
            p%active = .false.
            p(1 : coreNp)%active = .true.
            p%mass = 0.0_pDki
            p(1 : coreNp)%mass = 0.0_pDki

            ! Heaviside initial condition
            where (p(1 : coreNp)%loc(1) >= xmidpt) 
                p(1 : coreNp)%mass = 1.0_pDki
            endwhere


            ! point source initial condition near center

            ! if (my_rank == (floor(nnx/real(2,pDki))*nny + floor(nny/real(2,pDki)))) then
            !     p(coreNp/2)%mass = 1.0_pDki
            ! endif

            ! if (any(p(1 : coreNp)%mass == 1.0_pDki) .eqv. .true.) then
            !     Write(*,*) "Core number",my_rank,"has the intial mass."
            ! endif

            ! if (any(p(1 : coreNp)%mass == 1.0_pDki) .eqv. .true.) then
            !     spillX = p(coreNp/2)%loc(1)
            !     spillY = p(coreNp/2)%loc(2)
            ! endif       


end subroutine InitialSpacingAndIC







! Potentially use these subroutines in the future

! subroutine RandomlySpacedPointSource(masterP,xleftlim,xrightlim,ybottomlim,ytoplim,spillX,spillY,idx,numspills,Np)

!     type(ParticleType), intent(inout) :: masterP(:) ! particle array
!     real(pDki),         intent(inout) :: spillX(:), spillY(:)
!     integer,            intent(in   ) :: idx(:), numspills, Np
!     real(pDki),         intent(in   ) :: xleftlim, xrightlim, ybottomlim, ytoplim 

!     integer                           :: idxSpill(numspills), idxSmall(size(idx))

!     if (my_rank == master) then

!                call random_number(masterP%loc(1))
!                call random_number(masterP%loc(2))
!                masterP(1 : Np)%loc(1) = xleftlim + (xrightlim - xleftlim) * masterP(1 : Np)%loc(1)
!                masterP(1 : Np)%loc(2) = ybottomlim + (ytoplim-ybottomlim) * masterP(1 : Np)%loc(2)

!                masterP%active = .false.
!                masterP(1 : Np)%mass = 0.0_pDki

!                     where (masterP(1 : Np)%loc(1) >= (0.35*xrightlim) .and. masterP(1 : Np)%loc(1) <= (0.65*xrightlim) .and. & 
!                         masterP(1 : Np)%loc(2) >=  (0.35*ytoplim) .and. masterP(1 : Np)%loc(2) <= (0.65*ytoplim))

!                         masterP%active = .true. 

!                     endwhere

!                     idxSmall = pack(idx, masterP%active)
!                     call shuffleInt(idxSmall)

!                     idxSpill = idxSmall(1 : numspills)

!                     masterP(idxSpill)%mass = 1.0_pDki/real(numspills,pDki)

!                     spillX = masterP(idxSpill)%loc(1)
!                     spillY = masterP(idxSpill)%loc(2) 
                    
!                     masterP%active = .false.

!     endif

! end subroutine RandomlyspacedPointSource



! subroutine GetMyParticles(DDC_WhichWay, tmp_master, idxActive, idx, Dlims, nny, num_cores, Np, p)

!     type(ParticleType), intent(in   ) :: tmp_master(:)
!     type(ParticleType), intent(  out) :: p(:)
!     real(pDki),         intent(in   ) :: Dlims(4)
!     integer,            intent(in   ) :: DDC_WhichWay, nny, num_cores, Np, idx(:)

!     integer                           :: Nactive, idxActive(size(idx))
!     type(ParticleType)                :: masterP(Np)

!     masterP = tmp_master

!     select case (DDC_WhichWay)

!                     case (1) ! Vertical Slices

!                         if (my_rank /= (num_cores - 1)) then 
!                             where (masterP(1 : Np)%loc(1) >= Dlims(1) .and. masterP(1 : Np)%loc(1) < Dlims(2))

!                                 masterP(1 : Np)%active = .true.

!                             endwhere     

!                                 idxActive = 0
!                                 Nactive = count(masterP%active)
!             !                     ! Write(*,*) "My core has",Nactive,"particles!"
!                                 idxActive(1 : Nactive) = pack(idx, masterP%active)
!                                 p(1 : Nactive) = masterP(idxActive(1 : Nactive)) 
!                         else
!                             where (masterP(1 : Np)%loc(1) >= Dlims(1) .and. masterP(1 : Np)%loc(1) <= Dlims(2))

!                                 masterP(1 : Np)%active = .true.

!                             endwhere     
!                                 idxActive = 0
!                                 Nactive = count(masterP%active)
!                                 ! Write(*,*) "My core has",Nactive,"particles!"
!                                 idxActive(1 : Nactive) = pack(idx, masterP%active)
!                                 p(1 : Nactive) = masterP(idxActive(1 : Nactive)) 

!                         endif


!                     case (2) ! Checkerboard

!                         ! All cores not along top row or right column 
!                         ! Assign particles that are inside boundaries and equal to left and bottom thresholds

!                         if (my_rank /= 0) then

!                             if (nny > 1) then
!                                 if (mod(my_rank,nny) /= nny-1 .and. my_rank < num_cores - nny) then
!                                     where (masterP(1 : Np)%loc(1) >= Dlims(1) .and. masterP(1 : Np)%loc(2) >= Dlims(3) .and. &
!                                         masterP(1 : Np)%loc(1) < Dlims(2) .and. masterP(1 : Np)%loc(2) < Dlims(4))

!                                             masterP(1 : Np)%active = .true.

!                                     endwhere
!                                             idxActive = 0
!                                             Nactive = count(masterP%active)
!                                             idxActive(1 : Nactive) = pack(idx, masterP%active)
!                                             p(1 : Nactive) = masterP(idxActive(1 : Nactive))

!                                 endif
!                             endif
!                         endif

!                         ! Own case for Master since mod(my_rank,nny-1)=0, but we want to treat this area as the same
!                         ! as above case

!                         if (my_rank == 0) then

!                             if (num_cores > 2) then
!                                 where (masterP(1 : Np)%loc(1) >= Dlims(1) .and. masterP(1 : Np)%loc(2) >= Dlims(3) .and. &
!                                     masterP(1 : Np)%loc(1) < Dlims(2) .and. masterP(1 : Np)%loc(2) < Dlims(4))

!                                         masterP(1 : Np)%active = .true.

!                                 endwhere
!                                         idxActive = 0
!                                         Nactive = count(masterP%active)
!                                         idxActive(1 : Nactive) = pack(idx, masterP%active)
!                                         p(1 : Nactive) = masterP(idxActive(1 : Nactive))

!                             ! If only 2 cores, they will be side by side.. Need to include top boundary as well.              
!                             else if (num_cores == 2) then
!                                  where (masterP(1 : Np)%loc(1) >= Dlims(1) .and. masterP(1 : Np)%loc(2) >= Dlims(3) .and. &
!                                     masterP(1 : Np)%loc(1) < Dlims(2) .and. masterP(1 : Np)%loc(2) <= Dlims(4))

!                                         masterP(1 : Np)%active = .true.

!                                 endwhere     
!                                         idxActive = 0
!                                         Nactive = count(masterP%active)
!                                         ! Write(*,*) "Number of active particles on core",my_rank, "is",Nactive
!                                         idxActive(1 : Nactive) = pack(idx, masterP%active)
!                                         p(1 : Nactive) = masterP(idxActive(1 : Nactive))           

!                             ! If only 1 core, include all initial particles and activate them. 
!                             else
!                                         idxActive = 0
!                                         masterP(1 : Np)%active = .true.
!                                         Nactive = count(masterP%active)
!                                         ! Write(*,*) "My core has",Nactive,"particles!"
!                                         idxActive(1 : Nactive) = pack(idx, masterP%active)
!                                         p(1 : Nactive) = masterP(idxActive(1 : Nactive))

!                             endif 

!                         endif

!                         ! Cores along top row of responsibilities (exluding the top right corner)
!                         ! Include top boundary of particles as well

!                         if (my_rank /= 0) then
!                             if (nny > 1) then 
!                                 if(mod(my_rank,nny) == nny-1 .and. my_rank /= num_cores - 1) then
!                                     where (masterP(1 : Np)%loc(1) >= Dlims(1) .and. masterP(1 : Np)%loc(2) >= Dlims(3) .and. &
!                                         masterP(1 : Np)%loc(1) < Dlims(2) .and. masterP(1 : Np)%loc(2) <= Dlims(4))

!                                             masterP(1 : Np)%active = .true.

!                                     endwhere
!                                             idxActive = 0
!                                             Nactive = count(masterP%active)
!                                             idxActive(1 : Nactive) = pack(idx, masterP%active)
!                                             p(1 : Nactive) = masterP(idxActive(1 : Nactive))

!                                 endif
!                             endif
!                         endif

!                         ! Top right core alone
!                         ! Include top and right boundary
!                         if (my_rank /= 0) then
!                             if (my_rank == num_cores - 1) then
!                                 where (masterP(1 : Np)%loc(1) >= Dlims(1) .and. masterP(1 : Np)%loc(2) >= Dlims(3) .and. &
!                                     masterP(1 : Np)%loc(1) <= Dlims(2) .and. masterP(1 : Np)%loc(2) <= Dlims(4))

!                                         masterP(1 : Np)%active = .true.            

!                                 endwhere  
!                                         idxActive = 0
!                                         Nactive = count(masterP%active)
!                                         idxActive(1 : Nactive) = pack(idx, masterP%active)
!                                         p(1 : Nactive) = masterP(idxActive(1 : Nactive))                 
!                             endif
!                         endif


!                         ! Right column of cores (excluding top right core)
!                         ! Bottom and left boundaries included along with the right boundary

!                         if (my_rank >= num_cores - nny .and. my_rank /= num_cores - 1) then
!                             where (masterP(1 : Np)%loc(1) >= Dlims(1) .and. masterP(1 : Np)%loc(2) >= Dlims(3) .and. &
!                                 masterP(1 : Np)%loc(1) <= Dlims(2) .and. masterP(1 : Np)%loc(2) < Dlims(4))

!                                     masterP(1 : Np)%active = .true.          

!                             endwhere  
!                                     idxActive = 0
!                                     Nactive = count(masterP%active)
!                                     idxActive(1 : Nactive) = pack(idx, masterP%active)
!                                     p(1 : Nactive) = masterP(idxActive(1 : Nactive))     
!                         endif
                    
!                     case default
!                         print *, '**** ERROR: Invalid transfer mode ****'

!             end select



! end subroutine GetMyParticles




end module mod_DDC_2D_mpi_paper