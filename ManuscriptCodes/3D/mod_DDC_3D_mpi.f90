module mod_DDC_2D_mpi
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
integer               :: dim = 3 ! NOTE: this is the hard-coded 3 spatial dimension
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
    real(pDki) :: loc(3) ! real-valued spatial location
    real(pDki) :: mass ! real-valued mass
    logical    :: active ! whether a given particle is currently active
    logical    :: ghost(26) ! indicates whether a particle is a ghost particle and which direction it will be sent.
    logical    :: jumped(26) ! indicates whether the particle jumped subdomains and which way it went.

                                ! These two arrays correspond to each other directionally 
                                ! In order, the logical entries for array ghost() and jumped() will correspond to these directions:
                                ! Same plane 1-8: L, R, D, U, DL, DR, UL, UR
                                ! Plane above 9-17: M (directly above subdomain cuboid), L, R, D, U, DL, DR, UL, UR
                                ! Plane below 18-26: M (directly below subdomain cuboid), L, R, D, U, DL, DR, UL, UR
end type ParticleType

! sparse matrix type in coordinate format
type SparseMatType
    integer    :: row
    integer    :: col
    real(pDki) :: val
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
    integer, parameter        :: number = 5 ! number of fields in the ParticleType
    integer                   :: block_lengths(number)
    integer(mpi_address_kind) :: displacements(number)
    integer                   :: typelist(number)
    real(pDki)                :: r
    logical                   :: log
    integer                   :: errcode, ierr, disp1, disp2, disp3, disp4, disp5

    typelist = (/ realType, realType, mpi_logical, mpi_logical, mpi_logical /)
    block_lengths = (/ 3, 1, 1, 26, 26 /)

    disp1 = 0
    disp2 = int(3,mpi_address_kind)*sizeof(r)
    disp3 = disp2 + int(sizeof(r),mpi_address_kind)
    disp4 = disp3 + sizeof(log)
    disp5 = disp4 + int(26,mpi_address_kind)*sizeof(log)

    displacements = (/ disp1, disp2, disp3, disp4, disp5 /)

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
    real(pDki)                        :: normvec1(np), normvec2(np), normvec3(np) ! vector which will hold Normal(0, 1) values

    ! call N(0, 1) generator
    ! call init_random_seed()

    call box_mullerp(np, normvec1)
    call box_mullerp(np, normvec2)
    call box_mullerp(np, normvec3)


    p(idxActive)%loc(1) = p(idxActive)%loc(1) + sqrt(2.0_pDki * D * dt) * normvec1
    p(idxActive)%loc(2) = p(idxActive)%loc(2) + sqrt(2.0_pDki * D * dt) * normvec2
    p(idxActive)%loc(3) = p(idxActive)%loc(3) + sqrt(2.0_pDki * D * dt) * normvec3
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
subroutine swapDomains(pVecSize, Nactive, idxActive, Dlims, paddist, p, nnx, nny, DDC_WhichWay)
    integer,            intent(in   ) :: pVecSize ! total size of the local particle array
    integer,            intent(in   ) :: DDC_WhichWay ! which DDC method is being used
    integer,            intent(in   ) :: nny, nnx ! bring in the vert/horiz nodes
    integer,            intent(inout) :: Nactive ! number of active particles
    integer,            intent(inout) :: idxActive(:) ! indices of the active particles
    real(pDki),         intent(in   ) :: Dlims(6) ! subdomain boundary limits
    real(pDki),         intent(in   ) :: paddist ! ghost particle padding distance
    type(ParticleType), intent(inout) :: p(:) ! particle array

    ! ========================== LOCAL VARIABLES ===============================
    integer, allocatable :: idx(:) ! these are all used for bookkeeping in sending and receiving (hopefully what they are should be apparent by usage)
    integer              :: active(Nactive)
    integer              :: idxBig(5*Nactive)
    integer              :: numSendLeft, numSendRight, numSendUp, numSendDown, numSendDR, &
                                numSendDL, numSendUR, numSendUL, numSendAbove, numSendAboveU, numSendAboveD, &
                                numSendAboveL, numSendAboveR, numSendAboveDL, numSendAboveDR, numSendAboveUL, &
                                numSendAboveUR, numSendBelow, numSendBelowU, numSendBelowD, numSendBelowL, &
                                numSendBelowR, numSendBelowDL, numSendBelowDR, numSendBelowUL, numSendBelowUR, &
                                recvSizeBelow, recvSizeBelowU, recvSizeBelowD, recvSizeBelowL, recvSizeBelowR, &
                                recvSizeBelowDR, recvSizeBelowDL, recvSizeBelowUR, recvSizeBelowUL, recvSizeAbove, &
                                recvSizeAboveU, recvSizeAboveD, recvSizeAboveL, recvSizeAboveR, recvSizeAboveDR, &
                                recvSizeAboveDL, recvSizeAboveUL, recvSizeAboveUR, recvSizeRight, recvSizeLeft, &
                                recvSizeDown, recvSizeUp, recvSizeDL, recvSizeDR, recvSizeUL, recvSizeUR
    integer              :: nJumpedRight, nJumpedLeft, nJumpedDown, nJumpedUp, nJumpedUL, nJumpedUR, nJumpedDL, nJumpedDR, &
                                nJumpedAbove, nJumpedAboveU, nJumpedAboveD, nJumpedAboveL, nJumpedAboveR, nJumpedAboveDR, &
                                nJumpedAboveDL, nJumpedAboveUL, nJumpedAboveUR, nJumpedBelow, nJumpedBelowU, nJumpedBelowD, &
                                nJumpedBelowL, nJumpedBelowR, nJumpedBelowDR, nJumpedBelowDL, nJumpedBelowUL, nJumpedBelowUR, &
                                totRecv, nGhostRight, nGhostLeft, nGhostDown, nGhostUp, nGhostDR, nGhostDL, nGhostUR, nGhostUL, &
                                nGhostAbove, nGhostAboveU, nGhostAboveD, nGhostAboveL, nGhostAboveR, nGhostAboveDR, &
                                nGhostAboveDL, nGhostAboveUL, nGhostAboveUR, nGhostBelow, nGhostBelowU, nGhostBelowD, &
                                nGhostBelowL, nGhostBelowR, nGhostBelowDL, nGhostBelowDR, nGhostBelowUL, nGhostBelowUR, &
                                totRecvLR, totRecvLRU, totRecvCard, totRecvDiag1, totRecvDiag2, totRecvDiag3, totRecvDiag4, &
                                totRecvBelow1, totRecvBelow2, totRecvBelow3, totRecvBelow4, totRecvBelow5, totRecvBelow6, &
                                totRecvBelow7, totRecvBelow8, totRecvAbove1, totRecvAbove2, totRecvAbove3, totRecvAbove4, &
                                totRecvAbove5, totRecvAbove6, totRecvAbove7, totRecvAbove8, totRecvAbove9
    integer              :: i, notActive, xysheet

    type(ParticleType), allocatable :: tmp_pLeft(:), tmp_pRight(:), tmp_pDown(:), tmp_pUp(:),&
                                tmp_pDR(:), tmp_pDL(:), tmp_pUR(:), tmp_pUL(:), &
                                tmp_pAbove(:), tmp_pAboveU(:), tmp_pAboveD(:), tmp_pAboveL(:), &
                                tmp_pAboveR(:), tmp_pAboveDL(:), tmp_pAboveDR(:), tmp_pAboveUL(:), &
                                tmp_pAboveUR(:), tmp_pBelow(:), tmp_pBelowU(:), tmp_pBelowD(:), &
                                tmp_pBelowL(:), tmp_pBelowR(:), tmp_pBelowDL(:), tmp_pBelowDR(:), &
                                tmp_pBelowUL(:), tmp_pBelowUR(:) ! particle arrays for sending

    type(ParticleType), allocatable :: recv_pLeft(:), recv_pRight(:), recv_pDown(:),&
                            recv_pUp(:), recv_pUR(:), recv_pUL(:), recv_pDR(:),&
                            recv_pDL(:), recv_pAbove(:), recv_pAboveU(:), recv_pAboveD(:), &
                            recv_pAboveL(:), recv_pAboveR(:), recv_pAboveDR(:), recv_pAboveDL(:), &
                            recv_pAboveUL(:), recv_pAboveUR(:), recv_pBelow(:), recv_pBelowU(:), &
                            recv_pBelowD(:), recv_pBelowL(:), recv_pBelowR(:), recv_pBelowDL(:), &
                            recv_pBelowDR(:), recv_pBelowUL(:), recv_pBelowUR(:), recv_pTotal(:) ! particle arrays for receiving

    xysheet = nny*nnx

    ! initialize all of the relevant working info
    p%jumped(1) = .false.
    p%jumped(2) = .false.
    p%jumped(3) = .false.
    p%jumped(4) = .false.
    p%jumped(5) = .false.
    p%jumped(6) = .false.
    p%jumped(7) = .false.
    p%jumped(8) = .false.
    p%jumped(9) = .false.
    p%jumped(10) = .false.
    p%jumped(11) = .false.
    p%jumped(12) = .false.
    p%jumped(13) = .false.
    p%jumped(14) = .false.
    p%jumped(15) = .false.
    p%jumped(16) = .false.
    p%jumped(17) = .false.
    p%jumped(18) = .false.
    p%jumped(19) = .false.
    p%jumped(20) = .false.
    p%jumped(21) = .false.
    p%jumped(22) = .false.
    p%jumped(23) = .false.
    p%jumped(24) = .false.
    p%jumped(25) = .false.
    p%jumped(26) = .false.

    p%ghost(1) = .false.
    p%ghost(2) = .false.
    p%ghost(3) = .false.
    p%ghost(4) = .false.
    p%ghost(6) = .false.
    p%ghost(7) = .false.
    p%ghost(8) = .false.
    p%ghost(9) = .false.
    p%ghost(10) = .false.
    p%ghost(11) = .false.
    p%ghost(12) = .false.
    p%ghost(13) = .false.
    p%ghost(14) = .false.
    p%ghost(15) = .false.
    p%ghost(16) = .false.
    p%ghost(17) = .false.
    p%ghost(18) = .false.
    p%ghost(19) = .false.
    p%ghost(20) = .false.
    p%ghost(21) = .false.
    p%ghost(22) = .false.
    p%ghost(23) = .false.
    p%ghost(24) = .false.
    p%ghost(25) = .false.
    p%ghost(26) = .false.

    numSendLeft = 0
    numSendRight = 0
    numSendUp = 0
    numSendDown = 0
    numSendUL = 0
    numSendUR = 0
    numSendDL = 0
    numSendDR = 0
    numSendAbove = 0
    numSendAboveU = 0
    numSendAboveD = 0
    numSendAboveL = 0
    numSendAboveR  = 0
    numSendAboveDR = 0
    numSendAboveDL = 0
    numSendAboveUL = 0
    numSendAboveUR = 0
    numSendBelow = 0
    numSendBelowU = 0
    numSendBelowD = 0
    numSendBelowL = 0
    numSendBelowR = 0
    numSendBelowDL = 0
    numSendBelowDR = 0
    numSendBelowUR = 0
    numSendBelowUL = 0
    recvSizeRight = 0
    recvSizeLeft = 0
    recvSizeDown = 0
    recvSizeUp = 0
    recvSizeUR = 0
    recvSizeUL = 0
    recvSizeDL = 0
    recvSizeDR = 0
    recvSizeAbove = 0
    recvSizeAboveU = 0
    recvSizeAboveD = 0
    recvSizeAboveL = 0
    recvSizeAboveR = 0
    recvSizeAboveDL = 0
    recvSizeAboveDR = 0
    recvSizeAboveUL = 0
    recvSizeAboveUR = 0
    recvSizeBelow = 0
    recvSizeBelowU = 0
    recvSizeBelowD = 0
    recvSizeBelowL = 0
    recvSizeBelowR = 0
    recvSizeBelowDL = 0
    recvSizeBelowDR = 0
    recvSizeBelowUL = 0
    recvSizeBelowUR = 0

    nJumpedRight = 0
    nJumpedLeft = 0
    nJumpedDown = 0
    nJumpedUp = 0
    nJumpedDL = 0
    nJumpedDR = 0
    nJumpedUL = 0
    nJumpedUR = 0
    nJumpedAbove = 0
    nJumpedAboveU = 0
    nJumpedAboveD = 0
    nJumpedAboveL = 0
    nJumpedAboveR = 0
    nJumpedAboveDL = 0
    nJumpedAboveDR = 0
    nJumpedAboveUL = 0
    nJumpedAboveUR = 0
    nJumpedBelow = 0
    nJumpedBelowU = 0
    nJumpedBelowD = 0
    nJumpedBelowL = 0
    nJumpedBelowR = 0
    nJumpedBelowDL = 0
    nJumpedBelowDR = 0
    nJumpedBelowUL = 0
    nJumpedBelowUR = 0
    nGhostRight = 0
    nGhostLeft = 0
    nGhostDown = 0
    nGhostUp = 0
    nGhostUL = 0
    nGhostUR = 0
    nGhostDL = 0
    nGhostDR = 0
    nGhostAbove = 0
    nGhostAboveU = 0
    nGhostAboveD = 0
    nGhostAboveL = 0
    nGhostAboveR = 0
    nGhostAboveDL = 0
    nGhostAboveDR = 0
    nGhostAboveUL = 0
    nGhostAboveUR = 0
    nGhostBelow = 0
    nGhostBelowU = 0
    nGhostBelowD = 0
    nGhostBelowL = 0
    nGhostBelowR = 0
    nGhostBelowDL = 0
    nGhostBelowDR = 0
    nGhostBelowUL = 0
    nGhostBelowUR = 0

    ! this array will be used as an argument in the particle array
    allocate(idx(Nactive))
    active = idxActive(1 : Nactive)


! In order, the logical entries for array ghost() and jumped() will correspond to these directions:
        ! Same plane 1-8: L, R, D, U, DL, DR, UL, UR
        ! Plane above 9-17: M (directly above subdomain cuboid), L, R, D, U, DL, DR, UL, UR
        ! Plane below 18-26: M (directly below subdomain cuboid), L, R, D, U, DL, DR, UL, UR

    select case(DDC_WhichWay)


    !vertical slices
    case(1)

    ! first send particles to the left
    if (my_rank > 0) then
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
    endif
    ! if there are > 0 particles to be received, receive them
    if (my_rank < num_cores - 1 .and. recvSizeRight > 0) then
        call mpi_recv(recv_pRight(1 : recvSizeRight), recvSizeRight, mpi_pType,&
                      my_rank + 1, tag + my_rank, mpi_comm_world, status, ierror)

        ! Write(*,*) "Receiving on the right on core",my_rank,"is",recv_pRight(1 : recvSizeRight)
    endif

    ! now send particles to the right (this works the same as above, and, thus, is
    ! not commented)
    if (my_rank < num_cores - 1) then
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
    endif
    if (my_rank > 0 .and. recvSizeLeft > 0) then
        call mpi_recv(recv_pLeft(1 : recvSizeLeft), recvSizeLeft, mpi_pType,&
                      my_rank - 1, tag + my_rank, mpi_comm_world, status, ierror)
    endif

    !checkerboard method
    case(2)
    ! first send particles to the left (tag for all other left sends too)

    ! we will use dynamic allocation throughout this routine, as 3-d required a lot of memory. 
    allocate(tmp_pLeft(Nactive),tmp_pDL(Nactive),tmp_pUL(Nactive),tmp_pAboveL(Nactive),&
                    tmp_pAboveUL(Nactive),tmp_pAboveDL(Nactive),tmp_pBelowL(Nactive),&
                    tmp_pBelowUL(Nactive),tmp_pBelowDL(Nactive))

    if (mod(my_rank,xysheet) > (nny - 1)) then

        ! pair down to which core (on which plane) to send these particles 
        where (p(active)%loc(1) < Dlims(1))
            p(active)%jumped(1) = .true.
        endwhere

        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) > Dlims(4))
            p(active)%jumped(1) = .false.
            p(active)%jumped(7) = .true.
        endwhere

        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) > Dlims(4) &
                .and. p(active)%loc(3) > Dlims(6))
            p(active)%jumped(7) = .false.
            p(active)%jumped(16) = .true.
        endwhere

        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) > Dlims(4) &
                .and. p(active)%loc(3) < Dlims(5))
            p(active)%jumped(7) = .false.
            p(active)%jumped(25) = .true.
        endwhere

        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) < Dlims(3))
            p(active)%jumped(1) = .false.
            p(active)%jumped(5) = .true.
        endwhere

        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) < Dlims(3) &
                .and. p(active)%loc(3) > Dlims(6))
            p(active)%jumped(5) = .false.
            p(active)%jumped(14) = .true.
        endwhere

        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) < Dlims(3) &
                .and. p(active)%loc(3) < Dlims(5))
            p(active)%jumped(5) = .false.
            p(active)%jumped(23) = .true.
        endwhere

        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) > Dlims(3) &
                .and. p(active)%loc(2) < Dlims(4) .and. p(active)%loc(3) > Dlims(6))
            p(active)%jumped(1) = .false.
            p(active)%jumped(10) = .true.
        endwhere

        where (p(active)%loc(1) < Dlims(1) .and. p(active)%loc(2) > Dlims(3) &
                .and. p(active)%loc(2) < Dlims(4) .and. p(active)%loc(3) < Dlims(5))
            p(active)%jumped(1) = .false.
            p(active)%jumped(19) = .true.
        endwhere

        idx = 0
        nJumpedLeft = count(p(active)%jumped(1))
        idx(1 : nJumpedLeft) = pack(active, p(active)%jumped(1))
        ! put it in a temporary array for sending
        tmp_pLeft(1 : nJumpedLeft) = p(idx(1 : nJumpedLeft))

        idx = 0
        nJumpedDL = count(p(active)%jumped(5))
        idx(1 : nJumpedDL) = pack(active, p(active)%jumped(5))
        ! put it in a temporary array for sending
        tmp_pDL(1 : nJumpedDL) = p(idx(1 : nJumpedDL))

        idx = 0
        nJumpedUL = count(p(active)%jumped(7))
        idx(1 : nJumpedUL) = pack(active, p(active)%jumped(7))
        ! put it in a temporary array for sending
        tmp_pUL(1 : nJumpedUL) = p(idx(1 : nJumpedUL))

        idx = 0
        nJumpedAboveL = count(p(active)%jumped(10))
        idx(1 : nJumpedAboveL) = pack(active, p(active)%jumped(10))
        tmp_pAboveL(1 : nJumpedAboveL) = p(idx(1 : nJumpedAboveL))

        idx = 0
        nJumpedAboveUL = count(p(active)%jumped(16))
        idx(1 : nJumpedAboveUL) = pack(active, p(active)%jumped(16))
        tmp_pAboveUL(1 : nJumpedAboveUL) = p(idx(1 : nJumpedAboveUL))

        idx = 0
        nJumpedAboveDL = count(p(active)%jumped(14))
        idx(1 : nJumpedAboveDL) = pack(active, p(active)%jumped(14))
        tmp_pAboveDL(1 : nJumpedAboveDL) = p(idx(1 : nJumpedAboveDL))

        idx = 0
        nJumpedBelowL = count(p(active)%jumped(19))
        idx(1 : nJumpedBelowL) = pack(active, p(active)%jumped(19))
        tmp_pBelowL(1 : nJumpedBelowL) = p(idx(1 : nJumpedBelowL))

        idx = 0
        nJumpedBelowUL = count(p(active)%jumped(25))
        idx(1 : nJumpedBelowUL) = pack(active, p(active)%jumped(25))
        tmp_pBelowUL(1 : nJumpedBelowUL) = p(idx(1 : nJumpedBelowUL))

        idx = 0
        nJumpedBelowDL = count(p(active)%jumped(23))
        idx(1 : nJumpedBelowDL) = pack(active, p(active)%jumped(23))
        tmp_pBelowDL(1 : nJumpedBelowDL) = p(idx(1 : nJumpedBelowDL))



        ! tag the particles that will be sent as ghost particles
        where (p(active)%loc(1) >= Dlims(1) .and. p(active)%loc(1) < (Dlims(1) + paddist))
            p(active)%ghost(1) = .true.
        endwhere

        where((p(active)%ghost(1) .eqv. .true.) .and. (p(active)%loc(3) > (Dlims(6) - paddist)) .and. &
                (p(active)%loc(3) <= Dlims(6)))
            p(active)%ghost(10) = .true.
        endwhere

        where((p(active)%ghost(1) .eqv. .true.) .and. p(active)%loc(3) < (Dlims(5) + paddist) .and. &
                p(active)%loc(3) >= Dlims(5))
            p(active)%ghost(19) = .true.
        endwhere

        where((p(active)%ghost(1) .eqv. .true.) .and. p(active)%loc(2) > (Dlims(4) - paddist) .and. &
                p(active)%loc(2) <= Dlims(4))
            p(active)%ghost(7) = .true.
        endwhere

        where((p(active)%ghost(7) .eqv. .true.) .and. p(active)%loc(3) > (Dlims(6) - paddist) .and. &
                p(active)%loc(3) <= Dlims(6))
            p(active)%ghost(16) = .true.
        endwhere

        where((p(active)%ghost(7) .eqv. .true.) .and. p(active)%loc(3) < (Dlims(5) + paddist) .and. &
                p(active)%loc(3) >= Dlims(5))
            p(active)%ghost(25) = .true.
        endwhere

        where((p(active)%ghost(1) .eqv. .true.) .and. p(active)%loc(2) < (Dlims(3) + paddist) .and. &
                p(active)%loc(2) >= Dlims(3))
            p(active)%ghost(5) = .true.
        endwhere

        where((p(active)%ghost(5) .eqv. .true.) .and. p(active)%loc(3) > (Dlims(6) - paddist) .and. &
                p(active)%loc(3) <= Dlims(6))
            p(active)%ghost(14) = .true.
        endwhere

        where((p(active)%ghost(5) .eqv. .true.) .and. p(active)%loc(3) < (Dlims(5) + paddist) .and. &
                p(active)%loc(3) >= Dlims(5))
            p(active)%ghost(23) = .true.
        endwhere

        idx = 0
        nGhostLeft = count(p(active)%ghost(1))
        idx(1 : nGhostLeft) = pack(active, p(active)%ghost(1))
        numSendLeft = nJumpedLeft + nGhostLeft
        ! add the ghost particles to the temporary particles array
        tmp_pLeft(nJumpedLeft + 1 : numSendLeft) = p(idx(1 : nGhostLeft))

        idx = 0
        nGhostDL = count(p(active)%ghost(5))
        idx(1 : nGhostDL) = pack(active, p(active)%ghost(5))
        numSendDL = nJumpedDL + nGhostDL
        tmp_pDL(nJumpedDL + 1 : numSendDL) = p(idx(1 : nGhostDL))

        idx = 0
        nGhostUL = count(p(active)%ghost(7))
        idx(1 : nGhostUL) = pack(active, p(active)%ghost(7))
        numSendUL = nJumpedUL + nGhostUL
        tmp_pUL(nJumpedUL + 1 : numSendUL) = p(idx(1 : nGhostUL))

        idx = 0
        nGhostAboveL = count(p(active)%ghost(10))
        idx(1 : nGhostAboveL) = pack(active, p(active)%ghost(10))
        numSendAboveL = nJumpedAboveL + nGhostAboveL
        tmp_pAboveL(nJumpedAboveL + 1 : numSendAboveL) = p(idx(1 : nGhostAboveL))

        idx = 0
        nGhostAboveUL = count(p(active)%ghost(16))
        idx(1 : nGhostAboveUL) = pack(active, p(active)%ghost(16))
        numSendAboveUL = nJumpedAboveUL + nGhostAboveUL
        tmp_pAboveUL(nJumpedAboveUL + 1 : numSendAboveUL) = p(idx(1 : nGhostAboveUL))

        idx = 0
        nGhostAboveDL = count(p(active)%ghost(14))
        idx(1 : nGhostAboveDL) = pack(active, p(active)%ghost(14))
        numSendAboveDL = nJumpedAboveDL + nGhostAboveDL
        tmp_pAboveDL(nJumpedAboveDL + 1 : numSendAboveDL) = p(idx(1 : nGhostAboveDL))

        idx = 0
        nGhostBelowL = count(p(active)%ghost(19))
        idx(1 : nGhostBelowL) = pack(active, p(active)%ghost(19))
        numSendBelowL = nJumpedBelowL + nGhostBelowL
        tmp_pBelowL(nJumpedBelowL + 1 : numSendBelowL) = p(idx(1 : nGhostBelowL))

        idx = 0
        nGhostBelowUL = count(p(active)%ghost(25))
        idx(1 : nGhostBelowUL) = pack(active, p(active)%ghost(25))
        numSendBelowUL = nJumpedBelowUL + nGhostBelowUL
        tmp_pBelowUL(nJumpedBelowUL + 1 : numSendBelowUL) = p(idx(1 : nGhostBelowUL))

        idx = 0
        nGhostBelowDL = count(p(active)%ghost(23))
        idx(1 : nGhostBelowDL) = pack(active, p(active)%ghost(23))
        numSendBelowDL = nJumpedBelowDL + nGhostBelowDL
        tmp_pBelowDL(nJumpedBelowDL + 1 : numSendBelowDL) = p(idx(1 : nGhostBelowDL))


        ! turn off the ghost particle indicator for the particles that are staying
        p%ghost(1) = .false.
        p%ghost(5) = .false.
        p%ghost(7) = .false.
        p%ghost(10) = .false.
        p%ghost(14) = .false.
        p%ghost(16) = .false.
        p%ghost(19) = .false.
        p%ghost(23) = .false.
        p%ghost(25) = .false.

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

        where (p(active)%jumped(10)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(14)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(16)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(19)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(23)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(25)) 
            p(active)%active = .false.
        endwhere

        ! send the number of particles to be sent to the left
        call mpi_send(numSendLeft, 1, mpi_integer, my_rank - nny, tag + my_rank - nny,&
                      mpi_comm_world, ierror)
    endif

    allocate(recv_pRight(Nactive))

    if (mod(my_rank, xysheet) < xysheet - nny) then
        ! receive the number of particles to be received from the right
        call mpi_recv(recvSizeRight, 1, mpi_integer, my_rank + nny, tag + my_rank, &
                      mpi_comm_world, status, ierror)
    endif

    ! if there are > 0 particles to be sent, send them
    if (mod(my_rank, xysheet) > nny - 1 .and. numSendLeft > 0) then
        call mpi_send(tmp_pLeft(1 : numSendLeft), numSendLeft, mpi_pType, my_rank - nny,&
                      tag + my_rank - nny, mpi_comm_world, ierror)
    endif
    ! if there are > 0 particles to be received, receive them
    if (mod(my_rank, xysheet) < (xysheet - nny) .and. recvSizeRight > 0) then
        call mpi_recv(recv_pRight(1 : recvSizeRight), recvSizeRight, mpi_pType,&
                      my_rank + nny, tag + my_rank, mpi_comm_world, status, ierror)
    endif

    deallocate(tmp_pLeft)


    ! now send particles to the right (this works the same as above, and thus is
    ! less commented)

    allocate(tmp_pRight(Nactive),tmp_pDR(Nactive),tmp_pUR(Nactive),tmp_pAboveR(Nactive),&
                    tmp_pAboveUR(Nactive),tmp_pAboveDR(Nactive),tmp_pBelowR(Nactive),&
                    tmp_pBelowUR(Nactive),tmp_pBelowDR(Nactive))

    if (mod(my_rank, xysheet) < xysheet - nny) then

        where (p(active)%loc(1) > Dlims(2))
            p(active)%jumped(2) = .true.
        endwhere

        where ((p(active)%jumped(2) .eqv. .true.) .and. p(active)%loc(2) > Dlims(4))
            p(active)%jumped(2) = .false.
            p(active)%jumped(8) = .true.
        endwhere

        where ((p(active)%jumped(8) .eqv. .true.) .and. p(active)%loc(3) > Dlims(6))
            p(active)%jumped(8) = .false.
            p(active)%jumped(17) = .true.
        endwhere

        where ((p(active)%jumped(8) .eqv. .true.) .and. p(active)%loc(3) < Dlims(5))
            p(active)%jumped(8) = .false.
            p(active)%jumped(26) = .true.
        endwhere

        where ((p(active)%jumped(2) .eqv. .true.) .and. p(active)%loc(2) < Dlims(3))
            p(active)%jumped(2) = .false.
            p(active)%jumped(6) = .true.
        endwhere

        where ((p(active)%jumped(6) .eqv. .true.) .and. p(active)%loc(3) > Dlims(6))
            p(active)%jumped(6) = .false.
            p(active)%jumped(15) = .true.
        endwhere

        where ((p(active)%jumped(6) .eqv. .true.) .and. p(active)%loc(3) < Dlims(5))
            p(active)%jumped(6) = .false.
            p(active)%jumped(24) = .true.
        endwhere

        where ((p(active)%jumped(2) .eqv. .true.) .and. p(active)%loc(3) > Dlims(6))
            p(active)%jumped(2) = .false.
            p(active)%jumped(11) = .true.
        endwhere

        where ((p(active)%jumped(2) .eqv. .true.) .and. p(active)%loc(3) < Dlims(5))
            p(active)%jumped(2) = .false.
            p(active)%jumped(20) = .true.
        endwhere

        idx = 0
        nJumpedRight = count(p(active)%jumped(2))
        idx(1 : nJumpedRight) = pack(active, p(active)%jumped(2))
        ! put it in a temporary array for sending
        tmp_pRight(1 : nJumpedRight) = p(idx(1 : nJumpedRight))

        idx = 0
        nJumpedDR = count(p(active)%jumped(6))
        idx(1 : nJumpedDR) = pack(active, p(active)%jumped(6))
        ! put it in a temporary array for sending
        tmp_pDR(1 : nJumpedDR) = p(idx(1 : nJumpedDR))

        idx = 0
        nJumpedUR = count(p(active)%jumped(8))
        idx(1 : nJumpedUR) = pack(active, p(active)%jumped(8))
        ! put it in a temporary array for sending
        tmp_pUR(1 : nJumpedUR) = p(idx(1 : nJumpedUR))

        idx = 0
        nJumpedAboveR = count(p(active)%jumped(11))
        idx(1 : nJumpedAboveR) = pack(active, p(active)%jumped(11))
        tmp_pAboveR(1 : nJumpedAboveR) = p(idx(1 : nJumpedAboveR))

        idx = 0
        nJumpedAboveUR = count(p(active)%jumped(17))
        idx(1 : nJumpedAboveUR) = pack(active, p(active)%jumped(17))
        tmp_pAboveUR(1 : nJumpedAboveUR) = p(idx(1 : nJumpedAboveUR))

        idx = 0
        nJumpedAboveDR = count(p(active)%jumped(15))
        idx(1 : nJumpedAboveDR) = pack(active, p(active)%jumped(15))
        tmp_pAboveDR(1 : nJumpedAboveDR) = p(idx(1 : nJumpedAboveDR))

        idx = 0
        nJumpedBelowR = count(p(active)%jumped(20))
        idx(1 : nJumpedBelowR) = pack(active, p(active)%jumped(20))
        tmp_pBelowR(1 : nJumpedBelowR) = p(idx(1 : nJumpedBelowR))

        idx = 0
        nJumpedBelowUR = count(p(active)%jumped(26))
        idx(1 : nJumpedBelowUR) = pack(active, p(active)%jumped(26))
        tmp_pBelowUR(1 : nJumpedBelowUR) = p(idx(1 : nJumpedBelowUR))

        idx = 0
        nJumpedBelowDR = count(p(active)%jumped(24))
        idx(1 : nJumpedBelowDR) = pack(active, p(active)%jumped(24))
        tmp_pBelowDR(1 : nJumpedBelowDR) = p(idx(1 : nJumpedBelowDR))

        ! tag the particles that will be sent as ghost particles
        where (p(active)%loc(1) <= Dlims(2) .and. p(active)%loc(1) > (Dlims(2) - paddist))
            p(active)%ghost(2) = .true.
        endwhere

        where((p(active)%ghost(2) .eqv. .true.) .and. p(active)%loc(3) > (Dlims(6) - paddist) .and. &
                p(active)%loc(3) <= Dlims(6))
            p(active)%ghost(11) = .true.
        endwhere

        where((p(active)%ghost(2) .eqv. .true.) .and. p(active)%loc(3) < (Dlims(5) + paddist) .and. &
                p(active)%loc(3) >= Dlims(5))
            p(active)%ghost(20) = .true.
        endwhere

        where((p(active)%ghost(2) .eqv. .true.) .and. p(active)%loc(2) > (Dlims(4) - paddist) .and. &
                p(active)%loc(2) <= Dlims(4))
            p(active)%ghost(8) = .true.
        endwhere

        where((p(active)%ghost(8) .eqv. .true.) .and. p(active)%loc(3) > (Dlims(6) - paddist) .and. &
                p(active)%loc(3) <= Dlims(6))
            p(active)%ghost(17) = .true.
        endwhere

        where((p(active)%ghost(8) .eqv. .true.) .and. p(active)%loc(3) < (Dlims(5) + paddist) .and. &
                p(active)%loc(3) >= Dlims(5))
            p(active)%ghost(26) = .true.
        endwhere

        where((p(active)%ghost(2) .eqv. .true.) .and. p(active)%loc(2) < (Dlims(3) + paddist) .and. &
                p(active)%loc(2) >= Dlims(3))
            p(active)%ghost(6) = .true.
        endwhere

        where((p(active)%ghost(6) .eqv. .true.) .and. p(active)%loc(3) > (Dlims(6) - paddist) .and. &
                p(active)%loc(3) <= Dlims(6))
            p(active)%ghost(15) = .true.
        endwhere

        where((p(active)%ghost(6) .eqv. .true.) .and. p(active)%loc(3) < (Dlims(5) + paddist) .and. &
                p(active)%loc(3) >= Dlims(5))
            p(active)%ghost(24) = .true.
        endwhere

        idx = 0
        nGhostRight = count(p(active)%ghost(2))
        idx(1 : nGhostRight) = pack(active, p(active)%ghost(2))
        numSendRight = nJumpedRight + nGhostRight
        tmp_pRight(nJumpedRight + 1 : numSendRight) = p(idx(1 : nGhostRight))

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

        idx = 0
        nGhostAboveR = count(p(active)%ghost(11))
        idx(1 : nGhostAboveR) = pack(active, p(active)%ghost(11))
        numSendAboveR = nJumpedAboveR + nGhostAboveR
        tmp_pAboveR(nJumpedAboveR + 1 : numSendAboveR) = p(idx(1 : nGhostAboveR))

        idx = 0
        nGhostAboveUR = count(p(active)%ghost(17))
        idx(1 : nGhostAboveUR) = pack(active, p(active)%ghost(17))
        numSendAboveUR = nJumpedAboveUR + nGhostAboveUR
        tmp_pAboveUR(nJumpedAboveUR + 1 : numSendAboveUR) = p(idx(1 : nGhostAboveUR))

        idx = 0
        nGhostAboveDR = count(p(active)%ghost(15))
        idx(1 : nGhostAboveDR) = pack(active, p(active)%ghost(15))
        numSendAboveDR = nJumpedAboveDR + nGhostAboveDR
        tmp_pAboveDR(nJumpedAboveDR + 1 : numSendAboveDR) = p(idx(1 : nGhostAboveDR))

        idx = 0
        nGhostBelowR = count(p(active)%ghost(20))
        idx(1 : nGhostBelowR) = pack(active, p(active)%ghost(20))
        numSendBelowR = nJumpedBelowR + nGhostBelowR
        tmp_pBelowR(nJumpedBelowR + 1 : numSendBelowR) = p(idx(1 : nGhostBelowR))

        idx = 0
        nGhostBelowUR = count(p(active)%ghost(26))
        idx(1 : nGhostBelowUR) = pack(active, p(active)%ghost(26))
        numSendBelowUR = nJumpedBelowUR + nGhostBelowUR
        tmp_pBelowUR(nJumpedBelowUR + 1 : numSendBelowUR) = p(idx(1 : nGhostBelowUR))

        idx = 0
        nGhostBelowDR = count(p(active)%ghost(24))
        idx(1 : nGhostBelowDR) = pack(active, p(active)%ghost(24))
        numSendBelowDR = nJumpedBelowDR + nGhostBelowDR
        tmp_pBelowDR(nJumpedBelowDR + 1 : numSendBelowDR) = p(idx(1 : nGhostBelowDR))


        p%ghost(2) = .false.
        p%ghost(6) = .false.
        p%ghost(8) = .false.
        p%ghost(11) = .false.
        p%ghost(15) = .false.
        p%ghost(17) = .false.
        p%ghost(20) = .false.
        p%ghost(24) = .false.
        p%ghost(26) = .false.

        where (p(active)%jumped(2)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(6)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(8)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(11)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(15)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(17)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(20)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(24)) 
            p(active)%active = .false.
        endwhere

        where (p(active)%jumped(26)) 
            p(active)%active = .false.
        endwhere

        call mpi_send(numSendRight, 1, mpi_integer, my_rank + nny, tag + my_rank + nny,&
                      mpi_comm_world, ierror)
    endif

    allocate(recv_pLeft(3*Nactive))

    if (mod(my_rank, xysheet) > nny - 1) then
        call mpi_recv(recvSizeLeft, 1, mpi_integer, my_rank - nny, tag + my_rank, &
                      mpi_comm_world, status, ierror)
    endif

    if (mod(my_rank, xysheet) < (xysheet - nny) .and. numSendRight > 0) then
        call mpi_send(tmp_pRight(1 : numSendRight), numSendRight, mpi_pType, my_rank + nny,&
                      tag + my_rank + nny, mpi_comm_world, ierror)
    endif
    if (mod(my_rank, xysheet) > nny - 1 .and. recvSizeLeft > 0) then
        call mpi_recv(recv_pLeft(1 : recvSizeLeft), recvSizeLeft, mpi_pType,&
                      my_rank - nny, tag + my_rank, mpi_comm_world, status, ierror)
    endif

    deallocate(tmp_pRight)


    ! ! now send particles to up/down

        ! ! first we send the particles down

        allocate(tmp_pDown(Nactive),tmp_pAboveD(Nactive),tmp_pBelowD(Nactive))

        if (mod(my_rank,nny) > 0) then

            where (p(active)%loc(2) < Dlims(3) .and. p(active)%loc(1) > Dlims(1) .and. p(active)%loc(1) < Dlims(2))
                p(active)%jumped(3) = .true.
            endwhere

            where ((p(active)%jumped(3) .eqv. .true.) .and. p(active)%loc(3) > Dlims(6))
                p(active)%jumped(3) = .false.
                p(active)%jumped(12) = .true.
            endwhere

            where ((p(active)%jumped(3) .eqv. .true.) .and. p(active)%loc(3) < Dlims(5))
                p(active)%jumped(3) = .false.
                p(active)%jumped(21) = .true.
            endwhere

            idx = 0
            nJumpedDown = count(p(active)%jumped(3))
            idx(1 : nJumpedDown) = pack(active, p(active)%jumped(3))
            tmp_pDown(1 : nJumpedDown) = p(idx(1 : nJumpedDown))

            idx = 0
            nJumpedAboveD = count(p(active)%jumped(12))
            idx(1 : nJumpedAboveD) = pack(active, p(active)%jumped(12))
            tmp_pAboveD(1 : nJumpedAboveD) = p(idx(1 : nJumpedAboveD))

            idx = 0
            nJumpedBelowD = count(p(active)%jumped(21))
            idx(1 : nJumpedBelowD) = pack(active, p(active)%jumped(21))
            tmp_pBelowD(1 : nJumpedBelowD) = p(idx(1 : nJumpedBelowD))

            where (p(active)%loc(2) >= Dlims(3) .and. p(active)%loc(2) < (Dlims(3) + paddist))
                p(active)%ghost(3) = .true.
            endwhere

            where((p(active)%ghost(3) .eqv. .true.) .and. p(active)%loc(3) > (Dlims(6) - paddist) .and. &
                    p(active)%loc(3) <= Dlims(6))
                p(active)%ghost(12) = .true.
            endwhere

            where((p(active)%ghost(3) .eqv. .true.) .and. p(active)%loc(3) < (Dlims(5) + paddist) .and. &
                    p(active)%loc(3) >= Dlims(5))
                p(active)%ghost(21) = .true.
             endwhere

            idx = 0
            nGhostDown = count(p(active)%ghost(3))
            idx(1 : nGhostDown) = pack(active, p(active)%ghost(3))
            numSendDown = nJumpedDown + nGhostDown
            tmp_pDown(nJumpedDown + 1 : numSendDown) = p(idx(1 : nGhostDown))

            idx = 0
            nGhostAboveD = count(p(active)%ghost(12))
            idx(1 : nGhostAboveD) = pack(active, p(active)%ghost(12))
            numSendAboveD = nJumpedAboveD + nGhostAboveD
            tmp_pAboveD(nJumpedAboveD + 1 : numSendAboveD) = p(idx(1 : nGhostAboveD))

            idx = 0
            nGhostBelowD = count(p(active)%ghost(21))
            idx(1 : nGhostBelowD) = pack(active, p(active)%ghost(21))
            numSendBelowD = nJumpedBelowD + nGhostBelowD
            tmp_pBelowD(nJumpedBelowD + 1 : numSendBelowD) = p(idx(1 : nGhostBelowD))


            p%ghost(3) = .false.
            p%ghost(12) = .false.
            p%ghost(21) = .false.

            where (p(active)%jumped(3)) 
                p(active)%active = .false.
            endwhere

            where (p(active)%jumped(12)) 
                p(active)%active = .false.
            endwhere

            where (p(active)%jumped(21)) 
                p(active)%active = .false.
            endwhere

            call mpi_send(numSendDown, 1, mpi_integer, my_rank - 1, tag + my_rank - 1,&
                          mpi_comm_world, ierror)

            ! send diagonal particles if necessary

            ! down and to the right
            if (mod(my_rank, xysheet) < (xysheet - nny)) then
                call mpi_send(numSendDR, 1, mpi_integer, my_rank + nny - 1, tag + my_rank + nny - 1,&
                          mpi_comm_world, ierror)
            endif

            ! down and to the left
            if (mod(my_rank, xysheet) > (nny - 1)) then
                call mpi_send(numSendDL, 1, mpi_integer, my_rank - nny - 1, tag + my_rank - nny - 1,&
                          mpi_comm_world, ierror)
            endif


        endif

        allocate(recv_pUp(Nactive),recv_pUL(Nactive),recv_pUR(Nactive))

        ! recv number from above
        if (mod(my_rank,nny) /= (nny - 1)) then
            call mpi_recv(recvSizeUp, 1, mpi_integer, my_rank + 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif


        !recv number from UR
        if (mod(my_rank,nny) /= (nny - 1) .and. mod(my_rank, xysheet) < (xysheet - nny)) then
            call mpi_recv(recvSizeUR, 1, mpi_integer, my_rank + nny + 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif

        !recv number from UL
        if (mod(my_rank,nny) /= (nny - 1) .and. (mod(my_rank, xysheet) > (nny - 1))) then
            call mpi_recv(recvSizeUL, 1, mpi_integer, my_rank - nny + 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif

        !send array down
        if (mod(my_rank,nny) > 0 .and. numSendDown > 0) then
            call mpi_send(tmp_pDown(1 : numSendDown), numSendDown, mpi_pType, my_rank - 1,&
                          tag + my_rank - 1, mpi_comm_world, ierror)
        endif

        !recv array from above
        if (mod(my_rank,nny) /= (nny - 1) .and. recvSizeUp > 0) then
            call mpi_recv(recv_pUp(1 : recvSizeUp), recvSizeUp, mpi_pType,&
                          my_rank + 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        !send array DR
        if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) < (xysheet - nny) .and. numSendDR > 0) then
            call mpi_send(tmp_pDR(1 : numSendDR), numSendDR, mpi_pType, my_rank + nny - 1, tag + my_rank + nny - 1,&
                          mpi_comm_world, ierror)
        endif 

        !recv array from UL
        if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) > (nny - 1) .and. recvSizeUL > 0) then
            call mpi_recv(recv_pUL(1 : recvSizeUL), recvSizeUL, mpi_pType,&
                          my_rank - nny + 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        !send array DL
        if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) > (nny - 1) .and. numSendDL > 0) then
            call mpi_send(tmp_pDL(1 : numSendDL), numSendDL, mpi_pType, my_rank - nny - 1, tag + my_rank - nny - 1,&
                          mpi_comm_world, ierror)
        endif

        !recv array UR
        if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) < (xysheet - nny) .and. recvSizeUR > 0) then
            call mpi_recv(recv_pUR(1 : recvSizeUR), recvSizeUR, mpi_pType,&
                          my_rank + nny + 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        deallocate(tmp_pDL,tmp_pDR,tmp_pDown)


        ! ! now send particles up

        allocate(tmp_pUp(Nactive),tmp_pAboveU(Nactive),tmp_pBelowU(Nactive))


        if (mod(my_rank,nny) /= (nny - 1)) then

            where (p(active)%loc(2) > Dlims(4) .and. p(active)%loc(1) > Dlims(1) .and. p(active)%loc(1) < Dlims(2))
                p(active)%jumped(4) = .true.
            endwhere

            where ((p(active)%jumped(4) .eqv. .true.) .and. p(active)%loc(3) > Dlims(6))
                p(active)%jumped(4) = .false.
                p(active)%jumped(13) = .true.
            endwhere

            where ((p(active)%jumped(4) .eqv. .true.) .and. p(active)%loc(3) < Dlims(5))
                p(active)%jumped(4) = .false.
                p(active)%jumped(22) = .true.
            endwhere

            idx = 0
            nJumpedUp = count(p(active)%jumped(4))
            idx(1 : nJumpedUp) = pack(active, p(active)%jumped(4))
            tmp_pUp(1 : nJumpedUp) = p(idx(1 : nJumpedUp))

            idx = 0
            nJumpedAboveU = count(p(active)%jumped(13))
            idx(1 : nJumpedAboveU) = pack(active, p(active)%jumped(13))
            tmp_pAboveU(1 : nJumpedAboveU) = p(idx(1 : nJumpedAboveU))

            idx = 0
            nJumpedBelowU = count(p(active)%jumped(22))
            idx(1 : nJumpedBelowU) = pack(active, p(active)%jumped(22))
            tmp_pBelowU(1 : nJumpedBelowU) = p(idx(1 : nJumpedBelowU))


            where (p(active)%loc(2) <= Dlims(4) .and. p(active)%loc(2) > (Dlims(4) - paddist))
                p(active)%ghost(4) = .true.
            endwhere

            where((p(active)%ghost(4) .eqv. .true.) .and. p(active)%loc(3) > (Dlims(6) - paddist) .and. &
                    p(active)%loc(3) <= Dlims(6))
                p(active)%ghost(13) = .true.
            endwhere

            where((p(active)%ghost(4) .eqv. .true.) .and. p(active)%loc(3) < (Dlims(5) + paddist) .and. &
                    p(active)%loc(3) >= Dlims(5))
                p(active)%ghost(22) = .true.
             endwhere

            idx = 0
            nGhostUp = count(p(active)%ghost(4))
            idx(1 : nGhostUp) = pack(active, p(active)%ghost(4))
            numSendUp = nJumpedUp + nGhostUp
            tmp_pUp(nJumpedUp + 1 : numSendUp) = p(idx(1 : nGhostUp))

            idx = 0
            nGhostAboveU = count(p(active)%ghost(13))
            idx(1 : nGhostAboveD) = pack(active, p(active)%ghost(13))
            numSendAboveU = nJumpedAboveU + nGhostAboveU
            tmp_pAboveU(nJumpedAboveU + 1 : numSendAboveU) = p(idx(1 : nGhostAboveU))

            idx = 0
            nGhostBelowU = count(p(active)%ghost(22))
            idx(1 : nGhostBelowU) = pack(active, p(active)%ghost(22))
            numSendBelowU = nJumpedBelowU + nGhostBelowU
            tmp_pBelowU(nJumpedBelowU + 1 : numSendBelowU) = p(idx(1 : nGhostBelowU))

            p%ghost(4) = .false.
            p%ghost(13) = .false.
            p%ghost(22) = .false.

            where (p(active)%jumped(4)) 
                p(active)%active = .false.
            endwhere

            where (p(active)%jumped(13)) 
                p(active)%active = .false.
            endwhere

            where (p(active)%jumped(22)) 
                p(active)%active = .false.
            endwhere

            call mpi_send(numSendUp, 1, mpi_integer, my_rank + 1, tag + my_rank + 1,&
                          mpi_comm_world, ierror)

            !send size UR
            if (mod(my_rank, xysheet) < (xysheet - nny)) then
                call mpi_send(numSendUR, 1, mpi_integer, my_rank + nny + 1, tag + my_rank + nny + 1,&
                          mpi_comm_world, ierror)
            endif

            !send size UL
            if (mod(my_rank, xysheet) > (nny - 1)) then
                call mpi_send(numSendUL, 1, mpi_integer, my_rank - nny + 1, tag + my_rank - nny + 1,&
                          mpi_comm_world, ierror)
            endif
        endif

        allocate(recv_pDown(Nactive),recv_pDL(Nactive),recv_pDR(Nactive))

        !recv number from below
        if (mod(my_rank,nny) /= 0) then
            call mpi_recv(recvSizeDown, 1, mpi_integer, my_rank - 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif

        !send array up
        if (mod(my_rank,nny) /= nny-1 .and. numSendUp > 0) then
            call mpi_send(tmp_pUp(1 : numSendUp), numSendUp, mpi_pType, my_rank + 1,&
                          tag + my_rank + 1, mpi_comm_world, ierror)
        endif

        !recv array from below
        if (mod(my_rank,nny) /= 0 .and. recvSizeDown > 0) then
            call mpi_recv(recv_pDown(1 : recvSizeDown), recvSizeDown, mpi_pType,&
                          my_rank - 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        !recv size DL
        if (mod(my_rank,nny) /= 0 .and. mod(my_rank, xysheet) > nny - 1) then
            call mpi_recv(recvSizeDL, 1, mpi_integer, my_rank - nny - 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif

        !recv size DR
        if (mod(my_rank,nny) /= 0 .and. mod(my_rank, xysheet) < (xysheet - nny)) then
            call mpi_recv(recvSizeDR, 1, mpi_integer, my_rank + nny - 1, tag + my_rank, &
                          mpi_comm_world, status, ierror)
        endif

        !send array UR
        if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) < (xysheet - nny) .and. numSendUR > 0) then
            call mpi_send(tmp_pUR(1 : numSendUR), numSendUR, mpi_pType, my_rank + nny + 1,&
                          tag + my_rank + nny + 1, mpi_comm_world, ierror)
        endif

        !recv array DL
        if (mod(my_rank,nny) /= 0 .and. mod(my_rank, xysheet) > (nny - 1) .and. recvSizeDL > 0) then
            call mpi_recv(recv_pDL(1 : recvSizeDL), recvSizeDL, mpi_pType,&
                          my_rank - nny - 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        !send array UL
        if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) > (nny - 1) .and. numSendUL > 0) then
            call mpi_send(tmp_pUL(1 : numSendUL), numSendUL, mpi_pType, my_rank - nny + 1,&
                          tag + my_rank - nny + 1, mpi_comm_world, ierror)
        endif

        !recv array DR
        if (mod(my_rank,nny) /= 0 .and. mod(my_rank, xysheet) < (xysheet - nny) .and. recvSizeDR > 0) then
            call mpi_recv(recv_pDR(1 : recvSizeDR), recvSizeDR, mpi_pType,&
                          my_rank + nny - 1, tag + my_rank, mpi_comm_world, status, ierror)
        endif

        deallocate(tmp_pUp,tmp_pUR,tmp_pUL)

        !now have to send to above planes

        allocate(tmp_pAbove(Nactive))

        if (my_rank < num_cores - xysheet) then

            where ((p(active)%loc(1) >= Dlims(1)) .and. (p(active)%loc(1) <= Dlims(2)) .and. (p(active)%loc(2) >= Dlims(3)) &
                    .and. (p(active)%loc(2) <= Dlims(4)) .and. (p(active)%loc(3) > Dlims(6)))
                p(active)%jumped(9) = .true.
            endwhere


            idx = 0
            nJumpedAbove = count(p(active)%jumped(9))
            idx(1 : nJumpedAbove) = pack(active, p(active)%jumped(9))
            tmp_pAbove(1 : nJumpedAbove) = p(idx(1 : nJumpedAbove))

            where ((p(active)%loc(3) > (Dlims(6) - paddist)) .and. (p(active)%loc(3) <= Dlims(6)))
                p(active)%ghost(9) = .true.
            endwhere

            idx = 0
            nGhostAbove = count(p(active)%ghost(9))
            idx(1 : nGhostAbove) = pack(active, p(active)%ghost(9))
            numSendAbove = nJumpedAbove + nGhostAbove
            tmp_pAbove(nJumpedAbove + 1 : numSendAbove) = p(idx(1 : nGhostAbove))

            p%ghost(9) = .false.

            where (p(active)%jumped(9)) 
                p(active)%active = .false.
            endwhere

            !send size directly above (z-direction)
            call mpi_send(numSendAbove, 1, mpi_integer, my_rank + xysheet, tag + my_rank + xysheet,&
                          mpi_comm_world, ierror)

            !send size above (z) and right
            if (mod(my_rank, xysheet) < xysheet - nny) then
                call mpi_send(numSendAboveR, 1, mpi_integer, my_rank + xysheet + nny, tag + my_rank + xysheet + nny,&
                          mpi_comm_world, ierror)
            endif

            !send size above (z) and UR
            if (mod(my_rank, xysheet) < (xysheet - nny) .and. mod(my_rank, nny) /= (nny - 1)) then
                call mpi_send(numSendAboveUR, 1, mpi_integer, my_rank + xysheet + nny + 1, tag + my_rank + xysheet + nny + 1,&
                          mpi_comm_world, ierror)
            endif

            !send size above (z) and DR
            if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) < (xysheet - nny)) then
                    call mpi_send(numSendAboveDR, 1, mpi_integer, my_rank + xysheet + nny - 1, tag + my_rank + xysheet + nny - 1,&
                          mpi_comm_world, ierror)
            endif

            !send size above (z) and D
            if (mod(my_rank,nny) > 0) then
                call mpi_send(numSendAboveD, 1, mpi_integer, my_rank + xysheet - 1, tag + my_rank + xysheet - 1,&
                          mpi_comm_world, ierror)
            endif

            !send size above (z) and DL
            if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) > (nny - 1)) then
                call mpi_send(numSendAboveDL, 1, mpi_integer, my_rank + xysheet - nny - 1, tag + my_rank + xysheet - nny - 1,&
                          mpi_comm_world, ierror)
            endif

            !send size above (z) and L
            if (mod(my_rank,xysheet) > (nny - 1)) then
                call mpi_send(numSendAboveL, 1, mpi_integer, my_rank + xysheet - nny, tag + my_rank + xysheet - nny,&
                          mpi_comm_world, ierror)
            endif

            !send size above (z) and UL
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) > (nny - 1)) then
                call mpi_send(numSendAboveUL, 1, mpi_integer, my_rank + xysheet - nny + 1, tag + my_rank + xysheet - nny + 1,&
                          mpi_comm_world, ierror)
            endif

            !send size above (z) and U
            if (mod(my_rank,nny) /= (nny - 1)) then
                call mpi_send(numSendAboveU, 1, mpi_integer, my_rank + xysheet + 1, tag + my_rank + xysheet + 1,&
                          mpi_comm_world, ierror)
            endif

        endif


        !recv sizes from below
        if (my_rank >= xysheet) then

            !recv number from below (z)
            call mpi_recv(recvSizeBelow, 1, mpi_integer, my_rank - xysheet, tag + my_rank, &
                    mpi_comm_world, status, ierror)

            !recv size from below (z) and L
            if (mod(my_rank, xysheet) > (nny - 1)) then
                call mpi_recv(recvSizeBelowL, 1, mpi_integer, my_rank - xysheet - nny, tag + my_rank, &
                             mpi_comm_world, status, ierror)
            endif

            !recv size from below (z) and DL
            if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) > nny - 1) then
                call mpi_recv(recvSizeBelowDL, 1, mpi_integer, my_rank - xysheet - nny - 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from below (z) and UL
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) > (nny - 1)) then
                call mpi_recv(recvSizeBelowUL, 1, mpi_integer, my_rank - xysheet - nny + 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from below (z) and U
            if (mod(my_rank,nny) /= (nny - 1)) then
                call mpi_recv(recvSizeBelowU, 1, mpi_integer, my_rank - xysheet + 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from below (z) and UR
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) < (xysheet - nny)) then
                call mpi_recv(recvSizeBelowUR, 1, mpi_integer, my_rank - xysheet + nny + 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from below (z) and R
            if (mod(my_rank, xysheet) < xysheet - nny) then
                call mpi_recv(recvSizeBelowR, 1, mpi_integer, my_rank - xysheet + nny, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from below (z) and DR
            if (mod(my_rank,nny) /= 0 .and. mod(my_rank, xysheet) < (xysheet - nny)) then
                call mpi_recv(recvSizeBelowDR, 1, mpi_integer, my_rank - xysheet + nny - 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from below (z) and D
            if (mod(my_rank,nny) /= 0) then
                call mpi_recv(recvSizeBelowD, 1, mpi_integer, my_rank - xysheet - 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

        endif


        !send arrays above (z)
        if (my_rank < num_cores - xysheet) then

            !send array above (z)
            if (numSendAbove > 0) then
                call mpi_send(tmp_pAbove(1 : numSendAbove), numSendAbove, mpi_pType, my_rank + xysheet,&
                              tag + my_rank + xysheet, mpi_comm_world, ierror)
            endif

            !send array above (z) and right
            if (mod(my_rank, xysheet) < xysheet - nny .and. numSendAboveR > 0) then
                call mpi_send(tmp_pAboveR(1 : numSendAboveR), numSendAboveR, mpi_pType, my_rank + xysheet + nny,&
                              tag + my_rank + xysheet + nny, mpi_comm_world, ierror)
            endif

            !send array above (z) and UR
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) < (xysheet - nny) .and. numSendAboveUR > 0) then
                call mpi_send(tmp_pAboveUR(1 : numSendAboveUR), numSendAboveUR, mpi_pType, my_rank + xysheet + nny + 1,&
                            tag + my_rank + xysheet + nny + 1, mpi_comm_world, ierror)
            endif

            !send array above (z) and DR
            if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) < (xysheet - nny) .and. numSendAboveDR > 0) then
                call mpi_send(tmp_pAboveDR(1 : numSendAboveDR), numSendAboveDR, mpi_pType, my_rank + xysheet + nny - 1,&
                            tag + my_rank + xysheet + nny - 1, mpi_comm_world, ierror)
            endif

            !send array above (z) and D
            if (mod(my_rank,nny) > 0 .and. numSendAboveD > 0) then
                call mpi_send(tmp_pAboveD(1 : numSendAboveD), numSendAboveD, mpi_pType, my_rank + xysheet - 1,&
                            tag + my_rank + xysheet - 1, mpi_comm_world, ierror)
            endif

            !send array above (z) and DL
            if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) > (nny - 1) .and. numSendAboveDL > 0) then
                call mpi_send(tmp_pAboveDL(1 : numSendAboveDL), numSendAboveDL, mpi_pType, my_rank + xysheet - nny - 1,&
                            tag + my_rank + xysheet - nny - 1, mpi_comm_world, ierror)
            endif

            !send array above (z) and L
            if (mod(my_rank,xysheet) > (nny - 1) .and. numSendAboveL > 0) then
                call mpi_send(tmp_pAboveL(1 : numSendAboveL), numSendAboveL, mpi_pType, my_rank + xysheet - nny,&
                            tag + my_rank + xysheet - nny, mpi_comm_world, ierror)
            endif

            !send array above (z) and UL
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) > (nny - 1) .and. numSendAboveUL > 0) then
                call mpi_send(tmp_pAboveUL(1 : numSendAboveUL), numSendAboveUL, mpi_pType, my_rank + xysheet - nny + 1,&
                            tag + my_rank + xysheet - nny + 1, mpi_comm_world, ierror)
            endif

            !send array above (z) and U
            if (mod(my_rank,nny) /= (nny - 1) .and. numSendAboveU > 0) then
                call mpi_send(tmp_pAboveU(1 : numSendAboveU), numSendAboveU, mpi_pType, my_rank + xysheet + 1,&
                            tag + my_rank + xysheet + 1, mpi_comm_world, ierror)
            endif


        endif

        deallocate(tmp_pAbove,tmp_pAboveR,tmp_pAboveUR,tmp_pAboveDR,tmp_pAboveD,tmp_pAboveDL,&
                    tmp_pAboveL,tmp_pAboveUL,tmp_pAboveU)

        allocate(recv_pBelow(Nactive),recv_pBelowL(Nactive),recv_pBelowDL(Nactive),recv_pBelowUL(Nactive),&
                    recv_pBelowU(Nactive),recv_pBelowUR(Nactive),recv_pBelowR(Nactive),recv_pBelowDR(Nactive),&
                    recv_pBelowD(Nactive))


        !recv arrays from below (z)
        if (my_rank >= xysheet) then

            !recv array from below (z)
            if (recvSizeBelow > 0) then
                call mpi_recv(recv_pBelow(1 : recvSizeBelow), recvSizeBelow, mpi_pType,&
                              my_rank - xysheet, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from below (z) and left
            if (mod(my_rank, xysheet) > (nny - 1) .and. recvSizeBelowL > 0) then
                call mpi_recv(recv_pBelowL(1 : recvSizeBelowL), recvSizeBelowL, mpi_pType,&
                              my_rank - xysheet - nny, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from below (z) and DL
            if (mod(my_rank,nny) /= 0 .and. mod(my_rank, xysheet) > (nny - 1) .and. recvSizeBelowDL > 0) then
                call mpi_recv(recv_pBelowDL(1 : recvSizeBelowDL), recvSizeBelowDL, mpi_pType,&
                            my_rank - xysheet - nny - 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from below (z) and UL
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) > (nny - 1) .and. recvSizeBelowUL > 0) then
                call mpi_recv(recv_pBelowUL(1 : recvSizeBelowUL), recvSizeBelowUL, mpi_pType,&
                            my_rank - xysheet - nny + 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from below (z) and U
            if (mod(my_rank,nny) /= (nny - 1) .and. recvSizeBelowU > 0) then
                call mpi_recv(recv_pBelowU(1 : recvSizeBelowU), recvSizeBelowU, mpi_pType,&
                            my_rank - xysheet + 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from below (z) and UR
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) < (xysheet - nny) .and. recvSizeBelowUR > 0) then
                call mpi_recv(recv_pBelowUR(1 : recvSizeBelowUR), recvSizeBelowUR, mpi_pType,&
                            my_rank - xysheet + nny + 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from below (z) and R
            if (mod(my_rank, xysheet) < xysheet - nny .and. recvSizeBelowR > 0) then
                call mpi_recv(recv_pBelowR(1 : recvSizeBelowR), recvSizeBelowR, mpi_pType,&
                            my_rank - xysheet + nny, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from below (z) and DR
            if (mod(my_rank,nny) /= 0 .and. mod(my_rank, xysheet) < (xysheet - nny) .and. recvSizeBelowDR > 0) then
                call mpi_recv(recv_pBelowDR(1 : recvSizeBelowDR), recvSizeBelowDR, mpi_pType,&
                            my_rank - xysheet + nny - 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv size from below (z) and D
            if (mod(my_rank,nny) /= 0 .and. recvSizeBelowD > 0) then
                call mpi_recv(recv_pBelowD(1 : recvSizeBelowD), recvSizeBelowD, mpi_pType,&
                            my_rank - xysheet - 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif


        endif

        allocate(tmp_pBelow(Nactive))

        if (my_rank >= xysheet) then

            where ((p(active)%loc(1) >= Dlims(1)) .and. (p(active)%loc(1) <= Dlims(2)) .and. (p(active)%loc(2) >= Dlims(3)) &
                    .and. (p(active)%loc(2) <= Dlims(4)) .and. (p(active)%loc(3) < Dlims(5)))
                p(active)%jumped(18) = .true.
            endwhere

            idx = 0
            nJumpedBelow = count(p(active)%jumped(18))
            idx(1 : nJumpedBelow) = pack(active, p(active)%jumped(18))
            tmp_pBelow(1 : nJumpedBelow) = p(idx(1 : nJumpedBelow))

            where ((p(active)%loc(3) < (Dlims(5) + paddist)) .and. (p(active)%loc(3) >= Dlims(5)))
                p(active)%ghost(18) = .true.
            endwhere

            idx = 0
            nGhostBelow = count(p(active)%ghost(18))
            idx(1 : nGhostBelow) = pack(active, p(active)%ghost(18))
            numSendBelow = nJumpedBelow + nGhostBelow
            tmp_pBelow(nJumpedBelow + 1 : numSendBelow) = p(idx(1 : nGhostBelow))

            p%ghost(18) = .false.

            where (p(active)%jumped(18)) 
                p(active)%active = .false.
            endwhere

            !send size directly below (z-direction)
            call mpi_send(numSendBelow, 1, mpi_integer, my_rank - xysheet, tag + my_rank - xysheet,&
                          mpi_comm_world, ierror)

            !send size below (z) and right
            if (mod(my_rank, xysheet) < xysheet - nny) then
                call mpi_send(numSendBelowR, 1, mpi_integer, my_rank - xysheet + nny, tag + my_rank - xysheet + nny,&
                          mpi_comm_world, ierror)
            endif

            !send size below (z) and UR
            if (mod(my_rank, xysheet) < (xysheet - nny) .and. mod(my_rank, nny) /= (nny - 1)) then
                call mpi_send(numSendBelowUR, 1, mpi_integer, my_rank - xysheet + nny + 1, tag + my_rank - xysheet + nny + 1,&
                          mpi_comm_world, ierror)
            endif

            !send size below (z) and DR
            if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) < (xysheet - nny)) then
                    call mpi_send(numSendBelowDR, 1, mpi_integer, my_rank - xysheet + nny - 1, tag + my_rank - xysheet + nny - 1,&
                          mpi_comm_world, ierror)
            endif

            !send size below (z) and D
            if (mod(my_rank,nny) > 0) then
                call mpi_send(numSendBelowD, 1, mpi_integer, my_rank - xysheet - 1, tag + my_rank - xysheet - 1,&
                          mpi_comm_world, ierror)
            endif

            !send size below (z) and DL
            if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) > (nny - 1)) then
                call mpi_send(numSendBelowDL, 1, mpi_integer, my_rank - xysheet - nny - 1, tag + my_rank - xysheet - nny - 1,&
                          mpi_comm_world, ierror)
            endif

            !send size below (z) and L
            if (mod(my_rank,xysheet) > (nny - 1)) then
                call mpi_send(numSendBelowL, 1, mpi_integer, my_rank - xysheet - nny, tag + my_rank - xysheet - nny,&
                          mpi_comm_world, ierror)
            endif

            !send size below (z) and UL
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) > (nny - 1)) then
                call mpi_send(numSendBelowUL, 1, mpi_integer, my_rank - xysheet - nny + 1, tag + my_rank - xysheet - nny + 1,&
                          mpi_comm_world, ierror)
            endif

            !send size below (z) and U
            if (mod(my_rank,nny) /= (nny - 1)) then
                call mpi_send(numSendBelowU, 1, mpi_integer, my_rank - xysheet + 1, tag + my_rank - xysheet + 1,&
                          mpi_comm_world, ierror)
            endif

        endif


        !recv sizes from above
        if (my_rank < num_cores - xysheet) then

            !recv number from above (z)
            call mpi_recv(recvSizeAbove, 1, mpi_integer, my_rank + xysheet, tag + my_rank, &
                    mpi_comm_world, status, ierror)

            !recv size from above (z) and L
            if (mod(my_rank, xysheet) > (nny - 1)) then
                call mpi_recv(recvSizeAboveL, 1, mpi_integer, my_rank + xysheet - nny, tag + my_rank, &
                             mpi_comm_world, status, ierror)
            endif

            !recv size from above (z) and DL
            if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) > nny - 1) then
                call mpi_recv(recvSizeAboveDL, 1, mpi_integer, my_rank + xysheet - nny - 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from above (z) and UL
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) > (nny - 1)) then
                call mpi_recv(recvSizeAboveUL, 1, mpi_integer, my_rank + xysheet - nny + 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from above (z) and U
            if (mod(my_rank,nny) /= (nny - 1)) then
                call mpi_recv(recvSizeAboveU, 1, mpi_integer, my_rank + xysheet + 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from above (z) and UR
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) < (xysheet - nny)) then
                call mpi_recv(recvSizeAboveUR, 1, mpi_integer, my_rank + xysheet + nny + 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from above (z) and R
            if (mod(my_rank, xysheet) < xysheet - nny) then
                call mpi_recv(recvSizeAboveR, 1, mpi_integer, my_rank + xysheet + nny, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from above (z) and DR
            if (mod(my_rank,nny) /= 0 .and. mod(my_rank, xysheet) < (xysheet - nny)) then
                call mpi_recv(recvSizeAboveDR, 1, mpi_integer, my_rank + xysheet + nny - 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif

            !recv size from above (z) and D
            if (mod(my_rank,nny) /= 0) then
                call mpi_recv(recvSizeAboveD, 1, mpi_integer, my_rank + xysheet - 1, tag + my_rank, &
                              mpi_comm_world, status, ierror)
            endif


        endif


        !send arrays down (z)
        if (my_rank >= xysheet) then

            !send array below (z)
            if (numSendBelow > 0) then
                call mpi_send(tmp_pBelow(1 : numSendBelow), numSendBelow, mpi_pType, my_rank - xysheet,&
                              tag + my_rank - xysheet, mpi_comm_world, ierror)
            endif

            !send array below (z) and right
            if (mod(my_rank, xysheet) < xysheet - nny .and. numSendBelowR > 0) then
                call mpi_send(tmp_pBelowR(1 : numSendBelowR), numSendBelowR, mpi_pType, my_rank - xysheet + nny,&
                              tag + my_rank - xysheet + nny, mpi_comm_world, ierror)
            endif

            !send array below (z) and UR
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) < (xysheet - nny) .and. numSendBelowUR > 0) then
                call mpi_send(tmp_pBelowUR(1 : numSendBelowUR), numSendBelowUR, mpi_pType, my_rank - xysheet + nny + 1,&
                            tag + my_rank - xysheet + nny + 1, mpi_comm_world, ierror)
            endif

            !send array below (z) and DR
            if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) < (xysheet - nny) .and. numSendBelowDR > 0) then
                call mpi_send(tmp_pBelowDR(1 : numSendBelowDR), numSendBelowDR, mpi_pType, my_rank - xysheet + nny - 1,&
                            tag + my_rank - xysheet + nny - 1, mpi_comm_world, ierror)
            endif

            !send array below (z) and D
            if (mod(my_rank,nny) > 0 .and. numSendBelowD > 0) then
                call mpi_send(tmp_pBelowD(1 : numSendBelowD), numSendBelowD, mpi_pType, my_rank - xysheet - 1,&
                            tag + my_rank - xysheet - 1, mpi_comm_world, ierror)
            endif

            !send array below (z) and DL
            if (mod(my_rank,nny) > 0 .and. mod(my_rank, xysheet) > (nny - 1) .and. numSendBelowDL > 0) then
                call mpi_send(tmp_pBelowDL(1 : numSendBelowDL), numSendBelowDL, mpi_pType, my_rank - xysheet - nny - 1,&
                            tag + my_rank - xysheet - nny - 1, mpi_comm_world, ierror)
            endif

            !send array below (z) and L
            if (mod(my_rank,xysheet) > (nny - 1) .and. numSendBelowL > 0) then
                call mpi_send(tmp_pBelowL(1 : numSendBelowL), numSendBelowL, mpi_pType, my_rank - xysheet - nny,&
                            tag + my_rank - xysheet - nny, mpi_comm_world, ierror)
            endif

            !send size below (z) and UL
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) > (nny - 1) .and. numSendBelowUL > 0) then
                call mpi_send(tmp_pBelowUL(1 : numSendBelowUL), numSendBelowUL, mpi_pType, my_rank - xysheet - nny + 1,&
                            tag + my_rank - xysheet - nny + 1, mpi_comm_world, ierror)
            endif

            !send size below (z) and U
            if (mod(my_rank,nny) /= (nny - 1) .and. numSendBelowU > 0) then
                call mpi_send(tmp_pBelowU(1 : numSendBelowU), numSendBelowU, mpi_pType, my_rank - xysheet + 1,&
                            tag + my_rank - xysheet + 1, mpi_comm_world, ierror)
            endif

        endif

        deallocate(tmp_pBelow,tmp_pBelowR,tmp_pBelowUR,tmp_pBelowDR,tmp_pBelowD,tmp_pBelowDL,tmp_pBelowL,&
                    tmp_pBelowUL,tmp_pBelowU)

        allocate(recv_pAbove(Nactive),recv_pAboveL(Nactive),recv_pAboveDL(Nactive),recv_pAboveUL(Nactive),&
                    recv_pAboveU(Nactive),recv_pAboveUR(Nactive),recv_pAboveR(Nactive),recv_pAboveDR(Nactive),&
                    recv_pAboveD(Nactive))


        !recv arrays from above (z)
        if (my_rank < num_cores - xysheet) then

            !recv array from above (z)
            if (recvSizeAbove > 0) then
                call mpi_recv(recv_pAbove(1 : recvSizeAbove), recvSizeAbove, mpi_pType,&
                              my_rank + xysheet, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from above (z) and left
            if (mod(my_rank, xysheet) > (nny - 1) .and. recvSizeAboveL > 0) then
                call mpi_recv(recv_pAboveL(1 : recvSizeAboveL), recvSizeAboveL, mpi_pType,&
                              my_rank + xysheet - nny, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from above (z) and DL
            if (mod(my_rank,nny) /= 0 .and. mod(my_rank, xysheet) > (nny - 1) .and. recvSizeAboveDL > 0) then
                call mpi_recv(recv_pAboveDL(1 : recvSizeAboveDL), recvSizeAboveDL, mpi_pType,&
                            my_rank + xysheet - nny - 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from above (z) and UL
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) > (nny - 1) .and. recvSizeAboveUL > 0) then
                call mpi_recv(recv_pAboveUL(1 : recvSizeAboveUL), recvSizeAboveUL, mpi_pType,&
                            my_rank + xysheet - nny + 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from above (z) and U
            if (mod(my_rank,nny) /= (nny - 1) .and. recvSizeAboveU > 0) then
                call mpi_recv(recv_pAboveU(1 : recvSizeAboveU), recvSizeAboveU, mpi_pType,&
                            my_rank + xysheet + 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from above (z) and UR
            if (mod(my_rank,nny) /= (nny-1) .and. mod(my_rank, xysheet) < (xysheet - nny) .and. recvSizeAboveUR > 0) then
                call mpi_recv(recv_pAboveUR(1 : recvSizeAboveUR), recvSizeAboveUR, mpi_pType,&
                            my_rank + xysheet + nny + 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from above (z) and R
            if (mod(my_rank, xysheet) < xysheet - nny .and. recvSizeAboveR > 0) then
                call mpi_recv(recv_pAboveR(1 : recvSizeAboveR), recvSizeAboveR, mpi_pType,&
                            my_rank + xysheet + nny, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv array from above (z) and DR
            if (mod(my_rank,nny) /= 0 .and. mod(my_rank, xysheet) < (xysheet - nny) .and. recvSizeAboveDR > 0) then
                call mpi_recv(recv_pAboveDR(1 : recvSizeAboveDR), recvSizeAboveDR, mpi_pType,&
                            my_rank + xysheet + nny - 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif

            !recv size from above (z) and D
            if (mod(my_rank,nny) /= 0 .and. recvSizeAboveD > 0) then
                call mpi_recv(recv_pAboveD(1 : recvSizeAboveD), recvSizeAboveD, mpi_pType,&
                            my_rank + xysheet - 1, tag + my_rank, mpi_comm_world, status, ierror)
            endif


        endif
        


    end select

    deallocate(idx)
    allocate(idx(5*Nactive))
    idxBig = (/(i, i = 1, 5*Nactive)/)

    ! idx will show where to fill the newly received particles into the non-active spaces in "p"
    notActive = count(.not. p%active)
    ! allocate(idx(2*notActive), idxBig(2*notActive))
    idx(1 : notActive) = pack(idxBig, .not. p%active)


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
        totRecvDiag4 = totRecvDiag3 + recvSizeUR
        totRecvAbove1 = totRecvDiag4 + recvSizeAbove
        totRecvAbove2 = totRecvAbove1 + recvSizeAboveR
        totRecvAbove3 = totRecvAbove2 + recvSizeAboveUR
        totRecvAbove4 = totRecvAbove3 + recvSizeAboveDR
        totRecvAbove5 = totRecvAbove4 + recvSizeAboveD
        totRecvAbove6 = totRecvAbove5 + recvSizeAboveDL
        totRecvAbove7 = totRecvAbove6 + recvSizeAboveL
        totRecvAbove8 = totRecvAbove7 + recvSizeAboveUL
        totRecvAbove9 = totRecvAbove8 + recvSizeAboveU
        totRecvBelow1 = totRecvAbove9 + recvSizeBelow
        totRecvBelow2 = totRecvBelow1 + recvSizeBelowR
        totRecvBelow3 = totRecvBelow2 + recvSizeBelowUR
        totRecvBelow4 = totRecvBelow3 + recvSizeBelowDR
        totRecvBelow5 = totRecvBelow4 + recvSizeBelowD
        totRecvBelow6 = totRecvBelow5 + recvSizeBelowDL
        totRecvBelow7 = totRecvBelow6 + recvSizeBelowL
        totRecvBelow8 = totRecvBelow7 + recvSizeBelowUL
        totRecv       = totRecvBelow8 + recvSizeBelowU

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

    case(2)

        if (totRecv > 0) then
            recv_pLeft(recvSizeLeft + 1 : totRecvLR) = recv_pRight(1 : recvSizeRight)
            if (nny > 1) then
                recv_pLeft(totRecvLR + 1 : totRecvLRU) = recv_pUp(1 : recvSizeUp)
                recv_pLeft(totRecvLRU + 1 : totRecvCard) = recv_pDown(1 : recvSizeDown)
                recv_pLeft(totRecvCard + 1 : totRecvDiag1) = recv_pDL(1 : recvSizeDL)
                recv_pLeft(totRecvDiag1 + 1 : totRecvDiag2) = recv_pDR(1 : recvSizeDR)
                recv_pLeft(totRecvDiag2 + 1 : totRecvDiag3) = recv_pUL(1 : recvSizeUL)
                recv_pLeft(totRecvDiag3 + 1 : totRecvDiag4) = recv_pUR(1 : recvSizeUR)
                recv_pLeft(totRecvDiag4 + 1 : totRecvAbove1) = recv_pAbove(1 : recvSizeAbove)
                recv_pLeft(totRecvAbove1 + 1 : totRecvAbove2) = recv_pAboveR(1 : recvSizeAboveR)
                recv_pLeft(totRecvAbove2 + 1 : totRecvAbove3) = recv_pAboveUR(1 : recvSizeAboveUR)
                recv_pLeft(totRecvAbove3 + 1 : totRecvAbove4) = recv_pAboveDR(1 : recvSizeAboveDR)
                recv_pLeft(totRecvAbove4 + 1 : totRecvAbove5) = recv_pAboveD(1 : recvSizeAboveD)
                recv_pLeft(totRecvAbove5 + 1 : totRecvAbove6) = recv_pAboveDL(1 : recvSizeAboveDL)
                recv_pLeft(totRecvAbove6 + 1 : totRecvAbove7) = recv_pAboveL(1 : recvSizeAboveL)
                recv_pLeft(totRecvAbove7 + 1 : totRecvAbove8) = recv_pAboveUL(1 : recvSizeAboveUL)
                recv_pLeft(totRecvAbove8 + 1 : totRecvAbove9) = recv_pAboveU(1 : recvSizeAboveU)
                recv_pLeft(totRecvAbove9 + 1 : totRecvBelow1) = recv_pBelow(1 : recvSizeBelow)
                recv_pLeft(totRecvBelow1 + 1 : totRecvBelow2) = recv_pBelowR(1 : recvSizeBelowR)
                recv_pLeft(totRecvBelow2 + 1 : totRecvBelow3) = recv_pBelowUR(1 : recvSizeBelowUR)
                recv_pLeft(totRecvBelow3 + 1 : totRecvBelow4) = recv_pBelowDR(1 : recvSizeBelowDR)
                recv_pLeft(totRecvBelow4 + 1 : totRecvBelow5) = recv_pBelowD(1 : recvSizeBelowD)
                recv_pLeft(totRecvBelow5 + 1 : totRecvBelow6) = recv_pBelowDL(1 : recvSizeBelowDL)
                recv_pLeft(totRecvBelow6 + 1 : totRecvBelow7) = recv_pBelowL(1 : recvSizeBelowL)
                recv_pLeft(totRecvBelow7 + 1 : totRecvBelow8) = recv_pBelowUL(1 : recvSizeBelowUL)
                recv_pLeft(totRecvBelow8 + 1 : totRecv)       = recv_pBelowU(1 : recvSizeBelowU)
            endif
            do i = 1, totRecv
                p(idx(i)) = recv_pLeft(i)
                p(idx(i))%active = .true.
            enddo
        endif


    end select

    deallocate(idx,recv_pLeft,recv_pRight,recv_pUp,recv_pDown,recv_pDL,recv_pDR,recv_pUL,&
                recv_pUR,recv_pAbove,recv_pAboveR,recv_pAboveUR,recv_pAboveDR,recv_pAboveD,&
                recv_pAboveDL,recv_pAboveL,recv_pAboveUL,recv_pAboveU,recv_pBelow,recv_pBelowR,&
                recv_pBelowUR,recv_pBelowDR,recv_pBelowD,recv_pBelowDL,recv_pBelowL,recv_pBelowUL,&
                recv_pBelowU)
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
    integer,             intent(in   ) :: n ! number of entries in A
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
    integer                          :: start(n), finish(n), Nclose ! used for building the mass transfer matrix
    real(pDki)                       :: tmpmass(n) ! temporary array for holding particle masses
    Type(ParticleType)               :: tmp_p(n) ! temporary particle array for dealing with ghost particles
    integer                          :: idx(n) ! indexing array
    integer                          :: i
    integer                          :: nNotGhost, idxNotGhost(n), idxNotGhostTmp(n) ! indexing array for non-ghost particles
    logical                          :: logNotGhost(n) ! logical array for non-ghost particles

    tmp_p = p(idxActive(1 : n))
    
   ! Write(*,*) "Before DistmatSparse"

    ! build the pairwise distance matrix
    call build_DistmatSparse(n, cutdist, tmp_p, Emat, start, finish, Nclose)


    ! if (my_rank == master) then
    !     Write(*,*) "After DistmasSparse.. Before PmatSparse"
    ! endif
    ! build the matrix of co-location probabilities
    call build_PmatSparse(n, denom, start, finish, Emat)

   !  if (my_rank == master) then
   !     Write(*,*) "After PmatSparse.. Before EmatSparse"
   ! endif
    ! build the mass transfer matrix
    call build_EmatSparse(n, Nclose, Emat, beta)

   ! Write(*,*) "After EmatSparse"

    tmpmass = tmp_p%mass

    ! conduct the mass transfers with sparse matrix-vector multiplication
    call SP_matVecMult(Emat, tmpmass, Nclose, tmp_p%mass)

    ! only change the masses of non-ghost particles
    idx = (/(i, i = 1, n)/)

    logNotGhost = tmp_p%ghost(1) .or. tmp_p%ghost(2) .or. tmp_p%ghost(3) &
                .or. tmp_p%ghost(4) .or. tmp_p%ghost(5) .or. tmp_p%ghost(6) &
                .or. tmp_p%ghost(7) .or. tmp_p%ghost(8) .or. tmp_p%ghost(9) &
                .or. tmp_p%ghost(10) .or. tmp_p%ghost(11) .or. tmp_p%ghost(12) &
                .or. tmp_p%ghost(13) .or. tmp_p%ghost(14) .or. tmp_p%ghost(15) &
                .or. tmp_p%ghost(16) .or. tmp_p%ghost(17) .or. tmp_p%ghost(18) &
                .or. tmp_p%ghost(19) .or. tmp_p%ghost(20) .or. tmp_p%ghost(21) &
                .or. tmp_p%ghost(22) .or. tmp_p%ghost(23) .or. tmp_p%ghost(24) &
                .or. tmp_p%ghost(25) .or. tmp_p%ghost(26)

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
    logNotGhost = tmp_p%ghost(1) .or. tmp_p%ghost(2) .or. tmp_p%ghost(3) .or. tmp_p%ghost(4) .or. tmp_p%ghost(5) .or. tmp_p%ghost(6) &
                    .or. tmp_p%ghost(7) .or. tmp_p%ghost(8) .or. tmp_p%ghost(9) .or. tmp_p%ghost(10) .or. tmp_p%ghost(11) &
                    .or. tmp_p%ghost(12) .or. tmp_p%ghost(13) .or. tmp_p%ghost(14) .or. tmp_p%ghost(15) .or. tmp_p%ghost(16) &
                    .or. tmp_p%ghost(17) .or. tmp_p%ghost(18) .or. tmp_p%ghost(19) .or. tmp_p%ghost(20) .or. tmp_p%ghost(21) &
                    .or. tmp_p%ghost(22) .or. tmp_p%ghost(23) .or. tmp_p%ghost(24) .or. tmp_p%ghost(25) .or. tmp_p%ghost(26)
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
    integer,                          intent(  out) :: start(n), finish(n) ! indices (in the Distmat vectors) for the start and finish of each column of Distmat
    integer,                          intent(  out) :: Nclose ! total number of neighbor particles found by the kD search (i.e, the length of the vectors in Distmat)

    ! ========================== LOCAL VARIABLES ===============================
    type(kDRS_ResultType), allocatable :: neighbors(:) ! holds the results of the kD tree fixed radius search
    integer                            :: i ! loop iterator
    integer                            :: tmpstart ! temporary variable to handle the n+1^th calculation of start in the loop

    ! conduct the kD tree fixed-radius search
    ! NOTE: this returns squared distances between particles
!    Write(*,*) "Before f-r search"
    call fixedRadSearch(n, cutdist, p, neighbors)
   ! Write(*,*) "After f-r search"

    ! allocate Distmat to have length = total number of neighbors found by the kD search
    Nclose = 0
    Nclose = sum(neighbors%num)

    ! Write(*,*) "Nclose is",Nclose
    ! Write(*,*) "Rank",my_rank,"before allocation of Nclose =",Nclose

    allocate(Distmat(Nclose))

!    Write(*,*) "Before Distmat loop and n=",n
   ! Write(*,*) "Before Distmat loop and Nclose=",Nclose
    ! fill in Distmat
    tmpstart = 1
    do i = 1, n
        start(i) = tmpstart
        finish(i) = start(i) - 1 + neighbors(i)%num
!        Write(*,*) "start = ",start(i)
!        Write(*,*) "neighbors%num = ", neighbors(i)%num
!        Write(*,*) "finish = ", finish(i)
        Distmat(start(i) : finish(i))%col = i
        Distmat(start(i) : finish(i))%row = neighbors(i)%idx
        Distmat(start(i) : finish(i))%val = real(neighbors(i)%rad, pDki)
        tmpstart = finish(i) + 1
    enddo
    
!    Write(*,*) "After Distmat loop"
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
    real(kdkind)                    :: locs(dim,n), locs1(n), locs2(n), locs3(n) ! locations array to be passed to kD tree
    real(kdkind)                    :: r2 ! value of squared search radius for kD tree

    ! convert particle locations to kdkind
    ! locs = real(p%loc(1), kdkind)
    locs1 = real(p%loc(1), kdkind)
    locs2 = real(p%loc(2), kdkind)
    locs3 = real(p%loc(3), kdkind)
    locs(1,:) = locs1
    locs(2,:) = locs2
    locs(3,:) = locs3

    ! squared search cutoff distance in kdkind
    r2   = real(cutdist, kdkind)**2

    ! build the KD tree and search it
    call maketree(tree, 3, n, locs)

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
    integer,             intent(in   ) :: start(n), finish(n) ! indices (in the Distmat vectors) for the start and finish of each column of Distmat
    type(SparseMatType), intent(inout) :: Pmat(:) ! sparse matrix of co-location probabilities

    ! ========================== LOCAL VARIABLES ===============================
    integer :: i

    ! calculate the encounter probabilities
    ! NOTE: the leading constant term is neglected here, due to the normalization
    ! NOTE: Pmat%val is already squared separation distance because of the kD search
    Pmat%val = exp(-1.0_pDki*(Pmat%val / denom)) / (denom * sqrt(pi)**3)

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
    integer,             intent(in   ) :: Nclose ! total number of neighbor particles found by the kD search (i.e, the length of the vectors in Distmat)
    type(SparseMatType), intent(inout) :: Emat(:) ! sparse mass-transfer matrix
    real(pDki),          intent(in   ) :: beta ! beta parameter
    ! ========================== LOCAL VARIABLES ===============================
    integer    :: i, j
    integer    :: diag(n) ! linear indices of the diagonal elements of Pmat
    real(pDki) :: rowsum(n), colsum(n) ! arrays holding the row/col sums of Pmat

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

    real(pDki),         intent(in   ) :: Dlims(6) ! subdomain boundary limits
    integer,            intent(in   ) :: coreNp ! number of particles in local subdomain
    integer,            intent(in   ) :: nnx, nny ! grid specifications
    real(pDki),         intent(in   ) :: xmidpt, ymidpt ! midpoint of domain for IC
    real(pDki),         intent(  out) :: spillX, spillY ! location of spill
    type(ParticleType), intent(inout) :: p(:) ! particle array


            call random_number(p%loc(1))
            call random_number(p%loc(2))
            call random_number(p%loc(3))
            p(1 : coreNp)%loc(1) = Dlims(1) + (Dlims(2) - Dlims(1)) * p(1 : coreNp)%loc(1)
            p(1 : coreNp)%loc(2) = Dlims(3) + (Dlims(4) - Dlims(3)) * p(1 : coreNp)%loc(2)
            p(1 : coreNp)%loc(3) = Dlims(5) + (Dlims(6) - Dlims(5)) * p(1 : coreNp)%loc(3)            

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
            !     p(1)%mass = 1.0_pDki
            ! endif

            ! if (any(p(1 : coreNp)%mass == 1.0_pDki) .eqv. .true.) then
            !     Write(*,*) "Core number",my_rank,"has the intial mass."
            ! endif

            ! if (any(p(1 : coreNp)%mass == 1.0_pDki) .eqv. .true.) then
            !     spillX = p(1)%loc(1)
            !     spillY = p(1)%loc(2)
            ! endif       


end subroutine InitialSpacingAndIC

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



! ! ! subroutine EvenlySpaced()

! ! ! end subroutine EvenlySpaced


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




end module mod_DDC_2D_mpi
