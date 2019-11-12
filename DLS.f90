
MODULE DLS      

  implicit none
  private
  public :: infoDLS
  public :: DLS_MethodCount, DLS_LongName, DLS_ShortName
  public :: DLS_Partition, DLS_GroupSetup, DLS_Setup, DLS_Terminated
  public :: DLS_StartLoop, DLS_StartChunk, DLS_EndChunk, DLS_EndLoop
  public :: DLS_Finalize

  include 'mpif.h'

  ! include 'sched.h'
  integer, parameter :: maxProcs=128, maxChunks=500
  integer, parameter :: maxGrpSize=4, minGrpSize=2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DLS - a Fortran 90+MPI module for dynamic loop scheduling
!       on general-purpose clusters
! by RL Carino (rlc@erc.msstate.edu)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified by Ali Mohammed <ali.mohammed@unibas.ch>
!              (Nov. 2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This program is free software; you can redistribute it and/or modify it       
!under the terms of the license (GNU LGPL) which comes with this package.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!DLS4LB is designed for load balancing Fortran 90+MPI applications that
! contain Type I or II loops described below. The iterates are
! assumed be CPU-intensive, to be worth the load balancing overhead.
!
! Type I
!   ! begin i-loop
!   do i=1,N
!     ... i-iterate
!   end do
!   ! end i-loop
!
!
! Type II
!   ! begin j-loop
!   do j=1,M
!     ... part of j-iterate
!     ! begin i-loop
!     do i=1,N(j)
!     ... i-iterate of j-iterate
!     end do
!     ! end i-loop
!     ... part of j-iterate
!   end do
!   ! end j-loop
!
! Dynamic load balancing is achieved in the application
! by invoking DLS routines to manipulate the loop indices.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! For Type I loops, the application must be modified as:
!
!   use DLS
!   include 'mpif.h'
!   ...
!   type (infoDLS) info
!   integer method, iStart, iSize
!   integer iIters
!   double precision iTime
!   ...
!   method = ... (choice of loop scheduling technique)
!   call DLS_Setup (MPI_COMM_WORLD, info)
!   call DLS_StartLoop (info, 1, N, method)
!   do while ( .not. DLS_Terminated(info) )
!     call DLS_StartChunk (info, iStart, iSize)
!     ! begin i-loop, new extents
!     do i=iStart, iStart+iSize-1    ! 1,N
!       ... i-iterate
!     end do
!     ! end i-loop
!     call DLS_EndChunk (info)
!   end do ! while ( .not. DLS_Terminated(info) )
!   call DLS_EndLoop (info, iIters, iTime)
!
!
! where:
!
!   DLS_Setup (MPI_COMM_WORLD, info)
!      Initializes a dynamic load balancing environment on
!      MPI_COMM_WORLD. Information about this environment is
!      stored in "info" (see type declaration infoDLS below)
!
!   DLS_StartLoop (info, 1, N, method)
!      The synchronization point to start loop execution.
!      (1, N) is the loop range, and "method" is a user-
!      specified index [1..9] for the loop scheduling method
!      (see array DLS_LongName() below)
!
!   DLS_Terminated(info) 
!      returns .TRUE. if all loop iterates have been
!      executed
!
!   DLS_StartChunk (info, iStart, iSize)
!      Returns a range for a chunk of iterates, starting
!      at "iStart", with size "iSize"
!
!   do i=iStart, iStart+iSize-1    ! 1,N
!     ... i-iterate
!   end do
!      The original "serial" code, but now executing
!      a chunk of iterates instead of "i=1,N"
!
!   DLS_EndChunk (info)
!      Indicates the end of execution of a chunk of iterates
!
!   DLS_EndLoop (info, iIters, iTime)
!      The synchronization point to end loop execution. 
!      "iIters" is the number of iterates done by the calling
!      processor, and "iTime" is its cost (in seconds)
!      
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! For Type II loops, the application must be modified as:
!
!   use DLS
!   include 'mpif.h'
!   ...
!   type (infoDLS) iInfo, jInfo
!   integer iComm, jComm    ! communicators
!   integer iMethod, iStart, iSize, jMethod, jStart, jSize
!   integer iIters, jIters
!   double precision iTime, jTime
!   integer coordinator
!   ...
!   iMethod = ... (choice of loop scheduling technique)
!   jMethod = ... (choice of loop scheduling technique)
!   coordinator = 0
!   call DLS_GroupSetup (MPI_COMM_WORLD, coordinator, jInfo, iInfo)
!   call DLS_StartLoop (jInfo, 1, M, jMethod)
!   do while ( .not. DLS_Terminated(jInfo) )
!     call DLS_StartChunk (jInfo, jStart, jSize)
!
!     ! begin j-loop code, new extents
!     do j=jStart, jStart+jSjze-1    ! 1,M
!       ... part of j-iterate 
!
!       call DLS_StartLoop (iInfo, 1, N(j), iMethod)
!       do while ( .not. DLS_Terminated(iInfo) )
!         call DLS_StartChunk (iInfo, iStart, iSize)
!         ! begin i-loop, new extents
!         do i=iStart, iStart+iSize-1    ! 1,N
!           ... i-iterate of j-iterate
!         end do
!         ! end i-loop
!         call DLS_EndChunk (iInfo)
!       end do ! while ( .not. DLS_Terminated(iInfo) )
!       call DLS_EndLoop (iInfo, iIters, iTime)
!
!       ... part of j-iterate 
!     end do
!     ! end j-loop
!
!     call DLS_EndChunk (jInfo)
!   end do ! while ( .not. DLS_Terminated(jInfo) )
!   call DLS_EndLoop (jInfo, jIters, jTime)
!
! where:
!
!   DLS_GroupSetup (MPI_COMM_WORLD, coordinator, jInfo, iInfo)
!      Splits MPI_COMM_WORLD into non-overlapping "iComm"s and
!      a "jComm" comprised of the "coordinator" and the
!      foremen (rank 0 of the iComms).  Similar to DLS_Setup(), 
!      DLS_GroupSetup() initializes load balancing environments 
!      in "jComm" and "iComm", keeping the data in "jInfo" and 
!      "iInfo" respectively. Graphically:
!
!
!  ==============================================
!  = MPI_COMM_WORLD        .-------------.      =           
!  =    .-----------.       \  x  x   x /       =      
!  =     \ iComm x /         \ iComm   /        =      
!  =      \x  x x /           \     x /         =      
!  =       \ x   /             \x x  /          =      
!  =   .----\---/---------------\---/------.    =
!  =   |     \x/     coordinator \x/       |    =
!  =   |      .          x        .        |    =
!  =   |           .                 .     |    =
!  =   | jComm    /x\               /x\    |    =
!  =   .---------/---\-------------/---\---.    =
!  =            /iComm\           / x   \       =           
!  =           / x  x  \         / x  x  \      =           
!  =          -----------       / iComm x \     =           
!  =                           .-----------.    =           
!  ==============================================
!
!      The symbol "x" denotes a processor, "x"s inside
!      a triangle make up an "iComm", while "x"s inside
!      the rectangle make up the jComm. The j-iterates 
!      scheduled in "jComm" while the i-iterates for a
!      given j-iterate are scheduled in "iComm"
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, parameter :: STAT =  0
  integer, parameter :: SS   =  1
  integer, parameter :: FSC  =  2
  integer, parameter :: MFSC =  3
  integer, parameter :: GSS  =  4
  integer, parameter :: TSS  =  5
  integer, parameter :: FAC  =  6
  integer, parameter :: WF   =  7
  integer, parameter :: AWF  =  8
  integer, parameter :: AWFB =  9
  integer, parameter :: AWFC =  10
  integer, parameter :: AWFD =  11
  integer, parameter :: AWFE =  12
  integer, parameter :: AF   =  13
  integer, parameter :: SimAS =  14


! method names (long)
  integer, parameter :: DLS_MethodCount = 15 ! 10
  character (len=25), dimension(0:DLS_MethodCount-1) :: DLS_LongName = (/ &
     'STATIC SCHEDULING        ', 'SELF-SCHEDULING          ', 'FIXED SIZE SCHEDULING    ', &
     'MODIFIED FSC             ', 'GUIDED SELF SCHEDULING   ', 'TRAPIZIOD SELF-SCHEDULING', &
     'FACTORING                ', 'WEIGHTED FACTORING       ', 'ADAPTIVE WEIGHTED FAC    ', &
     'BATCH AWF                ', 'CHUNK AWF                ', 'CHUNK AWF (chunk time)   ', &
     'BATCH AWF (chunk time)   ', 'ADAPTIVE FACTORING       ', 'Simulation-assisted      ' &
     /)

! method names (short)
  character (len=5), dimension(0:DLS_MethodCount-1) :: DLS_ShortName = (/ &
     'STAT ', 'SS   ', 'FSC  ', 'MFSC ', 'GSS  ', 'TSS  ', 'FAC  ', 'WF   ', &
     'AWF  ', 'AWFB ', 'AWFC ', 'AWFD ', 'AWFE ' ,'AF   ', 'SimAS'/)

! message tags
  integer, parameter :: HNM_TAG =  9990
  integer, parameter :: CLR_TAG =  9980
  integer, parameter :: WRK_TAG =  9970
  integer, parameter :: REQ_TAG =  9960
  integer, parameter :: TRM_TAG =  9950
  integer, parameter :: END_TAG =  9940
  integer, parameter :: AWF_TAG =  9801

! loop scheduling variables
  type infoDLS
    integer :: comm, crew              ! communicators
    integer :: Foreman                 ! rank of foreman
    integer :: myRank                  ! rank of process 
    integer :: firstRank, lastRank     ! range of foreman's loops
    integer :: commSize, crewSize     ! size of my work group
    integer :: method
    integer :: FirstIter, LastIter, N  ! total iterates
    integer :: itersScheduled          ! total iterates scheduled
    integer :: batchSize               ! iterations in batch
    integer :: batchRem                ! remaining in batch
    integer :: minChunkSize            ! minimum chunk size
    integer :: maxChunkSize            ! maximum chunk size
    integer :: minChunk=-1
    integer :: breakAfter=-1
    integer :: requestWhen=-1
    integer :: chunkMFSC               ! assume # of chunks is that of FAC
    integer :: chunkFSC                ! FSC chunk size
    integer :: chunkStart              ! running start of chunks 
    integer :: numChunks               ! no. of chunks generated 
    integer :: probeFreq               ! iterates to do before message probe 
    integer :: sendRequest             ! iterates left when to send request 
    integer :: numENDed                ! no. of processes terminated 
    integer :: finishedOne             ! no. of processes that have finished first chunk 
    integer :: myExecs                 ! no. of chunks executed by this process 
    integer :: timeStep                ! no. of time-steps executed by this process 
    integer :: rStart, rSize           ! start, size of current chunk 
    integer :: wStart, wSize           ! start, size of remaining subchunk 
    integer :: nextStart, nextSize     ! foreman's response to advance request 
    integer :: subChunkSize
    integer :: myIters                 ! no. of iters executed by this process 
    logical :: gotWork                 ! termination flag 
    logical :: req4WRKsent             ! advance request sent toggle
    logical :: nextWRKrcvd             ! response recv toggle
    double precision :: sigma=0.0001         !sigma for FSC
    double precision :: h=0.0002              ! scheduling overhead h for FSC
    double precision,ALLOCATABLE,DIMENSION(:) :: weights        ! weights array for WF
    integer :: TSSchunk
    integer :: TSSdelta   
    double precision :: kopt0 ! , gN, gchunks, gsumt1, gsumt2, gsumovrhd, gsigma, gh ! FSC-A 
    double precision :: t0, t1, sumt1, sumt2 
    double precision :: mySumTimes, mySumSizes
    double precision :: workTime             ! useful work time 
    double precision,ALLOCATABLE,DIMENSION(:) :: stats 
!    double precision :: stats(0:3*maxProcs-1) ! performance stats
    !integer :: chunkMap(0:3*maxChunks-1) ! chunk map
  end type infoDLS

  ! scratch variables
  integer :: ierror
  integer :: chunkInfo(0:1)              ! Chunk info buffer (start, size)
  integer, dimension(MPI_STATUS_SIZE) :: mStatus, tStatus
  double precision :: perfInfo(0:3)      ! Performance info buffer (mu/time, sigma/size)
  
!.... DLS4LB variables ...Ali

      type (infoDLS) info
      integer method, iStart, iSize
      integer iIters
      double precision iTime
      INTEGER npini_BAK,npend_BAK
      DOUBLE PRECISION tgrav1, tgrav2

      public :: info
      public :: method, iStart, iSize
      public :: iIters
      public :: iTime
      public :: npini_BAK,npend_BAK
      public :: tgrav1, tgrav2


contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DLS_GroupSetup (world, coordinator, GS, LS)
    integer, intent (in) :: world, coordinator
    type (infoDLS), intent (out) :: GS, LS
    integer :: jcomm, icomm, tP

    call DLS_Partition(world, coordinator, jcomm, icomm)

    GS%comm = jcomm
    GS%crew = icomm
    LS%comm = icomm
    LS%crew = icomm
    if (icomm/=MPI_COMM_NULL) then ! workers
      call DLS_Setup (icomm, LS)
      call MPI_Comm_size(icomm, GS%crewSize, ierror) 
    else ! the scheduler
      LS%myRank = -1
      LS%commSize = 0
      LS%crewSize  = 0
      GS%crewSize  = 0
    end if
    ! foreman and scheduler
    if (jcomm/=MPI_COMM_NULL) then
      call MPI_Comm_size(jcomm, tP, ierror) 
      call MPI_Comm_rank(jcomm, GS%myRank, ierror) 
      GS%commSize = tP-1
      GS%firstRank = 1
      GS%lastRank = tP-1
    else ! workers
      GS%myRank = -1
      GS%commSize = 0
    end if
    return
  end subroutine DLS_GroupSetup


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DLS_Setup (icomm, info)
    integer, intent (in) :: icomm
    type (infoDLS), intent (out) :: info
    integer :: tP, worker
    call MPI_Comm_size(icomm, tP, ierror) 
    call MPI_Comm_rank(icomm, info%myRank, ierror) 
    info%comm = icomm
    info%crew = MPI_COMM_NULL
    info%commSize = tP
    info%firstRank = 0
    info%lastRank = tP-1
    info%Foreman = 0
    info%timeStep = 0

    ! .....Allocate data
    ALLOCATE ( info%stats(0:3*info%commSize) )
    ALLOCATE ( info%weights(0:info%commSize) )

      do worker=info%firstRank,info%lastRank
        info%weights(worker) = 1.0 !PE weight ...change if working on a heterogeneous system
        info%stats(3*worker) = -1.0    ! mu 
        info%stats(3*worker+1) = -1.0  ! sigma 
        info%stats(3*worker+2) = 0.0   ! performance data count 
      end do

    return
  end subroutine DLS_Setup


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DLS_StartLoop (info, firstIter, lastIter, imeth) 
    integer, intent (in) :: firstIter, lastIter, imeth
    type (infoDLS), intent (in out) :: info
    integer tSize, worker
    DOUBLE PRECISION n, K
    integer ::  NULLloc
    DOUBLE PRECISION awap, trw, weight
    integer i


    

    NULLloc = 0
    info%wSize = 0  ! remaining iterates in current chunk 
    info%gotWork = .true.
    info%workTime = 0.0 
    info%myIters = 0
    info%timeStep = info%timeStep + 1
    

!    ! .....Allocate data
!    ALLOCATE ( info%stats(0:3*info%commSize) ) 
!    ALLOCATE ( info%weights(0:info%commSize) )

    info%N = lastIter - firstIter + 1
    if ( (info%comm==MPI_COMM_NULL) .or. (info%N<=0) ) return

    info%method = imeth
    if (imeth>=DLS_MethodCount .or. imeth<0) then
      info%method = 0
    else
      info%method = imeth
    end if

! calculate AWF weights
    if ( (info%method == AWF) .and. (info%myRank == info%Foreman)) then
        if (info%timeStep == 1) then
            do i = info%firstRank,info%lastRank
                 info%weights(i) = 1.0d0
             end do
         else ! all ranks have wap
          awap = 0.0d0  ! average weighted performance
          do i=info%firstRank,info%lastRank
            !write(*,*) "rank", i, "info%stat", info%stats(3*i),info%stats(3*i+2)
            awap = awap + info%stats(3*i)
          end do
          awap = awap/info%commSize

          trw = 0.0d0  ! total ref weight (refwt(i) = awap/info%stats(3*i)
          do i=info%firstRank,info%lastRank
            trw = trw + awap/info%stats(3*i)
          end do

           do i=info%firstRank,info%lastRank
               info%weights(i) = ((awap/info%stats(3*i))*info%commSize)/trw
          end do

         end if
    end if

    info%FirstIter = firstIter 
    info%LastIter = lastIter 

    !info%chunkMap(0) = firstIter     ! start of data 
    !info%chunkMap(1) = info%N     ! size of data 
    !info%chunkMap(2) = 0  ! chunks in this rank 
    info%numChunks = 0 

    info%myExecs = 0 
    info%mySumTimes = 0.0 
    info%mySumSizes = 0.0 

    ! TSS
    info%TSSchunk = CEILING( DBLE(info%N) / DBLE( 2*info%commSize)) ! should be ceiled
    n = CEILING(2*DBLE(info%N)/DBLE(info%TSSchunk+1)) !n=2N/f+l ! should be ceiled
    info%TSSdelta = DBLE(info%TSSchunk - 1)/ DBLE(n-1)



    tSize = (info%N+info%commSize-1)/info%commSize
    info%chunkMFSC = (0.55+tSize*log(2.0d0)/log( dble (tSize) ) )
    info%kopt0 = sqrt(2.0d0)*info%N/( info%commSize*sqrt(log(1.0d0*info%commSize)) )


    K=(SQRT(2.0)*info%N*info%h)/(info%sigma*info%commSize*SQRT(LOG(DBLE(info%commSize))))
    K= K **(2.0/3.0)
    info%chunkFSC =  CEILING(K)

!   info%gN = 0.0
!   info%gchunks = 0.0
!   info%gsumt1 = 0.0
!   info%gsumt2 = 0.0
!   info%gsumovrhd = 0.0

    info%nextWRKrcvd = .false.
    info%req4WRKsent = .false.
    info%finishedOne = 0 

    info%probeFreq = max(1, info%breakAfter) 
    info%sendRequest = max(1, info%requestWhen) 

    if (info%myRank == info%Foreman) then
      info%chunkStart = firstIter
      info%itersScheduled = 0 
      info%batchSize = 0 
      info%batchRem = 0 
      info%numENDed = 0 
      info%numChunks = 0 

!      do worker=info%firstRank,info%lastRank
!        info%weights(worker) = 1.0 !PE weight ...change if working on heterogeneous system
!        info%stats(3*worker) = -1.0    ! mu 
!        info%stats(3*worker+1) = -1.0  ! sigma 
!        info%stats(3*worker+2) = 0.0   ! performance data count 
!      end do

      if (info%minChunk>0) then
        info%minChunkSize = info%minChunk 
      else
        info%minChunkSize = max(1,info%chunkMFSC/2)        ! default min chunk size 
      end if
      info%maxChunkSize = (info%N+2*info%commSize-1)/(2*info%commSize)

      ! send initial work to each processor
      do worker=info%firstRank,info%lastRank
        if (info%chunkStart < info%lastIter) then 
          call SendChunk (info, worker)
        else
            call MPI_Send (NULLloc, 0, MPI_INTEGER, worker, END_TAG, info%comm, ierror)  !end worker
            info%numENDed = info%numENDed + 1 ! increment ended workers
        end if
      end do
    end if

    return
  end subroutine DLS_StartLoop

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DLS_EndLoop (info, niters, worktime)
    type (infoDLS), intent (in out) :: info
    integer, intent (out) :: niters
    double precision, intent (out) :: worktime
    integer i, loc
    if (info%comm==MPI_COMM_NULL) return
    niters = info%myIters
    worktime = info%workTime 

!     call MPI_Barrier(info%comm, ierror) 
 
!... Communicate time-step performance data for AWF
     if (info%method==AWF) then
        ! timestepping adaptive weighted factoring 
        ! mu = (chunk work time)/(chunk size) 
        PerfInfo(0) = ( info%mySumTimes + (info%timeStep)*info%workTime ) / &
            ( info%mySumSizes + 1.0*(info%timeStep)*(info%myIters) )
        PerfInfo(1) = 0.0
        PerfInfo(2) = 1.0*info%timeStep
!        write(*,*) "time step", info%timeStep
         call MPI_Gather(PerfInfo, 3, MPI_DOUBLE_PRECISION, info%stats, 3,MPI_DOUBLE_PRECISION, info%Foreman, info%comm, ierror)

!        write(*,*) "rank ", info%myRank, "send its performance",  PerfInfo(0)
!          if (info%myRank == info%Foreman) then ! if foreman ...recieve all performance data
!   
!            do i=info%firstRank,info%lastRank
!               write(*,*) "performance data", i, info%stats(3*i), info%stats(3*i+1), info%stats(3*i+2)
!            end do
!          end if
       end if 

    !    call MPI_Barrier(info%comm, ierror) 
    return
  end subroutine DLS_EndLoop

 subroutine DLS_Finalize (info)
    type (infoDLS), intent (in out) :: info
   
    if (info%comm==MPI_COMM_NULL) return

    DEALLOCATE (info%stats)
    DEALLOCATE (info%weights)
    return
  end subroutine DLS_Finalize


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function DLS_Terminated (info)
    type (infoDLS), intent (in) :: info
    logical done
    integer i
    integer, dimension(MPI_STATUS_SIZE) :: tStatus

    if (info%N<=0) then
      done = .true.
    else
      if (info%comm==MPI_COMM_NULL .and. info%crew/=MPI_COMM_NULL) then
        call MPI_Recv (done, 1, MPI_LOGICAL, 0, TRM_TAG, info%crew, tStatus, ierror)
      else if (info%comm/=MPI_COMM_NULL .and. info%crew/=MPI_COMM_NULL) then
        done = (.not. info%gotWork) .and. (info%wSize==0)
        do i=1,info%crewSize-1
          call MPI_Send (done, 1, MPI_LOGICAL, i, TRM_TAG, info%crew, ierror)
        end do
      else
        done = (.not. info%gotWork) .and. (info%wSize==0) 
      end if 
    end if 
    DLS_Terminated = done
    return
  end function DLS_Terminated 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DLS_StartChunk (info, chunkStart, chunkSize)
    type (infoDLS), intent (in out) :: info
    integer, intent (out) :: chunkStart, chunkSize

    integer tSize, tStart, worker
    logical :: MsgInQueue                ! message came in 

    integer :: loc, maxRemaining         ! source of chunk to be migrated 
    integer :: i, j, NULLloc

if (info%comm==MPI_COMM_NULL) then ! I'm just a simple worker
  call MPI_Recv (chunkInfo, 2, MPI_INTEGER, 0, WRK_TAG, info%crew, tStatus, ierror)
  chunkStart = chunkInfo(0)
  chunkSize = chunkInfo(1)
                     
else ! I'm the coordinator, or a foreman
                 
    if (info%wSize == 0) then
      call MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, info%comm, mStatus, ierror) 
      MsgInQueue = .true.
    else
      call MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, info%comm, MsgInQueue, mStatus, ierror) 
    end if

    do while (MsgInQueue)  

      select case ( mStatus(MPI_TAG) ) 
  
        case (WRK_TAG) 
          call MPI_Recv (ChunkInfo, 2, MPI_INTEGER, mStatus(MPI_SOURCE), WRK_TAG, &
            info%comm, tStatus, ierror) 

          if (info%wSize == 0) then ! no pending chunk 
            info%t0 = MPI_Wtime() ! elapsed time for chunk starts here 
            info%wStart = ChunkInfo(0) 
            info%wSize = ChunkInfo(1) 
            info%rStart = info%wStart 
            info%rSize = info%wSize 
            info%req4WRKsent = .false.  ! cancel request for work 

            call SetBreaks (info)

!           write(*,*) 'WRK_TAG recv by', info%myRank, info%wStart, info%wSize

            info%sumt1 = 0.0  ! for mu/wap 
            info%sumt2 = 0.0  ! for sigma 
          else  ! current chunk is not finished  save as next chunk 
            info%nextStart = ChunkInfo(0) 
            info%nextSize = ChunkInfo(1) 
            info%nextWRKrcvd = .true.

!           write(*,*) 'WRK_TAG (adv) recv by', info%myRank, info%nextStart, info%nextSize

          end if

        case (REQ_TAG)  ! received by foreman only  
          worker = mStatus(MPI_SOURCE)
          call MPI_Recv (PerfInfo, 4, MPI_DOUBLE_PRECISION, worker, &
            REQ_TAG, info%comm, tStatus, ierror) 

          if (info%method==AWFB .or. info%method==AWFC .or. info%method==AF .or. &
              info%method==AWFD .or. info%method==AWFE) then
            loc = int (PerfInfo(2)) 
            info%stats(3*loc+2) = info%stats(3*loc+2)+1.0 
            ! adaptive methods 
            info%stats(3*loc) = PerfInfo(0) 
            info%stats(3*loc+1) = PerfInfo(1) 

            if (info%finishedOne /= info%commSize) then
              ! workers that have not finished a first chunk  
              !  assume the lowest performance 
              j = loc 
              do i=info%firstRank,info%lastRank
                if ( (info%stats(3*i+2) > 0.0) .and. &
                     (info%stats(3*i) < info%stats(3*j)) ) j = i 
              end do
              info%finishedOne = 0 
              do i=info%firstRank,info%lastRank
                if (info%stats(3*i+2) == 0.0) then
                  info%stats(3*i) = info%stats(3*j) 
                  info%stats(3*i+1) = info%stats(3*j+1) 
                else 
                   info%finishedOne = info%finishedOne + 1
                end if
              end do
            end if
!         else if (info%method==10) then
!           info%gchunks = info%gchunks + 1.0
!           info%gsumt1 = info%gsumt1 + PerfInfo(0)
!           info%gsumt2 = info%gsumt2 + PerfInfo(1)
!           info%gN = info%gN + PerfInfo(2)
!           info%gsumovrhd = info%gsumovrhd + PerfInfo(3)
!           info%gsigma = sqrt( (info%gsumt2 - info%gsumt1/info%gN) / (info%gN-1.0) )
!           info%gh = info%gsumovrhd/info%gchunks  ! average h
          end if

!         write(*,*) 'REQ_TAG from', worker, PerfInfo(0), PerfInfo(1), int(PerfInfo(2))

          ! any remaining unscheduled iterates ?
          if (info%chunkStart <= info%lastIter) then
            call SendChunk (info, worker) 
          else ! all iterates scheduled
            info%numENDed = info%numENDed + 1
            if (worker /= info%myRank) then
!             write(*,*) 'END_TAG to', worker
!             call MPI_Send (info%chunkMap, 3*(info%numChunks+1), MPI_INTEGER, worker, &
              call MPI_Send (NULLloc, 0, MPI_INTEGER, worker, &
                     END_TAG, info%comm, ierror) 
            end if
            info%gotWork = info%NumENDed/=info%commSize ! foreman exits?
          end if

        case (END_TAG)  ! received by workers only 
!         call MPI_Recv (info%chunkMap, 3*info%maxChunks-1, MPI_INTEGER, mStatus(MPI_SOURCE), &
          call MPI_Recv (NULLloc, 0, MPI_INTEGER, mStatus(MPI_SOURCE), &
              mStatus(MPI_TAG), info%comm, tStatus, ierror) 
          info%gotWork = .false.

      end select 
      call MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, info%comm, &
              MsgInQueue, mStatus, ierror) 

    end do ! while (MsgInQueue) 

    chunkStart = info%wStart
    chunkSize = min (info%wSize, info%probeFreq) 
    if (info%method==AF) then ! .or. info%method == 10) then 
      chunkSize = min(1, chunkSize)
    end if
    info%subChunkSize = chunkSize
    if (info%subChunkSize/=0) info%t1 = MPI_Wtime() 

   !write (*,*) "[start chunk] probe frequency, ", info%probeFreq
    ! relay chunkStart, chunkSize to icomm
  if (info%crew/=MPI_COMM_NULL) then
    chunkInfo(0) = chunkStart
    chunkInfo(1) = chunkSize
    do i=1,info%crewSize-1
      call MPI_Send (chunkInfo, 2, MPI_INTEGER, i, WRK_TAG, info%crew, ierror)
    end do      
  end if      
                   
end if ! (info%comm/=MPI_COMM_NULL) then
              

    return
  end subroutine DLS_StartChunk 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine SetBreaks (info)
    type (infoDLS), intent (in out) :: info
    if (info%myRank == info%Foreman) then
      ! when to check for messages 
      if (info%breakAfter<0) then
        info%probeFreq = max( 1, (info%wSize+info%commSize-1)/info%commSize/4 ) 
      else
        info%probeFreq = max( 1, info%breakAfter)
      end if

      ! how many iterates left before requesting next chunk 
      if (info%requestWhen<0) then
        info%sendRequest = info%probeFreq
      else
        info%sendRequest = info%requestWhen
      end if

    else ! not the foreman

      ! how many iterates left before requesting next chunk 
      if (info%requestWhen<0) then
        info%sendRequest = max( 1, (15*info%wSize)/100 )
      else
        info%sendRequest = info%requestWhen
      end if

      ! when to check for messages 
      info%probeFreq = max( 1, info%wSize-info%sendRequest)

    end if
    return
  end subroutine SetBreaks 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine GetChunkSize (info, rank, chunkSize)
    type (infoDLS), intent (in out) :: info
    integer, intent (in) :: rank
    integer, intent (out) :: chunkSize

    integer :: i, tChunk, rem
    double precision :: bigD, bigT               ! AF
    double precision :: awap, trw, weight        ! AWF
    

    rem = info%N-info%itersScheduled                       ! remaining
    select case ( info%method )

    case (STAT)
         tChunk = CEILING(DBLE (info%N)/ DBLE(info%commSize))
         info%batchSize = tChunk
         info%batchRem = min( info%batchSize, rem)

      case (SS)
         tChunk = 1
         info%batchSize = tChunk
         info%batchRem = min( info%batchSize, rem)
        
      case (FSC)
         tChunk = min( info%chunkFSC, rem)
         info%batchSize = tChunk
         info%batchRem = min( info%batchSize, rem)        
         

      case (MFSC) ! fixed size scheduling
        tChunk = min( info%chunkMFSC, rem)
        info%batchSize = tChunk
        info%batchRem = min( info%batchSize, rem)
     
      case (GSS) ! guided self scheduling
        tChunk = max( (rem+info%commSize-1)/info%commSize, info%minChunkSize )
        tChunk = min ( rem, tChunk )
        info%batchSize = tChunk
        info%batchRem = min( info%batchSize, rem)


      case (TSS) 
         tChunk = info%TSSchunk
        tChunk = min(rem, tChunk)
        tChunk = max( info%minChunkSize, tChunk)
        info%TSSchunk = tChunk - info%TSSdelta
        info%batchSize = tChunk
        info%batchRem = min( info%batchSize, rem)
        


      case (FAC) ! factoring
        if (info%batchRem == 0) then
          tChunk = max ( info%minChunkSize, (rem+2*info%commSize-1)/(2*info%commSize) )
          info%batchSize = info%commSize*tChunk
          info%batchRem = min (info%batchSize, rem)
        end if
        ! else use current batchSize
        tChunk = max( info%minChunkSize, info%batchSize/info%commSize )
        tChunk = min ( rem, tChunk )

      case (WF)
         if (info%batchRem == 0) then
            tChunk = max ( info%minChunkSize,(rem+2*info%commSize-1)/(2*info%commSize) )
            info%batchSize = info%commSize*tChunk
            info%batchRem = min (info%batchSize, rem);
         end if
         ! else use current batchSize */
         tChunk = max( info%minChunkSize, INT(info%batchSize/info%commSize*info%weights(rank)))
         tChunk = min ( rem, tChunk )
     
    case (AWF)
         if (info%batchRem == 0) then
            tChunk = max (info%minChunkSize,(rem+2*info%commSize-1)/(2*info%commSize) )
            info%batchSize = info%commSize*tChunk
            info%batchRem = min (info%batchSize, rem);
         end if
         ! else use current batchSize */
         tChunk = max( info%minChunkSize, INT(info%batchSize/info%commSize*info%weights(rank)))
         tChunk = min ( rem, tChunk )



      case (AWFB, AWFD) ! batch adaptive weighted factoring
        if (info%stats(3*rank) < 0.0) then
          tChunk = info%minChunkSize
          info%batchSize = min(rem, tChunk)
          info%batchRem = info%batchSize
        else ! all ranks have wap
          awap = 0.0d0  ! average weighted performance
          do i=info%firstRank,info%lastRank
            awap = awap + info%stats(3*i)
          end do
          awap = awap/info%commSize
     
          trw = 0.0d0  ! total ref weight (refwt(i) = awap/info%stats(3*i)
          do i=info%firstRank,info%lastRank
            trw = trw + awap/info%stats(3*i)
          end do
     
          ! normalized weight for rank
          weight = ((awap/info%stats(3*rank))*info%commSize)/trw
          info%weights(rank) = weight
     
          if (info%batchRem == 0) then
            tChunk = max( info%minChunkSize, (rem+2*info%commSize-1)/(2*info%commSize) )
            info%batchSize = info%commSize*tChunk
            info%batchRem = min (info%batchSize, rem)
          end if
          ! else use current batchSize
          tChunk = weight*(info%batchSize/info%commSize) + 0.55d0
          tChunk = max( info%minChunkSize, tChunk)
          tChunk = min ( rem, tChunk )
        end if

      case (AWFC, AWFE) ! chunk adaptive weighted factoring
        if (info%stats(3*rank) < 0.0) then
          tChunk = info%minChunkSize
        else ! all ranks have wap
          awap = 0.0d0  ! average weighted performance
          do i=info%firstRank,info%lastRank
            awap = awap + info%stats(3*i)
          end do
          awap = awap/info%commSize

          trw = 0.0d0  ! total ref weight (refwt(i) = awap/info%stats(3*i)
          do i=info%firstRank,info%lastRank
            trw = trw + awap/info%stats(3*i)
          end do

          ! normalized weight for rank
          weight = ((awap/info%stats(3*rank))*info%commSize)/trw
          info%weights(rank) = weight
          tChunk = weight*((rem+2*info%commSize-1)/(2*info%commSize)) + 0.55d0
        end if
        tChunk = max( info%minChunkSize, tChunk)
        info%batchSize = tChunk
        info%batchRem = min(rem, tChunk)

      case (AF) ! adaptive factoring
        if (info%stats(3*rank) < 0.0) then
          tChunk = info%minChunkSize
        else
          bigD = 0.0d0
          bigT = 0.0d0
          do i=info%firstRank,info%lastRank
            bigD = bigD + info%stats(3*i+1)/info%stats(3*i)
            bigT = bigT + 1.0d0/info%stats(3*i)
          end do
          bigT = 1.0d0/bigT
          ! compute chunk size for rank
          tChunk = 0.55d0 + (0.5d0*(bigD + 2.0d0*bigT*rem -  &
                         sqrt(bigD*(bigD + 4.0d0*bigT*rem)))/info%stats(3*rank))
          tChunk = min( info%maxChunkSize, tChunk)
        end if
        tChunk = max( info%minChunkSize, tChunk)
        info%batchSize = tChunk
        info%batchRem = min( info%batchSize, rem)
     
!     case (10) ! experimental FSC
!       if (info%gN <= 0.0) then
!         tChunk = info%minChunkSize
!       else
!         tChunk = ((info%kopt0*info%gh/info%gsigma)**2)**(1.0d0/3.0d0)
!       end if
!       tChunk = max( info%minChunkSize, tChunk)
!       info%batchSize = tChunk
!       info%batchRem = min( info%batchSize, rem)

      case default ! Unsupported fall back to STATIC
        write (*,*) 'Unsupported DLS technique, fall back to STATIC'
        tChunk = (info%N+info%commSize-1)/info%commSize
        i = mod(info%N, info%commSize)
        if ((i>0) .and. (rank>=i) ) tChunk=tChunk-1
        tChunk = min( tChunk, rem)
        info%batchSize = tChunk
        info%batchRem = min( info%batchSize, rem)

    end select

    !write (*,*) rank, 'chunk size ', tChunk
    chunkSize = min(info%batchRem, tChunk)
!     write(*,*) '[GetChunkSize] rank: ',rank,' chunk',tChunk
!    write (*,*) '[Getchunk] tchunk = ', tChunk
!    write (*,*) '[getchunk] FSC chunk',  info%chunkMFSC

    ! adjust remaining in batch
    info%batchRem = info%batchRem - chunkSize
    if ( (info%batchRem > 0) .and. (info%batchRem <= info%minChunkSize) ) then
      chunkSize = chunkSize + info%batchRem
      info%batchRem = 0
    end if

    return
  end subroutine GetChunkSize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine SendChunk (info, worker)
    type (infoDLS), intent (in out) :: info
    integer, intent (in) :: worker

    integer :: chunkSize

    call GetChunkSize (info, worker, chunkSize) 

    ChunkInfo(0) = info%chunkStart      
    ChunkInfo(1) = chunkSize 

    !write(*,*) '[SendChunk] rank: ',worker,' start', ChunkInfo(0), 'chunk size', ChunkInfo(1)

    if(worker == info%Foreman)  then

          if (info%wSize == 0) then ! no pending chunk
            info%t0 = MPI_Wtime() ! elapsed time for chunk starts here
            info%wStart = ChunkInfo(0)
            info%wSize = ChunkInfo(1)
            info%rStart = info%wStart
            info%rSize = info%wSize
            info%req4WRKsent = .false.  ! cancel request for work

            call SetBreaks (info)

           !write(*,*) 'WRK_TAG recv by', info%myRank, info%wStart, info%wSize

            info%sumt1 = 0.0  ! for mu/wap
            info%sumt2 = 0.0  ! for sigma
          else  ! current chunk is not finished  save as next chunk
            info%nextStart = ChunkInfo(0)
            info%nextSize = ChunkInfo(1)
            info%nextWRKrcvd = .true.

           !write(*,*) 'WRK_TAG (adv) recv by', info%myRank,
           !info%nextStart,info%nextSize

          end if

      else
       call MPI_Send (ChunkInfo, 2, MPI_INTEGER, worker, WRK_TAG, info%comm,ierror)
       !call MPI_ISEND(ChunkInfo, 2, MPI_INTEGER, worker, WRK_TAG,
       !info%comm,!request,ierror)

     end if

!
!   write(*,*) 'WRK_TAG to', worker, ChunkInfo(0), ChunkInfo(1) !, &
!       info%batchSize, info%batchRem
!
    ! update counters 
    info%chunkStart     = info%chunkStart     + chunkSize 
    info%itersScheduled = info%itersScheduled + chunkSize 

    info%numChunks = info%numChunks + 1
    !info%chunkMap(2) = info%numChunks
    !info%chunkMap(3*info%numChunks  ) = ChunkInfo(0) 
    !info%chunkMap(3*info%numChunks+1) = ChunkInfo(1) 
    !info%chunkMap(3*info%numChunks+2) = worker 
    return
  end subroutine SendChunk 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DLS_EndChunk (info)
    type (infoDLS), intent (in out) :: info
    double precision :: tk
    integer :: loc, i, j 
  
    if (info%comm==MPI_COMM_NULL) return

    if (info%subChunkSize==0) return
    
    tk = MPI_Wtime()
    info%t1 = tk - info%t1 
    info%wStart = info%wStart + info%subChunkSize 
    info%wSize = info%wSize - info%subChunkSize 
    info%sumt1 = info%sumt1 + info%t1 
    info%workTime = info%workTime + info%t1 
    if (info%method == AF) then !  .or. info%method == 10) then 
      info%sumt2 = info%sumt2 + info%t1**2
    end if


    if (info%wSize == 0) then ! chunk finished 

      if (info%method==AWFB .or. info%method==AWFC) then
        ! adaptive weighted factoring, work time 
        info%mySumTimes = info%mySumTimes + (1+info%myExecs)*info%sumt1
        info%mySumSizes = info%mySumSizes + (1.0+info%myExecs)*info%rSize
      else if (info%method==AWFD .or. info%method==AWFE) then
        ! adaptive weighted factoring, elapsed time 
        info%mySumTimes = info%mySumTimes + (1+info%myExecs)*(tk-info%t0)
        info%mySumSizes = info%mySumSizes + (1.0+info%myExecs)*info%rSize
      end if

      if(info%method /= AF) then
      ! reset accumulators
      info%myIters = info%myIters + info%rSize
      info%myExecs = info%myExecs + 1
      info%sumt1 = 0.0  ! for mu 
      info%sumt2 = 0.0  ! for sigma 
      info%rSize = 0
      end if

      if (info%nextWRKrcvd) then ! foreman already responded to advance request 
        info%t0 = MPI_Wtime() ! elapsed time for chunk starts here 
        info%wStart = info%nextStart
        info%wSize = info%nextSize
        info%rStart = info%wStart
        info%rSize = info%wSize

        call SetBreaks (info)

        info%nextSize = 0
        info%nextWRKrcvd = .false.
        info%req4WRKsent = .false.
      end if
    end if ! if (info%wSize == 0) 

    ! send request ?
    if ( (info%wSize<=info%sendRequest) .and. (.not. info%req4WRKsent) ) then 
      if (info%method==STAT .or. info%method==MFSC .or. &
          info%method==SS .or. info%method==FSC .or. &
          info%method==TSS .or. info%method==WF .or. info%method==AWF .or. &
          info%method==GSS .or. info%method==FAC) then
        ! non-adaptive methods 
        PerfInfo(0) = 0.0 
        PerfInfo(1) = 0.0 
        PerfInfo(2) = 1.0*info%myRank  
      else if (info%method==AWFB .or. info%method==AWFC) then
        ! non-timestepping adaptive weighted factoring 
        ! mu = (chunk work time)/(chunk size) 
        PerfInfo(0) = ( info%mySumTimes + (info%myExecs+1)*info%sumt1 ) / &
            ( info%mySumSizes + 1.0*(info%myExecs+1)*(info%rSize-info%wSize) ) 
        PerfInfo(1) = 0.0 
        PerfInfo(2) = 1.0*info%myRank  
      else if (info%method==AF) then
        ! adaptive factoring 
        PerfInfo(0) = info%sumt1/(info%rSize-info%wSize)         ! mu 
        if ((info%rSize-info%wSize) > 1) then 
          PerfInfo(1) = (info%sumt2 - PerfInfo(0)*PerfInfo(0)*(info%rSize-info%wSize))/ &
            (info%rSize-info%wSize-1)  ! sigma 
          if (PerfInfo(1) < 0.0) PerfInfo(1) = 0.0 
          PerfInfo(1) = sqrt(PerfInfo(1))
        else 
          PerfInfo(1) = 0.0 
        end if
          if (info%wSize == 0) then
              ! reset accumulators
             info%myIters = info%myIters + info%rSize
             info%myExecs = info%myExecs + 1
             info%sumt1 = 0.0  ! for mu 
             info%sumt2 = 0.0  ! for sigma 
             info%rSize = 0
          end if
        PerfInfo(2) = 1.0*info%myRank  
      else if (info%method==AWFD .or. info%method==AWFE) then
        ! non-timestepping adaptive weighted factoring 
        ! mu = (chunk elapsed time)/(chunk size) 
        PerfInfo(0) = ( info%mySumTimes + (info%myExecs+1)*(tk-info%t0)) / &
            ( info%mySumSizes + 1.0*(info%myExecs+1)*(info%rSize-info%wSize) ) 
        PerfInfo(1) = 0.0 
        PerfInfo(2) = 1.0*info%myRank  
!     else if (info%method==10) then
!       PerfInfo(0) = info%sumt1 
!       PerfInfo(1) = info%sumt2 
!       PerfInfo(2) = info%rSize-info%wSize
!       PerfInfo(3) = (tk-info%t0) - info%sumt1  
      end if
!
!     write(*,*) 'REQ_TAG by', info%myRank, PerfInfo(0)
!
   
       if(info%myRank == info%Foreman)  then
          !update performance data
       if (info%method==AWFB .or. info%method==AWFC .or. info%method==AF .or.&
            info%method==AWFD .or. info%method==AWFE) then
            loc = int (PerfInfo(2))
            info%stats(3*loc+2) = info%stats(3*loc+2)+1.0
             ! adaptive methods
             info%stats(3*loc) = PerfInfo(0)
             info%stats(3*loc+1) = PerfInfo(1)
    !
             if (info%finishedOne /= info%commSize) then
              ! workers that have not finished a first chunk
              !  assume the lowest performance
              j = loc
              do i=info%firstRank,info%lastRank
                if ( (info%stats(3*i+2) > 0.0) .and. &
                     (info%stats(3*i) < info%stats(3*j)) ) j = i
              end do
              info%finishedOne = 0
              do i=info%firstRank,info%lastRank
                if (info%stats(3*i+2) == 0.0) then
                  info%stats(3*i) = info%stats(3*j)
                  info%stats(3*i+1) = info%stats(3*j+1)
                else
                   info%finishedOne = info%finishedOne + 1
                end if
              end do
             end if
           end if

! get more work to myself
          if (info%chunkStart < info%lastIter) then
!            print *, '[end chunk] obtaining work ... current chunk start',
!            info%chunkStart
             info%req4WRKsent = .true.
             info%nextWRKrcvd = .false.
            call SendChunk (info, info%myRank)
!            print *, "[end chunk] wsize: ", info%wSize

          else if (info%wsize == 0 .and. info%chunkStart >= info%lastIter) then
! all iterates scheduled
            info%numENDed = info%numENDed + 1
 !           print *, '[end chunk] ended ranks ', info%numENDed, 'chunk',
 !           info%subChunkSize
            info%gotWork = info%NumENDed/=info%commSize ! foreman exits?
           end if

       else
          call MPI_Send (PerfInfo, 4, MPI_DOUBLE_PRECISION, info%Foreman, &
          REQ_TAG, info%comm, ierror)
          info%req4WRKsent = .true.
          info%nextWRKrcvd = .false.
       end if !...if foreman

    end if ! if (...info%sendRequest...) 
    

    return
  end subroutine DLS_EndChunk


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DLS_Partition(world, coordinator, jcomm, icomm)
    implicit none
  
    integer, intent (in) :: world
    integer, intent (in) :: coordinator
    integer, intent (out) :: jcomm, icomm
  
    character (len = MPI_MAX_PROCESSOR_NAME), dimension(0:maxProcs-1) :: procName
    integer, dimension(MPI_STATUS_SIZE) :: tStatus
  
    integer, dimension(0:maxProcs) :: grpSize, grpIdx
    integer, dimension (0:maxProcs-1) :: foremanOf, localRank
    integer ierr, idx, i, j, mrg1, mrg2, minGrpRank, nGrps
    integer gRank, gSize
  
    call MPI_Comm_rank (world, gRank, ierr)
    call MPI_Comm_size (world, gSize, ierr)
    call MPI_Get_Processor_name (procName(gRank), i, ierr)
  
    if (gRank==coordinator) then
  !   write(*,*) coordinator, 'is', trim(procName(coordinator))
      do i = 1,gSize-1
        call MPI_Probe(MPI_ANY_SOURCE, HNM_TAG, world, tStatus, ierr)
        call MPI_Get_count (tStatus, MPI_CHARACTER, mrg1, ierr)
        mrg2 = tStatus(MPI_SOURCE)
        call MPI_Recv (procName(mrg2), mrg1, MPI_CHARACTER, &
          mrg2, HNM_TAG, world, tStatus, ierr)
  !     write(*,*) mrg2, 'is', trim(procName(mrg2))
      end do
      !
      foremanOf = -1
      nGrps = 0
      grpSize = 0
      do i = 0,gSize-1
        if (foremanOf(i)==-1) then
          nGrps = nGrps+1
          foremanOf(i) = nGrps
          grpSize(nGrps) = 1
  !       write (*,*) i, nGrps, grpSize(nGrps), trim(procName(i))
          idx = index(procName(i),"-")+3
          do j = i+1,gSize-1
            if (procName(j)(:idx)==procName(i)(:idx)) then
              ! impose max group size
              if (grpSize(nGrps)==maxGrpSize) then
                nGrps = nGrps+1
                grpSize(nGrps) = 1
              else
                grpSize(nGrps) = 1+grpSize(nGrps)
              end if
              foremanOf(j) = nGrps
  !           write (*,*) j, nGrps, grpSize(nGrps), trim(procName(j))
            end if
          end do
        end if
      end do
  !   impose min group size without violating max group size
      do
    !   sort groups according to size
        do i = 1,maxProcs
          grpIdx(i) = i
        end do
        do i = 1,nGrps-1
          do j = i+1,nGrps
            if (grpSize(grpIdx(i))<grpSize(grpIdx(j))) then
              ierr = grpIdx(i)
              grpIdx(i) = grpIdx(j)
              grpIdx(j) = ierr
            end if 
          end do
        end do
  !  
  !     write(*,*) nGrps, 'groupings...'
  !     do i = 1,nGrps
  !       write(*,*) 'Group', grpIdx(i), '; size = ', grpSize(grpIdx(i))
  !       do j = 0,gSize-1
  !         if (foremanOf(j)==grpIdx(i)) then
  !           write (*,*) j, trim(procName(j))
  !         end if
  !       end do
  !     end do
  !  
      !   find two smallest group
        mrg1 = grpIdx(nGrps)
        mrg2 = grpIdx(nGrps-1)
        if (grpSize(mrg1)>=minGrpSize) exit
        if ((grpSize(mrg1)+grpSize(mrg2))>maxGrpSize) exit
  !     write(*,*) 'Merge grp1:', mrg1, grpSize(mrg1), 'with grp2:', mrg2, grpSize(mrg2)
      !   merge groups
        do j = 0,gSize-1
          if (foremanOf(j)==mrg1) then
            foremanOf(j) = mrg2
          end if
        end do
        grpSize(mrg2) = grpSize(mrg1)+grpSize(mrg2)
      !   renumber groups
        do i = mrg1,nGrps-1
          do j = 0,gSize-1
            if (foremanOf(j)==i+1) foremanOf(j)=i
          end do
          grpSize(i) = grpSize(i+1)
        end do
        nGrps = nGrps-1
      end do

      write(*,*) nGrps, 'groupings...'
      do i = 1,nGrps
        write(*,*) 'Group', grpIdx(i), '; size = ', grpSize(grpIdx(i))
        do j = 0,gSize-1
          if (foremanOf(j)==grpIdx(i)) then
            write (*,*) j, trim(procName(j))
          end if
        end do
      end do

  ! assign lowest ranks in groups as foremen
      foremanOf(gRank) = 0
      do i = 1,nGrps
        minGrpRank = gSize
        idx = grpIdx(i)
        ierr = 0
        do j = 0,gSize-1
          if (foremanOf(j)==idx) then
            ierr = ierr+1
            if (j<minGrpRank) minGrpRank = j
          end if
        end do
  !     write(*,*) 'Group foremanOf is', minGrpRank, 'size is', ierr
        do j = 0,gSize-1
          if (foremanOf(j)==idx) then
            foremanOf(j) = gSize+minGrpRank
          end if
        end do
      end do
      foremanOf(0:gSize-1) = foremanOf(0:gSize-1)-gSize
      foremanOf(gRank) = 0
      do i = 0,gSize-1 ! send colors
        if (i==gRank) cycle
        call MPI_Send (foremanOf, gSize, MPI_INTEGER, i, CLR_TAG, world, ierr)
      end do
    else ! send proc name to coordinator, wait for group color
      idx = len_trim(procName(gRank))
      call MPI_Sendrecv (procName(gRank), idx, MPI_CHARACTER, &
           coordinator, HNM_TAG, foremanOf, gSize, MPI_INTEGER, &
           coordinator, CLR_TAG, world, tStatus, ierr)
    end if
  
    localRank = -1 ! map global ranks to local ranks
    do i = 0,gSize-1
      idx = foremanOf(i) ! increment no. of workers of foremanOf idx
      localRank(idx) = localRank(idx)+1 ! temporary counter; reset in next loop 
      localRank(i) = localRank(idx) ! this worker's rank
    end do
    do i = 0,gSize-1 ! reset local rank of a foremanOf to 0 
      if (foremanOf(i)==i) localRank(i)=0 
    end do
  
    ! split processors into work groups 
    if (gRank==coordinator) then
      idx = MPI_UNDEFINED
    else
      idx = foremanOf(gRank)
    end if
    call MPI_Comm_split(world, idx, gRank, icomm, ierr)
  
    ! make control group for the foremen and coordinator
    if (gRank==coordinator .or. foremanOf(gRank)==grank) then
      idx = 0
    else
      idx = MPI_UNDEFINED
    end if
    call MPI_Comm_split(world, idx, gRank, jcomm, ierr)
  
  ! if (icomm==MPI_COMM_NULL) then 
  !   call MPI_Comm_size (jcomm, mrg1, ierr)
  !   call MPI_Comm_rank (jcomm, mrg2, ierr) 
  !   write(*,*) gRank, ' is the coordinator,  ', mrg2, '/', mrg1, ' of jcomm'
  ! else if (jcomm==MPI_COMM_NULL) then 
  !   call MPI_Comm_size (icomm, j, ierr)
  !   call MPI_Comm_rank (icomm, i, ierr) 
  !   write(*,*) gRank, ' is ', i, '/', j, ' of local icomm'
  ! else
  !   call MPI_Comm_size (icomm, j, ierr)
  !   call MPI_Comm_rank (icomm, i, ierr) 
  !   call MPI_Comm_size (jcomm, mrg1, ierr)
  !   call MPI_Comm_rank (jcomm, mrg2, ierr) 
  !   write(*,*) gRank, ' is ', i, '/', j, ' of local icomm, ', mrg2, '/', mrg1, ' of jcomm'
  ! end if
  
    return
  end subroutine DLS_Partition

end module DLS      
