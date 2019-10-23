# DLS4LB
Dynamic Loop Self-scheduling For Load Balancing (DLS4LB) is an MPI-Based load balancing library. It is implemented in C and FORTRAN (F90) programming languages to support scientific applications executed on High Performance Computing (HPC) systems.

DLS4LB library is based on the DLB_tool developed by Ricolindo L. Carino (rlc@erc.msstate.edu) and Ioana Banicescu (ioana@cse.msstate.edu), see publication [1].
It is modified and extended by Ali Mohammed (ali.mohammed@unibas.ch) to support more scheduling techniques and some bug fixes, see publication [2].
Also, DLS4LB is extended to support Simulation assisted scheduling Algorithm Selection (SimAS) as the fifteen DLS technique. 
That is, the most efficient DLS technique will be selected dynamically during execution based on simulation results. Please read publication [2] for more details.

The DLS4LB parallelizes and load balances scientific applications that contain simple parallel loops (1D loops) or nested parallel loops (2D loops). The tool employs a master-worker model where workers request work from the master whenever they become free. The master serves work requests and assigns workers chunks of loop iterations according to the selected DLS technique. The master also doubles as a worker and executes chunks of loop iterations when it is not serving any requests.


The DLS4LB library supports fourteen scheduling techniques:
  1.  Straightforward static scheduluing (STATIC)
  2.  Self Scheduling (SS)
  3.  Fixed size chunk (FSC)
  4.  Modified fixed size chhunk (mFSC)
  5.  Guided self-scheduling (GSS) 
  6.  Trapezoid self-scheduling (TSS)
  7.  Factoring (FAC)
  8.  Weighted factoring (WF)
  9.  Adaptive weighted factoring (AWF)
  10. Adaptive weighted factoring - Batch (AWF-B)
  11. Adaptive weighted factoring - Chunk (AWF-C)
  12. Adaptive weighted factoring - Batch with scheduling overhead (AWF-D)
  13. Adaptive weighted factoring - Chunk with scheduling overhead (AWF-E)
  14. Adaptive factoring (AF)
  
 
The DLS4LB is designed to load balance scientific applications with minimum code changes. Below is a simple example on how to use the DLS4LB tool for simple 1D loops (Type I) and nested 2D loops (Type II).
  
Type I
   ! begin i-loop
   do i=1,N
     ... i-iterate
   end do
   ! end i-loop


 Type II
   ! begin j-loop
   do j=1,M
     ... part of j-iterate
     ! begin i-loop
     do i=1,N(j)
     ... i-iterate of j-iterate
     end do
     ! end i-loop
     ... part of j-iterate
   end do
   ! end j-loop

 Dynamic load balancing is achieved in the application
 by invoking DLS routines to manipulate the loop indices.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 For Type I loops, the application must be modified as:

   use DLS
   include 'mpif.h'
   ...
   type (infoDLS) info
   integer method, iStart, iSize
   integer iIters
   double precision iTime
   ...
   method = ... (choice of loop scheduling technique)
   call DLS_Setup (MPI_COMM_WORLD, info)
   call DLS_StartLoop (info, 1, N, method)
   do while ( .not. DLS_Terminated(info) )
     call DLS_StartChunk (info, iStart, iSize)
     ! begin i-loop, new extents
     do i=iStart, iStart+iSize-1    ! 1,N
       ... i-iterate
     end do
     ! end i-loop
     call DLS_EndChunk (info)
   end do ! while ( .not. DLS_Terminated(info) )
   call DLS_EndLoop (info, iIters, iTime)


 where:

   DLS_Setup (MPI_COMM_WORLD, info)
      Initializes a dynamic load balancing environment on
      MPI_COMM_WORLD. Information about this environment is
      stored in "info" (see type declaration infoDLS below)

   DLS_StartLoop (info, 1, N, method)
      The synchronization point to start loop execution.
      (1, N) is the loop range, and "method" is a user-
      specified index [1..9] for the loop scheduling method
      (see array DLS_LongName() below)

   DLS_Terminated(info) 
      returns .TRUE. if all loop iterates have been
      executed

   DLS_StartChunk (info, iStart, iSize)
      Returns a range for a chunk of iterates, starting
      at "iStart", with size "iSize"

   do i=iStart, iStart+iSize-1    ! 1,N
     ... i-iterate
   end do
      The original "serial" code, but now executing
      a chunk of iterates instead of "i=1,N"

   DLS_EndChunk (info)
      Indicates the end of execution of a chunk of iterates

   DLS_EndLoop (info, iIters, iTime)
      The synchronization point to end loop execution. 
      "iIters" is the number of iterates done by the calling
      processor, and "iTime" is its cost (in seconds)
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 For Type II loops, the application must be modified as:

   use DLS
   include 'mpif.h'
   ...
   type (infoDLS) iInfo, jInfo
   integer iComm, jComm    ! communicators
   integer iMethod, iStart, iSize, jMethod, jStart, jSize
   integer iIters, jIters
   double precision iTime, jTime
   integer coordinator
   ...
   iMethod = ... (choice of loop scheduling technique)
   jMethod = ... (choice of loop scheduling technique)
   coordinator = 0
   call DLS_GroupSetup (MPI_COMM_WORLD, coordinator, jInfo, iInfo)
   call DLS_StartLoop (jInfo, 1, M, jMethod)
   do while ( .not. DLS_Terminated(jInfo) )
     call DLS_StartChunk (jInfo, jStart, jSize)

   ! begin j-loop code, new extents
     do j=jStart, jStart+jSjze-1    ! 1,M
       ... part of j-iterate 

   call DLS_StartLoop (iInfo, 1, N(j), iMethod)
       do while ( .not. DLS_Terminated(iInfo) )
         call DLS_StartChunk (iInfo, iStart, iSize)
         ! begin i-loop, new extents
         do i=iStart, iStart+iSize-1    ! 1,N
           ... i-iterate of j-iterate
         end do
       end i-loop
       call DLS_EndChunk (iInfo)
      end do ! while ( .not. DLS_Terminated(iInfo) )
      call DLS_EndLoop (iInfo, iIters, iTime)

   ... part of j-iterate 
     end do
   ! end j-loop

  call DLS_EndChunk (jInfo)
   end do ! while ( .not. DLS_Terminated(jInfo) )
   call DLS_EndLoop (jInfo, jIters, jTime)

 where:

   DLS_GroupSetup (MPI_COMM_WORLD, coordinator, jInfo, iInfo)
      Splits MPI_COMM_WORLD into non-overlapping "iComm"s and
      a "jComm" comprised of the "coordinator" and the
      foremen (rank 0 of the iComms).  Similar to DLS_Setup(), 
      DLS_GroupSetup() initializes load balancing environments 
      in "jComm" and "iComm", keeping the data in "jInfo" and 
      "iInfo" respectively. Graphically:


  ==============================================
  = MPI_COMM_WORLD        .-------------.      =           
  =    .-----------.       \  x  x   x /       =      
  =     \ iComm x /         \ iComm   /        =      
  =      \x  x x /           \     x /         =      
  =       \ x   /             \x x  /          =      
  =   .----\---/---------------\---/------.    =
  =   |     \x/     coordinator \x/       |    =
  =   |      .          x        .        |    =
  =   |           .                 .     |    =
  =   | jComm    /x\               /x\    |    =
  =   .---------/---\-------------/---\---.    =
  =            /iComm\           / x   \       =           
  =           / x  x  \         / x  x  \      =           
  =          -----------       / iComm x \     =           
  =                           .-----------.    =           
  ==============================================

  The symbol "x" denotes a processor, "x"s inside
  a triangle make up an "iComm", while "x"s inside
  the rectangle make up the jComm. The j-iterates 
  scheduled in "jComm" while the i-iterates for a
  given j-iterate are scheduled in "iComm"

Publications:
============
[1] Ricolindo L. Cariño and Ioana Banicescu. A Tool for a Two- level Dynamic Load Balancing Strategy in Scientific Applica- tions. Scalable Computing: Practice and Experience, 8(3), 2007.
[2] Ali Mohammed and Florina M. Ciorba. Research Report - University of Basel, Switzerland. https://drive.switch.ch/index.php/s/imo0kyLe3PmETWL. [Online; accessed 20 May 2019]. October 2018.
