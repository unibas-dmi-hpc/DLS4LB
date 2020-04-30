#pragma once
/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

#include "SchedulingActor.hpp"
#include "mpi.h"

class Master : SchedulingActor
{

int rank;
int nWorkers; //number of workers
MPI_Comm commWorld; // MPI comm 
int probFrequency; //how often master check for requests
DLS *schMethod;
Loop *parLoop;


public:
   int totExeIters; //total executed iterations
   int totExeTime; //total execution time
   Master(int rank=0, int commSize=1, MPI_Comm comm=MPI_COMM_WORLD,int probFrequency=10);
   /* Should be called before the start of the loop, master needs to know loop properities, DLS method, 
     minChunkSize, TSSChunkSize, and Scheduling overhead */ 
   void startLoop(DLS *method, int probFrequency, Loop *cLoop);
   Chunk startChunk();
   void endChunk();
   void endLoop();
   void finalize();

   ~Master();
};

