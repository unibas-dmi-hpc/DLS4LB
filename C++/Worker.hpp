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


class Worker: SchedulingActor
{

int rank;
double mu; //mean execution time
double sigma; //standard deviation of execution time
double weight; //relative weight
MPI_Comm commWorld; // MPI comm 
int requestWhen;  //when to request new chunk, before the end of the current chunk
DLS *schMethod;

public:
   int totExeIters; //total executed iterations
   int totExeTime; //total execution time
   Worker(int rank=1, MPI_Comm comm=MPI_COMM_WORLD);
   virtual void startLoop(DLS *method, int requestWhen=10, Loop *cLoop = nullptr);
   virtual Chunk startChunk();
   virtual void endChunk();
   virtual void endLoop();
   virtual void finalize();
   
   ~Worker();

};
