/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

#include "Master.hpp"

 
Master::Master(int rank, int commSize, MPI_Comm comm,int probFrequency)
{
    this->rank = rank;
    nWorkers = commSize;
    commWorld = comm;
    this->probFrequency = probFrequency; 
}

Master::~Master()
{}


void Master::startLoop(DLS *method, int probFrequency, Loop *cLoop)
{
    parLoop = cLoop;
    schMethod = method;
    this->probFrequency = probFrequency;
}

Chunk Master::startChunk()
{}

void Master::endChunk()
{}
	
void Master::endLoop()
{}
	
void Master::finalize()
{}
  
