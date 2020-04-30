/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

#include "Worker.hpp"


 
Worker::Worker(int rank, MPI_Comm comm)
{
    this->rank = rank;
    commWorld = comm;
    
}


void Worker::startLoop(DLS *method, int requestWhen, Loop *cLoop)
{
 schMethod = method;
this->requestWhen = requestWhen;
}

Chunk Worker::startChunk()
{}

void Worker::endChunk()
{}
	
void Worker::endLoop()
{}
	
void Worker::finalize()
{}

Worker::~Worker()
{
}  
