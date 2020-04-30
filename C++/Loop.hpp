#pragma once

/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

#include "Chunk.hpp"

struct Loop
{

int N; //Number of iterations
// Available iterations for each rank ...the ones that ranks have its data .
//In replicated or centralized data, all iterations will be in rank 0, the master.
int *availableIters;  
int remainingIters; //Number of remaining iterations
double meanItersExeTime;
double stdItersExeTime;

void ** data; // pointer to data related to this chunk

int * dataSources; // list of ranks who has each element of data

int * dataOffsets; // list of offsets where to start reading data in each data element

int * dataSizes; //list of how many elements to read from each data element for this chunk

Loop(int numIters);
Loop(int numIters, double mu, double sigma);

//allocate a chunk with "size" for the requesting worker with rank - takes into account the data distribution, i.e., allocates a chunk
//of the required size that the requesting worker preferably has its data
Chunk allocateChunk(int size, int rank);
~Loop();

};
