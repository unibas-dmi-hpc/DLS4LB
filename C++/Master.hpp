#pragma once

#include "SchedulingActor.hpp"
#include "Chunk.hpp"

class Master: SchedulingActor
{

int rank;
double mu; //mean execution time
double sigma; //standard deviation of execution time
double weight; //relative weight


int calculateChunk();

public:
   Master(int rank);
   Master(int rank, double weight);
};
