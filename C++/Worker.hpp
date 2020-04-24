#pragma once

#include "SchedulingActor.hpp"
#include "Chunk.hpp"

class Worker: SchedulingActor
{

int rank;
double mu; //mean execution time
double sigma; //standard deviation of execution time
double weight; //relative weight

public:
   Worker(int rank);
   Worker(int rank, double weight);
};
