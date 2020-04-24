#pragma once

/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/


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
