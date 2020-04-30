#include "Loop.hpp"

/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

Loop::Loop(int numIters)
{
    N = numIters;
    remainingIters = N;
    availableIters = new int[N]; 
}

Loop::Loop(int numIters, double mu, double sigma)
{
    meanItersExeTime = mu;
    stdItersExeTime = sigma;
    N= numIters;
}
         

Loop::~Loop()
{
   delete [] availableIters;  
   availableIters = nullptr;
}
  
