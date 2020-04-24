/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

#include "Worker"


class Worker : SchedulingActor
{

    public: 
        Worker(int id)
        {
          rank = id; //set rank own id
        }
      
        Worker(int id, double relWeight)
        {
          rank = id;
          weight = relWeight;
        }

        void startLoop(int Nitrs, int SchMethod)
        {
        }

	Chunk startChunk()
        {
        }
	void endChunk()
        {
        }
	
        void endLoop()
        {
        }
	
        void finalize()
        {
        }
  
};
