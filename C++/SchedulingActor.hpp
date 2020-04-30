#pragma once

/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

#include "Chunk.hpp"
#include "Loop.hpp"
#include "DLS.hpp"

class SchedulingActor
{
    public:
        virtual void startLoop(DLS *, int, Loop *) =0;
	virtual Chunk startChunk() = 0;
	virtual void endChunk() = 0;
	virtual void endLoop() = 0;
	virtual void finalize()  = 0;
        virtual ~SchedulingActor() {}
};
