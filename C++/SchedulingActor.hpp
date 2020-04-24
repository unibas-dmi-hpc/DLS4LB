#pragma once

#include "Chunk.hpp"

class SchedulingActor
{
    public:
        virtual void startLoop(int Nitrs, int SchMethod) =0;
	virtual Chunk startChunk() = 0;
	virtual void endChunk() = 0;
	virtual void endLoop() = 0;
	virtual void finalize()  = 0;
};
