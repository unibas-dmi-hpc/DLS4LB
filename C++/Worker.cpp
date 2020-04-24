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
