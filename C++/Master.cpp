#include "Master"


class Master : SchedulingActor
{

    int calculateChunk()
    {
    }

    public: 
        Master(int id)
        {
          rank = id; //set rank own id
        }
      
        Master(int id, double relWeight)
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
