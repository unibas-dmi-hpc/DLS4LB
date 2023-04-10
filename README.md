This branch provides the python interface of the library.  The interface relies on mpi4py. So you needd to install it before using this interface. 
Example of how one can use the library can be found in parallelpyth.py. This code shows an example of the Mandelbrot kernel written with the python interface of the library
For more details one can read the following report
A. Mendhe, "Dynamic loop self-scheduling on distributed-memory systems with Python", University of Basel, October 2021
https://hpc.dmi.unibas.ch/en/scientific_output/completed_theses_and_projects/
In general one can do the following
```
from mpi4py import MPI
import DLS as cdls
cdls.infoDLS(..) #here one should pass all the configuration parameters to the library, e.g., comm_world, size, rank, startloop_index, number of iterations, master, scheduling method, ...
cdls.StartLoop(info)
while(not cdls.DLS_Terminated(info)):
    start, chunks = cdls.DLS_StartChunk(info)
    if start < (start+chunks):
      \\do something
    cdls.DLS_EndChunk(info)    
cdls.DLS_EndLoop(info) 
}
```
Please be aware that this branch is still under development and the authors of this extension are Abhilash Mendhe, Ahmed Eleliemy, and Florina M. Ciorba.
For bugs, we encourage everyone to create issues on the repository.
For further questions, please contact ahmed.eleliemy@unibas.ch or florina.ciorba@unibas.ch
