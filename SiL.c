#include "SiL.h"




void SiL_update(infoDLS *info, infoSiL *sil_data)
{

int myid;
int selected_dls = 1;
MPI_Comm_rank(MPI_COMM_WORLD, &myid);

if(myid == 0) // I am the master
{
 //printf("master entering update\n");
	if(sil_data->SiL_flag)
        {    
            if((info->N - info->itersScheduled) <= info->commSize) // almost finished .. do not simulate anymore
            {
               return ;
            }
            sil_data->SiL_t2 = MPI_Wtime();
            if((sil_data->SiL_t2 - sil_data->last_check) < 5)
            {
              return;
            }
            if(!sil_data->sim_finished) // if last time I checked the simulation was not finished
            {
               
               selected_dls = SiL_select(sil_data->platform,info->commSize,info->N,sil_data->SiL_last_t);
               if(selected_dls == -1) // simulation is not finsihed
               {
                  sil_data->sim_finished = 0; //simulation is not finsihed
               }
               else
               {
                  sil_data->sim_finished = 1;
                  info->method = selected_dls;
                  printf("selected dls, %d \n", selected_dls);
               }
               sil_data->last_check = sil_data->SiL_t2;
            }
            else if((sil_data->SiL_t2 - sil_data->SiL_t1) >= sil_data->period)
             {
                 sil_data->SiL_t1 = sil_data->SiL_t2;
                 sil_data->SiL_last_t++;
              	 SiL_sim(sil_data->sim_path, sil_data->platform,info->commSize,info->N,sil_data->tasks_file,info->h_overhead,info->sigma,info->itersScheduled,sil_data->SiL_t2 - sil_data->SiL_t1,sil_data->period,sil_data->SiL_last_t);
                 sil_data->sim_finished = 0;
                 sil_data->last_check = sil_data->SiL_t2;
	     }


        }
}
}




// returns the selected method (SS is default)
int SiL_setup(infoSiL *sil_data, char *sim_path, char *plaform_file, char *tasks_file, int method, int P, int N, double h, double sigma)
{

int selected_method = method;
int myid;

MPI_Comm_rank(MPI_COMM_WORLD, &myid);

 if(method == SiL)
 {
     sil_data->SiL_flag = 1;
     sil_data->sim_path = sim_path;
     sil_data->default_DLS = 6;
     selected_method = sil_data->default_DLS; //default
     sil_data->period = 50; //simulation period
     sil_data->platform = plaform_file;
     sil_data->tasks_file = tasks_file;

     if(myid== 0) // I am the master
     {
       // start SiL timer and start simulations
       sil_data->SiL_t1 =  MPI_Wtime();
       sil_data->SiL_last_t = 0;
       sil_data->last_check =  sil_data->SiL_t1;
       SiL_sim(sim_path, sil_data->platform,P,N,sil_data->tasks_file,h,sigma,0,0,sil_data->period, sil_data->SiL_last_t); // simulate at time zero, no tasks scheduled to the end of the execution
       sil_data->sim_finished = 0;
       printf("selected dls, %d \n", selected_method);
     }
 }
 else
 {
    sil_data->SiL_flag = 0; // disable SiL
 }

return selected_method;
}




int SiL_setup_full(infoSiL *sil_data, char *sim_path, char *plaform_file, char *tasks_file, int method, int P, int N, double h, double sigma, int default_DLS, double period)
{

int selected_method = method;
int myid;

MPI_Comm_rank(MPI_COMM_WORLD, &myid);

 if(method == SiL)
 {
     sil_data->SiL_flag = 1;
     sil_data->sim_path = sim_path;
     sil_data->default_DLS = default_DLS;
     selected_method = sil_data->default_DLS; //default
     sil_data->period = period; //simulation period
     sil_data->platform = plaform_file;
     sil_data->tasks_file = tasks_file;

     if(myid== 0) // I am the master
     {
       // start SiL timer and start simulations
       sil_data->SiL_t1 =  MPI_Wtime();
       sil_data->SiL_last_t = 0;
       sil_data->last_check =  sil_data->SiL_t1;
       SiL_sim(sim_path, sil_data->platform,P,N,sil_data->tasks_file,h,sigma,0,0,sil_data->period, sil_data->SiL_last_t); // simulate at time zero, no tasks scheduled to the end of the execution
       sil_data->sim_finished = 0;
       printf("selected dls, %d \n", selected_method);
     }
 }
 else
 {
    sil_data->SiL_flag = 0; // disable SiL
 }

return selected_method;
}




void SiL_sim(char *sim_path, char *platform, int nhosts,int ntasks,char *flops_file,double h, double sigma,int start_task, double start_time, double period, int id)
{
   int selected_dls = -1;
   int i;
   char cmd[1000];
   FILE *inputFile;
   int bufferSize = 500;
   char bufferLine[bufferSize];
   int method;
   int max_finished = 0;
   double min_time = 0;
   int DLS_portfolio[8];
    DLS_portfolio[0] = 6; //WF
    DLS_portfolio[1] = 2; //FSC
    DLS_portfolio[2] = 7; //AWF-B
    DLS_portfolio[3] = 8; //AWF-C
    DLS_portfolio[4] = 9; //AWF-D
    DLS_portfolio[5] = 10; //AWF-E
    DLS_portfolio[6] = 11; //AF
    DLS_portfolio[7] = 12; //mFSC


  char filename[500];
   //char filename[]="tmp.XXXXXX";

   sprintf(filename,"%s_%d_%d_%d",platform,nhosts,ntasks,id);
   //tmpnam(filename);
   //inputFile = fdopen(mkstemp(filename),"w+");
   //printf("%s\n",filename);

   //call inloop simulator for all the DLS 
   for(i=0; i<8 ; i++)
   {
    sprintf(cmd, "taskset -c %s %s %s %d 1 %d %d %s %lf %lf %d %lf %lf --log=root.threshold:critical >> %s &", "19-20", sim_path,platform,nhosts,ntasks,DLS_portfolio[i],flops_file,h, sigma,start_task, start_time, period, filename);
    //printf("%s \n", cmd);
    system(cmd);
   }
/*
   //compare the results 
   inputFile = fopen( filename , "r");
   if (inputFile == NULL)
   {
      printf("Error: Can not open file: %s \n",filename);
      return 1;
   }

    for (i = 0; i < 10; i++)
    {
      fgets(bufferLine, 500, inputFile);
      //printf("%d. I read: %s\n", i, bufferLine);
      sscanf(bufferLine, "Method, %d,  Parallel time, %lf, scheduled tasks, %d, finished tasks, %d",&method,&parallel_time[i],&scheduled[i],&finsihed_tasks[i]);
      //printf("%d. I read: %s, finish percent: %lf \n", i, bufferLine,(double) finsihed_tasks[i]/ (double) scheduled[i]);
    }
        
*/
return;

}

int SiL_select(char *platform, int nhosts,int ntasks, int id)
{

int num_DLS = 8;
FILE *inputFile;
char filename[500];
int bufferSize = 500;
char bufferLine[bufferSize];
int finsihed_tasks[num_DLS];
double parallel_time[num_DLS];
int method[num_DLS];
int scheduled[num_DLS];
int max_finished = 0;
int i;
double min_completion_ratio = 0.8; // minimum acceptable ratio
double max_ratio  = 0.0;

double min_time = 0;
double completion_ratio[num_DLS];
int selected_dls = -1;
char cmd[1000];

sprintf(filename,"%s_%d_%d_%d",platform,nhosts,ntasks,id);
//printf("trying to read file %s \n", filename);
inputFile = fopen( filename , "r");


   if (inputFile == NULL)
   {
      printf("Error: Can not open file: %s \n",filename);
      printf("returned 1 on error \n");
      return -1;
   }

    for (i = 0; i < num_DLS; i++)
    {
      fgets(bufferLine, 500, inputFile);
      //printf("%d. I read: %s\n", i, bufferLine);
      sscanf(bufferLine, "Method, %d,  Parallel time, %lf, scheduled tasks, %d, finished tasks, %d",&method[i],&parallel_time[i],&scheduled[i],&finsihed_tasks[i]);
//      printf("%d. I read: %s, finish percent: %lf \n", i, bufferLine,(double) finsihed_tasks[i]/ (double) scheduled[i]);
    }
 
 //calculate ratios
   for (i =0; i< num_DLS; i++)
   {
     completion_ratio[i] = (double) finsihed_tasks[i]/ (double) scheduled[i];
     if((method[i] > 12) || (method[i] < 1))
     {return -1;} // simulation is not finished ... you are reading garbage
     if(completion_ratio[i] >  max_ratio)
     {
       max_ratio = completion_ratio[i];
     }
   }

   if( min_completion_ratio > max_ratio)
   {
     min_completion_ratio = max_ratio;
   }

   //compare
   for (i = 0; i < num_DLS; i++)
   {
     if(completion_ratio[i] < min_completion_ratio)
     {
        continue; // do not consider this technique
     }
     if(finsihed_tasks[i] > max_finished)
     {
        max_finished = finsihed_tasks[i];
        min_time = parallel_time[i];
        selected_dls = i;
     }
     if(finsihed_tasks[i] ==  max_finished)
     {
        if(parallel_time[i] < min_time)
        {
          min_time = parallel_time[i];
          selected_dls = i;
        }
     }

   }
  // printf("selected dls, %d \n", method[selected_dls]);

    // delete the file to clean up
    sprintf(cmd, "rm  %s", filename);
    system(cmd);

return method[selected_dls];
}

