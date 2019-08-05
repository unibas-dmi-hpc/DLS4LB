
#ifndef SiL_H
#define SiL_H

#include <stdio.h>
#include "mpi.h"
#include "DLS.h"



typedef struct
{

 int SiL_flag;
 double SiL_t1, SiL_t2;
 char *tasks_file;
 char *platform;
 int SiL_last_t;
 int sim_finished;
 char *sim_path;
 double last_check; 
 int default_DLS; // default DLS
 double period; // simulation period
} infoSiL;


void SiL_sim(char *sim_path, char *platform, int nhosts,int ntasks,char *flops_file,double h, double sigma,int start_task, double start_time, double period, int id);
int SiL_select(char *platform, int nhosts,int ntasks, int id);
//void SiL_setup(infoDLS *info, infoSiL *sil_data, char *plaformfile, char *task_file);
int SiL_setup(infoSiL *sil_data, char *sim_path, char *plaformfile, char *task_file, int method, int P, int N, double h, double sigma);
int SiL_setup_full(infoSiL *sil_data, char *sim_path, char *plaform_file, char *tasks_file, int method, int P, int N, double h, double sigma, int default_DLS, double period);
void SiL_update(infoDLS *iInfo, infoSiL *sil_data);


#endif
