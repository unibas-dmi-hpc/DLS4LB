/*********************************************************************************
 * Copyright (c) 2007                                                            *
 * Ricolindo L. Carino and Ioana Banicescu                                       *
 * Mississippi State University, USA                                             *
 * All rights reserved.                                                          *
 *                                                                               *
 * Modified - June 2018                                                          *
 * Ali Mohammed <ali.mohammed@unibas.ch>                                         *
 * University of Basel, Switzerland                                              *
 *                                                                               *
 * This program is free software; you can redistribute it and/or modify it       *
 * under the terms of the license (GNU LGPL) which comes with this package.      *
 *********************************************************************************/
#ifndef DLB_H
#define DLB_H

#include <string.h>
#include <stdlib.h>
#define maxProcs    128
#define maxChunks   500
#define maxGrpSize    4
#define minGrpSize    2

#define HNM_TAG   9990
#define CLR_TAG   9980
#define WRK_TAG   9970
#define REQ_TAG   9960
#define TRM_TAG   9950
#define END_TAG   9940

#define min(a, b)       ((a) < (b) ? (a) : (b))
#define max(a, b)       ((a) > (b) ? (a) : (b))

#define STATIC  0
#define SS      1
#define FSC     2
#define mFSC    3
#define GSS     4
#define TSS     5
#define FAC     6
#define WF      7
#define AWF     8
#define AWF_B   9
#define AWF_C   10
#define AWF_D   11
#define AWF_E   12
#define AF      13
#define SimAS   14



typedef struct
{
  MPI_Comm comm, crew;
  int commSize, crewSize;
  int foreman, myRank, firstRank, lastRank;
  int method;
  int firstIter, lastIter, N, itersScheduled;
  int batchSize, batchRem, minChunkSize, maxChunkSize;
  int minChunk, breakAfter, requestWhen, chunkFSC, chunkMFSC;
  int chunkStart, probeFreq, sendRequest, subChunkSize;
  int numChunks, numENDed, finishedOne;
  int myExecs, myIters;
  int rStart, rSize, wStart, wSize, nextStart, nextSize;
  int gotWork, req4WRKsent, nextWRKrcvd;
  double kopt0, workTime;
  double t0, t1, sumt1, sumt2, mySumTimes, mySumSizes;
  double *stats;
  double h_overhead;
  double sigma;
  double *weights;
  int TSSchunk;
  int TSSdelta;
  int timeStep;
//  int chunkMap[3*maxChunks];
  } infoDLS;


char * get_DLS_ShortName(int DLS);
void DLS_Setup ( MPI_Comm, infoDLS * );
//void DLS_Parameters_Setup ( MPI_Comm, infoDLS *, int numProcs, int requestWhen, int breakAfter, int minChunk );
void DLS_Parameters_Setup( MPI_Comm icomm, infoDLS *info, int numProcs, int requestWhen, int breakAfter, int minChunk, double h_overhead, double sigma, int nKNL, double Xeon_speed, double KNL_speed);
void DLS_GroupSetup ( MPI_Comm, int, infoDLS *, infoDLS *);
void DLS_StartLoop ( infoDLS *, int, int, int);
int  DLS_Terminated ( infoDLS *);
void DLS_StartChunk (infoDLS *, int *, int *);
void DLS_EndChunk ( infoDLS *);
void DLS_EndLoop (infoDLS *, int *, double *);
void DLS_Finalize (infoDLS *info);
void SetBreaks(infoDLS *info);
#endif 
