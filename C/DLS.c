/*********************************************************************************
 * Copyright (c) 2007                                                            *
 * Ricolindo L. Carino and Ioana Banicescu                                       *
 * Mississippi State University, USA                                             *
 * All rights reserved.                                                          *
 *                                                                               *
 * Modified - Oct. 2018                                                          *
 * Ali Mohammed <ali.mohammed@unibas.ch>                                         *
 * University of Basel, Switzerland                                              *
 *                                                                               *
 * This program is free software; you can redistribute it and/or modify it       *
 * under the terms of the license (GNU LGPL) which comes with this package.      *
 *********************************************************************************/

#include <stdio.h>
#include <math.h>
#include "mpi.h"

#include "DLS.h"



int DLS_MethodCount = 15;

char *DLS_ShortName [] = {
     "STATIC", "SS", "FSC", "mFSC", "GSS", "TSS","FAC", "WF", "AWF",
     "AWF-B", "AWF-C", "AWF-D", "AWF-E", "AF", "SimAS"
     };

char *DLS_LongName [] = {
     "STATIC SCHEDULING",
     "SELF-SCHEDULING",
     "FIXED SIZE-SCHEDULING", "MODIFIED FSC", "GUIDED SELF-SCHEDULING",
     "Trapezoid Self-Scheduling","FACTORING",              "WEIGHTED FACTORING", "ADAPTIVE WEIGHTED FACTORING",
     "BATCH AWF",              "CHUNK AWF",
     "BATCH AWF (chunk times)" ,"CHUNK AWF (chunk times)", "ADAPTIVE FACTORING",  "Simulation-assisted"
     };


char * get_DLS_ShortName(int DLS)
{

return DLS_ShortName[DLS];
}

void Partition (MPI_Comm world, int coordinator,
MPI_Comm *jcomm, MPI_Comm *icomm)
{
    char procName[maxProcs][MPI_MAX_PROCESSOR_NAME];
    MPI_Status tStatus;

    int grpSize[maxProcs+1], grpIdx[maxProcs+1];
    int foremanOf[maxProcs], localRank[maxProcs];
    int ierr, idx, i, j, mrg1, mrg2, minGrpRank, nGrps;
    int gRank, gSize;

    MPI_Comm_rank (world, &gRank);
    MPI_Comm_size (world, &gSize);
    MPI_Get_processor_name (procName[gRank], &i);

    if (gRank==coordinator) {
        /* printf("%d - %s is coordinator\n", gRank, procName[coordinator]); */
        for (i=1;i<gSize;i++) {
            MPI_Probe(MPI_ANY_SOURCE, HNM_TAG, world, &tStatus);
            MPI_Get_count (&tStatus, MPI_CHAR, &mrg1);
            mrg2 = tStatus.MPI_SOURCE;
            MPI_Recv (procName[mrg2], mrg1, MPI_CHAR,
                mrg2, HNM_TAG, world, &tStatus);
            /* printf("%d - %s\n", mrg2, procName[mrg2]); */
        }

        for (i=0;i<gSize;i++) {
            foremanOf[i] = -1;
            grpSize[i] = 0;
        }
        nGrps = 0;
        for (i=0;i<gSize;i++) {
            if (foremanOf[i]==-1) {
                nGrps = nGrps+1;
                foremanOf[i] = nGrps;
                grpSize[nGrps] = 1;
                /* printf("i=%d nGrps=%d size=%d name=%s\n", 
                       i, nGrps, grpSize[nGrps], procName[i]); */
                /* idx = strchr(procName[i],"-")+3; */
                idx = strstr(procName[i],".");
                for (j=i+1;j<gSize;j++) {
                    if (0 == strncmp(procName[j], procName[i], idx)) {
                        /* impose max group size */
                        if (grpSize[nGrps]==maxGrpSize) {
                            nGrps = nGrps+1;
                            grpSize[nGrps] = 1;
                        }
                        else grpSize[nGrps] = 1+grpSize[nGrps];
                        foremanOf[j] = nGrps;
                        /* printf("i=%d nGrps=%d size=%d name=%s\n",
                             j, nGrps, grpSize[nGrps], procName[j]); */
                    }
                }
            }
        }

        /* impose min group size without violating max group size */
        while (1) {
            /* sort groups according to size */
            for (i=1;i<maxProcs;i++) grpIdx[i] = i;
            for (i=1;i<nGrps;i++) {
                for (j=i+1;j<=nGrps;j++)
                    if (grpSize[grpIdx[i]]<grpSize[grpIdx[j]]) {
                        ierr = grpIdx[i];
                        grpIdx[i] = grpIdx[j];
                        grpIdx[j] = ierr;
                    }
            }
            /*
            printf("%d groupings\n",nGrps);
            for (i=1;i<=nGrps;i++) {
                printf("Group %d has size %d\n", grpIdx[i], grpSize[grpIdx[i]]);
                for (j=0;j<gSize;j++)
                    if (foremanOf[j]==grpIdx[i])
                        printf("%d - %s\n", j, procName[j]);
            }
            */

            /* find two smallest group */
            mrg1 = grpIdx[nGrps];
            mrg2 = grpIdx[nGrps-1];
            if (grpSize[mrg1]>=minGrpSize) break;
            if ((grpSize[mrg1]+grpSize[mrg2])>maxGrpSize) break;

            /*
            printf("Merge grp1: (%d,%d) with grp2: (%d,%d)\n",
                mrg1, grpSize[mrg1], mrg2, grpSize[mrg2]); */
            /* merge groups */
            for (j=0;j<gSize;j++)
                if (foremanOf[j]==mrg1) foremanOf[j] = mrg2;
            grpSize[mrg2] = grpSize[mrg1]+grpSize[mrg2];
            /* renumber groups */
            for (i=mrg1;i<nGrps;i++) {
                for (j=0;j<gSize;j++)
                    if (foremanOf[j]==(i+1)) foremanOf[j]=i;
                grpSize[i] = grpSize[i+1];
            }
            nGrps = nGrps-1;
        }

        /* assign lowest ranks in groups as foremen */
        foremanOf[gRank] = 0;
        for (i=1;i<=nGrps;i++) {
            minGrpRank = gSize;
            idx = grpIdx[i];
            ierr = 0;
            for (j=0;j<gSize;j++) {
                if (foremanOf[j]==idx) {
                    ierr = ierr+1;
                    if (j<minGrpRank) minGrpRank = j;
                }
            }
            /* printf("Group %d, foreman=%d, size=%d\n", idx, minGrpRank, ierr); */
            for (j=0;j<gSize;j++) {
                if (foremanOf[j]==idx) foremanOf[j] = gSize+minGrpRank;
            }
        }
        for (j=0;j<gSize;j++)
            foremanOf[j] = foremanOf[j]-gSize;
        foremanOf[gRank] = 0;
        for (i=0;i<gSize;i++) { /*  send colors */
            if (i!=gRank) MPI_Send (foremanOf, gSize, MPI_INT, i, CLR_TAG, world);
        }


    }

    else { /*  send proc name to coordinator, wait for group color */
        idx = strlen(procName[gRank]);
        MPI_Sendrecv (procName[gRank], idx, MPI_CHAR, coordinator, HNM_TAG,
            foremanOf, gSize, MPI_INT, coordinator, CLR_TAG, world, &tStatus);
    }

    /*  map global ranks to local ranks */
    for (i=0;i<gSize;i++) localRank[i] = -1;
    for (i=0;i<gSize;i++)  {
        idx = foremanOf[i]; /*  increment no. of workers of foremanOf idx */
        localRank[idx] = localRank[idx]+1; /*  temporary counter; reset in next loop  */
        localRank[i] = localRank[idx]; /*  this worker's rank */
    }
    /*  reset local rank of a foremanOf to 0  */
    for (i=0;i<gSize;i++)
        if (foremanOf[i]==i) localRank[i]=0;

    /*  split processors into work groups  */
    if (gRank==coordinator)
        idx = MPI_UNDEFINED;
    else
        idx = foremanOf[gRank];
    MPI_Comm_split(world, idx, gRank, icomm);

    /*  make control group for the foremen and coordinator */
    if ((gRank==coordinator) || (foremanOf[gRank]==gRank))
        idx = 0;
    else
        idx = MPI_UNDEFINED;
    MPI_Comm_split(world, idx, gRank, jcomm);
}

void GetChunkSize ( infoDLS *info, int rank, int *chunkSize )
{
    int i, tChunk, rem;
    double bigD, bigT, awap, trw, weight, K;

    rem = info->N-info->itersScheduled;
    switch ( info->method ) {

    case STATIC:
         tChunk = ceil((double) info->N/ (double) info->commSize);
         info->batchSize = tChunk;
         info->batchRem = min( info->batchSize, rem);
         break;
    case SS:
         tChunk = 1;
         info->batchSize = tChunk;
         info->batchRem = min( info->batchSize, rem);
         break;
    case FSC:
        tChunk = min( info->chunkFSC, rem);
        info->batchSize = tChunk;
        info->batchRem = min( info->batchSize, rem);
        break;

    case mFSC:
        tChunk = min( info->chunkMFSC, rem);
        info->batchSize = tChunk;
        info->batchRem = min( info->batchSize, rem);
        break;

    case GSS:
        tChunk = max( (rem+info->commSize-1)/info->commSize, info->minChunkSize );
        tChunk = min ( rem, tChunk );
        info->batchSize = tChunk;
        info->batchRem = min( info->batchSize, rem);
        break;

    case TSS:
	tChunk = info->TSSchunk;
        tChunk = min(rem, tChunk);
        tChunk = max( info->minChunkSize, tChunk);
        info->TSSchunk = tChunk - info->TSSdelta;
        info->batchSize = tChunk;
        info->batchRem = min( info->batchSize, rem);
        break;

    case FAC:
        if (info->batchRem == 0) {
            tChunk = max ( info->minChunkSize, (rem+2*info->commSize-1)/(2*info->commSize) );
            info->batchSize = info->commSize*tChunk;
            info->batchRem = min (info->batchSize, rem);
        }
        /* else use current batchSize */
        tChunk = max( info->minChunkSize, info->batchSize/info->commSize );
        tChunk = min ( rem, tChunk );
        break;

    case WF:
    case AWF:
         if (info->batchRem == 0) {
            tChunk = max ( info->minChunkSize, (rem+2*info->commSize-1)/(2*info->commSize) );
            info->batchSize = info->commSize*tChunk;
            info->batchRem = min (info->batchSize, rem);
         }
         /* else use current batchSize */
         tChunk = max( info->minChunkSize, info->batchSize/info->commSize*info->weights[rank]);
         tChunk = min ( rem, tChunk );
         break;

    case AWF_B: 
    case AWF_D:
        if (info->stats[3*rank] < 0.0) {
            tChunk = info->minChunkSize;
            info->batchSize = min(rem, tChunk);
            info->batchRem = info->batchSize;
        }
        else { /* all ranks have wap */
            awap = 0.0;  /* average weighted performance */
            for (i=info->firstRank;i<=info->lastRank;i++)
                awap = awap + info->stats[3*i];
            awap = awap/info->commSize;

            trw = 0.0;  /* total ref weight (refwt(i) = awap/info->stats[3*i] */
            for (i=info->firstRank;i<=info->lastRank;i++)
                trw = trw + awap/info->stats[3*i];

            /* normalized weight for rank */
            weight = ((awap/info->stats[3*rank])*info->commSize)/trw;

            if (info->batchRem == 0) {
                tChunk = max( info->minChunkSize, (rem+2*info->commSize-1)/(2*info->commSize) );
                info->batchSize = info->commSize*tChunk;
                info->batchRem = min (info->batchSize, rem);
            }
            /* else use current batchSize */
            tChunk = weight*(info->batchSize/info->commSize) + 0.55;
            tChunk = max( info->minChunkSize, tChunk);
            tChunk = min ( rem, tChunk );
        }
        break;

    case AWF_C: 
    case AWF_E:
        if (info->stats[3*rank] < 0.0)
            tChunk = info->minChunkSize;
        else { /* all ranks have wap */
            awap = 0.0;  /* average weighted performance */
            for (i=info->firstRank;i<=info->lastRank;i++)
                awap = awap + info->stats[3*i];
            awap = awap/info->commSize;

            trw = 0.0;  /* total ref weight (refwt(i) = awap/info->stats[3*i) */
            for (i=info->firstRank;i<=info->lastRank;i++)
                trw = trw + awap/info->stats[3*i];

            /* normalized weight for rank */
            weight = ((awap/info->stats[3*rank])*info->commSize)/trw;
            tChunk = weight*((rem+2*info->commSize-1)/(2*info->commSize)) + 0.55;
        }
        tChunk = max( info->minChunkSize, tChunk);
        info->batchSize = tChunk;
        info->batchRem = min(rem, tChunk);
        break;

    case AF:
        if (info->stats[3*rank] < 0.0)
            tChunk = info->minChunkSize;
        else {
            bigD = 0.0;
            bigT = 0.0;
            for (i=info->firstRank;i<=info->lastRank;i++) {
                bigD = bigD + info->stats[3*i+1]/info->stats[3*i];
                bigT = bigT + 1.0/info->stats[3*i];
            }
            bigT = 1.0/bigT;
            /* compute chunk size for rank */
            tChunk = 0.55 + (0.5*(bigD + 2.0*bigT*rem -
                sqrt(bigD*(bigD + 4.0*bigT*rem)))/info->stats[3*rank]);
            tChunk = min( info->maxChunkSize, tChunk);
        }
        tChunk = max( info->minChunkSize, tChunk);
        info->batchSize = tChunk;
        info->batchRem = min( info->batchSize, rem);
        break;

    default:
        printf("Unsupported DLS technique, fall back to STATIC\n");
        tChunk = (info->N+info->commSize-1)/info->commSize;
        i = info->N % info->commSize;
        if ((i>0) && (rank>=i) ) tChunk=tChunk-1;
        tChunk = min( tChunk, rem);
        info->batchSize = tChunk;
        info->batchRem = min( info->batchSize, rem);

    }

    //printf("Rank: %d, Method: %d, chunk size is %d \n", rank, info->method, tChunk);
    *chunkSize = min(info->batchRem, tChunk);

    /* adjust remaining in batch */
    info->batchRem = info->batchRem - *chunkSize;
    if ( (info->batchRem > 0) && (info->batchRem <= info->minChunkSize) ) {
        *chunkSize = *chunkSize + info->batchRem;
        info->batchRem = 0;
    }
}


void SendChunk ( infoDLS *info, int worker )
{
    int chunkSize, chunkInfo[2];

    GetChunkSize (info, worker, &chunkSize);

    chunkInfo[0] = info->chunkStart;
    chunkInfo[1] = chunkSize;
    //printf("send to worker %d \n", worker);
    //
    if(worker == info->foreman)
    {

       if (info->wSize == 0)  // no pending chunk
       {
            info->t0 = MPI_Wtime(); // elapsed time for chunk starts here
            info->wStart = chunkInfo[0];
            info->wSize = chunkInfo[1];
            info->rStart = info->wStart;
            info->rSize = info->wSize;
            info->req4WRKsent = 0;  // cancel request for work

            SetBreaks(info);

           //printf("WRK_TAG recv by %d %d %d", info->myRank, info->wStart, info->wSize);

            info->sumt1 = 0.0;  //for mu/wap
            info->sumt2 = 0.0;  // for sigma
         }
         else  //current chunk is not finished  save as next chunk
         {
            info->nextStart = chunkInfo[0];
            info->nextSize = chunkInfo[1];
            info->nextWRKrcvd = 1;

           //printf('WRK_TAG (adv) recv by %d, %d, %d \n", info->myRank,info->nextStart,info->nextSize);
         }
      }
      else
      {
          MPI_Send (chunkInfo, 2, MPI_INT, worker, WRK_TAG, info->comm);
      }
    //printf("message sent \n");
    info->chunkStart     = info->chunkStart     + chunkSize;
    info->itersScheduled = info->itersScheduled + chunkSize;

    info->numChunks = info->numChunks + 1;
//    info->chunkMap[2] = info->numChunks;
//    info->chunkMap[3*info->numChunks  ] = chunkInfo[0];
//    info->chunkMap[3*info->numChunks+1] = chunkInfo[1];
//    info->chunkMap[3*info->numChunks+2] = worker;
    // printf("... sent start=%d, size=%d to %d\n", chunkInfo[0], chunkInfo[1], worker); 
}


void  SetBreaks ( infoDLS *info )
{
    if (info->myRank == info->foreman) {
        /* when to check for messages */
        if (info->breakAfter<0)
            info->probeFreq = max( 1, (info->wSize+info->commSize-1)/info->commSize/4 );
        else
            info->probeFreq = max( 1, info->breakAfter);

        /* how many iterates left before requesting next chunk */
        if (info->requestWhen<0)
            info->sendRequest = info->probeFreq;
        else
            info->sendRequest = info->requestWhen;
    }
    else { /* not the foreman */

        /* how many iterates left before requesting next chunk */
        if (info->requestWhen<0)
            info->sendRequest = max( 1, (15*info->wSize)/100 );
        else
            info->sendRequest = info->requestWhen;

        /* when to check for messages */
        info->probeFreq = max( 1, info->wSize-info->sendRequest);

    }
}


void DLS_Setup ( MPI_Comm icomm, infoDLS *info )
{
    int tP;

    MPI_Comm_size(icomm, &tP);
    MPI_Comm_rank(icomm, &(info->myRank));
    info->comm = icomm;
    info->crew = MPI_COMM_NULL;
    info->commSize = tP;
    info->firstRank = 0;
    info->lastRank = tP-1;
    info->foreman = 0;
    info->breakAfter = -1;
    info->probeFreq -1;
    info->minChunk = tP;
    info->stats = malloc(3*tP*sizeof(double));   
    info->weights = malloc(tP*sizeof(double));
    info->timeStep = 0;
}
void DLS_Parameters_Setup( MPI_Comm icomm, infoDLS *info, int numProcs, int requestWhen, int breakAfter, int minChunk, double h_overhead, double sigma, int nKNL, double Xeon_speed, double KNL_speed )
{
    int tP;
    double total_sum = 0.0;
    double core_speed = 0.0;
    int i;

    MPI_Comm_size(icomm, &tP);
    MPI_Comm_rank(icomm, &(info->myRank));
    info->comm = icomm;
    info->crew = MPI_COMM_NULL;
    info->commSize = tP;
    info->firstRank = 0;
    info->lastRank = tP-1;
    info->foreman = 0;
    info->breakAfter = breakAfter;
    info->requestWhen= requestWhen;
    info->probeFreq -1;
    info->minChunk = minChunk;
    info->stats = malloc(3*numProcs*sizeof(double));
    // h and sigma for FSC orginial equation
    info->h_overhead = h_overhead;
    info->sigma = sigma;
    // calculate weights for WF .. assume two types of processors, Xeon and KNL
    info->weights = malloc(numProcs*sizeof(double));
    info->timeStep = 0;
    total_sum = nKNL * KNL_speed + (numProcs - nKNL) * Xeon_speed ;
    for (i = 0; i < numProcs; i++)
    {
       //initialize mu, sigma, and performance data count
       info->stats[3*i]   = -1; //mu
       info->stats[3*i+1] = -1; //sigma
       info->stats[3*i+2] = 0; //performance data count
 	
       // initialize weights
       if(i<(numProcs - nKNL))
       {core_speed = Xeon_speed;}
       else
       {core_speed = KNL_speed; }
       info->weights[i] =  core_speed/total_sum * numProcs;
    }
  
 
}


void DLS_GroupSetup ( MPI_Comm world, int coordinator,
infoDLS *GS, infoDLS *LS )
{
    MPI_Comm jcomm, icomm;
    int tP;

    Partition(world, coordinator, &jcomm, &icomm);

    GS->comm = jcomm;
    GS->crew = icomm;
    LS->comm = icomm;
    LS->crew = icomm;
    if (icomm!=MPI_COMM_NULL) {
        DLS_Setup (icomm, LS);
        MPI_Comm_size(icomm, &(GS->crewSize));
    }
    else {
        LS->myRank = -1;
        LS->commSize = 0;
        LS->crewSize  = 0;
        GS->crewSize  = 0;
    }
    if (jcomm!=MPI_COMM_NULL) {
        MPI_Comm_size(jcomm, &tP);
        MPI_Comm_rank(jcomm, &(GS->myRank));
        GS->commSize = tP-1;
        GS->firstRank = 1;
        GS->lastRank = tP-1;
        GS->breakAfter = -1;
        GS->probeFreq -1;
        GS->minChunk = tP-1;
    }
    else {
        GS->myRank = -1;
        GS->commSize = 0;
    }
}



void DLS_StartLoop ( infoDLS *info, int firstIter, int lastIter, int imeth )
{
    int tSize, worker;
    int NULLloc;
    double K;
    info->wSize = 0;  /* remaining iterates in current chunk */
    info->gotWork = 1; /*.true.; */
    info->workTime = 0.0;
    info->myIters = 0;
    info->N = lastIter - firstIter + 1;
    info->timeStep = info->timeStep + 1;
    int i;
    double awap, trw;

    if ( (info->comm==MPI_COMM_NULL) || (info->N<=0) ) return;

    info->method = imeth;
    if ( (imeth>DLS_MethodCount-1) || (imeth<0) )
        info->method = 0;
    else
        info->method = imeth;

    info->firstIter = firstIter;
    info->lastIter = lastIter;
    info->N = lastIter - firstIter + 1;
    // TSS
    info->TSSchunk = ceil((double) info->N / ((double) 2*info->commSize)); 
    int n = ceil(2*info->N/(info->TSSchunk+1)); //n=2N/f+l
    info->TSSdelta = (double) (info->TSSchunk - 1)/(double) (n-1);

 //  printf("ID: %d, TSS chunk %d, delta: %d\n",info->myRank,info->TSSchunk,info->TSSdelta);
 //  info->chunkMap[0] = firstIter;    /* start of data */
 //  info->chunkMap[1] = info->N;    /* size of data */
 //  info->chunkMap[2] = 0; /* chunks in this rank */


//  calculate AWF weights
    if ( (info->method == AWF) && (info->myRank == info->foreman))
    {
        if (info->timeStep == 1) //first timeStep
        {
            for( i = info->firstRank; i<=info->lastRank; i++)
            {
                 info->weights[i] = 1.0;
             }
         }
         else // all ranks have wap
         {
          awap = 0.0;  // average weighted performance
          for(i=info->firstRank; i<= info->lastRank; i++)
          {
            //printf("rank %d: %lf", i, info->stat[3*i]);
            awap = awap + info->stats[3*i];
          }
          awap = awap/info->commSize;

          trw = 0.0;  // total ref weight (refwt(i) = awap/info%stats(3*i)
          for( i=info->firstRank; i<=info->lastRank; i++)
          {
            trw = trw + awap/info->stats[3*i];
          }

           for( i=info->firstRank; i<=info->lastRank; i++)
           {
               info->weights[i] = ((awap/info->stats[3*i])*info->commSize)/trw;
           }

         }
    }


    info->numChunks = 0;

    info->myExecs = 0;
    info->mySumTimes = 0.0;
    info->mySumSizes = 0.0;
    tSize = (info->N+info->commSize-1)/info->commSize;
    info->chunkMFSC = (0.55+tSize*log(2.0)/log( (1.0*tSize) ) );
    info->kopt0 = sqrt(2.0)*info->N/( info->commSize*sqrt(log(1.0*info->commSize)) );

   //calculate FSC chunk
    K=(sqrt(2)*info->N*info->h_overhead)/(info->sigma*info->commSize*sqrt(log(info->commSize)));
    K=pow(K, 2.0/3.0);
    info->chunkFSC = (int) ceil(K);
 

  
    info->nextWRKrcvd = 0;
    info->req4WRKsent = 0;
    info->finishedOne = 0;

    info->probeFreq = max(1, info->breakAfter);
    info->sendRequest = max(1, info->requestWhen);

    if (info->myRank == info->foreman) {
        info->chunkStart = firstIter;
        info->itersScheduled = 0;
        info->batchSize = 0;
        info->batchRem = 0;
        info->numENDed = 0;
        info->numChunks = 0;
   
   
        if (info->minChunk>0)
            info->minChunkSize = info->minChunk;
        else
            info->minChunkSize = max(1,info->chunkMFSC/2);       /* default min chunk size */
        info->maxChunkSize = (info->N+2*info->commSize-1)/(2*info->commSize);

        /* send initial work to each processor */
        for (worker=info->firstRank;worker<=info->lastRank;worker++)
        { 
            if (info->chunkStart < info->lastIter)
            { 
                 SendChunk (info, worker);
            }
            else
            {
                 MPI_Send (&NULLloc, 0, MPI_INT, worker, END_TAG, info->comm);  //end worker
                 info->numENDed++; // increment ended workers

            }
         }
    }
}

int DLS_Terminated ( infoDLS *info ) {
    int i, done;
    MPI_Status tStatus;

    if (info->N<=0) done=1;
    else {
      if ( (info->comm==MPI_COMM_NULL) &&  (info->crew!=MPI_COMM_NULL) )
        MPI_Recv (&done, 1, MPI_INT, 0, TRM_TAG, info->crew, &tStatus);
      else
        if ( (info->comm!=MPI_COMM_NULL) && (info->crew!=MPI_COMM_NULL) ) {
            done = (info->gotWork==0) && (info->wSize==0);
            for (i=1;i<info->crewSize;i++)
                MPI_Send (&done, 1, MPI_INT, i, TRM_TAG, info->crew);
        }
        else
            done = (info->gotWork==0) && (info->wSize==0);
    }
    return (done);
}



void DLS_StartChunk (infoDLS *info, int *chunkStart, int *chunkSize)
{
    int tSize, tStart, worker;
    int MsgInQueue;              /* message came in */
    int loc, maxRemaining;        /* source of chunk to be migrated */
    int i, j, NULLloc, chunkInfo[2];
    double perfInfo[4];
    MPI_Status mStatus, tStatus;

    if (info->comm==MPI_COMM_NULL) { /* I'm just a simple worker */
        MPI_Recv (chunkInfo, 2, MPI_INT, 0, WRK_TAG, info->crew, &tStatus);
        *chunkStart = chunkInfo[0];
        *chunkSize = chunkInfo[1];
    }
    else { /* I'm the coordinator, or a foreman */

        if (info->wSize == 0) {
            MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, info->comm, &mStatus);
            MsgInQueue = 1; /*.true. */
        }
        else
            MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, info->comm, &MsgInQueue, &mStatus);

        while (MsgInQueue) {

            switch ( mStatus.MPI_TAG ) {

            case (WRK_TAG):
                MPI_Recv (chunkInfo, 2, MPI_INT, mStatus.MPI_SOURCE, WRK_TAG, info->comm, &tStatus);

                if (info->wSize == 0) { /* no pending chunk */
                    info->t0 = MPI_Wtime(); /* elapsed time for chunk starts here */
                    info->wStart = chunkInfo[0];
                    info->wSize = chunkInfo[1];
                    info->rStart = info->wStart;
                    info->rSize = info->wSize;
                    info->req4WRKsent = 0; /* cancel request for work */
                    SetBreaks (info);
                    info->sumt1 = 0.0; /* for mu/wap */
                    info->sumt2 = 0.0; /* for sigma */
                }
                else { /* current chunk is not finished  save as next chunk */
                    info->nextStart = chunkInfo[0];
                    info->nextSize = chunkInfo[1];
                    info->nextWRKrcvd = 1; /*.true. */
                }
                break;

            case (REQ_TAG): /* received by foreman only */
                worker = mStatus.MPI_SOURCE;
                MPI_Recv (perfInfo, 4, MPI_DOUBLE, worker, REQ_TAG, info->comm, &tStatus);
                /* printf("... recv REQ_TAG from %d\n", worker); */

                if ( (info->method==AF) || (info->method==AWF_B) || (info->method==AWF_C) ||
                    (info->method==AWF_D) || (info->method==AWF_E) ) {
                    loc = perfInfo[2];
                    info->stats[3*loc+2] = info->stats[3*loc+2]+1.0;
                    /* adaptive methods */
                    info->stats[3*loc] = perfInfo[0];
                    info->stats[3*loc+1] = perfInfo[1];

                    if (info->finishedOne != info->commSize) {
                        /* workers that have not finished a first chunk */
                        /*  assume the lowest performance */
                        j = loc;
                        for (i=info->firstRank;i<=info->lastRank;i++)
                            if ( (info->stats[3*i+2] > 0.0) &&
                                (info->stats[3*i] < info->stats[3*j]) ) j = i;
                        info->finishedOne = 0;
                        for (i=info->firstRank;i<=info->lastRank;i++)
                            if (info->stats[3*i+2] == 0.0) {
                                info->stats[3*i] = info->stats[3*j];
                                info->stats[3*i+1] = info->stats[3*j+1];
                            }
                            else
                                info->finishedOne = info->finishedOne + 1;
                    }
                }

                /* any remaining unscheduled iterates ? */
                if (info->chunkStart <= info->lastIter)
                    SendChunk (info, worker);
                else { /* all iterates scheduled */
                    info->numENDed = info->numENDed + 1;
                    if (worker != info->myRank) {
                        MPI_Send (&NULLloc, 0, MPI_INT, worker, END_TAG, info->comm);
                        /* printf("... sent END_TAG to %d\n", worker); */
                    }
                    info->gotWork = (info->numENDed!=info->commSize); /* foreman exits? */
                }
                break;

            case (END_TAG): /* received by workers only */
                MPI_Recv (&NULLloc, 0, MPI_INT, mStatus.MPI_SOURCE,
                    mStatus.MPI_TAG, info->comm, &tStatus);
                info->gotWork = 0;
                break;

            } /* switch */
            MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, info->comm, &MsgInQueue, &mStatus);
        } /* while (MsgInQueue) */

        *chunkStart = info->wStart;
        *chunkSize = min (info->wSize, info->probeFreq);
        if (info->method == AF) *chunkSize = min(1, *chunkSize);
        info->subChunkSize = *chunkSize;
        if (info->subChunkSize!=0) info->t1 = MPI_Wtime();

        /* relay chunkStart, chunkSize to icomm */
        if (info->crew!=MPI_COMM_NULL) {
            chunkInfo[0] = *chunkStart;
            chunkInfo[1] = *chunkSize;
            for (i=1;i<info->crewSize;i++)
                MPI_Send (chunkInfo, 2, MPI_INT, i, WRK_TAG, info->crew);
        }
    } /* (info->comm!=MPI_COMM_NULL) { */
}


void  DLS_EndChunk (infoDLS *info)
{
    double tk, perfInfo[4];
    int loc;
    int i, j;

    if (info->comm==MPI_COMM_NULL) return;

    if (info->subChunkSize==0) return;

    tk = MPI_Wtime();
    info->t1 = tk - info->t1;
    info->wStart = info->wStart + info->subChunkSize;
    info->wSize = info->wSize - info->subChunkSize;
    info->sumt1 = info->sumt1 + info->t1;
    info->workTime = info->workTime + info->t1;
    if (info->method == AF) info->sumt2 = info->sumt2 + info->t1*info->t1;



 if (info->wSize == 0) { /* chunk finished */

        if ( (info->method==AWF_B) || (info->method==AWF_C) ) {
            /* adaptive weighted factoring, work time */
            info->mySumTimes = info->mySumTimes + (1+info->myExecs)*info->sumt1;
            info->mySumSizes = info->mySumSizes + (1.0+info->myExecs)*info->rSize;
        }
        else if ( (info->method==AWF_D) || (info->method==AWF_E) ) {
            /* adaptive weighted factoring, elapsed time */
            info->mySumTimes = info->mySumTimes + (1+info->myExecs)*(tk-info->t0);
            info->mySumSizes = info->mySumSizes + (1.0+info->myExecs)*info->rSize;
        }
      
        if(info->method !=AF)
        {
            /* reset accumulators */
            info->myIters = info->myIters + info->rSize;
            info->myExecs = info->myExecs + 1;
            info->sumt1 = 0.0; /* for mu */
            info->sumt2 = 0.0; /* for sigma */
            info->rSize = 0;
        }

        if (info->nextWRKrcvd) { /* foreman already responded to advance request */
            info->t0 = MPI_Wtime(); /* elapsed time for chunk starts here */
            info->wStart = info->nextStart;
            info->wSize = info->nextSize;
            info->rStart = info->wStart;
            info->rSize = info->wSize;

            SetBreaks (info);

            info->nextSize = 0;
            info->nextWRKrcvd = 0;
            info->req4WRKsent = 0;
        }
    } /* if (info->wSize == 0) */




    /* send request ? */
    if ( (info->wSize<=info->sendRequest) && (info->req4WRKsent==0) ) {
        switch (info->method) {
        case STATIC:
        case SS: 
        case FSC:
        case mFSC: 
        case GSS:
        case TSS: 
        case FAC:
        case WF:
        case AWF: 
            perfInfo[0] = 0.0;
            perfInfo[1] = 0.0;
            perfInfo[2] = 1.0*info->myRank;
            break;

        case AWF_B: 
        case AWF_C: /* mu = (chunk work time)/(chunk size) */
            perfInfo[0] = ( info->mySumTimes + (info->myExecs+1)*info->sumt1 ) /
                ( info->mySumSizes + 1.0*(info->myExecs+1)*(info->rSize-info->wSize) );
            perfInfo[1] = 0.0;
            perfInfo[2] = 1.0*info->myRank;
            break;

        case AF:
            perfInfo[0] = info->sumt1/(info->rSize-info->wSize);       /* mu */
            if ((info->rSize-info->wSize) > 1) {
                perfInfo[1] = (info->sumt2 - perfInfo[0]*perfInfo[0]*(info->rSize-info->wSize)) /
                    (info->rSize-info->wSize-1); /* sigma */
                if (perfInfo[1] < 0.0) perfInfo[1] = 0.0;
                perfInfo[1] = sqrt(perfInfo[1]);
            }
            else perfInfo[1] = 0.0;
            perfInfo[2] = 1.0*info->myRank;
           
            if (info->wSize == 0)
            {  // reset accumulators
             info->myIters = info->myIters + info->rSize;
             info->myExecs = info->myExecs + 1;
             info->sumt1 = 0.0;  // for mu 
             info->sumt2 = 0.0;  // for sigma 
             info->rSize = 0;
            }

            break;

        case AWF_D: 
        case AWF_E: /* mu = (chunk elapsed time)/(chunk size) */
            perfInfo[0] = ( info->mySumTimes + (info->myExecs+1)*(tk-info->t0)) /
                ( info->mySumSizes + 1.0*(info->myExecs+1)*(info->rSize-info->wSize) );
            perfInfo[1] = 0.0;
            perfInfo[2] = 1.0*info->myRank;
            break;
        }
        
       if(info->myRank == info->foreman)
       {   //update performance data
       if (info->method==AWF_B || info->method==AWF_C || info->method==AF || info->method==AWF_D || info->method==AWF_E)
        {
            loc = (int) perfInfo[2];
            info->stats[3*loc+2] = info->stats[3*loc+2]+1.0;
             //adaptive methods
             info->stats[3*loc] = perfInfo[0];
             info->stats[3*loc+1] = perfInfo[1];

             if (info->finishedOne != info->commSize)
             {
              // workers that have not finished a first chunk
              //  assume the lowest performance
              j = loc;
              for(i=info->firstRank; i<=info->lastRank; i++)
              {
                if ( (info->stats[3*i+2] > 0.0) && (info->stats[3*i] < info->stats[3*j]) )
			 j = i;
              }
              info->finishedOne = 0;
              for( i=info->firstRank; i<=info->lastRank;i++)
              {
                if (info->stats[3*i+2] == 0.0)
                {
                  info->stats[3*i] = info->stats[3*j];
                  info->stats[3*i+1] = info->stats[3*j+1];
                }
                else
                {
                   info->finishedOne = info->finishedOne + 1;
                 }
               }
              }
           }

	  // get more work to myself
          if (info->chunkStart < info->lastIter)
          {
//           printf("[end chunk] obtaining work ... current chunk start, %d \n", info->chunkStart);
             info->req4WRKsent = 1;
             info->nextWRKrcvd = 0;
             SendChunk (info, info->myRank);
//           printf("[end chunk] wsize: %d \n", info->wSize);
          }
          else if (info->wSize == 0 && info->chunkStart >= info->lastIter)
          {
	    // all iterates scheduled
            info->numENDed = info->numENDed + 1;
//          printf( "[end chunk] ended ranks %d chunk %d \n, info->numENDed,info->subChunkSize);
            info->gotWork = info->numENDed!=info->commSize; // foreman exits?
           }
       }
       else
       {
        MPI_Send (perfInfo, 4, MPI_DOUBLE, info->foreman, REQ_TAG, info->comm);
        info->req4WRKsent = 1; /*.true. */
        info->nextWRKrcvd = 0;
       }
    } /* if (...info->sendRequest...) */

}


void  DLS_EndLoop (infoDLS *info, int *niters, double *worktime)
{
    double perfInfo[4];

    if (info->comm==MPI_COMM_NULL) return;
    *niters = info->myIters;
    *worktime = info->workTime;
    //MPI_Barrier(info->comm);

//... Communicate time-step performance data for AWF
     if (info->method==AWF)
     {
        // timestepping adaptive weighted factoring 
        // mu = (chunk work time)/(chunk size) 
        perfInfo[0] = ( info->mySumTimes + (info->timeStep)*info->workTime ) / ( info->mySumSizes + 1.0*(info->timeStep)*(info->myIters));
        perfInfo[1] = 0.0;
        perfInfo[2] = 1.0*info->timeStep;
//      printf("time step %d\n", info->timeStep);
        MPI_Gather(perfInfo, 3, MPI_DOUBLE_PRECISION, info->stats, 3,MPI_DOUBLE_PRECISION, info->foreman, info->comm);

//      printf("rank %d send its performance %lf \n", info->myRank, PerfInfo[0]);
//          if (info->myRank == info->foreman)  // if foreman ...recieve all performance data
//          {
//            for(i=info->firstRank; i<=info->lastRank, i++)
//               printf( "performance data %d  %lf %lf %lf", i, info->stats[3*i], info->stats[3*i+1], info->stats[3*i+2]);
//          }
       } 

//    MPI_Barrier(info->comm) 

}


void  DLS_Finalize (infoDLS *info)
{   
    if (info->comm==MPI_COMM_NULL) 
        return;

    free(info->stats);
    free (info->weights);
    return;
}

