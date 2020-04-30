#pragma once
/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

#include <iostream> 
#include <string>

#define dlsMethodCount 15

#define STATIC  0
#define SS      1
#define FSC     2
#define mFSC    3
#define GSS     4
#define TSS     5
#define FAC     6
#define WF      7
#define AWF     8
#define AWFB    9
#define AWFC   10
#define AWFD   11
#define AWFE   12
#define AF     13
#define SimAS  14

using namespace std;

class DLS
{
double *mu; //mean execution time
double *sigma; //standard deviation of execution time
double *weights; //relative weight
int method;
int minChunkSize;
double schOverhead; //h for FSC
int TSSDelta; // TSS delta 
int TSSFinalChunk;
int nWorkers;


static string dlsShortName[dlsMethodCount];

static string dlsLongName [dlsMethodCount];

public:
	DLS(const int meth=0,const int P=1, const int minCSize=1, double w[]=nullptr); //For all DLS, exept FSC and TSS, see below
	DLS(const int meth, const int P, const int minCSize, const double schH); //for FSC
        DLS(const int meth,const int P, const int minCSize, const int TSSFChunk); //for TSS
        int calculateChunk();
        static string getDLSShortName(int method);
        static string getDLSLongName(int method);
        ~DLS();
};
