#pragma once

/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

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

#include "Master.hpp"
#include "Worker.hpp"


int dlsMethodCount = 15;

char *dlsShortName [] = {
     "STATIC", "SS", "FSC", "mFSC", "GSS", "TSS","FAC", "WF", "AWF",
     "AWF-B", "AWF-C", "AWF-D", "AWF-E", "AF", "SimAS"
     };

char *dlsLongName [] = {
     "STATIC SCHEDULING",
     "SELF-SCHEDULING",
     "FIXED SIZE-SCHEDULING", "MODIFIED FSC", "GUIDED SELF-SCHEDULING",
     "Trapezoid Self-Scheduling","FACTORING",              "WEIGHTED FACTORING", "ADAPTIVE WEIGHTED FACTORING",
     "BATCH AWF",              "CHUNK AWF",
     "BATCH AWF (chunk times)" ,"CHUNK AWF (chunk times)", "ADAPTIVE FACTORING",  "Simulation-assisted"
     };


char * getDLSShortName(int DLS)
{
return dlsShortName[DLS];
}

char * getDLSLongName(int DLS)
{
return dlsLongName[DLS];
}
