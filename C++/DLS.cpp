
/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

#include "DLS.hpp"

string DLS::dlsShortName[dlsMethodCount] = {
     "STATIC", "SS", "FSC", "mFSC", "GSS", "TSS","FAC", "WF", "AWF",
     "AWF-B", "AWF-C", "AWF-D", "AWF-E", "AF", "SimAS"
     };

string DLS::dlsLongName [dlsMethodCount]= {
     "STATIC SCHEDULING",
     "SELF-SCHEDULING",
     "FIXED SIZE-SCHEDULING", "MODIFIED FSC", "GUIDED SELF-SCHEDULING",
     "Trapezoid Self-Scheduling","FACTORING",              "WEIGHTED FACTORING", "ADAPTIVE WEIGHTED FACTORING",
     "BATCH AWF",              "CHUNK AWF",
     "BATCH AWF (chunk times)" ,"CHUNK AWF (chunk times)", "ADAPTIVE FACTORING",  "Simulation-assisted"
     };

DLS::DLS(const int meth,const int P, const int minCSize, double w[])
{
    //Check is method is FSC or TSS, you need to specify some more info
    if(meth == FSC)
    {
        std::cout << "Scheduling overhead is required for FSC!!" << std::endl;
        return;
    }
    if(meth == TSS)
    {
        std::cout << "Final chunk size is required for TSS!!" << std::endl;
        return;
    }
    method = meth;
    minChunkSize = minCSize;
    nWorkers = P;
    //if method needs weights
    if((method == WF) || (method == AWF) || (method == AWFB) || (method == AWFC) || (method ==AWFD) ||(method == AWFE))
    {
        //if weights are NULL, create it and set it to one
	if(w == nullptr)
        {
            weights = new double[nWorkers];
        }
        //else make local weights equal to weights
        else
        {
            weights = w;
        }
    }
}

DLS::DLS(const int meth, const int P, const int minCSize, const double schH) // for FSC
{
        
    method = meth;
    minChunkSize = minCSize;
    nWorkers = P;

    schOverhead = schH;
}


DLS::DLS(const int meth,const int P, const int minCSize, const int TSSFChunk) //for TSS
{

    method = meth;
    minChunkSize = minCSize;
    nWorkers = P;

    TSSFinalChunk = TSSFChunk;
}
 

string DLS::getDLSShortName(int method)
{
    return dlsShortName[method];
}

string DLS::getDLSLongName(int method)
{
    return dlsLongName[method];
}

DLS::~DLS()
{
    if(weights != nullptr)
    {
        delete [] weights;
        weights = nullptr;
    }
    if(mu != nullptr)
    {
        delete [] mu;
        mu = nullptr;
    }
    if(sigma != nullptr)
    {
        delete [] sigma;
        sigma = nullptr;
    }
}

//calculate the chunk size according to the selected DLS technique
int DLS::calculateChunk()
{
}



