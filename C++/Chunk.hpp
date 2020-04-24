#pragma once

/********************************************************************************
* Ali Mohammed <ali.mohammed@unibas.ch>                                         *
* University of Basel, Switzerland                                              *
*                                                                               *
* This program is free software; you can redistribute it and/or modify it       *
* under the terms of the license (GNU LGPL) which comes with this package.      *
********************************************************************************/

class Chunk
{
int start; //start element of the chunk
int size; //chunk size
void ** data; // pointer to data related to this chunk

int * dataSources; // list of ranks who has each element of data

int * dataOffsets; // list of offsets where to start reading data in each data element

int * dataSizes; //list of how many elements to read from each data element for this chunk

};
