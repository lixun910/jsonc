//
//  main.cpp
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#include <iostream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include "shapefil.h"

using namespace std;

void Extract(SHPHandle hSHP, int nShapeType)
{
    if (nShapeType == SHPT_NULL || nShapeType == SHPT_POINT ||
        nShapeType == SHPT_POINTM || nShapeType == SHPT_POINTZ ) {
        return;
    }
    int index = -1;
    vector<int> lines;
    vector<int> rings;
    vector<int> coordinates;
    

    for (int i=0; i < hSHP->nRecords; i++) {
        SHPObject *psShape;
        psShape = SHPReadObject(hSHP, i);
        
        if (nShapeType == SHPT_ARC || nShapeType == SHPT_ARCM ||
            nShapeType == SHPT_ARCZ) {
            //extract line
            for(int j = 0, iPart = 1; j < psShape->nVertices; j++ ) {
                
            }
        }
    }
    
}
void Join(SHPHandle hSHP)
{
    int n = hSHP->nRecords;
   
    int* visitedByIndex = new int[n];
    int* leftByIndex    = new int[n];
    int* rightByIndex   = new int[n];

    for (int i=0; i < n; ++i) {
        visitedByIndex[i] = leftByIndex[i] = rightByIndex[i] = -1;
    }
    
    for (int i=0; i < n; ++i) {
        
    }
}




int main(int argc, const char * argv[])
{
    // insert code here...
    std::cout << "Hello, World!\n";
    
    SHPHandle	hSHP;
    int		nShapeType, nEntities, i, nInvalidCount=0;
    double 	adfMinBound[4], adfMaxBound[4];
    
    hSHP = SHPOpen( "", "rb" );
	
    if( hSHP == NULL )
    {
        printf( "Unable to open:%s\n", argv[1] );
        exit( 1 );
    }
    
    SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

    printf("ShapeType %d", nShapeType);
    
    SHPClose( hSHP );

    return 0;
}

