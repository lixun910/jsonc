//
//  main.cpp
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#include <iostream>
#include <array>
#include <string.h>
#include <stdlib.h>
#include "shapefil.h"


int main(int argc, const char * argv[])
{
    // insert code here...
    std::cout << "Hello, World!\n";
    
    SHPHandle	hSHP;
    int		nShapeType, nEntities, i, nInvalidCount=0;
    double 	adfMinBound[4], adfMaxBound[4];
    
    hSHP = SHPOpen( argv[1], "rb" );

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

