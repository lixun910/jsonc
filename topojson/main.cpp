//
//  main.cpp
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#include <iostream>
#include <vector>
#include <map>
#include <string.h>
#include <stdlib.h>

#include "shapefil.h"
#include "topojson.h"

using namespace std;

int main(int argc, const char * argv[])
{
    // insert code here...
    std::cout << "Hello, World!\n";
    
    Topojson topo( "C:\\Users\\Xun\\Documents\\GitHub\\PySAL-Viz\\test_data\\NAT.shp" );
	topo.Extract();

    return 0;
}

