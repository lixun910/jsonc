//
//  hashtable.h
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#ifndef __topojson__topojson__
#define __topojson__topojson__

#include <vector>
#include <map>
#include "point.h"
#include "shapefil.h"

using namespace std;

class Topojson {
    int index;
    vector<map<int,int> > lines;
    vector<map<int,int> > rings;
    vector<point> coordinates;

	SHPHandle hSHP;
	int nShapeType;
	int nEntities;
	double adfMinBound[4];
	double adfMaxBound[4];

public:
	Topojson(const char* shp_path);
	~Topojson();

	void Extract();
	void Join();    
};

#endif /* defined(__topojson__) */
