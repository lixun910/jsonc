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
#include "hashtable.h"

using namespace std;

class Arc {
public:
    Arc* next;
    int first;
    int second;
public:
    Arc();
    Arc(int _first, int _second);
    ~Arc();
    
    Arc & operator=(const Arc& arc);
};

class Topojson {
    int index;
    vector<Arc> lines;
    vector<Arc> rings;
    vector<point> coordinates;

	SHPHandle hSHP;
	int nShapeType;
	int nEntities;
	double adfMinBound[4];
	double adfMaxBound[4];

    vector<int> indexes;
	vector<int> visitedByIndex;
	vector<int> leftByIndex;
	vector<int> rightByIndex;
	vector<int> junctionByIndex;
	int junctionCount; // upper bound on number of junctions
    
public:
	Topojson(const char* shp_path);
	~Topojson();

	void Extract();
	Hashset* Join();
	void Cut();
    void Dedup();
    
private:
    void sequence(int i, int previousIndex, int currentIndex, int nextIndex);
    
    void rotateArray(vector<point>& array, int start, int end, int offset);
    
    void reverse(vector<point>& array, int start, int end);
};

#endif /* defined(__topojson__) */
