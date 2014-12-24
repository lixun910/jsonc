//
//  hashtable.h
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#ifndef __topojson__hashtable__
#define __topojson__hashtable__

#include <vector>
#include <exception>

#include "point.h"

using namespace std;

class FullHashmapException: public exception
{
    virtual const char* what() const throw()
    {
        return "Full hashmap";
    }
};

class Hashmap {
	vector<int> keystore;
	vector<int> valstore;
	int size;
	int hm_mask;
	int hm_free;
    int keyEmpty;
    int valEmpty;

    
public:
    Hashmap(int _size);
    
    int set(int key, int value, vector<point>& coords);
    
    int maybeSet(int key, int value, vector<point>& coords);
   
    int get(int key, int missingValue, vector<point>& coords);
    
    vector<int> keys();
};

class FullHashsetException: public exception
{
    virtual const char* what() const throw()
    {
        return "Full hashset";
    }
};

class Hashset {
    vector<int> store;
    int size;
	int hs_mask;
	int hs_free;
    int empty;
    
public:
    Hashset(int _size);
   
    bool add(int value, vector<point>& coords);
    
    bool has(int value, vector<point>& coords);
    
    vector<int> values();
};
#endif /* defined(__topojson__hashtable__) */
