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

class FullHashmapException: public exception
{
    virtual const char* what() const throw()
    {
        return "Full hashmap";
    }
};

class Hashmap {
	vector<point*> keystore;
	vector<vector<Arc*>> valstore;
	int size;
	int hm_mask;
	int hm_free;
    point* keyEmpty;
    Arc* valEmpty;

    
public:
    Hashmap(int _size);
    
    void set(point* key, vector<Arc*>& value);
    
    vector<Arc*> maybeSet(point* key, vector<Arc*>& value);
   
    vector<Arc*> get(point* key);
    
    vector<point*> keys();
};

class FullHashsetException: public exception
{
    virtual const char* what() const throw()
    {
        return "Full hashset";
    }
};

class Hashset {
    vector<point*> store;
    int size;
	int hs_mask;
	int hs_free;
    point* empty;
    
public:
    Hashset(int _size);
   
    bool add(point* value);
    
    bool has(point* value);
    
    vector<point*> values();
};
#endif /* defined(__topojson__hashtable__) */
