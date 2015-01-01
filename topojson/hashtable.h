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

class FullHashsetException: public exception
{
    virtual const char* what() const throw()
    {
        return "Full hashset";
    }
};

class NoMatchKeyFoundException: public exception
{
    virtual const char* what() const throw()
    {
        return "No match key found";
    }
};
class FullHashmapException: public exception
{
    virtual const char* what() const throw()
    {
        return "Full hashmap";
    }
};

class Arc {
    int hashCode;
public:
    Arc* next;
    int first; // index of coordinates
    int second;
    
public:
    Arc();
    Arc(int _first, int _second);
    ~Arc();
    
    Arc & operator=(const Arc& arc);
    
    int hash();
    bool equals(Arc* arcB);
};

class Hashmap {
	int size;
	int hm_mask;
	int hm_free;
    
	vector<int> keystore;
	vector<int> valstore;
    int keyEmpty;
    int valEmpty;

public:
    Hashmap(int _size);
    
    void set(int key, int value, vector<point*>& coords);
    
    int maybeSet(int key, int value, vector<point*>& coords);
   
    int get(int key, vector<point*>& coords);
    
    vector<int> keys();
};

class PointHashmap {
	int size;
	int hm_mask;
	int hm_free;
    
	vector<point*> keystore;
	vector<vector<Arc*>> valstore;
    point* keyEmpty;
    vector<Arc*> valEmpty;

public:
    PointHashmap(int _size);
    
    void set(point* key, vector<Arc*>& value);
    
    vector<Arc*> maybeSet(point* key, vector<Arc*>& value);
   
    vector<Arc*> get(point* key);
    
    vector<point*> keys();
};



class ArcHashmap {
    
	int size;
	int hm_mask;
	int hm_free;
    
	vector<Arc*> keystore;
	vector<int> valstore;
    Arc* keyEmpty;
    int valEmpty;

public:
    ArcHashmap(int _size);
    
    void set(Arc* key, int value);
   
    int get(Arc* key);
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
