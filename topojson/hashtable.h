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

using namespace std;

class Hashmap {
	vector<int> keystore;
	vector<int> valstore;
	int mask;
	int size;
	int free;
    
public:
    int set(int key, int value);
    
    int maybeSet(int key, int value);
    
    Hashmap(int size);
    
    void* peak(int key);
    
    void* get(int key);
    
    
};
#endif /* defined(__topojson__hashtable__) */
