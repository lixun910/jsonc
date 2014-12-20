//
//  hashtable.h
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#ifndef __topojson__hashtable__
#define __topojson__hashtable__

#include <iostream>

class Hashtable {
    int size;
    int mask;
    
private:
    int retFunc(float x, float y);
    
    bool equal(int* keyA, int* keyB);
    
public:
    Hashtable(int size);
    
    void* peak(int key);
    
    void* get(int key);
    
    
};
#endif /* defined(__topojson__hashtable__) */
