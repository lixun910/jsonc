//
//  hashtable.cpp
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//
#include <math.h>

#include "hashtable.h"
#include "point.h"

Hashmap::Hashmap(int _size)
{
    size = _size;
    keyEmpty = -1;
    valEmpty = -1;
    
    size = 1 << max(4, int(ceil(log(size)/log(2.0))));
    hm_mask= size -1;
    hm_free = size;
    
    for (int i=0; i < size; i++) {
        keystore.push_back(keyEmpty);
        valstore.push_back(valEmpty);
    }
}

int Hashmap::set(int key, int value, vector<point>& coords)
{
    int index = coords[key].hash() & hm_mask;
    int matchKey = keystore[index];
    int collisions = 0;
    
    while (matchKey != keyEmpty) {
        if ( coords[matchKey] == coords[key] ) {
            return valstore[index] = value;
        }
        if ( ++collisions >= size) {
            throw new FullHashmapException();
        }
        matchKey = keystore[index = (index+1) & hm_mask];
    }
    keystore[index] = key;
    valstore[index] = value;
    --hm_free;
    return value;
}

int Hashmap::maybeSet(int key, int value, vector<point>& coords)
{
    int index = coords[key].hash() & hm_mask;
    int matchKey = keystore[index];
    int collisions = 0;
    
    while (matchKey != keyEmpty) {
        if (coords[matchKey] == coords[key]) {
            return valstore[index];
        }
        if (++collisions >= size) {
            throw new FullHashmapException();
        }
        matchKey = keystore[index = (index+1) & hm_mask];
    }
    keystore[index] = key;
    valstore[index] = value;
    --hm_free;
    return value;
}

int Hashmap::get(int key, int missingValue, vector<point> &coords)
{
    int index = coords[key].hash() & hm_mask;
    int matchKey = keystore[index];
    int collisions = 0;
    
    while (matchKey != keyEmpty) {
        if (coords[matchKey] == coords[key]) {
            return valstore[index];
        }
        if (++collisions >= size) {
            break;
        }
        matchKey = keystore[index = (index + 1) & hm_mask];
    }
    return missingValue;
}

vector<int> Hashmap::keys()
{
    vector<int> keys;
    for (int i=0; i < keystore.size(); i++) {
        int matchKey = keystore[i];
        if (matchKey != keyEmpty) {
            keys.push_back(matchKey);
        }
    }
    return keys;
}


Hashset::Hashset(int _size)
{
    size = _size;
    empty = -1;
    
    size = 1 << max(4, int(ceil(log(size)/log(2.0))));
    hs_mask= size -1;
    hs_free = size;
    
    for (int i=0; i < size; i++) {
        store.push_back(empty);
    }
}

bool Hashset::add(int value, vector<point>& coords)
{
    int index = coords[value].hash() & hs_mask;
    int match = store[index];
    int collisions = 0;
   
    while (match != empty) {
        if (coords[match] == coords[value]) {
            return true;
        }
        if (++collisions >= size) {
            throw new FullHashsetException();
        }
        match = store[index = (index + 1) & hs_mask];
    }
    store[index] = value;
    --hs_free;
    return true;
}

bool Hashset::has(int value, vector<point>& coords)
{
    int index = coords[value].hash() & hs_mask;
    int match = store[index];
    int collisions = 0;

    while (match != empty) {
        if (coords[match] == coords[value]) {
            return true;
        }
        if (++collisions >= size) {
            break;
        }
        match = store[index = (index + 1) & hs_mask];
    }
    return false;
}

vector<int> Hashset::values()
{
    vector<int> values;
    for (int i=0; i < store.size(); i++) {
        int match = store[i];
        if (match != empty) {
            values.push_back(match);
        }
    }
    return values;
}