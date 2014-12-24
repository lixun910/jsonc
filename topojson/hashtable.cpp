//
//  hashtable.cpp
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//
#include <algorithm>

#include "hashtable.h"
#include "point.h"

using namespace std;

Arc::Arc()
{
    next = NULL;
}

Arc::Arc(int _first, int _second)
{
    first = _first;
    second = _second;
    next = NULL;
}

Arc::~Arc()
{
    if (next) {
        delete next;
        next = NULL;
    }
}

Arc& Arc::operator=(const Arc &arc)
{
    first = arc.first;
    second = arc.second;
    next = arc.next;
    
    return *this;
}

Hashmap::Hashmap(int _size)
{
    size = _size;
    keyEmpty = NULL;
    valEmpty = NULL;
    
    size = 1 << max(4, int(ceil(log(size)/log(2.0))));
    hm_mask= size -1;
    hm_free = size;
    
    for (int i=0; i < size; i++) {
        keystore.push_back(keyEmpty);
        valstore.push_back(valEmpty);
    }
}

void Hashmap::set(point* key, vector<Arc*>& value)
{
    int index = key->hash() & hm_mask;
    point* matchKey = keystore[index];
    int collisions = 0;
    
    while (matchKey != keyEmpty) {
        if ( *matchKey == *key ) {
            valstore[index] = value;
        }
        if ( ++collisions >= size) {
            throw new FullHashmapException();
        }
		index = (index+1) & hm_mask;
        matchKey = keystore[index];
    }
    keystore[index] = key;
    valstore[index] = value;
    --hm_free;
}

vector<Arc*> Hashmap::maybeSet(point* key, vector<Arc*>& value)
{
    int index = key->hash() & hm_mask;
    point* matchKey = keystore[index];
    int collisions = 0;
    
    while (matchKey != keyEmpty) {
        if (*matchKey == *key) {
            return valstore[index];
        }
        if (++collisions >= size) {
            throw new FullHashmapException();
        }
		index = (index+1) & hm_mask;
        matchKey = keystore[index];
    }
    keystore[index] = key;
    valstore[index] = value;
    --hm_free;
    return value;
}

vector<Arc*> Hashmap::get(point* key)
{
    int index = key->hash() & hm_mask;
    point* matchKey = keystore[index];
    int collisions = 0;
    
    while (matchKey != keyEmpty) {
        if (*matchKey == *key) {
            return valstore[index];
        }
        if (++collisions >= size) {
            break;
        }
		index = (index + 1) & hm_mask;
        matchKey = keystore[index];
    }
    return vector<Arc*>();
}

vector<point*> Hashmap::keys()
{
    vector<point*> keys;
    for (int i=0; i < keystore.size(); i++) {
        point* matchKey = keystore[i];
        if (matchKey != keyEmpty) {
            keys.push_back(matchKey);
        }
    }
    return keys;
}


Hashset::Hashset(int _size)
{
    size = _size;
    empty = NULL;
    
    size = 1 << max(4, int(ceil(log(size)/log(2.0))));
    hs_mask= size -1;
    hs_free = size;
    
    for (int i=0; i < size; i++) {
        store.push_back(empty);
    }
}

bool Hashset::add(point* value)
{
    int index = value->hash() & hs_mask;
    point* match = store[index];
    int collisions = 0;
   
    while (match != empty) {
        if (*match == *value) {
            return true;
        }
        if (++collisions >= size) {
            throw new FullHashsetException();
        }
		index = (index + 1) & hs_mask;
        match = store[index];
    }
    store[index] = value;
    --hs_free;
    return true;
}

bool Hashset::has(point* value)
{
    int index = value->hash() & hs_mask;
    point* match = store[index];
    int collisions = 0;

    while (match != empty) {
        if (*match == *value) {
            return true;
        }
        if (++collisions >= size) {
            break;
        }
		index = (index + 1) & hs_mask;
        match = store[index];
    }
    return false;
}

vector<point*> Hashset::values()
{
    vector<point*> values;
    for (int i=0; i < store.size(); i++) {
        point* match = store[i];
        if (match != empty) {
            values.push_back(match);
        }
    }
    return values;
}