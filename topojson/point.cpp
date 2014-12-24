//
//  point.cpp
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#include "point.h"

point::point(const point& pt)
{
    x = pt.x;
    y = pt.y;
    hashCode = -1;
}

point::point(double _x, double _y)
{
	x = _x;
	y = _y;
    hashCode = -1;
}

int point::hash()
{
    if (hashCode != -1) {
        byte buffer[16] = {0};
        memcpy(buffer, &x, sizeof(double));
        memcpy(buffer + sizeof(double), &y, sizeof(double));
        
        unsigned int uint0 = *((int*)buffer);
        unsigned int uint1 = *((int*)(buffer + 4));
        unsigned int uint2 = *((int*)(buffer + 8));
        unsigned int uint3 = *((int*)(buffer + 12));
        
        int hash = uint0 ^ uint1;
        hash = hash << 5 ^ hash >> 7 ^ uint2 ^ uint3;
       
        hashCode = hash & 0x7fffffff;
    }
    return hashCode;
}

bool point::operator==(const point & pt)
{
    return x == pt.x && y == pt.y;
}

bool point::equal(const point & pt)
{
    return x == pt.x && y == pt.y;
}

point& point::operator=(const point & pt)
{
    if (this == &pt) return *this;
    
    if (x != pt.x || y != pt.y) {
        hashCode = -1;
    }
    x = pt.x;
    y = pt.y;
    
    return *this;
}