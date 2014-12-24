//
//  point.h
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#ifndef __topojson__point__
#define __topojson__point__

#include <iostream>

typedef unsigned char byte;

class point {
    double x;
    double y;
    int hashCode;
    
public:
	point(double _x, double _y);
	point(const point& pt);
   
    bool operator==(const point& pt);
    
    point & operator=(const point& pt);
    
    int hash();
};
#endif /* defined(__topojson__point__) */
