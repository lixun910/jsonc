//
//  hashtable.h
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#include <vector>
#include <map>
#include <string.h>
#include <stdlib.h>

#include "shapefil.h"
#include "point.h"
#include "topojson.h"
#include "hashtable.h"

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

Arc & Arc::operator=(const Arc &arc)
{
    first = arc.first;
    second = arc.second;
    next = arc.next;
    
    return *this;
}

Topojson::Topojson(const char* shp_path)
{
	index = -1;
	junctionCount = 0;
    
    hSHP = SHPOpen( shp_path, "rb" );
    if(hSHP == NULL)
    {
        printf( "Unable to open:%s\n", shp_path );
		return;
    }
    SHPGetInfo(hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

	Extract();
    
    Hashmap indexByPoint(coordinates.size() * 1.4);
    for (int i=0; i < coordinates.size(); i++) {
        indexes[i] = indexByPoint.maybeSet(i, i, coordinates);
    }
}

Topojson::~Topojson()
{
	SHPClose(this->hSHP);
}

void Topojson::sequence(int i, int previousIndex, int currentIndex, int nextIndex)
{
    if (visitedByIndex[currentIndex] == i) {
        // ignore self-intersection
        return;
    }
    visitedByIndex[currentIndex] = i;
    int leftIndex = leftByIndex[currentIndex];
    
    if (leftIndex >= 0) {
        int rightIndex = rightByIndex[currentIndex];
        if ((leftIndex != previousIndex || rightIndex != nextIndex)
            && (leftIndex != nextIndex || rightIndex != previousIndex)) {
            ++junctionCount;
            junctionByIndex[currentIndex] = 1;
        }
    } else {
        leftByIndex[currentIndex] = previousIndex;
        rightByIndex[currentIndex] = nextIndex;
    }
}

void Topojson::rotateArray(vector<point>& array, int start, int end, int offset)
{
    reverse(array, start, end);
    reverse(array, start, start + offset);
    reverse(array, start + offset, end);
}

void Topojson::reverse(vector<point> &array, int start, int end)
{
    for (int mid = start + ((end-- - start) >> 1); start < mid; ++start, --end) {
        point t(array[start]);
        array[start] = array[end];
        array[end] = t;
    }
}

void Topojson::Extract()
{
    if (nShapeType == SHPT_NULL || nShapeType == SHPT_POINT ||
        nShapeType == SHPT_POINTM || nShapeType == SHPT_POINTZ ) {
        return;
    }
    for (int i=0; i < hSHP->nRecords; i++) {
        SHPObject *psShape;
        psShape = SHPReadObject(hSHP, i);
		int nParts = psShape->nParts > 0 ? psShape->nParts : 1;
		for (int j=0; j < nParts; j++ ) {
			int start = psShape->panPartStart[j];
			int end = j < nParts-1 ? psShape->panPartStart[j+1] : psShape->nVertices;

            Arc arc;
			arc.first = index + 1;
			index += end - start;
			arc.second = index;

			if (nShapeType == SHPT_ARC || nShapeType == SHPT_ARCM ||
				nShapeType == SHPT_ARCZ) {
				lines.push_back(arc);
			} else if (nShapeType == SHPT_POLYGON || nShapeType == SHPT_POLYGONZ ||
				nShapeType == SHPT_POLYGONM) {
				rings.push_back(arc);
			}
			for (int k=start; k < end; k++ ) {
				coordinates.push_back(point(psShape->padfX[k], psShape->padfY[k]));
			}
		}
    }
}

Hashset* Topojson::Join()
{

	for (int i=0; i < coordinates.size(); i++) {
		visitedByIndex.push_back(-1);
		leftByIndex.push_back(-1);
		rightByIndex.push_back(-1);
	}

	for (int i=0; i < lines.size(); i++) {
		int lineStart = lines[i].first;
		int lineEnd = lines[i].second;
		int currentIndex = indexes[lineStart];
        int nextIndex = indexes[++lineStart];
        
        ++junctionCount;
        junctionByIndex[currentIndex] = 1; // start
        
        while (++lineStart <= lineEnd) {
            int previousIndex = currentIndex;
            currentIndex = nextIndex;
            nextIndex = indexes[lineStart];
            sequence(i, previousIndex, currentIndex, nextIndex);
        }
        ++junctionCount;
        junctionByIndex[nextIndex] = 1; // end
	}
    
    for (int i=0; i < coordinates.size(); i++) {
        visitedByIndex[i] = -1;
    }
    
    for (int i=0; i < rings.size(); i++) {
        Arc ring = rings[i];
        int ringStart = ring.first + 1;
        int ringEnd = ring.second;
        int previousIndex = indexes[ringEnd - 1];
        int currentIndex = indexes[ringStart - 1];
        int nextIndex = indexes[ringStart];
        sequence(i, previousIndex, currentIndex, nextIndex);
        while (++ringStart <= ringEnd) {
            previousIndex = currentIndex;
            currentIndex = nextIndex;
            nextIndex = indexes[ringStart];
            sequence(i, previousIndex, currentIndex, nextIndex);
        }
    }

    Hashset* junctionByPoint = new Hashset(junctionCount * 1.4);
    
    // Convert back to a standard hashset by point for caller convenience.
    for (int i=0; i < coordinates.size(); i++) {
        int j = indexes[i];
        if (junctionByIndex[j]) {
            junctionByPoint->add(j, coordinates);
        }
    }
    
    return junctionByPoint;
}

void Topojson::Cut()
{
    Hashset* junctions = Join();
   
    for (int i=0; i < lines.size(); i++) {
        Arc line = lines[i];
        int lineMid = line.first;
        int lineEnd = line.second;
        while (++lineMid < lineEnd) {
            if (junctions->has(lineMid, coordinates)) {
                Arc* next = new Arc(lineMid, line.second);
                line.second = lineMid;
                line.next = next;
                //line = *next;
            }
        }
    }
    
    for (int i=0; i < rings.size(); i++) {
        Arc ring = rings[i];
        int ringStart = ring.first;
        int ringMid = ringStart;
        int ringEnd = ring.second;
        bool ringFixed = junctions->has(ringStart, coordinates);
        while (++ringMid < ringEnd) {
            if (junctions->has(ringMid, coordinates)) {
                if (ringFixed) {
                    Arc* next = new Arc(ringMid, ring.second);
                    ring.second = ringMid;
                    ring.next = next;
                    //ring = *next;
                } else { // For the first junction, we can rotate rather than cut.
                    rotateArray(coordinates, ringStart, ringEnd, ringEnd - ringMid);
                    coordinates[ringEnd] = coordinates[ringStart];
                    ringFixed = true;
                    ringMid = ringStart; // restart; we may have skipped junctions
                }
            }
        }
    }
}

void Topojson::Dedup()
{
    int arcCount = lines.size() + rings.size();
    
    // Count the number of (non-unique) arcs to initialize the hashmap safely.
    for (int i=0; i < lines.size(); i++) {
        Arc line = lines[i];
        Arc* nextline = line.next;
        while (nextline) {
            ++arcCount;
            nextline = nextline->next;
        }
    }
    for (int i=0; i < rings.size(); i++) {
        Arc ring = rings[i];
        Arc* nextring = ring.next;
        while (nextring) {
            ++arcCount;
            nextring = nextring->next;
        }
    }
    
    Hashmap arcsByEnd(arcCount * 2 * 1.4);
    
    for (int i=0; i < lines.size(); i++) {
        Arc line = lines[i];
        Arc* nextline = line.next;
        
        do {
            dedupLine(nextline);
            nextline = nextline->next;
        } while (nextline);
    }
    
    for (int i=0; i < rings.size(); i++) {
        Arc ring = rings[i];
        if (ring.next) {
            Arc* nextring = &ring;
            do {
                dedupLine(nextring);
                nextring = nextring->next;
            } while (nextring);
        } else {
            dedupRing(&ring);
        }
    }
}
