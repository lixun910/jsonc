//
//  hashtable.h
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#include <cmath>
#include <vector>
#include <map>
#include <string.h>
#include <stdlib.h>

#include "shapefil.h"
#include "point.h"
#include "topojson.h"
#include "hashtable.h"

using namespace std;

GeoObject::GeoObject(int nShapeType, SHPObject* psShape)
{
    int nParts = psShape->nParts > 0 ? psShape->nParts : 1;
    for (int j=0; j < nParts; j++ ) {
        int start = psShape->panPartStart[j];
        int end = j < nParts-1 ? psShape->panPartStart[j+1] : psShape->nVertices;
        int n = end - start;
        point* pt;
        for (int k=start; k < end; k++ ) {
            pt = new point(psShape->padfX[k], psShape->padfY[k]);
            coords.push_back(pt);
        }
        if (nShapeType == SHPT_ARC || nShapeType == SHPT_ARCM || nShapeType == SHPT_ARCZ) {
            if (n  < 2) { // must have 2+
                point* new_pt = new point(pt->x, pt->y);
                coords.push_back(new_pt);
            }
        } else if (nShapeType == SHPT_POLYGON || nShapeType == SHPT_POLYGONM
                   || nShapeType == SHPT_POLYGONZ) {
            if (n  < 4) { // must have 4+
                for (int k=0; k < 4 - n; k++) {
                    point* new_pt = new point(pt->x, pt->y);
                    coords.push_back(new_pt);
                }
            }
        }
    }
}

GeoObject::~GeoObject()
{
    for (int i=0; i < coords.size(); i++) {
        if (coords[i]) {
            delete coords[i];
            coords[i] = NULL;
        }
    }
    coords.empty();
}

void GeoObject::quantizePoint(point* pt)
{
    pt->x = round((pt->x + x0) * kx);
    pt->y = round((pt->y + y0) * ky);
}

void GeoObject::quantize(double _x0, double _y0, double _kx, double _ky)
{
    x0 = _x0;
    y0 = _y0;
    kx = _kx;
    ky = _ky;
    
    int i = 0;
    int j = 1;
    int n = coords.size();
    
    point* pi = coords[i];
    quantizePoint(pi);
    
    point* pj;
    double px = pi->x;
    double py = pi->y;
    double x;
    double y;
   
    while (++i < n) {
        pi = coords[i];
        quantizePoint(pi);
        x = pi->x;
        y = pi->y;
        if ( x != px || y != py ) {
            pj = coords[j++];
            pj->x = x;
            px = x;
            pj->y = y;
            py = y;
        }
    }
    
    length = j;
}

void GeoObject::normalizePoint(point *pt, double y0e, double y1e, double x0, double y0, double y1)
{
    if (pt->y <= y0e) {
        pt->x = 0;
        pt->y = y0;
    } else {
        if (pt->y >= y1e) {
            pt->x = 0;
            pt->y = y1;
        } else {
            pt->x = x0;
        }
    }
}

void GeoObject::stitch(int nShapeType, SHPObject *psShape, vector<vector<point*> >& fragments)
{
    
    if ( nShapeType != SHPT_POLYGON && nShapeType != SHPT_POLYGONM
        && nShapeType != SHPT_POLYGONZ ) {
        return;
    }

    int nParts = psShape->nParts > 0 ? psShape->nParts : 1;
    
    for (int p=0; p < nParts; p++ ) {
        int start = psShape->panPartStart[p];
        int end = p < nParts-1 ? psShape->panPartStart[p+1] : psShape->nVertices;
        int n = end - start;
        bool a = false;
        bool b = false;
        bool c = false;
        int i0 = -1;
        int i = 0;
        while (i < n ) {
            double x = psShape->padfX[i];
            double y = psShape->padfY[i];
            bool antimeridian = abs(abs(x) - 180.0) < 1e-2;
            bool polar = abs(abs(y) - 90.0) < 1e-2;
            if (antimeridian || polar) {
                if (! (a || b || c)) {
                    i0 = i;
                }
                if (antimeridian) {
                    if (a) {
                        c = true;
                    } else {
                        a = true;
                    }
                }
                if (polar) {
                    b = true;
                }
            }
            if ((!antimeridian && !polar) || i == n-1) {
                if ( a && b && c ) {
                    //del line[i0:i]
                    for (int j = i-1; j >= i0; j--) {
                        coords.erase(coords.begin() + start + j);
                    }
                    n -= i - i0;
                    i = i0;
                }
                a = false;
                b = false;
                c = false;
            }
            i += 1;
        }
    }
    /*
    double ε = 1e-2;
    double x0 = -180.0;
    double x0e = x0 + ε;
    double x1 = 180.0;
    double x1e = x1 - ε;
    double y0 = -90.0;
    double y0e = y0 + ε;
    double y1 = 90.0;
    double y1e = y1 - ε;
    
    double kx = this->kx;
    double ky = this->ky;
    double dx = 1 / this->x0;
    double dy = 1 / this->y0;
    
    x0 = round((x0 - dx) / kx);
    x1 = round((x1 - dx) / kx);
    y0 = round((y0 - dy) / ky);
    y1 = round((y1 - dy) / ky);
    x0e = round((x0e - dx) / kx);
    x1e = round((x1e - dx) / kx);
    y0e = round((y0e - dy) / ky);
    y1e = round((y1e - dy) / ky);
    
    int nParts = psShape->nParts > 0 ? psShape->nParts : 1;
    
    for (int p=0; p < nParts; p++ ) {
        int start = psShape->panPartStart[p];
        int end = p < nParts-1 ? psShape->panPartStart[p+1] : psShape->nVertices;
      
        vector<point*> ring;
        for (int i = start; i < end ; i++ ) ring.push_back(coords[i]);
        fragments.push_back(ring);
        
        int n = ring.size(); // n vertices in polygon
        
        for (int i = 0; i < n; i++ ) {
            point* pt = ring[i];
            double x = pt->x;
            double y = pt->y;
            
            // If this is an antimeridian or polar point…
            if (x <= x0e || x >= x1e || y <= y0e || y >= y1e) {
                
                // Advance through any antimeridian or polar points…
                int k = i + 1;
                for (; k < n; ++k) {
                    point* ptk = ring[k];
                    double xk = ptk->x;
                    double yk = ptk->y;
                    if (xk > x0e && xk < x1e && yk > y0e && yk < y1e)
                        break;
                }
                
                // If this was just a single antimeridian or polar point,
                // we don’t need to cut this ring into a fragment;
                // we can just leave it as-is.
                if (k == i + 1) continue;
                
                // Otherwise, if this is not the first point in the ring,
                // cut the current fragment so that it ends at the current point.
                // The current point is also normalized for later joining.
                if (i) {
                    vector<point*> fragmentBefore;
                    for (int m=0; m < i+1; m++) fragmentBefore.push_back(ring[m]);
                    normalizePoint(fragmentBefore[fragmentBefore.size() -1], y0e, y1e, x0, y0, y1);
                    fragments[fragments.size() -1] = fragmentBefore;
                    //fragments.push_back(
                }
                // If the ring started with an antimeridian fragment,
                // we can ignore that fragment entirely.
                else {
                    fragments.pop_back();
                }
                
                // If the remainder of the ring is an antimeridian fragment,
                // move on to the next ring.
                if (k >= n) break;
                
                // Otherwise, add the remaining ring fragment and continue.
                vector<point*> fragmentRemain;
                for (int m=k-1; m < n; m++) fragmentRemain.push_back(ring[m]);
                fragments.push_back(fragmentRemain);
                normalizePoint(ring[0], y0e, y1e, x0, y0, y1);
                //ring.polygon = polygon;
                i = -1;
                n = fragmentRemain.size();
            }
        }
    }
    */
}

Topojson::Topojson(const char* shp_path)
{
	junctionCount = 0;
    
    hSHP = SHPOpen( shp_path, "rb" );
    if(hSHP == NULL) {
        printf( "Unable to open:%s\n", shp_path );
		return;
    }
	double adfMinBound[4]; // The X, Y, Z and M minimum values
	double adfMaxBound[4];
    SHPGetInfo(hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

    bounds[0] = adfMinBound[0];
    bounds[1] = adfMinBound[1];
    bounds[2] = adfMaxBound[0];
    bounds[3] = adfMaxBound[1];
   
  
    if (Q0) {
        prequantize();
    }
    
    stich();

    computeTopology();
    
    if (Q1 && Q1 != Q0) {
        postquantize();
    }
    
    if (Q1) {
        Delta();
    }
}

Topojson::~Topojson()
{
	SHPClose(this->hSHP);
}

void Topojson::computeTopology()
{
    // extract
	Extract();
    
    // cut with join
    Cut();
    
    // dedup
    Dedup();
   
    // index: topology {
    //  coordinates: coordinates
    //  lines: lines
    //  rings: rings
    //  objects: geojson
    //  arcs: [] -->dedup
    indexByArc = new ArcHashmap(arcs.size() * 1.4);
    
    for (int i=0; i < arcs.size(); i++) {
        indexByArc->set(arcs[i], i);
        vector<point*> newArcs;
        for (int j=arcs[i]->first; j < arcs[i]->second + 1; j++) {
            newArcs.push_back(coordinates[j]);
        }
        topology_arcs.push_back(newArcs);
    }
   
    for (int i=0; i < objects.size(); i++) {
        for (int j=0; j < objects[i]->arcs.size(); j++) {
            vector<int> indexes = indexArcs(objects[i]->arcs[i]);
            objects[i]->index_arcs.push_back(indexes);
        }
    }
    
}

void Topojson::Save()
{
    // Convert to geometry objects.
    
}

void Topojson::Delta()
{
    int i = -1;
    int n = topology_arcs.size();
   
    while (++i < n ) {
        vector<point*> arc = topology_arcs[i];
        int j = 0;
        int m = arc.size();
        point* pt = arc[0];
        double x0 = pt->x;
        double y0 = pt->y;
        double x1;
        double y1;
        while (++j < m) {
            pt = arc[j];
            x1 = pt->x;
            y1 = pt->y;
            arc[j]->x = x1 -x0;
            arc[j]->y = y1 - y0;
            x0 = x1;
            y0 = y1;
        }
    }
}

vector<int> Topojson::indexArcs(Arc* arc)
{
    vector<int> indexes;
    do {
        int index = indexByArc->get(arc);
        indexes.push_back(arc->first < arc->second ? index : ~index);
        arc = arc->next;
    } while (arc);
    return indexes;
}

void Topojson::checkCoordSystem()
{
    bool oversize = bounds[0] < -180 - ε
    || bounds[1] < -90 - ε
    || bounds[2] > 180 + ε
    || bounds[3] > 90 + ε;
}

void Topojson::prequantize()
{
    double x0 = bounds[0];
    double y0 = bounds[1];
    double x1 = bounds[2];
    double y1 = bounds[3];
    
    double kx = x1 -x0 ? (Q1 - 1) / (x1 - x0) * Q0 / Q1 : 1;
    double ky = y1 -y0 ? (Q1 - 1) / (y1 - y0) * Q0 / Q1 : 1;
 
    scale[0] = 1 / kx;
    scale[1] = 1 / ky;
    translate[0] = x0;
    translate[1] = y0;
    
    
    for (int i=0; i < hSHP->nRecords; i++) {
        SHPObject *psShape = SHPReadObject(hSHP, i);
        GeoObject* object = new GeoObject(hSHP->nShapeType, psShape);
        object->quantize(-x0, -y0, kx, ky);
        objects.push_back(object);
        delete psShape;
    }
    
}

void Topojson::postquantize()
{
    if (Q0) {
        if (Q1 == Q0) {
            return;
        }
        double k = Q1 / Q0;
        
        for (int i=0; i < hSHP->nRecords; i++) {
            SHPObject *psShape = SHPReadObject(hSHP, i);
            GeoObject* object = new GeoObject(hSHP->nShapeType, psShape);
            object->quantize(0, 0, k, k);
            delete psShape;
        }
        scale[0] /= k;
        scale[1] /= k;
    } else {
        
        double x0 = bounds[0];
        double x1 = bounds[1];
        double y0 = bounds[2];
        double y1 = bounds[3];
        
        double kx = x1 -x0 ? (Q1 - 1) / (x1 - x0) * Q0 / Q1 : 1;
        double ky = y1 -y0 ? (Q1 - 1) / (y1 - y0) * Q0 / Q1 : 1;
        for (int i=0; i < hSHP->nRecords; i++) {
            SHPObject *psShape = SHPReadObject(hSHP, i);
            GeoObject* object = new GeoObject(hSHP->nShapeType, psShape);
            object->quantize(-x0, -y0, kx, ky);
            delete psShape;
        }
    }
}

void Topojson::stich()
{
    vector<vector<point*> > fragments;
    vector<vector<point*> > fragment_polygon;
    
    int nShapeType = hSHP->nShapeType;
    
    for (int r=0; r < hSHP->nRecords; r++) {
        SHPObject *psShape = SHPReadObject(hSHP, r);
        objects[r]->stitch(nShapeType, psShape, fragments);
        delete psShape;
    }
   
    /*
    // Now stitch the fragments back together into rings.
    // To connect the fragments start-to-end, create a simple index by end.
    map<point*, vector<point*> > fragmentByStart;
    map<point*, vector<point*> > fragmentByEnd;
    map<point*, vector<point*> >::iterator it;
  
    // For each fragment…
    for (int i=0; i < fragments.size(); i++) {
        vector<point*> fragment = fragments[i];
        point* start = fragment[0];
        point* end = fragment[fragment.size() - 1];
        
        // If this fragment is closed, add it as a standalone ring.
        if (start[0] == end[0] && start[1] == end[1]) {
            //fragment_polygon.push_back(fragment);
            fragments[i].empty();
            continue;
        }
        
        //fragment.index = i;
        fragmentByEnd[end] = fragment;
        fragmentByStart[start] = fragment;
    }
    
    // For each open fragment…
    for (int i=0; i < fragments.size(); i++) {
        vector<point*> fragment = fragments[i];
        
        if (fragment.size() > 0) {
            
            point* start = fragment[0];
            point* end = fragment[fragment.size()- 1];
            
            // If this fragment is closed, add it as a standalone ring.
            if (start[0] == end[0] && start[1] == end[1]) {
                //fragment.polygon.push_back(fragment);
                continue;
            }
            
            vector<point*> startFragment = fragmentByEnd[start];
            vector<point*> endFragment = fragmentByStart[end];
            if (startFragment.size() > 0) {
                startFragment.pop_back(); // drop the shared coordinate
                //fragments[startFragment.index] = null;
                fragment = startFragment.concat(fragment);
                fragment.polygon = startFragment.polygon;
                
                if (startFragment === endFragment) {
                    // Connect both ends to this single fragment to create a ring.
                    fragment.polygon.push(fragment);
                } else {
                    fragment.index = n++;
                    fragments.push(fragmentByStart[fragment[0]] = fragmentByEnd[fragment[fragment.length - 1]] = fragment);
                }
            } else if (endFragment.size() > 0) {
                delete fragmentByStart[end];
                delete fragmentByEnd[endFragment[endFragment.length - 1]];
                fragment.pop_back(); // drop the shared coordinate
                fragment = fragment.concat(endFragment);
                fragment.polygon = endFragment.polygon;
                fragment.index = n++;
                fragments[endFragment.index] = null;
                fragments.push(fragmentByStart[fragment[0]] = fragmentByEnd[fragment[fragment.length - 1]] = fragment);
            } else {
                fragment.push(fragment[0]); // close ring
                fragment.polygon.push(fragment);
            }
        }
    }
    */
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

void Topojson::rotateArray(vector<point*>& array, int start, int end, int offset)
{
    reverse(array, start, end);
    reverse(array, start, start + offset);
    reverse(array, start + offset, end);
}

void Topojson::reverse(vector<point*>& array, int start, int end)
{
    for (int mid = start + ((end-- - start) >> 1); start < mid; ++start, --end) {
        point* t = new point(*array[start]);
        array[start] = array[end];
        array[end] = t;
    }
}

bool Topojson::equalLine(Arc* arcA, Arc* arcB)
{
	int ia = arcA->first;
	int ib = arcB->first;
	int ja = arcA->second;
	int jb = arcB->second;
	if (ia - ja != ib - jb) return false;
	for (; ia <= ja; ++ia, ++ib) {
		if (!coordinates[ia]->equal(*coordinates[ib])) return false;
	}
    return true;
}

bool Topojson::reverseEqualLine(Arc* arcA, Arc* arcB)
{
	int ia = arcA->first;
	int ib = arcB->first;
	int ja = arcA->second;
	int jb = arcB->second;
	if (ia - ja != ib - jb) return false;
	for (; ia <= ja; ++ia, --jb) {
		if (!coordinates[ia]->equal(*coordinates[ib])) return false;
	}
    return true;
}

bool Topojson::equalRing(Arc* arcA, Arc* arcB)
{
	int ia = arcA->first;
	int ib = arcB->first;
	int ja = arcA->second;
	int jb = arcB->second;
	int n = ja - jb;
	if (n != jb - ib) return false;
	int ka = findMinimumOffset(arcA);
	int kb = findMinimumOffset(arcB);
	for (int i=0; i < n; ++i) {
		if (!coordinates[ia + (i + ka)%n]->equal(*coordinates[ib + (i+kb)%n]))
            return false;
	}
    return true;
}

int Topojson::findMinimumOffset(Arc* arc)
{
	int start = arc->first;
	int end = arc->second;
	int mid = start;
	int minimum = mid;
	point* minimumPoint = coordinates[mid];
	while (++mid < end) {
      point* pt = coordinates[mid];
      if (pt->x < minimumPoint->x || (pt->x == minimumPoint->x && pt->y < minimumPoint->y)) {
        minimum = mid;
        minimumPoint = pt;
      }
    }
    return minimum - start;
}

void Topojson::Extract()
{
    if (nShapeType == SHPT_NULL || nShapeType == SHPT_POINT ||
        nShapeType == SHPT_POINTM || nShapeType == SHPT_POINTZ ) {
        return;
    }
    int index = -1;
    
    for (int i=0; i < hSHP->nRecords; i++) {
        SHPObject *psShape = SHPReadObject(hSHP, i);
		int nParts = psShape->nParts > 0 ? psShape->nParts : 1;
		for (int j=0; j < nParts; j++ ) {
			int start = psShape->panPartStart[j];
			int end = j < nParts-1 ? psShape->panPartStart[j+1] : psShape->nVertices;
            int n = end - start;
            
			for (int k=start; k < end; k++ ) {
                point* pt = new point(psShape->padfX[k], psShape->padfY[k]);
				coordinates.push_back(pt);
                index++;
			}
           
            Arc* arc = new Arc();
			arc->first = index -n + 1;
			arc->second = index;

			if (nShapeType == SHPT_ARC || nShapeType == SHPT_ARCM ||
				nShapeType == SHPT_ARCZ) {
				lines.push_back(arc);
			} else if (nShapeType == SHPT_POLYGON || nShapeType == SHPT_POLYGONZ ||
				nShapeType == SHPT_POLYGONM) {
				rings.push_back(arc);
			}
            objects[i]->arcs.push_back(arc);
		}
    }
}

Hashset* Topojson::Join()
{
    // index()
    vector<int> indexes;
    //Hashmap indexByPoint(coordinates.size() * 1.4);
    
    for (int i=0; i < coordinates.size(); i++) {
        //int val = indexByPoint.maybeSet(i, i, coordinates);
        indexes.push_back(i);
    }
   
    // init visitedByIndex, leftByIndex, rightByIndex
	for (int i=0; i < coordinates.size(); i++) {
		visitedByIndex.push_back(-1);
		leftByIndex.push_back(-1);
		rightByIndex.push_back(-1);
        junctionByIndex.push_back(0);
	}

	for (int i=0; i < lines.size(); i++) {
		int lineStart = lines[i]->first;
		int lineEnd = lines[i]->second;
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
        Arc* ring = rings[i];
        int ringStart = ring->first + 1;
        int ringEnd = ring->second;
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
            junctionByPoint->add(coordinates[j]);
        }
    }
    
    return junctionByPoint;
}

void Topojson::Cut()
{
    Hashset* junctions = Join();
   
    for (int i=0; i < lines.size(); i++) {
        Arc* line = lines[i];
        int lineMid = line->first;
        int lineEnd = line->second;
        while (++lineMid < lineEnd) {
            if (junctions->has(coordinates[lineMid])) {
                Arc* next = new Arc(lineMid, line->second);
                line->second = lineMid;
                line->next = next;
                line = next;
            }
        }
    }
    
    for (int i=0; i < rings.size(); i++) {
        Arc* ring = rings[i];
        int ringStart = ring->first;
        int ringMid = ringStart;
        int ringEnd = ring->second;
        bool ringFixed = junctions->has(coordinates[ringStart]);
        
        while (++ringMid < ringEnd) {
            if (junctions->has(coordinates[ringMid])) {
                if (ringFixed) {
                    Arc* next = new Arc(ringMid, ring->second);
                    ring->second = ringMid;
                    ring->next = next;
                    ring = next;
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
        Arc* line = lines[i];
        Arc* nextline = line->next;
        while (nextline) {
            ++arcCount;
            nextline = nextline->next;
        }
    }
    for (int i=0; i < rings.size(); i++) {
        Arc* ring = rings[i];
        Arc* nextring = ring->next;
        while (nextring) {
            ++arcCount;
            nextring = nextring->next;
        }
    }
    
    arcsByEnd = new PointHashmap(arcCount * 2 * 1.4);
    
    for (int i=0; i < lines.size(); i++) {
        Arc* nextline = lines[i];
        do {
            dedupLine(nextline);
            nextline = nextline->next;
        } while (nextline);
    }
    
    for (int i=0; i < rings.size(); i++) {
        Arc* ring = rings[i];
        if (ring->next) {  // arc is no longer closed
            Arc* nextring = ring;
            do {
                dedupLine(nextring);
                nextring = nextring->next;
            } while (nextring);
        } else {
            dedupRing(ring);
        }
    }
}

void Topojson::dedupLine(Arc* arc)
{
	// Does this arc match an existing arc in order?
	point* startPoint = coordinates[arc->first];
	vector<Arc*> startArcs = arcsByEnd->get(startPoint);
	for (int i=0; i < startArcs.size(); i++) {
		Arc* startArc = startArcs[i];	
		if (equalLine(startArc, arc)) {
			arc->first = startArc->first;
			arc->second = startArc->second;
			return;
		}
	}

	// Does this arc match an existing arc in reverse order?
	point* endPoint = coordinates[arc->second];
	vector<Arc*> endArcs = arcsByEnd->get(endPoint);
	for (int i=0; i < endArcs.size(); i++) {
		Arc* endArc = endArcs[i];	
		if (reverseEqualLine(endArc, arc)) {
			arc->first = endArc->first;
			arc->second = endArc->second;
			return;
		}
	}

	if (startArcs.size() > 0) {
		startArcs.push_back(arc);
	} else {
        vector<Arc*> tempArcs;
        tempArcs.push_back(arc);
        arcsByEnd->set(startPoint, tempArcs);
    }
    
	if (endArcs.size() > 0) {
		endArcs.push_back(arc);
	} else {
        vector<Arc*> tempArcs;
        tempArcs.push_back(arc);
        arcsByEnd->set(endPoint, tempArcs);
    }
    
    arcs.push_back(arc);
}


void Topojson::dedupRing(Arc* arc)
{
    // Does this arc match an existing line in order, or reverse order?
    // Rings are closed, so their start point and end point is the same.
    point* endPoint = coordinates[arc->first];
    vector<Arc*> endArcs = arcsByEnd->get(endPoint);
    for (int i=0; i < endArcs.size(); i++) {
        Arc* endArc = endArcs[i];
        if (equalRing(endArc, arc)) {
            arc->first = endArc->first;
            arc->second = endArc->second;
            return;
        }
        if (reverseEqualLine(endArc, arc)) {
            arc->first = endArc->second;
            arc->second = endArc->first;
            return;
        }
    }
    
    // Otherwise, does this arc match an existing ring in order, or reverse order?
    endPoint = coordinates[arc->first + findMinimumOffset(arc)];
    endArcs.empty();
    endArcs = arcsByEnd->get(endPoint);
    for (int i=0; i < endArcs.size(); i++) {
        Arc* endArc = endArcs[i];
        if (equalRing(endArc, arc)) {
            arc->first = endArc->first;
            arc->second = endArc->second;
            return;
        }
        if (reverseEqualLine(endArc, arc)) {
            arc->first = endArc->second;
            arc->second = endArc->first;
            return;
        }
    }
    
    if (endArcs.size() > 0) {
        endArcs.push_back(arc);
    } else {
        vector<Arc*> tempArcs;
        tempArcs.push_back(arc);
        arcsByEnd->set(endPoint, tempArcs);
    }
    
    arcs.push_back(arc);
}

