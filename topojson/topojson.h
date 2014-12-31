//
//  hashtable.h
//  topojson
//
//  Created by Xun Li on 12/19/14.
//  Copyright (c) 2014 Xun Li. All rights reserved.
//

#ifndef __topojson__topojson__
#define __topojson__topojson__

#include <vector>
#include <map>
#include "point.h"
#include "shapefil.h"
#include "hashtable.h"

using namespace std;

class SphericalCoordSystemException: public exception
{
    virtual const char* what() const throw()
    {
        return "spherical coordinates outside of [±180°, ±90°]";
    }
};

class GeoObject {
    int length;
    double x0;
    double y0;
    double kx;
    double ky;
    vector<point*> coords;
public:
    vector<Arc*> arcs;
    vector<vector<int> > index_arcs;
    
public:
    GeoObject(int nShapeType, SHPObject* psShape);
    ~GeoObject();
    
    void quantizePoint(point* pt);
    void quantize(double _x0, double _y0, double _kx, double _ky);
    
    void normalizePoint(point *pt, double y0e, double y1e, double x0, double y0, double y1);
    void stitch(int nShapeType, SHPObject* psShape, vector<vector<point*> >& fragments);
};

class Topojson {
    //Shapefile
	SHPHandle hSHP;
	int nShapeType;
	int nEntities;
	double bounds[4];
    vector<GeoObject*> objects;
    
    //
    double scale[2];
    double translate[2];
    vector<Arc*> lines;
    vector<Arc*> rings;
    vector<point*> coordinates;
    vector<Arc*> arcs;
    vector<vector<point*> > topology_arcs;
    
    //vector<vector<int> > objects_indexes;
    
    ArcHashmap* indexByArc;
	PointHashmap* arcsByEnd;
    
    //Join
	vector<int> visitedByIndex;
	vector<int> leftByIndex;
	vector<int> rightByIndex;
	vector<int> junctionByIndex;
	int junctionCount; // upper bound on number of junctions
   
    double ε = 1e-6;
    double Q0 = 1e4;
    double Q1 = 1e4;
    double minimumArea = 0.0;
    bool stitchPoles = true;

public:
	Topojson(const char* shp_path);
	~Topojson();

    void Save();

private:
    void checkCoordSystem();
    void prequantize();
    void stich();
    void postquantize();
    void computeTopology();
private:
    void Delta();
	void Extract();
    
	Hashset* Join();
    
	void Cut();
    
    void Dedup();
    
    void sequence(int i, int previousIndex, int currentIndex, int nextIndex);
    
    void rotateArray(vector<point*>& array, int start, int end, int offset);
    
    void reverse(vector<point*>& array, int start, int end);

	void dedupLine(Arc* line);

	void dedupRing(Arc* ring);

	bool equalLine(Arc* arcA, Arc* arcB);
    
	bool reverseEqualLine(Arc* arcA, Arc* arcB);
    
	bool equalRing(Arc* arcA, Arc* arcB);

	int findMinimumOffset(Arc* arc);
    
    vector<int> indexArcs(Arc* arc);
};

#endif /* defined(__topojson__) */
