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

using namespace std;

Topojson::Topojson(const char* shp_path)
{
	this->index = -1;
    this->hSHP = SHPOpen( shp_path, "rb" );
    if( this->hSHP == NULL )
    {
        printf( "Unable to open:%s\n", shp_path );
		return;
    }
    SHPGetInfo( this->hSHP, &this->nEntities, &this->nShapeType, this->adfMinBound, this->adfMaxBound );

	this->Extract();
}

Topojson::~Topojson()
{
	SHPClose(this->hSHP);
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

			map<int, int> arc;
			arc[0] = index + 1;
			index += end - start;
			arc[1] = index;

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

void Topojson::Join()
{
	vector<int> visitedByIndex;
	vector<int> leftByIndex;
	vector<int> rightByIndex;
	vector<int> junctionByIndex;
	int junctionCount = 0;

	for (int i=0; i < this->coordinates.size(); i++) {
		visitedByIndex.push_back(-1);
		leftByIndex.push_back(-1);
		rightByIndex.push_back(-1);
	}

	for (int i=0; i < this->lines.size(); i++) {
		int lineStart = lines[i][0];
		int lineEnd = lines[i][1];
		int currentIndex =0; 
	}
}