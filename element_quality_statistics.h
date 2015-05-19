/*
 * element_quality_statistics.cpp
 *
 *  Created on: 17.04.2012
 *      Author: Martin Stepniewski
 */


#ifndef __ELEMENT_QUALITY_STATISTICS_H__
#define __ELEMENT_QUALITY_STATISTICS_H__

/* system includes */
#include <stddef.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "lib_grid/lib_grid.h"
#include "elem_stat_util.h"
//#include "volume_calculation.h"



using namespace std;


namespace ug {


////////////////////////////////////////////////////////////////////////////////////////////
//	MinAngleHistogram
template <class TIterator, class TAAPosVRT>
void MinAngleHistogram(Grid& grid, 	TIterator elementsBegin,
									TIterator elementsEnd,
									TAAPosVRT& aaPos,
									uint stepSize)
{
	//PROFILE_FUNC();
//	Initialization
	vector<number> minAngles;
	typename TIterator::value_type refElem = *elementsBegin;

//	Calculate the minAngle of every element
	for(TIterator iter = elementsBegin; iter != elementsEnd; ++iter)
	{
		number curMinAngle = CalculateMinAngle(grid, *iter, aaPos);
		minAngles.push_back(curMinAngle);
	}

//	Sort the calculated minAngles in an ascending way
	sort (minAngles.begin(), minAngles.end());


//	Evaluate the minimal and maximal degree rounding to 10
	int minDeg = round(number(minAngles.front()) / 10.0) * 10;
	int maxDeg = round(number(minAngles.back()) / 10.0) * 10;

//	Expand minDeg and maxDeg by plus minus 10 degrees or at least to 0 or 180 degress
	if((minDeg-10) > 0)
		minDeg = minDeg - 10;
	else
		minDeg = 0;

	if((maxDeg+10) < 180)
		maxDeg = maxDeg + 10;
	else
		maxDeg = 180;

//	Evaluate the number of ranges in respect to the specified step size
	uint numRanges = floor((maxDeg-minDeg) / stepSize);
	vector<uint> counter(numRanges, 0);

//	Count the elements in their corresponding minAngle range
	for(uint i = 0; i < minAngles.size(); ++i)
	{
		number minAngle = minAngles[i];
		for (uint range = 0; range < numRanges; range++)
		{
			if (minAngle < minDeg + (range+1)*stepSize)
			{
				++counter[range];
				break;
			}
		}
	}

//	----------------------------------------
//	Histogram table output section: (THIRDS)
//	----------------------------------------

//	Divide the output table into three thirds (columnwise)
	uint numRows = ceil(number(numRanges) / 3.0);

//	Create table object
	ug::Table<std::stringstream> minAngleTable(numRows, 6);

//	Specific element header
	//UG_LOG(endl << "MinAngle-Histogram for '" << refElem->reference_object_id() << "' elements");
	UG_LOG(endl << "(*) MinAngle-Histogram for '" << refElem->base_object_id() << "d' elements");
	UG_LOG(endl);

//	First third
	uint i = 0;
	for(; i < numRows; ++i)
	{
		minAngleTable(i, 0) << minDeg + i*stepSize << " - " << minDeg + (i+1)*stepSize << " deg : ";
		minAngleTable(i, 1) << counter[i];
	}

//	Second third
//	Check, if second third of table is needed
	if(i < counter.size())
	{
		for(; i < 2*numRows; ++i)
		{
			minAngleTable(i-numRows, 2) << minDeg + i*stepSize << " - " << minDeg + (i+1)*stepSize << " deg : ";
			minAngleTable(i-numRows, 3) << counter[i];
		}
	}

//	Third third
	if(i < counter.size())
//	Check, if third third of table is needed
	{
		for(; i < numRanges; ++i)
		{
			minAngleTable(i-2*numRows, 4) << minDeg + i*stepSize << " - " << minDeg + (i+1)*stepSize << " deg : ";
			minAngleTable(i-2*numRows, 5) << counter[i];
		}
	}

//	Output table
	UG_LOG(endl << minAngleTable);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	MaxAngleHistogram
template <class TIterator, class TAAPosVRT>
void MaxAngleHistogram(Grid& grid, 	TIterator elementsBegin,
                       TIterator elementsEnd,
                       TAAPosVRT& aaPos,
                       uint stepSize)
{
    //PROFILE_FUNC();
    //	Initialization
    vector<number> maxAngles;
    typename TIterator::value_type refElem = *elementsBegin;

    //	Calculate the minAngle of every element
    for(TIterator iter = elementsBegin; iter != elementsEnd; ++iter)
    {
        number curMaxAngle = CalculateMaxAngle(grid, *iter, aaPos);
        maxAngles.push_back(curMaxAngle);
    }

    //	Sort the calculated minAngles in an ascending way
    sort (maxAngles.begin(), maxAngles.end());


    //	Evaluate the minimal and maximal degree rounding to 10
    int minDeg = round(number(maxAngles.front()) / 10.0) * 10;
    int maxDeg = round(number(maxAngles.back()) / 10.0) * 10;

    //	Expand minDeg and maxDeg by plus minus 10 degrees or at least to 0 or 180 degress
    if((minDeg-10) > 0)
        minDeg = minDeg - 10;
    else
        minDeg = 0;

    if((maxDeg+10) < 180)
        maxDeg = maxDeg + 10;
    else
        maxDeg = 180;

    //	Evaluate the number of ranges in respect to the specified step size
    uint numRanges = floor((maxDeg-minDeg) / stepSize);
    vector<uint> counter(numRanges, 0);

    //	Count the elements in their corresponding minAngle range
    for(uint i = 0; i < maxAngles.size(); ++i)
    {
        number minAngle = maxAngles[i];
        for (uint range = 0; range < numRanges; range++)
        {
            if (minAngle < minDeg + (range+1)*stepSize)
            {
                ++counter[range];
                break;
            }
        }
    }

    //	----------------------------------------
    //	Histogram table output section: (THIRDS)
    //	----------------------------------------

    //	Divide the output table into three thirds (columnwise)
    uint numRows = ceil(number(numRanges) / 3.0);

    //	Create table object
    ug::Table<std::stringstream> maxAngleTable(numRows, 6);

    //	Specific element header
    //UG_LOG(endl << "MaxAngle-Histogram for '" << refElem->reference_object_id() << "' elements");
    UG_LOG(endl << "(*) MaxAngle-Histogram for '" << refElem->base_object_id() << "d' elements");
    UG_LOG(endl);

    //	First third
    uint i = 0;
    for(; i < numRows; ++i)
    {
        maxAngleTable(i, 0) << minDeg + i*stepSize << " - " << minDeg + (i+1)*stepSize << " deg : ";
        maxAngleTable(i, 1) << counter[i];
    }

    //	Second third
    //	Check, if second third of table is needed
    if(i < counter.size())
    {
        for(; i < 2*numRows; ++i)
        {
            maxAngleTable(i-numRows, 2) << minDeg + i*stepSize << " - " << minDeg + (i+1)*stepSize << " deg : ";
            maxAngleTable(i-numRows, 3) << counter[i];
        }
    }

    //	Third third
    if(i < counter.size())
        //	Check, if third third of table is needed
    {
        for(; i < numRanges; ++i)
        {
            maxAngleTable(i-2*numRows, 4) << minDeg + i*stepSize << " - " << minDeg + (i+1)*stepSize << " deg : ";
            maxAngleTable(i-2*numRows, 5) << counter[i];
        }
    }

    //	Output table
    UG_LOG(endl << maxAngleTable);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	PrintAngleStatistics2d
template <class TAAPosVRT>
void PrintAngleStatistics2d(Grid& grid, GridObjectCollection& goc, int level, TAAPosVRT& aaPos)
{
	int i = level;

//	Calculate and output standard deviation for triangular/quadrilateral angles
	if(goc.num<Triangle>(i) > 0 || goc.num<Quadrilateral>(i) > 0)
	{
		number sd_tri = 0.0;
		number sd_quad = 0.0;
		number mean_tri = 0.0;
		number mean_quad = 0.0;
		number regAngle = 90.0;

		vector<number> vAnglesOut;

		for(FaceIterator fIter = goc.begin<Face>(i); fIter != goc.end<Face>(i); ++fIter)
		{
			vAnglesOut.clear();
			Face* f = *fIter;

			if(f->reference_object_id() == ROID_TRIANGLE)
			{
				regAngle = 60.0;
				CalculateFaceAngles(vAnglesOut, grid, f, aaPos);

				for(size_t k = 0; k < vAnglesOut.size(); ++k)
				{
					sd_tri += (regAngle-vAnglesOut[k])*(regAngle-vAnglesOut[k]);
					mean_tri += vAnglesOut[k];
				}
			}

			if(f->reference_object_id() == ROID_QUADRILATERAL)
			{
				regAngle = 90.0;
				CalculateFaceAngles(vAnglesOut, grid, f, aaPos);

				for(size_t k = 0; k < vAnglesOut.size(); ++k)
				{
					sd_quad += (regAngle-vAnglesOut[k])*(regAngle-vAnglesOut[k]);
					mean_quad += vAnglesOut[k];
				}
			}
		}

		if(goc.num<Triangle>(i) > 0)
		{
			sd_tri *= (1.0/(3*goc.num<Triangle>(i)));
			sd_tri = sqrt(sd_tri);
			mean_tri *= (1.0/(3*goc.num<Triangle>(i)));
		}

		if(goc.num<Quadrilateral>(i) > 0)
		{
			sd_quad *= (1.0/(4*goc.num<Quadrilateral>(i)));
			sd_quad = sqrt(sd_quad);
			mean_quad *= (1.0/(4*goc.num<Quadrilateral>(i)));
		}

		UG_LOG("(*) Standard deviation of face angles to regular case" << endl);
		UG_LOG("	(60° for triangles, 90° for quadrilaterals)" << endl);
		UG_LOG(endl);
		UG_LOG("	Triangles (" << goc.num<Triangle>(i) << "):" << endl);
		if(goc.num<Triangle>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_tri << endl);
			UG_LOG("		mean = " << mean_tri << endl);
		}
		UG_LOG(endl);
		UG_LOG("	Quadrilaterals (" << goc.num<Quadrilateral>(i) << "):" << endl);
		if(goc.num<Quadrilateral>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_quad << endl);
			UG_LOG("		mean = " << mean_quad << endl);
		}
		UG_LOG(endl);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	PrintAngleStatistics3d
template <class TAAPosVRT>
void PrintAngleStatistics3d(Grid& grid, GridObjectCollection& goc, int level, TAAPosVRT& aaPos)
{
	int i = level;

//	Calculate and output standard deviation for triangular/quadrilateral angles
	if(goc.num<Triangle>(i) > 0 || goc.num<Quadrilateral>(i) > 0)
	{
		number sd_tri = 0.0;
		number sd_quad = 0.0;
		number mean_tri = 0.0;
		number mean_quad = 0.0;
		number regAngle = 90.0;

		vector<number> vAnglesOut;

		for(FaceIterator fIter = goc.begin<Face>(i); fIter != goc.end<Face>(i); ++fIter)
		{
			vAnglesOut.clear();
			Face* f = *fIter;

			if(f->reference_object_id() == ROID_TRIANGLE)
			{
				regAngle = 60.0;
				CalculateFaceAngles(vAnglesOut, grid, f, aaPos);

				for(size_t k = 0; k < vAnglesOut.size(); ++k)
				{
					sd_tri += (regAngle-vAnglesOut[k])*(regAngle-vAnglesOut[k]);
					mean_tri += vAnglesOut[k];
				}
			}

			if(f->reference_object_id() == ROID_QUADRILATERAL)
			{
				regAngle = 90.0;
				CalculateFaceAngles(vAnglesOut, grid, f, aaPos);

				for(size_t k = 0; k < vAnglesOut.size(); ++k)
				{
					sd_quad += (regAngle-vAnglesOut[k])*(regAngle-vAnglesOut[k]);
					mean_quad += vAnglesOut[k];
				}
			}
		}

		if(goc.num<Triangle>(i) > 0)
		{
			sd_tri *= (1.0/(3*goc.num<Triangle>(i)));
			sd_tri = sqrt(sd_tri);
			mean_tri *= (1.0/(3*goc.num<Triangle>(i)));
		}

		if(goc.num<Quadrilateral>(i) > 0)
		{
			sd_quad *= (1.0/(4*goc.num<Quadrilateral>(i)));
			sd_quad = sqrt(sd_quad);
			mean_quad *= (1.0/(4*goc.num<Quadrilateral>(i)));
		}

		UG_LOG("(*) Standard deviation of face angles to regular case" << endl);
		UG_LOG("	(60° for triangles, 90° for quadrilaterals)" << endl);
		UG_LOG(endl);
		UG_LOG("	Triangles (" << goc.num<Triangle>(i) << "):" << endl);
		if(goc.num<Triangle>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_tri << endl);
			UG_LOG("		mean = " << mean_tri << endl);
		}
		UG_LOG(endl);
		UG_LOG("	Quadrilaterals (" << goc.num<Quadrilateral>(i) << "):" << endl);
		if(goc.num<Quadrilateral>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_quad << endl);
			UG_LOG("		mean = " << mean_quad << endl);
		}
		UG_LOG(endl);
	}


//	Calculate and output standard deviation for tetrahedral/hexahedral angles
	if(goc.num<Tetrahedron>(i) > 0 || goc.num<Hexahedron>(i) > 0)
	{
		number sd_tet = 0.0;
		number sd_hex = 0.0;
		number mean_tet = 0.0;
		number mean_hex = 0.0;
		number regVolDihedral = 90.0;

		vector<number> vDihedralsOut;

		for(VolumeIterator vIter = goc.begin<Volume>(i); vIter != goc.end<Volume>(i); ++vIter)
		{
			vDihedralsOut.clear();
			Volume* vol = *vIter;

			if(vol->reference_object_id() == ROID_TETRAHEDRON)
			{
				regVolDihedral = 70.5288;
				CalculateVolumeDihedrals(vDihedralsOut, grid, vol, aaPos);

				for(size_t k = 0; k < vDihedralsOut.size(); ++k)
				{
					sd_tet += (regVolDihedral-vDihedralsOut[k])*(regVolDihedral-vDihedralsOut[k]);
					mean_tet += vDihedralsOut[k];
				}
			}

			if(vol->reference_object_id() == ROID_HEXAHEDRON)
			{
				regVolDihedral = 90.0;
				CalculateVolumeDihedrals(vDihedralsOut, grid, vol, aaPos);

				for(size_t k = 0; k < vDihedralsOut.size(); ++k)
				{
					sd_hex += (regVolDihedral-vDihedralsOut[k])*(regVolDihedral-vDihedralsOut[k]);
					mean_hex += vDihedralsOut[k];
				}
			}
		}

		if(goc.num<Tetrahedron>(i) > 0)
		{
			sd_tet *= (1.0/(6*goc.num<Tetrahedron>(i)));
			sd_tet = sqrt(sd_tet);
			mean_tet *= (1.0/(6*goc.num<Tetrahedron>(i)));
		}

		if(goc.num<Hexahedron>(i) > 0)
		{
			sd_hex *= (1.0/(12*goc.num<Hexahedron>(i)));
			sd_hex = sqrt(sd_hex);
			mean_hex *= (1.0/(12*goc.num<Hexahedron>(i)));
		}

		UG_LOG("(*) Standard deviation of dihedral angles to regular case" << endl);
		UG_LOG("	(70.5288° for tetrahedrons, 90° for hexahedrons)" << endl);
		UG_LOG(endl);
		UG_LOG("	Tetrahedrons (" << goc.num<Tetrahedron>(i) << "):" << endl);
		if(goc.num<Tetrahedron>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_tet << endl);
			UG_LOG("		mean = " << mean_tet << endl);
		}
		UG_LOG(endl);
		UG_LOG("	Hexahedrons (" << goc.num<Hexahedron>(i) << "):" << endl);
		if(goc.num<Hexahedron>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_hex << endl);
			UG_LOG("		mean = " << mean_hex << endl);
		}
		UG_LOG(endl);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	ElementQualityStatistics
////////////////////////////////////////////////////////////////////////////////////////////

//	Wrapper
void ElementQualityStatistics(MultiGrid& mg, int dim);

void ElementQualityStatistics(Grid& grid, int dim);

//	Actual procedures
void ElementQualityStatistics2d(Grid& grid, GridObjectCollection goc);
void ElementQualityStatistics3d(Grid& grid, GridObjectCollection goc);


}	 
#endif  //__ELEMENT_QUALITY_STATISTICS_H__

