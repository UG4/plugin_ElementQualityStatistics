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
#include "volume_calculation.h"



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
	UG_LOG(endl << "MinAngle-Histogram for '" << refElem->base_object_id() << "d' elements");
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
//	SubdivisionTetGridSmooth
////////////////////////////////////////////////////////////////////////////////////////////

void SubdivisionTetGridSmooth(MultiGrid& mg, MGSubsetHandler& sh);
void SubdivisionTetGridSmoothBasic(MultiGrid& mg, MGSubsetHandler& sh);
void MoveVertexToSmoothTetGridSubdivisionPosition(MultiGrid& mg, Vertex* vrt, 	Grid::VertexAttachmentAccessor<APosition>& aaPos,
																				Grid::VertexAttachmentAccessor<APosition>& aaSmoothPos);

////////////////////////////////////////////////////////////////////////////////////////////
//	ElementQualityStatistics
////////////////////////////////////////////////////////////////////////////////////////////

//	Wrapper
//void ElementQualityStatistics(MultiGrid& mg, int level);
void ElementQualityStatistics(MultiGrid& mg);

void ElementQualityStatistics(Grid& grid);

//	Actual procedures
void ElementQualityStatistics2d(Grid& grid, GridObjectCollection goc);
void ElementQualityStatistics3d(Grid& grid, GridObjectCollection goc);


}	 
#endif  //__ELEMENT_QUALITY_STATISTICS_H__

