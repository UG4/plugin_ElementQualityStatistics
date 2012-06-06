/*
 * element_quality_statistics.cpp
 *
 *  Created on: 17.04.2012
 *      Author: Martin Stepniewski
 */


#ifndef __ElementQualityStatistics_h__
#define __ElementQualityStatistics_h__

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
#include "volume_calculation.h"



using namespace std;


namespace ug {


////////////////////////////////////////////////////////////////////////////////////////////
//	CollectAssociatedSides
////////////////////////////////////////////////////////////////////////////////////////////

///	Collects all edges (= 2) which exist in the given face and which share the given vertex.
UG_API
inline void CollectAssociatedSides(EdgeBase* sidesOut[2], Grid& grid, Face* f, VertexBase* vrt);

///	Collects all faces (= 2) which exist in the given volume and which share the given edge.
UG_API
inline void CollectAssociatedSides(Face* sidesOut[2], Grid& grid, Volume* v, EdgeBase* e);






////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinAngle
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMinAngle(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Face* f, TAAPosVRT& aaPos);

//	Tetrahedron
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos);

//	Prism
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Prism* prism, TAAPosVRT& aaPos);

//	Pyramid
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos);

//	Volume
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Volume* v, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinDihedral
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Face* f, TAAPosVRT& aaPos);

//	Tetrahedron
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos);

//	Prism
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Prism* prism, TAAPosVRT& aaPos);

//	Pyramid
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos);

//	Volume
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Volume* v, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMaxAngle
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Face* f, TAAPosVRT& aaPos);

//	Tetrahedron
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos);

//	Prism
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Prism* prism, TAAPosVRT& aaPos);

//	Pyramid
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos);

//	Volume
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Volume* v, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMaxDihedral
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Face* f, TAAPosVRT& aaPos);

//	Tetrahedron
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos);

//	Prism
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Prism* prism, TAAPosVRT& aaPos);

//	Pyramid
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos);

//	Volume
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Volume* v, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinTriangleHeight
template <class TAAPosVRT>
number CalculateMinTriangleHeight(Face* face, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateAspectRatio
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Constrained Triangles)
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Face* face, TAAPosVRT& aaPos);

//	Tetrahedron
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos);

/*
//	Prism
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Prism* prism, TAAPosVRT& aaPos);

//	Pyramid
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos);
*/

//	Volume
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Volume* vol, TAAPosVRT& aaPos);




////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithSmallestMinAngle
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithSmallestMinAngle(Grid& grid, 	TIterator volumesBegin,
												TIterator volumesEnd, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	FindVolumeWithSmallestMinDihedral
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindVolumeWithSmallestMinDihedral(Grid& grid, 	TIterator elementsBegin,
												TIterator elementsEnd, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithLargestMaxAngle
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithLargestMaxAngle(Grid& grid, 	TIterator volumesBegin,
											TIterator volumesEnd, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	FindVolumeWithLargestMaxDihedral
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindVolumeWithLargestMaxDihedral(Grid& grid, 	TIterator elementsBegin,
												TIterator elementsEnd, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	FindLargestFace
template <class TIterator, class TAAPosVRT>
Face* FindLargestFace(TIterator facesBegin, TIterator facesEnd, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	FindSmallestVolume
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindSmallestVolume(TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	FindLargestVolume
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindLargestVolume(TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithSmallestAspectRatio
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithSmallestAspectRatio(Grid& grid, 	TIterator elemsBegin,
												TIterator elemsEnd, TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithLargestAspectRatio
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithLargestAspectRatio(Grid& grid, 	TIterator elemsBegin,
												TIterator elemsEnd, TAAPosVRT& aaPos);






////////////////////////////////////////////////////////////////////////////////////////////
//	MinAngleHistogram
template <class TIterator, class TAAPosVRT>
void MinAngleHistogram(Grid& grid, 	TIterator elementsBegin,
									TIterator elementsEnd,
									TAAPosVRT& aaPos,
									uint stepSize);


////////////////////////////////////////////////////////////////////////////////////////////
//	ElementQualityStatistics
////////////////////////////////////////////////////////////////////////////////////////////

//	Wrapper
//void ElementQualityStatistics(MultiGrid& mg, int level);
void ElementQualityStatistics(MultiGrid& mg);

void ElementQualityStatistics(Grid& grid);

//	Actual procedures
void ElementQualityStatistics2d(Grid& grid, GeometricObjectCollection goc);
void ElementQualityStatistics3d(Grid& grid, GeometricObjectCollection goc);




}	 
#endif

