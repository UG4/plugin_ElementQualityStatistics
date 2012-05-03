/*
 * element_quality_statistics.cpp
 *
 *  Created on: 17.04.2012
 *      Author: Martin Stepniewski
 */


#include "element_quality_statistics.h"
#include "common/util/table.h"

#ifdef __GXX_EXPERIMENTAL_CXX0X__
	#include <functional>
	template <class T> function<bool(T)> isInRange(T low, T high)
	{
		return [low,high](T value)
		{
			return std::fabs(value - low) >= std::numeric_limits<T>::epsilon() \
                   && std::fabs(value - high) <= std::numeric_limits<T>::epsilon();
		};
	}
#endif

namespace ug
{



////////////////////////////////////////////////////////////////////////////////////////////
//	CollectAssociatedSides
////////////////////////////////////////////////////////////////////////////////////////////

///	Collects all edges (= 2) which exist in the given face and which share the given vertex.
UG_API
inline void CollectAssociatedSides(EdgeBase* sidesOut[2], Grid& grid, Face* f, VertexBase* vrt)
{
	vector<EdgeBase*> vNeighbourEdgesToVertex;
	sidesOut[0] = NULL;
	sidesOut[1] = NULL;

	CollectEdges(vNeighbourEdgesToVertex, grid, vrt, true);
	for(uint i = 0; i < vNeighbourEdgesToVertex.size(); ++i)
	{
		if(FaceContains(f, vNeighbourEdgesToVertex[i]) == true)
		{
			UG_ASSERT(	sidesOut[1] == NULL,
						"Only two edges may be adjacent to a vertex in a face element.")

			if(sidesOut[0] == NULL)
				sidesOut[0] = vNeighbourEdgesToVertex[i];
			else
				sidesOut[1] = vNeighbourEdgesToVertex[i];
		}
	}

	UG_ASSERT(	sidesOut[1] != NULL,
				"Exactly two edges should be adjacent to a vertex in a face element.")
}

///	Collects all faces (= 2) which exist in the given volume and which share the given edge.
UG_API
inline void CollectAssociatedSides(Face* sidesOut[2], Grid& grid, Volume* v, EdgeBase* e)
{
	vector<Face*> vNeighbourFacesToEdge;
	sidesOut[0] = NULL;
	sidesOut[1] = NULL;

	CollectFaces(vNeighbourFacesToEdge, grid, e, true);
	for(uint i = 0; i < vNeighbourFacesToEdge.size(); ++i)
	{
		if(VolumeContains(v, vNeighbourFacesToEdge[i]) == true)
		{
			UG_ASSERT(	sidesOut[1] == NULL,
						"Only two faces may be adjacent to an edge in a volume element.")

			if(sidesOut[0] == NULL)
				sidesOut[0] = vNeighbourFacesToEdge[i];
			else
				sidesOut[1] = vNeighbourFacesToEdge[i];
		}
	}

	UG_ASSERT(	sidesOut[1] != NULL,
				"Exactly two faces should be adjacent to an edge in a volume element.")
}






////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinAngle

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMinAngle(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Face* f, TAAPosVRT& aaPos)
{
//	in the current implementation this method requires, that all edges
//	are created for all faces.
//TODO: improve this!
	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
	{
		LOG("WARNING: autoenabling GRIDOPT_AUTOGENERATE_SIDES in GetNeighbours(Face).\n");
		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
	}

//	Get type of vertex attachment in aaPos and define it as ValueType
	typedef typename TAAPosVRT::ValueType ValueType;

//	Initialization
	uint numFaceVrts = f->num_vertices();
	ValueType vNorm1, vNorm2;
	ValueType vDir1, vDir2;
	number minAngle = 180.0;
	number tmpAngle;
	EdgeBase* vNeighbourEdgesToVertex[2];

//	Iterate over all face vertices
	for(uint vrtIter = 0; vrtIter < numFaceVrts; ++vrtIter)
	{
		VertexBase* vrt = f->vertex(vrtIter);

	//	Get adjacent edges at the current vertex and calculate the angle between their normals
		CollectAssociatedSides(vNeighbourEdgesToVertex, grid, f, vrt);

	//	Calculate direction vectors of the current two adjacent edges
	//	!!!	Beware of the correct order of the vertices	to get correct angle value !!!
		VecSubtract(vDir1,
					aaPos[vNeighbourEdgesToVertex[0]->vertex(0)],
					aaPos[vNeighbourEdgesToVertex[0]->vertex(1)]);

		VecSubtract(vDir2,
					aaPos[vNeighbourEdgesToVertex[1]->vertex(1)],
					aaPos[vNeighbourEdgesToVertex[1]->vertex(0)]);

	//	Normalize
		VecNormalize(vDir1, vDir1);
		VecNormalize(vDir2, vDir2);

	//	Calculate current angle
		tmpAngle = acos(VecDot(vDir1, vDir2));

	//	Check for minimality
		if(tmpAngle < minAngle)
		{
			minAngle = tmpAngle;
		}
	}

//	Transform minAngle from RAD to DEG
	minAngle = 180/PI * minAngle;

	return minAngle;
}

//	Tetrahedron
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	return CalculateMinAngle(grid, static_cast<Volume*>(tet), aaPos);
}

//	Prism
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	return CalculateMinAngle(grid, static_cast<Volume*>(prism), aaPos);
}

//	Pyramid
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	return CalculateMinAngle(grid, static_cast<Volume*>(pyr), aaPos);
}

//	Volume (For volume elements the smallest dihedral will be calculated.)
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Volume* v, TAAPosVRT& aaPos)
{
//	in the current implementation this method requires, that all edges
//	are created for all faces.
//TODO: improve this!
	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
	{
		LOG("WARNING: autoenabling GRIDOPT_AUTOGENERATE_SIDES in GetNeighbours(Face).\n");
		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
	}

//	Initialization
	uint numElementEdges = v->num_edges();
	vector3 vNorm1, vNorm2;
	number minAngle = 360.0;
	number tmpAngle;
	Face* vNeighbourFacesToEdge[2];

//	Iterate over all element edges
	for(uint eIter = 0; eIter < numElementEdges; ++eIter)
	{
		EdgeBase* e = grid.get_edge(v, eIter);

	//	Get adjacent faces at the current edge and calculate the angle between their normals
		CollectAssociatedSides(vNeighbourFacesToEdge, grid, v, e);

		CalculateNormal(vNorm1, vNeighbourFacesToEdge[0], aaPos);
		CalculateNormal(vNorm2, vNeighbourFacesToEdge[1], aaPos);

	/*	!!!	Beware of the correct direction normals to get correct angle value !!!
		INFO:	Angles of a regular tetrahedron:
				- "Tetrahedron (dihedral) angle x" 	= 109,471...deg
				- "Face-to-face angle y"			= 180deg - x = 70,52...deg

		--> Old version:
		//VecScale(vNorm1, vNorm1, -1);
		//tmpAngle = acos(VecDot(vNorm1, vNorm2));

		--> New version:
			(s.	"Qualitäts-Metriken und Optimierung von Tetraeder-Netzen",
				 Seminararbeit von Johannes Ahlmann, Universität Karlsruhe)
	*/
		tmpAngle = acos(VecDot(vNorm1, vNorm2));
		tmpAngle = PI - tmpAngle;

	//	Check for minimality
		if(tmpAngle < minAngle)
		{
			minAngle = tmpAngle;
		}
	}

//	Transform minAngle from RAD to DEG
	minAngle = 180/PI * minAngle;

	return minAngle;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinTriangleHeight
template <class TAAPosVRT>
number CalculateMinTriangleHeight(Face* face, TAAPosVRT& aaPos)
{
	if(face->num_vertices() == 3)
	{
	//	Get type of vertex attachment in aaPos and define it as ValueType
		typedef typename TAAPosVRT::ValueType ValueType;

		number minHeight, tmpMinHeight;
		ValueType v = aaPos[face->vertex(2)];
		ValueType dir;

	//	Calculate start height and set to minHeight
		VecSubtract(dir, aaPos[face->vertex(1)], aaPos[face->vertex(0)]);
		minHeight = DistancePointToRay(	v, aaPos[face->vertex(0)], dir);

		for(uint i = 1; i < 3; ++i)
		{
			v = aaPos[face->vertex((i+2)%3)];
			VecSubtract(dir, aaPos[face->vertex((i+1)%3)], aaPos[face->vertex((i))]);
			tmpMinHeight = DistancePointToRay(v, aaPos[face->vertex(i )], dir);

			if(tmpMinHeight < minHeight)
			{
				minHeight = tmpMinHeight;
			}
		}

		return minHeight;
	}
	else
		UG_ASSERT(false, "Error. Face is not a triangle.");

	return NAN;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateAspectRatio

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Face* face, TAAPosVRT& aaPos)
{
	number AspectRatio;
	number maxEdgeLength;

//	Collect element edges, find longest edge and calculate its length
	vector<EdgeBase*> edges;
	CollectAssociated(edges, grid, face);
	EdgeBase* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
	maxEdgeLength = EdgeLength(longestEdge, aaPos);

	switch (face->reference_object_id())
	{
		case ROID_TRIANGLE:
		{
			/* MINHEIGHT / MAXEDGELENGTH
			 * optimal Aspect Ratio of a regular triangle
			 * Q = sqrt(3)/2 * a / a = 0.86602540378444...
			 */

		//	Calculate minimal triangle height
			number minTriangleHeight = CalculateMinTriangleHeight(face, aaPos);

		//	Calculate the aspect ratio
			AspectRatio = minTriangleHeight / maxEdgeLength;

			return AspectRatio;
		}

		case ROID_QUADRILATERAL:
		{
			/* AREA / MAXEDGELENGTH */
		//	Calculate the element area
			number area = FaceArea(face, aaPos);

		//	Calculate the aspect ratio
			AspectRatio = area / maxEdgeLength;

			return AspectRatio;
		}

		default:
		 	UG_ASSERT(false, "Error. Unknown element type for aspect ratio calculation.");
	}

	return NAN;
}


//	Tetrahedron
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	number AspectRatio;
	number maxEdgeLength;

//	Collect element edges, find longest edge and calculate its length
	vector<EdgeBase*> edges;
	CollectAssociated(edges, grid, tet);
	EdgeBase* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
	maxEdgeLength = EdgeLength(longestEdge, aaPos);

	/* MINHEIGHT / MAXEDGELENGTH
	 * optimal Aspect Ratio of a regular tetrahedron
	 * Q = sqrt(2/3) * a / a = 0.81...
	 */

//	Calculate the aspect ratio
	AspectRatio = CalculateTetrahedronAspectRatio(grid, tet, aaPos);

	return AspectRatio;
}

//	Prism
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	number AspectRatio;
	number maxEdgeLength;

//	Collect element edges, find longest edge and calculate its length
	vector<EdgeBase*> edges;
	CollectAssociated(edges, grid, prism);
	EdgeBase* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
	maxEdgeLength = EdgeLength(longestEdge, aaPos);

	/* VOLUME / MAXEDGELENGTH */
//	Calculate the element volume
	number volume = CalculateVolume(*prism, aaPos);

//	Calculate the aspect ratio
	AspectRatio = volume / maxEdgeLength;

	return AspectRatio;
}

//	Pyramid
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	number AspectRatio;
	number maxEdgeLength;

//	Collect element edges, find longest edge and calculate its length
	vector<EdgeBase*> edges;
	CollectAssociated(edges, grid, pyr);
	EdgeBase* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
	maxEdgeLength = EdgeLength(longestEdge, aaPos);

	/* VOLUME / MAXEDGELENGTH */
//	Calculate the element volume
	number volume = CalculateVolume(*pyr, aaPos);

//	Calculate the aspect ratio
	AspectRatio = volume / maxEdgeLength;

	return AspectRatio;
}

//	Volume
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Volume* vol, TAAPosVRT& aaPos)
{
	switch (vol->reference_object_id())
	{
		case ROID_TETRAHEDRON:
		{
			return CalculateAspectRatio(grid, static_cast<Tetrahedron*>(vol), aaPos);
		}

		case ROID_PRISM:
		{
			return CalculateAspectRatio(grid, static_cast<Prism*>(vol), aaPos);
		}

		case ROID_PYRAMID:
		{
			return CalculateAspectRatio(grid, static_cast<Pyramid*>(vol), aaPos);
		}

		default:
		 	UG_ASSERT(false, "Error. Unknown element type for aspect ratio calculation.");
	}

	return NAN;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithSmallestMinAngle
template <class TElem, class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithSmallestMinAngle(Grid& grid, TIterator elementsBegin, TIterator elementsEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	typename TIterator::value_type elementWithSmallestMinAngle = *elementsBegin;
	number smallestMinAngle = CalculateMinAngle(grid, elementWithSmallestMinAngle, aaPos);
	++elementsBegin;

//	compare all volumes and find that one with smallest minAngle
	for(; elementsBegin != elementsEnd; ++elementsBegin)
	{
		typename TIterator::value_type curElement = *elementsBegin;
		number curSmallestMinAngle = CalculateMinAngle(grid, curElement, aaPos);

		if(curSmallestMinAngle < smallestMinAngle)
		{
			elementWithSmallestMinAngle = curElement;
			smallestMinAngle = curSmallestMinAngle;
		}
	}

	return elementWithSmallestMinAngle;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindLargestFace
template <class TIterator, class TAAPosVRT>
Face* FindLargestFace(TIterator facesBegin, TIterator facesEnd, TAAPosVRT& aaPos)
{
	//	if facesBegin equals facesEnd, then the list is empty and we can
	//	immediately return NULL
		if(facesBegin == facesEnd)
			return NULL;

	//	the first face is the first candidate for the smallest face.
		Face* largestFace = *facesBegin;
		number largestArea = FaceArea(largestFace, aaPos);
		++facesBegin;

		for(; facesBegin != facesEnd; ++facesBegin){
			Face* curFace = *facesBegin;
			number curArea = FaceArea(curFace, aaPos);
			if(curArea > largestArea){
				largestFace = curFace;
				largestArea = curArea;
			}
		}

		return largestFace;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindSmallestVolumeElement
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindSmallestVolumeElement(TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	typename TIterator::value_type smallestVolume = *volumesBegin;
	number smallestVolumeVolume = CalculateVolume(*smallestVolume, aaPos);
	++volumesBegin;

//	compare all tetrahedrons and find minimal volume
	for(; volumesBegin != volumesEnd; ++volumesBegin)
	{
		Volume* curVolume = *volumesBegin;
		number curVolumeVolume = CalculateVolume(*curVolume, aaPos);

		if(curVolumeVolume < smallestVolumeVolume)
		{
			smallestVolume = curVolume;
			smallestVolumeVolume = curVolumeVolume;
		}
	}

	return smallestVolume;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindLargestVolumeElement
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindLargestVolumeElement(TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	typename TIterator::value_type largestVolume = *volumesBegin;
	number largestVolumeVolume = CalculateVolume(*largestVolume, aaPos);
	++volumesBegin;

//	compare all tetrahedrons and find minimal volume
	for(; volumesBegin != volumesEnd; ++volumesBegin)
	{
		Volume* curVolume = *volumesBegin;
		number curVolumeVolume = CalculateVolume(*curVolume, aaPos);

		if(curVolumeVolume > largestVolumeVolume)
		{
			largestVolume = curVolume;
			largestVolumeVolume = curVolumeVolume;
		}
	}

	return largestVolume;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithSmallestAspectRatio
template <class TElem, class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithSmallestAspectRatio(Grid& grid, 	TIterator elemsBegin,
												TIterator elemsEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	typename TIterator::value_type elementWithSmallestAspectRatio = *elemsBegin;
	number smallestAspectRatio = CalculateAspectRatio(grid, elementWithSmallestAspectRatio, aaPos);
	++elemsBegin;

//	compare all tetrahedrons and find that one with minimal aspect ratio
	for(; elemsBegin != elemsEnd; ++elemsBegin)
	{
		typename TIterator::value_type curElement = *elemsBegin;
		//TElem* curElement = *elemsBegin;
		number curSmallestAspectRatio = CalculateAspectRatio(grid, curElement, aaPos);

		if(curSmallestAspectRatio < smallestAspectRatio)
		{
			elementWithSmallestAspectRatio = curElement;
			smallestAspectRatio = curSmallestAspectRatio;
		}
	}

	return elementWithSmallestAspectRatio;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithLargestAspectRatio
template <class TElem, class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithLargestAspectRatio(Grid& grid,  	TIterator elemsBegin,
												TIterator elemsEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	typename TIterator::value_type elementWithLargestAspectRatio = *elemsBegin;
	number largestAspectRatio = CalculateAspectRatio(grid, elementWithLargestAspectRatio, aaPos);
	++elemsBegin;

//	compare all tetrahedrons and find that one with maximal aspect ratio
	for(; elemsBegin != elemsEnd; ++elemsBegin)
	{
		typename TIterator::value_type curElement = *elemsBegin;
		//TElem* curElement = *elemsBegin;
		number curSmallestAspectRatio = CalculateAspectRatio(grid, curElement, aaPos);

		if(curSmallestAspectRatio > largestAspectRatio)
		{
			elementWithLargestAspectRatio = curElement;
			largestAspectRatio = curSmallestAspectRatio;
		}
	}

	return elementWithLargestAspectRatio;
}






////////////////////////////////////////////////////////////////////////////////////////////
//	MinAngleHistogram
template <class TIterator, class TAAPosVRT>
void MinAngleHistogram(Grid& grid, 	TIterator elementsBegin,
									TIterator elementsEnd,
									TAAPosVRT& aaPos,
									uint stepSize)
{
//	Initialization
	vector<number> minAngles;
	uint numRanges = floor(180 / stepSize);
	vector<uint> counter(numRanges, 0);
	typename TIterator::value_type refElem = *elementsBegin;

//	Calculate the minAngle of every element
	for(TIterator iter = elementsBegin; iter != elementsEnd; ++iter)
	{
		number curMinAngle = CalculateMinAngle(grid, *iter, aaPos);
		minAngles.push_back(curMinAngle);
	}

//	Sort the calculated minAngles in an ascending way
	sort (minAngles.begin(), minAngles.end());

//	Count the elements in their corresponding minAngle range (180/stepSize ranges)
	for(uint i = 0; i < minAngles.size(); ++i)
	{
		number minAngle = minAngles[i];
		#ifdef __GXX_EXPERIMENTAL_CXX0X__
			for (uint range = 0; range < numRanges; range++)
			{
				if (isInRange(range*stepSize, (range+1)*stepSize)(minAngle))
				{
					++counter[range];
					break;
				}
			}
		#else
			for (uint range = 0; range < numRanges; range++)
			{
				if (minAngle < (range+1)*stepSize)
				{
					++counter[range];
					break;
				}
			}
		#endif
	}

/*
//	-------------------------------
//	Histogram table output section:

//	Divide table into two halfs
	uint numRows = ceil(number(numRanges) / 2.0);
	ug::Table<std::stringstream> minAngleTable(numRows, 4);

//	Specific element header
	UG_LOG(endl << "MinAngle-Histogram for '" << refElem->reference_object_id() << "' elements");
	UG_LOG(endl);

//	First half
	uint i = 0;
	for(; i < numRows; ++i)
	{
		minAngleTable(i, 0) << i*stepSize << " - " << (i+1)*stepSize << " degrees : ";
		minAngleTable(i, 1) << counter[i];
	}

//	Second half
	for(; i < numRanges; ++i)
	{
		if(i < numRanges-1)
		{
			minAngleTable(i-numRows, 2) << i*stepSize << " - " << (i+1)*stepSize << " degrees : ";
			minAngleTable(i-numRows, 3) << counter[i];
		}
	//	Last entry needs special treatment
		else
		{
			minAngleTable(i-numRows, 2) << i*stepSize << " - " << 180 << " degrees : ";
			minAngleTable(i-numRows, 3) << counter[i];
		}
	}
*/

//	-------------------------------
//	Histogram table output section:

//	Divide table into three thirds
	uint numRows = ceil(number(numRanges) / 3.0);
	ug::Table<std::stringstream> minAngleTable(numRows, 6);

//	Specific element header
	UG_LOG(endl << "MinAngle-Histogram for '" << refElem->reference_object_id() << "' elements");
	UG_LOG(endl);

//	First third
	uint i = 0;
	for(; i < numRows; ++i)
	{
		minAngleTable(i, 0) << i*stepSize << " - " << (i+1)*stepSize << " degrees : ";
		minAngleTable(i, 1) << counter[i];
	}

//	Second third
	for(; i < 2*numRows; ++i)
	{
		minAngleTable(i-numRows, 2) << i*stepSize << " - " << (i+1)*stepSize << " degrees : ";
		minAngleTable(i-numRows, 3) << counter[i];
	}

//	Third third
	for(; i < numRanges; ++i)
	{
		if(i < numRanges-1)
		{
			minAngleTable(i-2*numRows, 4) << i*stepSize << " - " << (i+1)*stepSize << " degrees : ";
			minAngleTable(i-2*numRows, 5) << counter[i];
		}
	//	Last entry needs special treatment
		else
		{
			minAngleTable(i-2*numRows, 4) << i*stepSize << " - " << 180 << " degrees : ";
			minAngleTable(i-2*numRows, 5) << counter[i];
		}
	}

//	Output table
	UG_LOG(endl << minAngleTable);


}


////////////////////////////////////////////////////////////////////////////////////////////
//	ElementQualityStatistics
////////////////////////////////////////////////////////////////////////////////////////////

//	Wrapper
void ElementQualityStatistics(MultiGrid& mg, int level)
{
	ElementQualityStatistics(mg, mg.get_geometric_objects(level));
}

void ElementQualityStatistics(Grid& grid)
{
	ElementQualityStatistics(grid, grid.get_geometric_objects());
}

//	Actual procedure
void ElementQualityStatistics(Grid& grid, GeometricObjectCollection goc)
{
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPosition);

//	Numbers
	number smallestVolumeVolume;
	number largestVolumeVolume;
	number minFaceAngle;
	number minVolumeAngle;
	number shortestEdgeLength;
	number longestEdgeLength;
	number minTetrahedronAspectRatio;
	number maxTetrahedronAspectRatio;

//	Elements
	Volume* smallestVolume;
	Volume* largestVolume;
	Volume* volumeWithSmallestMinAngle;
	Tetrahedron* tetrahedronWithSmallestAspectRatio;
	Tetrahedron* tetrahedronWithLargestAspectRatio;
	EdgeBase* shortestEdge;
	EdgeBase* longestEdge;
	Face* faceWithSmallestMinAngle;


//	Basic grid properties

//	Check dimension
	if(grid.num_volumes() > 0)
	{
	//	Volumes
		smallestVolume = FindSmallestVolumeElement(	grid.volumes_begin(),
												grid.volumes_end(),
												aaPos);

		largestVolume = FindLargestVolumeElement(	grid.volumes_begin(),
											grid.volumes_end(),
											aaPos);

		volumeWithSmallestMinAngle = FindElementWithSmallestMinAngle<Volume>(grid,
																	grid.volumes_begin(),
																	grid.volumes_end(),
																	aaPos);
	//	Numbers
		smallestVolumeVolume = CalculateVolume(*smallestVolume,	aaPos);
		largestVolumeVolume = CalculateVolume(*largestVolume, aaPos);
		minVolumeAngle = CalculateMinAngle(grid, volumeWithSmallestMinAngle, aaPos);
	}

//	Tetrahedron section
	if(grid.num<Tetrahedron>() > 0)
	{
	//	Tetrahedrons
		tetrahedronWithSmallestAspectRatio = FindElementWithSmallestAspectRatio<Tetrahedron>(	grid,
																								grid.begin<Tetrahedron>(),
																								grid.end<Tetrahedron>(),
																								aaPos);

		tetrahedronWithLargestAspectRatio = FindElementWithLargestAspectRatio<Tetrahedron>(	grid,
																							grid.begin<Tetrahedron>(),
																							grid.end<Tetrahedron>(),
																							aaPos);

	//	Numbers
		minTetrahedronAspectRatio = CalculateAspectRatio(grid, tetrahedronWithSmallestAspectRatio, aaPos);
		maxTetrahedronAspectRatio = CalculateAspectRatio(grid, tetrahedronWithLargestAspectRatio, aaPos);
	}


//	Elements
	shortestEdge = FindShortestEdge(grid.edges_begin(),
									grid.edges_end(),
									aaPos);

	longestEdge = FindLongestEdge(	grid.edges_begin(),
									grid.edges_end(),
									aaPos);

	faceWithSmallestMinAngle = FindElementWithSmallestMinAngle<Face>(grid,
															grid.begin<Face>(),
															grid.end<Face>(),
															aaPos);

//	Numbers
	shortestEdgeLength = EdgeLength(shortestEdge, aaPos);
	longestEdgeLength = EdgeLength(longestEdge, aaPos);
	minFaceAngle = CalculateMinAngle(grid, faceWithSmallestMinAngle, aaPos);


//	Table summary
	ug::Table<std::stringstream> table(10, 2);
	table(0, 0) << "Number of volumes"; 	table(0, 1) << grid.num_volumes();
	table(1, 0) << "Number of faces"; 		table(1, 1) << grid.num_faces();
	table(2, 0) << "Number of vertices";	table(2, 1) << grid.num_vertices();

	table(3, 0) << "Shortest edge length";	table(3, 1) << shortestEdgeLength;
	table(4, 0) << "Longest edge length";	table(4, 1) << longestEdgeLength;

	table(5, 0) << "Smallest face angle";	table(5, 1) << minFaceAngle;

	if(grid.num_volumes() > 0)
	{
		table(6, 0) << "Smallest volume";	table(6, 1) << smallestVolumeVolume;
		table(7, 0) << "Largest volume";	table(7, 1) << largestVolumeVolume;
	}

	if(grid.num<Tetrahedron>() > 0)
	{
		table(8, 0) << "Smallest tetrahedron AR";	table(8, 1) << minTetrahedronAspectRatio;
		table(9, 0) << "Largest tetrahedron AR";	table(9, 1) << maxTetrahedronAspectRatio;
	}


//	Output section
	UG_LOG(endl << "--------------------------------------------------------------------------" << endl);
	UG_LOG("Grid quality statistics:" << endl << endl);
	UG_LOG(table);

	MinAngleHistogram(grid, grid.begin<Face>(), grid.end<Face>(), aaPos, 10);

	if(grid.num_volumes() > 0)
	{
		MinAngleHistogram(grid, grid.volumes_begin(), grid.volumes_end(), aaPos, 10);
	}
	UG_LOG("--------------------------------------------------------------------------" << endl << endl);



}




















}




