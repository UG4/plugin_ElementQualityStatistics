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
//	CalculateMinFaceAngle
template <class TAAPosVRT>
number CalculateMinFaceAngle(Grid& grid, Face* f, TAAPosVRT& aaPos)
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


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinVolumeAngle
number CalculateMinVolumeAngle(	Grid& grid, Volume* v,
								Grid::VertexAttachmentAccessor<APosition>& aaPos)
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
number CalculateMinTriangleHeight(Triangle* tri, TAAPosVRT& aaPos)
{
//	Get type of vertex attachment in aaPos and define it as ValueType
	typedef typename TAAPosVRT::ValueType ValueType;

	number minHeight, tmpMinHeight;
	ValueType v = aaPos[tri->vertex(2)];

//	Calculate start height and set to minHeight
	minHeight = DistancePointToRay(	v,
									aaPos[tri->vertex(0)],
									aaPos[tri->vertex(1)] - aaPos[tri->vertex(0)]);

	for(uint i = 1; i < 3; ++i)
	{
		v = aaPos[tri->vertex((i+2)%3)];
		tmpMinHeight = DistancePointToRay(	v,
											aaPos[tri->vertex(i )],
											aaPos[tri->vertex((i+1)%3)] - aaPos[tri->vertex((i))]);

		if(tmpMinHeight < minHeight)
		{
			minHeight = tmpMinHeight;
		}
	}

	return minHeight;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateAspectRatio

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, TElem*, TAAPosVRT& aaPos);

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
			number minTriangleHeight = CalculateMinTriangleHeight(static_cast<Triangle*>(face), aaPos);

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
	number minTetrahedronHeight;

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
	number minTetrahedronHeight;

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
//	FindFaceWithSmallestMinAngle
template <class TIterator, class TAAPosVRT>
Face* FindFaceWithSmallestMinAngle(Grid& grid, TIterator facesBegin, TIterator facesEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	Face* faceWithSmallestMinAngle = *facesBegin;
	number smallestMinAngle = CalculateMinFaceAngle(grid, faceWithSmallestMinAngle, aaPos);
	++facesBegin;

//	compare all volumes and find that one with smallest minAngle
	for(; facesBegin != facesEnd; ++facesBegin)
	{
		Face* curFace = *facesBegin;
		number curSmallestMinAngle = CalculateMinFaceAngle(grid, curFace, aaPos);

		if(curSmallestMinAngle < smallestMinAngle)
		{
			faceWithSmallestMinAngle = curFace;
			smallestMinAngle = curSmallestMinAngle;
		}
	}

	return faceWithSmallestMinAngle;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindVolumeWithSmallestMinAngle
template <class TIterator, class TAAPosVRT>
Volume* FindVolumeWithSmallestMinAngle(Grid& grid, TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	Volume* volumeWithSmallestMinAngle = *volumesBegin;
	number smallestMinAngle = CalculateMinVolumeAngle(grid, volumeWithSmallestMinAngle, aaPos);
	++volumesBegin;

//	compare all volumes and find that one with smallest minAngle
	for(; volumesBegin != volumesEnd; ++volumesBegin)
	{
		Volume* curVolume = *volumesBegin;
		number curSmallestMinAngle = CalculateMinVolumeAngle(grid, curVolume, aaPos);

		if(curSmallestMinAngle < smallestMinAngle)
		{
			volumeWithSmallestMinAngle = curVolume;
			smallestMinAngle = curSmallestMinAngle;
		}
	}

	return volumeWithSmallestMinAngle;
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
//	FindSmallestVolume
template <class TIterator, class TAAPosVRT>
Volume* FindSmallestVolume(TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	Volume* smallestVolume = *volumesBegin;
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
//	FindLargestVolume
template <class TIterator, class TAAPosVRT>
Volume* FindLargestVolume(TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	Volume* largestVolume = *volumesBegin;
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
template <class TIterator, class TAAPosVRT>
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
//	MinFaceAngleHistogram
void MinFaceAngleHistogram(Grid& grid)
{
//	Quality histograms
	//ANumber aMinFaceAngle;
	//Grid::FaceAttachmentAccessor<ANumber> aaFaceMinAngle(grid, aFaceMinAngle);

//	Initialization
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPosition);
	vector<number> MinFaceAngles;
	uint counter[18] = {0};

//	Calculate the MinAngle of every face
	for(FaceIterator fIter = grid.faces_begin(); fIter != grid.faces_end(); ++fIter)
	{
		number curFaceMinAngle = CalculateMinFaceAngle(grid, *fIter, aaPos);
		MinFaceAngles.push_back(curFaceMinAngle);
	}

//	Sort the calculated MinAngles in an ascending way
	sort (MinFaceAngles.begin(), MinFaceAngles.end());

//	Count the faces in their corresponding MinAngle range (18 ranges)
	for(uint i = 0; i < MinFaceAngles.size(); ++i)
	{
		number MinFaceAngle = MinFaceAngles[i];
				#ifdef __GXX_EXPERIMENTAL_CXX0X__
					for (uint range = 0; range < 18; range++) {
						if (isInRange(range*10, (range+1)*10)(MinFaceAngle)) {
							++counter[range];
							break;
						}
					}
				#else
					for (uint range = 0; range < 18; range++) {
						if (MinFaceAngle < (range+1)*10) {
							++counter[range];
							break;
						}
					}
				#endif
	}

//	MinFaceAngle-Histogram output table
	UG_LOG(endl << "MinFaceAngle-Histogram\n" << endl);
	UG_LOG("   " << 	" 0"  << " - " << "10"  << " degrees   :   "  << counter[0] << "      |      " <<
						" 90" << " - " << "100" << " degrees   :   "  << counter[9] << endl);
	for(uint i = 1; i < 9; ++i)
		{
			UG_LOG("   " << 	i*10 		<< " - " << (i+1)*10 	<< " degrees   :   "  << counter[i] << "      |      " <<
								(i+9)*10 	<< " - " << (i+9+1)*10 	<< " degrees   :   "  << counter[i+9] << endl);
		}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	MinFaceAngleHistogram
void MinVolumeAngleHistogram(Grid& grid)
{
//	Quality histograms
	//ANumber aMinFaceAngle;
	//Grid::FaceAttachmentAccessor<ANumber> aaFaceMinAngle(grid, aFaceMinAngle);

//	Initialization
	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPosition);
	vector<number> MinVolumeAngles;
	uint counter[18] = {0};

//	Calculate the MinAngle of every volume
	for(VolumeIterator vIter = grid.volumes_begin(); vIter != grid.volumes_end(); ++vIter)
	{
		number curVolumeMinAngle = CalculateMinVolumeAngle(grid, *vIter, aaPos);
		MinVolumeAngles.push_back(curVolumeMinAngle);
	}

//	Sort the calculated MinAngles in an ascending way
	sort (MinVolumeAngles.begin(), MinVolumeAngles.end());

//	Count the volumes in their corresponding MinAngle range (18 ranges)
	for(uint i = 0; i < MinVolumeAngles.size(); ++i)
	{
		number MinVolumeAngle = MinVolumeAngles[i];

		#ifdef __GXX_EXPERIMENTAL_CXX0X__
			for (uint range = 0; range < 18; range++)
			{
				if (isInRange(range*10, (range+1)*10)(MinVolumeAngle))
				{
					++counter[range];
					break;
				}
			}
		#else
			for (uint range = 0; range < 18; range++)
			{
				if (MinVolumeAngle < (range+1)*10)
				{
					++counter[range];
					break;
				}
			}
		#endif

	}

//	MinVolumeAngle-Histogram output table
	UG_LOG(endl << "MinVolumeAngle-Histogram\n" << endl);
	UG_LOG("   " << 	" 0"  << " - " << "10"  << " degrees   :   "  << counter[0] << "      |      " <<
						" 90" << " - " << "100" << " degrees   :   "  << counter[9] << endl);
	for(uint i = 1; i < 9; ++i)
		{
			UG_LOG("   " << 	i*10 		<< " - " << (i+1)*10 	<< " degrees   :   "  << counter[i] << "      |      " <<
								(i+9)*10 	<< " - " << (i+9+1)*10 	<< " degrees   :   "  << counter[i+9] << endl);
		}
}






////////////////////////////////////////////////////////////////////////////////////////////
//	element_quality_statistics
void ElementQualityStatistics(Grid& grid)
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
	Tetrahedron* tetrahedronWithLargestAspectRatio = *grid.begin<Tetrahedron>();
	Tetrahedron* tetrahedronWithSmallestAspectRatio = *grid.begin<Tetrahedron>();
	EdgeBase* shortestEdge;
	EdgeBase* longestEdge;
	Face* faceWithSmallestMinAngle;


//	Basic grid properties

//	Check dimension
	if(grid.num_volumes() > 0)
	{
	//	Volumes
		smallestVolume = FindSmallestVolume(	grid.volumes_begin(),
												grid.volumes_end(),
												aaPos);

		largestVolume = FindLargestVolume(	grid.volumes_begin(),
											grid.volumes_end(),
											aaPos);

		volumeWithSmallestMinAngle = FindVolumeWithSmallestMinAngle(grid,
																	grid.volumes_begin(),
																	grid.volumes_end(),
																	aaPos);
	//	Numbers
		smallestVolumeVolume = CalculateVolume(*smallestVolume,	aaPos);
		largestVolumeVolume = CalculateVolume(*largestVolume, aaPos);
		minVolumeAngle = CalculateMinVolumeAngle(grid, volumeWithSmallestMinAngle, aaPos);
	}

//	Tetrahedron section
	if(grid.num<Tetrahedron>() > 0)
	{
	//	Tetrahedrons
		tetrahedronWithSmallestAspectRatio = FindElementWithSmallestAspectRatio(	grid,
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

	faceWithSmallestMinAngle = FindFaceWithSmallestMinAngle(grid,
															grid.faces_begin(),
															grid.faces_end(),
															aaPos);

//	Numbers
	shortestEdgeLength = EdgeLength(shortestEdge, aaPos);
	longestEdgeLength = EdgeLength(longestEdge, aaPos);
	minFaceAngle = CalculateMinFaceAngle(grid, faceWithSmallestMinAngle, aaPos);





//	Output section
	UG_LOG(endl << "--------------------------------------------------------------------------" << endl);
	UG_LOG("Grid quality statistics:" << endl << endl);

	UG_LOG("   " << "Number of volumes     = 	" << grid.num_volumes() << endl);
	UG_LOG("   " << "Number of faces       = 	" << grid.num_faces() << endl);
	UG_LOG("   " << "Number of vertices    = 	" << grid.num_vertices() << endl << endl);

	UG_LOG("   " << "Shortest edge length  = 	" << shortestEdgeLength << "   |   " <<
					"Longest edge length   = 	" << longestEdgeLength 	<< endl);

	if(grid.num_volumes() > 0)
	{
	UG_LOG("   " << "Smallest volume       = 	" << smallestVolumeVolume << "   |   " <<
					"Largest volume        = 	" << largestVolumeVolume  << endl);
	UG_LOG("   " << "Smallest volume angle = 	" << minVolumeAngle << endl);
	}

	UG_LOG("   " << "Smallest face angle   = 	" << minFaceAngle << endl);

	if(grid.num<Tetrahedron>() > 0)
	{
	UG_LOG(endl);
	UG_LOG("   " 	<< "Smallest tetrahedron" << endl <<
		   "   "	<< " aspect ratio         =      " << minTetrahedronAspectRatio << endl);
	UG_LOG("   "    << "Largest tetrahedron" << endl <<
		   "   "	<< " aspect ratio         =      " << maxTetrahedronAspectRatio  << endl);
	}

	MinFaceAngleHistogram(grid);

	if(grid.num_volumes() > 0)
	{
	MinVolumeAngleHistogram(grid);
	}

	UG_LOG("--------------------------------------------------------------------------" << endl << endl);

	ug::Table<std::stringstream> table(3, 2);
	table(0, 0) << "Number of volumes"; table(0, 1) << grid.num_volumes();
	table(1, 0) << "Number of faces"; table(1, 1) << grid.num_faces();
	table(2, 0) << "Number of vertices"; table(2, 1) << grid.num_vertices();

	std::cout << table;
}




















}




