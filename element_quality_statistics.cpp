/*
 * element_quality_statistics.cpp
 *
 *  Created on: 17.04.2012
 *      Author: Martin Stepniewski
 */


#include "element_quality_statistics.h"

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
			UG_ASSERT(sidesOut[1] == NULL, "Only two edges may be adjacent to a vertex in a face element.")
			if(sidesOut[0] == NULL)
				sidesOut[0] = vNeighbourEdgesToVertex[i];
			else
				sidesOut[1] = vNeighbourEdgesToVertex[i];
			}
		}

	UG_ASSERT(sidesOut[1] != NULL, "Exactly two edges should be adjacent to a vertex in a face element.")
}

///	Collects all faces (= 2) which exist in the given volume and which share the given vertex.
/*
UG_API
inline void CollectAssociatedSides(Face* sidesOut[2], Grid& grid, Volume* v, VertexBase* vrt)
{
	vector<Face*> vNeighbourFacesToVertex;
	sidesOut[0] = NULL;
	sidesOut[1] = NULL;

	CollectFaces(vNeighbourFacesToVertex, grid, vrt, true);
	for(uint i = 0; i < vNeighbourFacesToVertex.size(); ++i)
	{
		if(VolumeContains(v, vNeighbourFacesToVertex[i]) == true)
		{
			UG_ASSERT(sidesOut[1] == NULL, "Only two faces may be adjacent to vertex in a volume element.")
			if(sidesOut[0] == NULL)
				sidesOut[0] = vNeighbourFacesToVertex[i];
			else
				sidesOut[1] = vNeighbourFacesToVertex[i];
			}
		}

	UG_ASSERT(sidesOut[1] != NULL, "Exactly two faces should be adjacent to a vertex in a volume element.")
}
*/

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
			UG_ASSERT(sidesOut[1] == NULL, "Only two faces may be adjacent to an edge in a volume element.")
			if(sidesOut[0] == NULL)
				sidesOut[0] = vNeighbourFacesToEdge[i];
			else
				sidesOut[1] = vNeighbourFacesToEdge[i];
			}
		}

	UG_ASSERT(sidesOut[1] != NULL, "Exactly two faces should be adjacent to an edge in a volume element.")
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinFaceAngle
number CalculateMinFaceAngle(Grid& grid, Face* f, Grid::VertexAttachmentAccessor<APosition> aaPos)
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
	uint numElementVrts = f->num_vertices();
	vector3 vNorm1, vNorm2;
	vector3 vDir1, vDir2;
	number minAngle = 180.0;
	number tmpAngle;
	EdgeBase* vNeighbourEdgesToVertex[2];

//	Iterate over all face vertices
	for(uint vrtIter = 0; vrtIter < numElementVrts; ++vrtIter)
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
//	CalculateMinVolumeAngle
number CalculateMinVolumeAngle(Grid& grid, Volume* v, Grid::VertexAttachmentAccessor<APosition> aaPos)
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
//	CalculateVolumeMinHeight
number CalculateMinVolumeHeight(const Volume& vol,
								Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	switch (vol.reference_object_id())
	{
		case ROID_TETRAHEDRON:
			return CalculateMinVolumeHeight(static_cast<Tetrahedron>(vol), aaPos);
		/*
		case ROID_PRISM:
			return CalculateMinVolumeHeight(static_cast<Prism>(vol), aaPos);
		case ROID_PYRAMID:
			return CalculateMinVolumeHeight(static_cast<Pyramid>(vol), aaPos);
		case ROID_HEXAHEDRON:
		*/
			return CalculateMinVolumeHeight(static_cast<Hexahedron>(vol), aaPos);
		default:
			UG_ASSERT(false, "dont know how to calculate height for given volume.");
	}

	return NAN;
}

number CalculateMinVolumeHeight(const Tetrahedron& tet,
								Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	vector3& a = aaPos[tet.vertex(0)];
	vector3& b = aaPos[tet.vertex(1)];
	vector3& c = aaPos[tet.vertex(2)];
	vector3& d = aaPos[tet.vertex(3)];

	return CalculateMinTetrahedronHeight(a, b, c, d);
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
//	FindTetrahedronWithSmallestAspectRatio
template <class TIterator, class TAAPosVRT>
Volume* FindTetrahedronWithSmallestAspectRatio(Grid& grid, TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	Volume* volumeWithSmallestAspectRatio = *volumesBegin;
	number smallestAspectRatio = CalculateVolume(*volumeWithSmallestAspectRatio, aaPos);
	++volumesBegin;

//	compare all tetrahedrons and find that one with minimal aspect ratio
	for(; volumesBegin != volumesEnd; ++volumesBegin)
	{
		Volume* curVolume = *volumesBegin;
		number curSmallestAspectRatio = CalculateTetrahedronAspectRatio(grid, *curVolume);

		if(curSmallestAspectRatio < smallestAspectRatio)
		{
			volumeWithSmallestAspectRatio = curVolume;
			smallestAspectRatio = curSmallestAspectRatio;
		}
	}

	return volumeWithSmallestAspectRatio;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindTetrahedronWithLargestAspectRatio
template <class TIterator, class TAAPosVRT>
Volume* FindTetrahedronWithLargestAspectRatio(Grid& grid, TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	//	if(volumesBegin == volumesBegin)
	//		return NULL;

//	Initializations
	Volume* volumeWithLargestAspectRatio = *volumesBegin;
	number largestAspectRatio = CalculateVolume(*volumeWithLargestAspectRatio, aaPos);
	++volumesBegin;

//	compare all tetrahedrons and find that one with maximal aspect ratio
	for(; volumesBegin != volumesEnd; ++volumesBegin)
	{
		Volume* curVolume = *volumesBegin;
		number curSmallestAspectRatio = CalculateTetrahedronAspectRatio(grid, *curVolume);

		if(curSmallestAspectRatio > largestAspectRatio)
		{
			volumeWithLargestAspectRatio = curVolume;
			largestAspectRatio = curSmallestAspectRatio;
		}
	}

	return volumeWithLargestAspectRatio;
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
void element_quality_statistics(Grid& grid)
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
	Volume* tetrahedronWithLargestAspectRatio;
	Volume* tetrahedronWithSmallestAspectRatio;
	EdgeBase* shortestEdge;
	EdgeBase* longestEdge;
	Face* faceWithSmallestMinAngle;

	bool bTetrahedronCheck;



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

	//	Check, if grid contains tetrahedrons
		for(VolumeIterator vIter = grid.volumes_begin(); vIter != grid.volumes_end(); ++vIter)
		{
			Volume* v = *vIter;
			if(v->reference_object_id() == ROID_TETRAHEDRON)
			{
				bTetrahedronCheck = true;
				break;
			}
		}
	}

//	Tetrahedron section
	if(bTetrahedronCheck)
		{
		//	Tetrahedrons
			tetrahedronWithSmallestAspectRatio = FindTetrahedronWithSmallestAspectRatio(grid,
																						grid.volumes_begin(),
																						grid.volumes_end(),
																						aaPos);

			tetrahedronWithLargestAspectRatio = FindTetrahedronWithLargestAspectRatio(	grid,
																						grid.volumes_begin(),
																						grid.volumes_end(),
																						aaPos);
		//	Numbers
			minTetrahedronAspectRatio = CalculateTetrahedronAspectRatio(grid, *tetrahedronWithSmallestAspectRatio);
			maxTetrahedronAspectRatio = CalculateTetrahedronAspectRatio(grid, *tetrahedronWithLargestAspectRatio);
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

	if(bTetrahedronCheck)
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

}




















}




