/*
 * elem_stat_util.h
 *
 *  Created on: 17.04.2012
 *      Author: Martin Stepniewski
 */


#ifndef __ELEM_STAT_UTIL_h__
#define __ELEM_STAT_UTIL_h__

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/distributed_grid.h"
#endif

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
/**	This method uses Grid::mark **/
UG_API
inline void CollectAssociatedSides(EdgeBase* sidesOut[2], Grid& grid, Face* f, Vertex* vrt)
{
	//PROFILE_BEGIN(CollectAssociatedSides_VERTEX);
	sidesOut[0] = NULL;
	sidesOut[1] = NULL;

	grid.begin_marking();
	for(size_t i = 0; i < f->num_vertices(); ++i){
		grid.mark(f->vertex(i));
	}

	//vector<EdgeBase*> vNeighbourEdgesToVertex;
	//CollectAssociated(vNeighbourEdgesToVertex, grid, vrt, true);

	Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(vrt);
	//Grid::AssociatedEdgeIterator iterEnd = vNeighbourEdgesToVertex.end();
	for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(vrt); iter != iterEnd; ++iter)
	//for(Grid::AssociatedEdgeIterator iter = vNeighbourEdgesToVertex.begin(); iter != iterEnd; ++iter)
	{
		EdgeBase* e = *iter;
		if(grid.is_marked(e->vertex(0)) && grid.is_marked(e->vertex(1))){
			UG_ASSERT(	sidesOut[1] == NULL,
						"Only two edges may be adjacent to a vertex in a face element.");

			if(sidesOut[0] == NULL)
				sidesOut[0] = e;
			else
				sidesOut[1] = e;
		}
	}

	grid.end_marking();
	UG_ASSERT(	sidesOut[1] != NULL,
				"Exactly two edges should be adjacent to a vertex in a face element.")
}

///	Collects all faces (= 2) which exist in the given volume and which share the given edge.
/**	This method uses Grid::mark **/
UG_API
inline void CollectAssociatedSides(Face* sidesOut[2], Grid& grid, Volume* v, EdgeBase* e)
{
	//PROFILE_BEGIN(CollectAssociatedSides_EDGE);
	sidesOut[0] = NULL;
	sidesOut[1] = NULL;

	grid.begin_marking();

	for(size_t i = 0; i < v->num_vertices(); ++i)
		grid.mark(v->vertex(i));

	vector<Face*> vNeighbourFacesToEdge;
	CollectAssociated(vNeighbourFacesToEdge, grid, e, true);

	//Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(e);
	Grid::AssociatedFaceIterator iterEnd = vNeighbourFacesToEdge.end();
	//for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(e); iter != iterEnd; ++iter)
	for(Grid::AssociatedFaceIterator iter = vNeighbourFacesToEdge.begin(); iter != iterEnd; ++iter)
	{
		Face* f = *iter;

	//	check whether all vertices of f are marked
		bool allMarked = true;
		for(size_t i = 0; i < f->num_vertices(); ++i){
			if(!grid.is_marked(f->vertex(i))){
				allMarked = false;
				break;
			}
		}

		if(allMarked){
			if(FaceContains(f, e)){
				UG_ASSERT(	sidesOut[1] == NULL,
							"Only two faces may be adjacent to an edge in a volume element.")

				if(sidesOut[0] == NULL)
					sidesOut[0] = f;
				else
					sidesOut[1] = f;
			}
		}
	}

	grid.end_marking();

	UG_ASSERT(	sidesOut[1] != NULL,
				"Exactly two faces should be adjacent to an edge in a volume element.")
}



////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateSubsetSurfaceArea
number CalculateSubsetSurfaceArea(MultiGrid& mg, int subsetIndex, MGSubsetHandler& sh);

////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateSubsetVolume
number CalculateSubsetVolume(MultiGrid& mg, int subsetIndex, MGSubsetHandler& sh);



////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinAngle
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMinAngle(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Face* f, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
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
		Vertex* vrt = f->vertex(vrtIter);

	//	Get adjacent edges at the current vertex and calculate the angle between their normals
		CollectAssociatedSides(vNeighbourEdgesToVertex, grid, f, vrt);

	//	Calculate vExtrDir vectors of the current two adjacent edges
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

//	Volume (For volume elements the minimum of dihedral and edge angle will be returned.)
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Volume* v, TAAPosVRT& aaPos)

//	CORRECT VERSION BUT IN COMBINATION WITH SIMULATION ONLY UNDER DEBUG
{
	//PROFILE_FUNC();
	number minDihedral;
	number tmpMinEdgeAngle;
	number minEdgeAngle = 360.0;

//	Calculate the minimal dihedral
	minDihedral = CalculateMinDihedral(grid, v, aaPos);

//	Calculate the minimal edge angle
	for(uint i = 0; i < v->num_faces(); ++i)
	{
		tmpMinEdgeAngle = CalculateMinAngle(grid, grid.get_face(v, i), aaPos);
		if(tmpMinEdgeAngle < minEdgeAngle)
		{
			minEdgeAngle = tmpMinEdgeAngle;
		}
	}

//	return the minimum of minimal dihedral and minimal edge angle
	return min(minDihedral, minEdgeAngle);
}

//	INCORRECT VERSION
/*{
	UG_LOG("------CalculateMinAngle() ------" << endl);
	//PROFILE_FUNC();
	number minDihedral;
	number tmpMinEdgeAngle;
	number minEdgeAngle = 360;

//	Calculate the minimal dihedral
	minDihedral = CalculateMinDihedral(grid, v, aaPos);

//	Calculate the minimal edge angle
	for(uint i = 0; i < v->num_faces(); ++i)
	{
		minEdgeAngle = CalculateMinAngle(grid, grid.get_face(v, i), aaPos);
		UG_LOG("tmpMinEdgeAngle = " << tmpMinEdgeAngle << "minEdgeAngle = " << minEdgeAngle << endl);
		if(minEdgeAngle < tmpMinEdgeAngle)
		{
			minEdgeAngle = tmpMinEdgeAngle;
		}
	}

//	return the minimum of minimal dihedral and minimal edge angle
	return min(minDihedral, minEdgeAngle);
}*/

//	INCORRECT VERSION LOGGED
/*{
	//PROFILE_FUNC();
	number minDihedral;
	number tmpMinEdgeAngle;
	number minEdgeAngle = 360;

//	Calculate the minimal dihedral
	minDihedral = CalculateMinDihedral(grid, v, aaPos);

//	Calculate the minimal edge angle
	for(uint i = 0; i < v->num_faces(); ++i)
	{
		minEdgeAngle = CalculateMinAngle(grid, grid.get_face(v, i), aaPos);
		if(minEdgeAngle < tmpMinEdgeAngle)
		{
			minEdgeAngle = tmpMinEdgeAngle;
		}
	}

//	return the minimum of minimal dihedral and minimal edge angle
	return min(minDihedral, minEdgeAngle);
}*/


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinDihedral
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Tetrahedron
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(tet), aaPos);
}

//	Prism
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(prism), aaPos);
}

//	Pyramid
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(pyr), aaPos);
}

//	Volume
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Volume* v, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
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
	number minDihedral = 360.0;
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

	/*	!!!	Beware of the correct vExtrDir normals to get correct angle value !!!
		INFO:	Angles of a regular tetrahedron:
				- "Tetrahedron (dihedral) angle x" 	= 109,471...deg
				- "Face-to-face angle y"			= 180deg - x = 70,52...deg

		--> Old version:
		//VecScale(vNorm1, vNorm1, -1);
		//tmpAngle = acos(VecDot(vNorm1, vNorm2));

		--> New version:
			(s.	"Qualit�ts-Metriken und Optimierung von Tetraeder-Netzen",
				 Seminararbeit von Johannes Ahlmann, Universit�t Karlsruhe)
	*/
		tmpAngle = acos(VecDot(vNorm1, vNorm2));
		tmpAngle = PI - tmpAngle;

	//	Check for minimality
		if(tmpAngle < minDihedral)
		{
			minDihedral = tmpAngle;
		}
	}

//	Transform minAngle from RAD to DEG
	minDihedral = 180/PI * minDihedral;

	return minDihedral;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMaxAngle
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Face* f, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
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
	number maxAngle = 0;
	number tmpAngle;
	EdgeBase* vNeighbourEdgesToVertex[2];

//	Iterate over all face vertices
	for(uint vrtIter = 0; vrtIter < numFaceVrts; ++vrtIter)
	{
		Vertex* vrt = f->vertex(vrtIter);

	//	Get adjacent edges at the current vertex and calculate the angle between their normals
		CollectAssociatedSides(vNeighbourEdgesToVertex, grid, f, vrt);

	//	Calculate vExtrDir vectors of the current two adjacent edges
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

	//	Check for maximality
		if(tmpAngle > maxAngle)
		{
			maxAngle = tmpAngle;
		}
	}

//	Transform maxAngle from RAD to DEG
	maxAngle = 180/PI * maxAngle;

	return maxAngle;
}

//	Tetrahedron
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	return CalculateMaxAngle(grid, static_cast<Volume*>(tet), aaPos);
}

//	Prism
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	return CalculateMaxAngle(grid, static_cast<Volume*>(prism), aaPos);
}

//	Pyramid
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	return CalculateMaxAngle(grid, static_cast<Volume*>(pyr), aaPos);
}

//	Volume (For volume elements the maximum of dihedral and edge angle will be returned.)
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Volume* v, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
	number maxDihedral;
	number tmpMaxEdgeAngle;
	number maxEdgeAngle = 0.0;

//	Calculate the maximal dihedral
	maxDihedral = CalculateMaxDihedral(grid, v, aaPos);

//	Calculate the maximal edge angle
	for(uint i = 0; i < v->num_faces(); ++i)
	{
		tmpMaxEdgeAngle = CalculateMaxAngle(grid, grid.get_face(v, i), aaPos);

		if(tmpMaxEdgeAngle > maxEdgeAngle)
		{
			maxEdgeAngle = tmpMaxEdgeAngle;
		}
	}

//	return the maximum of maximal dihedral and maximal edge angle
	return max(maxDihedral, maxEdgeAngle);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMaxDihedral
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Tetrahedron
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(tet), aaPos);
}

//	Prism
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(prism), aaPos);
}

//	Pyramid
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(pyr), aaPos);
}

//	Volume
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Volume* v, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
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
	number maxDihedral = 0;
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

	/*	!!!	Beware of the correct vExtrDir normals to get correct angle value !!!
		INFO:	Angles of a regular tetrahedron:
				- "Tetrahedron (dihedral) angle x" 	= 109,471...deg
				- "Face-to-face angle y"			= 180deg - x = 70,52...deg

		--> Old version:
		//VecScale(vNorm1, vNorm1, -1);
		//tmpAngle = acos(VecDot(vNorm1, vNorm2));

		--> New version:
			(s.	"Qualit�ts-Metriken und Optimierung von Tetraeder-Netzen",
				 Seminararbeit von Johannes Ahlmann, Universit�t Karlsruhe)
	*/
		tmpAngle = acos(VecDot(vNorm1, vNorm2));
		tmpAngle = PI - tmpAngle;

	//	Check for maximality
		if(tmpAngle > maxDihedral)
		{
			maxDihedral = tmpAngle;
		}
	}

//	Transform maxDihedral from RAD to DEG
	maxDihedral = 180/PI * maxDihedral;

	return maxDihedral;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinTriangleHeight
template <class TAAPosVRT>
number CalculateMinTriangleHeight(Face* face, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
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


/*
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
		//	MINHEIGHT / MAXEDGELENGTH
		//  optimal Aspect Ratio of a regular triangle
		//  Q = sqrt(3)/2 * a / a = 0.86602540378444...

		//	Calculate minimal triangle height
			number minTriangleHeight = CalculateMinTriangleHeight(face, aaPos);

		//	Calculate the aspect ratio
			AspectRatio = minTriangleHeight / maxEdgeLength;

			return AspectRatio;
		}

		case ROID_QUADRILATERAL:
		{
		//  AREA / MAXEDGELENGTH

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

	//	MINHEIGHT / MAXEDGELENGTH
	// 	optimal Aspect Ratio of a regular tetrahedron
	//	 Q = sqrt(2/3) * a / a = 0.81...

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

//  VOLUME / MAXEDGELENGTH

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

//  VOLUME / MAXEDGELENGTH

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
*/


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateAspectRatio
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Constrained Triangles supported)
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Face* face, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
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
		//	MINHEIGHT / MAXEDGELENGTH
		//  optimal Aspect Ratio of a regular triangle
		//  Q = sqrt(3)/2 * a / a = 0.86602540378444...

		//	Calculate minimal triangle height
			number minTriangleHeight = CalculateMinTriangleHeight(face, aaPos);

		//	Calculate the aspect ratio
			AspectRatio = minTriangleHeight / maxEdgeLength;
			return AspectRatio;
		}

		default:
		{
		 	UG_THROW("Note: Currently only faces of type triangle supported in aspect ratio calculation.");
		 	break;
		}
	}

	return NAN;
}

//	Tetrahedron
template <class TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
	number AspectRatio;

	//	MINHEIGHT / MAXEDGELENGTH
	// 	optimal Aspect Ratio of a regular tetrahedron
	//	 Q = sqrt(2/3) * a / a = 0.81...

//	Calculate the aspect ratio
	AspectRatio = CalculateTetrahedronAspectRatio(grid, tet, aaPos);

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

		default:
		{
		 	UG_THROW("Note: Currently only volumes of type tetrahedron supported in aspect ratio calculation.");
		 	break;
		}
	}

	return NAN;
}






////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithSmallestMinAngle
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithSmallestMinAngle(Grid& grid, TIterator elementsBegin, TIterator elementsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	if(elementsBegin == elementsEnd)
		return NULL;

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
//	FindVolumeWithSmallestMinDihedral
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindVolumeWithSmallestMinDihedral(Grid& grid, TIterator elementsBegin, TIterator elementsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	if(elementsBegin == elementsEnd)
		return NULL;

//	Initializations
	typename TIterator::value_type elementWithSmallestMinDihedral = *elementsBegin;
	number smallestMinDihedral = CalculateMinDihedral(grid, elementWithSmallestMinDihedral, aaPos);
	++elementsBegin;

//	compare all volumes and find that one with smallest minAngle
	for(; elementsBegin != elementsEnd; ++elementsBegin)
	{
		typename TIterator::value_type curElement = *elementsBegin;
		number curSmallestMinDihedral = CalculateMinDihedral(grid, curElement, aaPos);

		if(curSmallestMinDihedral < smallestMinDihedral)
		{
			elementWithSmallestMinDihedral = curElement;
			smallestMinDihedral = curSmallestMinDihedral;
		}
	}

	return elementWithSmallestMinDihedral;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithLargestMaxAngle
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithLargestMaxAngle(Grid& grid, TIterator elementsBegin, TIterator elementsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	if(elementsBegin == elementsEnd)
		return NULL;

//	Initializations
	typename TIterator::value_type elementWithLargestMaxAngle = *elementsBegin;
	number largestMaxAngle = CalculateMaxAngle(grid, elementWithLargestMaxAngle, aaPos);
	++elementsBegin;

//	compare all volumes and find that one with largest maxAngle
	for(; elementsBegin != elementsEnd; ++elementsBegin)
	{
		typename TIterator::value_type curElement = *elementsBegin;
		number curLargestMaxAngle = CalculateMaxAngle(grid, curElement, aaPos);

		if(curLargestMaxAngle > largestMaxAngle)
		{
			elementWithLargestMaxAngle = curElement;
			largestMaxAngle = curLargestMaxAngle;
		}
	}

	return elementWithLargestMaxAngle;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindVolumeWithLargestMaxDihedral
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindVolumeWithLargestMaxDihedral(Grid& grid, TIterator elementsBegin, TIterator elementsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	if(elementsBegin == elementsEnd)
		return NULL;

//	Initializations
	typename TIterator::value_type elementWithLargestMaxDihedral = *elementsBegin;
	number largestMaxDihedral = CalculateMaxDihedral(grid, elementWithLargestMaxDihedral, aaPos);
	++elementsBegin;

//	compare all volumes and find that one with largest maxDihedral
	for(; elementsBegin != elementsEnd; ++elementsBegin)
	{
		typename TIterator::value_type curElement = *elementsBegin;
		number curLargestMaxDihedral = CalculateMaxDihedral(grid, curElement, aaPos);

		if(curLargestMaxDihedral > largestMaxDihedral)
		{
			elementWithLargestMaxDihedral = curElement;
			largestMaxDihedral = curLargestMaxDihedral;
		}
	}

	return elementWithLargestMaxDihedral;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindLargestFace
template <class TIterator, class TAAPosVRT>
Face* FindLargestFace(TIterator facesBegin, TIterator facesEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
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
FindSmallestVolume(TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	if(volumesBegin == volumesEnd)
		return NULL;

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
FindLargestVolume(TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	if(volumesBegin == volumesEnd)
		return NULL;

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
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithSmallestAspectRatio(Grid& grid, 	TIterator elemsBegin,
												TIterator elemsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	if(elemsBegin == elemsEnd)
		return NULL;

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
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithLargestAspectRatio(Grid& grid,  	TIterator elemsBegin,
												TIterator elemsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	if(elemsBegin == elemsEnd)
		return NULL;

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


}	 
#endif

