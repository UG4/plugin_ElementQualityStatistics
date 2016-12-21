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
//#include "volume_calculation.h"


using namespace std;


namespace ug {



////////////////////////////////////////////////////////////////////////////////////////////
//	CollectAssociatedSides
////////////////////////////////////////////////////////////////////////////////////////////

///	Collects all edges (= 2) which exist in the given face and which share the given vertex.
/**	This method uses Grid::mark **/
UG_API
inline void CollectAssociatedSides(Edge* sidesOut[2], Grid& grid, Face* f, Vertex* vrt)
{
	//PROFILE_BEGIN(CollectAssociatedSides_VERTEX);
	sidesOut[0] = NULL;
	sidesOut[1] = NULL;

	grid.begin_marking();
	for(size_t i = 0; i < f->num_vertices(); ++i){
		grid.mark(f->vertex(i));
	}

	//vector<Edge*> vNeighbourEdgesToVertex;
	//CollectAssociated(vNeighbourEdgesToVertex, grid, vrt, true);

	Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(vrt);
	//Grid::AssociatedEdgeIterator iterEnd = vNeighbourEdgesToVertex.end();
	for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(vrt); iter != iterEnd; ++iter)
	//for(Grid::AssociatedEdgeIterator iter = vNeighbourEdgesToVertex.begin(); iter != iterEnd; ++iter)
	{
		Edge* e = *iter;
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
inline void CollectAssociatedSides(Face* sidesOut[2], Grid& grid, Volume* v, Edge* e)
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
	Edge* vNeighbourEdgesToVertex[2];
	Vertex* adjacentVrt1;
	Vertex* adjacentVrt2;

//	Iterate over all face vertices
	for(uint vrtIter = 0; vrtIter < numFaceVrts; ++vrtIter)
	{
		Vertex* vrt = f->vertex(vrtIter);

	//	Get adjacent edges at the current vertex and calculate the angle between their normals
		CollectAssociatedSides(vNeighbourEdgesToVertex, grid, f, vrt);

	//	Calculate vDir vectors of the current two adjacent edges
	//	!!!	Beware of the correct order of the vertices	to get correct angle value !!!
		if(vrt != vNeighbourEdgesToVertex[0]->vertex(0))
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(0);
		else
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(1);

		if(vrt != vNeighbourEdgesToVertex[1]->vertex(0))
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(0);
		else
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(1);

		VecSubtract(vDir1, aaPos[adjacentVrt1], aaPos[vrt]);
		VecSubtract(vDir2, aaPos[adjacentVrt2], aaPos[vrt]);

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

//	INFO:
//	For volume elements the minAngle corresponds to the smallest dihedral

//	Tetrahedron
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(tet), aaPos);
}

//	Prism
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(prism), aaPos);
}

//	Pyramid
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(pyr), aaPos);
}

//	Volume
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Volume* vol, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, vol, aaPos);
}


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
		Edge* e = grid.get_edge(v, eIter);

	//	Get adjacent faces at the current edge and calculate the angle between their normals
	//	!!!	Beware of the correct vExtrDir normals to get correct angle value !!!
		CollectAssociatedSides(vNeighbourFacesToEdge, grid, v, e);

		Vertex* adjacentVrt1;
		Vertex* adjacentVrt2;
		vector3 vDir1;
		vector3 vDir2;
		vector3 vDir3;

	//	Get vertex of each associated face, which is not part of the given edge
		for(size_t i = 0; i < vNeighbourFacesToEdge[0]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(1))
				adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(i);
		}

		for(size_t i = 0; i < vNeighbourFacesToEdge[1]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(1))
				adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(i);
		}

		/*
		if(vNeighbourFacesToEdge[0]->vertex(0) != e->vertex(0) && vNeighbourFacesToEdge[0]->vertex(0) != e->vertex(1))
			adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(0);
		else if(vNeighbourFacesToEdge[0]->vertex(1) != e->vertex(0) && vNeighbourFacesToEdge[0]->vertex(1) != e->vertex(1))
			adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(1);
		else
			adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(2);

		if(vNeighbourFacesToEdge[1]->vertex(0) != e->vertex(0) && vNeighbourFacesToEdge[1]->vertex(0) != e->vertex(1))
			adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(0);
		else if(vNeighbourFacesToEdge[1]->vertex(1) != e->vertex(0) && vNeighbourFacesToEdge[1]->vertex(1) != e->vertex(1))
			adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(1);
		else
			adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(2);
		*/

		VecSubtract(vDir1, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		VecSubtract(vDir2, aaPos[adjacentVrt1], aaPos[e->vertex(0)]);
		VecSubtract(vDir3, aaPos[e->vertex(0)], aaPos[adjacentVrt2]);

		VecCross(vNorm1, vDir1, vDir2);
		VecCross(vNorm2, vDir1, vDir3);

		VecNormalize(vNorm1, vNorm1);
		VecNormalize(vNorm2, vNorm2);

	/*
		INFO:	Angles of a regular tetrahedron:
				- "Tetrahedron (dihedral) angle x" 	= 109,471...deg
				- "Face-to-face angle y"			= 180deg - x = 70,52...deg

		--> Old version:
		//VecScale(vNorm1, vNorm1, -1);
		//tmpAngle = acos(VecDot(vNorm1, vNorm2));

		--> New version:
			(s.	"Qualitaets-Metriken und Optimierung von Tetraeder-Netzen",
				 Seminararbeit von Johannes Ahlmann, Universitaet Karlsruhe)
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
	Edge* vNeighbourEdgesToVertex[2];
	Vertex* adjacentVrt1;
	Vertex* adjacentVrt2;

//	Iterate over all face vertices
	for(uint vrtIter = 0; vrtIter < numFaceVrts; ++vrtIter)
	{
		Vertex* vrt = f->vertex(vrtIter);

	//	Get adjacent edges at the current vertex and calculate the angle between their normals
		CollectAssociatedSides(vNeighbourEdgesToVertex, grid, f, vrt);

	//	Calculate vExtrDir vectors of the current two adjacent edges
	//	!!!	Beware of the correct order of the vertices	to get correct angle value !!!
		if(vrt != vNeighbourEdgesToVertex[0]->vertex(0))
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(0);
		else
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(1);

		if(vrt != vNeighbourEdgesToVertex[1]->vertex(0))
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(0);
		else
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(1);

		VecSubtract(vDir1, aaPos[adjacentVrt1], aaPos[vrt]);
		VecSubtract(vDir2, aaPos[adjacentVrt2], aaPos[vrt]);

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

//	INFO:
//	For volume elements the maxAngle corresponds to the largest dihedral

//	Tetrahedron
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(tet), aaPos);
}

//	Prism
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(prism), aaPos);
}

//	Pyramid
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(pyr), aaPos);
}

//	Volume
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Volume* vol, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, vol, aaPos);
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
		Edge* e = grid.get_edge(v, eIter);

	//	Get adjacent faces at the current edge and calculate the angle between their normals
	//	!!!	Beware of the correct vExtrDir normals to get correct angle value !!!
		CollectAssociatedSides(vNeighbourFacesToEdge, grid, v, e);

		Vertex* adjacentVrt1;
		Vertex* adjacentVrt2;
		vector3 vDir1;
		vector3 vDir2;
		vector3 vDir3;

		for(size_t i = 0; i < vNeighbourFacesToEdge[0]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(1))
				adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(i);
		}

		for(size_t i = 0; i < vNeighbourFacesToEdge[1]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(1))
				adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(i);
		}

		/*
		if(vNeighbourFacesToEdge[0]->vertex(0) != e->vertex(0) && vNeighbourFacesToEdge[0]->vertex(0) != e->vertex(1))
			adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(0);
		else if(vNeighbourFacesToEdge[0]->vertex(1) != e->vertex(0) && vNeighbourFacesToEdge[0]->vertex(1) != e->vertex(1))
			adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(1);
		else
			adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(2);

		if(vNeighbourFacesToEdge[1]->vertex(0) != e->vertex(0) && vNeighbourFacesToEdge[1]->vertex(0) != e->vertex(1))
			adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(0);
		else if(vNeighbourFacesToEdge[1]->vertex(1) != e->vertex(0) && vNeighbourFacesToEdge[1]->vertex(1) != e->vertex(1))
			adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(1);
		else
			adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(2);
		*/

		VecSubtract(vDir1, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		VecSubtract(vDir2, aaPos[adjacentVrt1], aaPos[e->vertex(0)]);
		VecSubtract(vDir3, aaPos[e->vertex(0)], aaPos[adjacentVrt2]);

		VecCross(vNorm1, vDir1, vDir2);
		VecCross(vNorm2, vDir1, vDir3);

		VecNormalize(vNorm1, vNorm1);
		VecNormalize(vNorm2, vNorm2);

	/*
		INFO:	Angles of a regular tetrahedron:
				- "Tetrahedron (dihedral) angle x" 	= 109,471...deg
				- "Face-to-face angle y"			= 180deg - x = 70,52...deg

		--> Old version:
		//VecScale(vNorm1, vNorm1, -1);
		//tmpAngle = acos(VecDot(vNorm1, vNorm2));

		--> New version:
			(s.	"Qualitaets-Metriken und Optimierung von Tetraeder-Netzen",
				 Seminararbeit von Johannes Ahlmann, Universitaet Karlsruhe)
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


//	Volume
template <class TAAPosVRT>
void CalculateVolumeDihedrals(vector<number>& vDihedralsOut, Grid& grid, Volume* v, TAAPosVRT& aaPos)
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
	number tmpAngle;
	Face* vNeighbourFacesToEdge[2];

//	Iterate over all element edges
	for(uint eIter = 0; eIter < numElementEdges; ++eIter)
	{
		Edge* e = grid.get_edge(v, eIter);

	//	Get adjacent faces at the current edge and calculate the angle between their normals
	//	!!!	Beware of the correct vExtrDir normals to get correct angle value !!!
		CollectAssociatedSides(vNeighbourFacesToEdge, grid, v, e);

		Vertex* adjacentVrt1;
		Vertex* adjacentVrt2;
		vector3 vDir1;
		vector3 vDir2;
		vector3 vDir3;

		for(size_t i = 0; i < vNeighbourFacesToEdge[0]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(1))
				adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(i);
		}

		for(size_t i = 0; i < vNeighbourFacesToEdge[1]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(1))
				adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(i);
		}

		VecSubtract(vDir1, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		VecSubtract(vDir2, aaPos[adjacentVrt1], aaPos[e->vertex(0)]);
		VecSubtract(vDir3, aaPos[e->vertex(0)], aaPos[adjacentVrt2]);

		VecCross(vNorm1, vDir1, vDir2);
		VecCross(vNorm2, vDir1, vDir3);

		VecNormalize(vNorm1, vNorm1);
		VecNormalize(vNorm2, vNorm2);

	/*
		INFO:	Angles of a regular tetrahedron:
				- "Tetrahedron (dihedral) angle x" 	= 109,471...deg
				- "Face-to-face angle y"			= 180deg - x = 70,52...deg

		--> Old version:
		//VecScale(vNorm1, vNorm1, -1);
		//tmpAngle = acos(VecDot(vNorm1, vNorm2));

		--> New version:
			(s.	"Qualitaets-Metriken und Optimierung von Tetraeder-Netzen",
				 Seminararbeit von Johannes Ahlmann, Universitaet Karlsruhe)
	*/

		tmpAngle = acos(VecDot(vNorm1, vNorm2));
		tmpAngle = PI - tmpAngle;

	//	Transform maxDihedral from RAD to DEG
		tmpAngle = 180.0/PI * tmpAngle;

		vDihedralsOut.push_back(tmpAngle);
	}
}


//	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
void CalculateFaceAngles(vector<number>& vAnglesOut, Grid& grid, Face* f, TAAPosVRT& aaPos)
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
	number tmpAngle;
	Edge* vNeighbourEdgesToVertex[2];
	Vertex* adjacentVrt1;
	Vertex* adjacentVrt2;

//	Iterate over all face vertices
	for(uint vrtIter = 0; vrtIter < numFaceVrts; ++vrtIter)
	{
		Vertex* vrt = f->vertex(vrtIter);

	//	Get adjacent edges at the current vertex and calculate the angle between their normals
		CollectAssociatedSides(vNeighbourEdgesToVertex, grid, f, vrt);

	//	Calculate vDir vectors of the current two adjacent edges
	//	!!!	Beware of the correct order of the vertices	to get correct angle value !!!
		if(vrt != vNeighbourEdgesToVertex[0]->vertex(0))
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(0);
		else
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(1);

		if(vrt != vNeighbourEdgesToVertex[1]->vertex(0))
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(0);
		else
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(1);

		VecSubtract(vDir1, aaPos[adjacentVrt1], aaPos[vrt]);
		VecSubtract(vDir2, aaPos[adjacentVrt2], aaPos[vrt]);

	//	Normalize
		VecNormalize(vDir1, vDir1);
		VecNormalize(vDir2, vDir2);

	//	Calculate current angle
		tmpAngle = acos(VecDot(vDir1, vDir2));

	//	Transform minAngle from RAD to DEG
		tmpAngle = 180/PI * tmpAngle;

		vAnglesOut.push_back(tmpAngle);
	}
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
	vector<Edge*> edges;
	CollectAssociated(edges, grid, face);
	Edge* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
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
	vector<Edge*> edges;
	CollectAssociated(edges, grid, tet);
	Edge* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
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
	vector<Edge*> edges;
	CollectAssociated(edges, grid, prism);
	Edge* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
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
	vector<Edge*> edges;
	CollectAssociated(edges, grid, pyr);
	Edge* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
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
	vector<Edge*> edges;
	CollectAssociated(edges, grid, face);
	Edge* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
	maxEdgeLength = EdgeLength(longestEdge, aaPos);

	switch (face->reference_object_id())
	{
		case ROID_TRIANGLE:
		{
		/*
		 * optimal Aspect Ratio of a regular tetrahedron with edge lengths a:
		 * Q = hmin/lmax = sqrt(3)/2*a / a = 0.866...
		 *
		 * Info: return value is normalized by factor 2/sqrt(3)
		 * (s. Shewchuk 2002)
		 */

		//	Calculate minimal triangle height
			number minTriangleHeight = CalculateMinTriangleHeight(face, aaPos);

		//	Calculate the aspect ratio
			AspectRatio = 2/std::sqrt(3.0) * minTriangleHeight / maxEdgeLength;
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

	/*
	 * optimal Aspect Ratio of a regular tetrahedron with edge lengths a:
	 * Q = hmin/lmax = sqrt(2/3)*a / a = 0.8164...
	 *
	 * Info: return value is normalized by factor sqrt(3/2)
	 * (s. Shewchuk 2002)
	 */

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
//	CalculateAspectRatio
////////////////////////////////////////////////////////////////////////////////////////////

//	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateVolToRMSFaceAreaRatio(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

//	Face (Triangles and Constrained Triangles supported)
template <class TAAPosVRT>
number CalculateVolToRMSFaceAreaRatio(Grid& grid, Face* face, TAAPosVRT& aaPos)
{
	UG_THROW("CalculateVolToRMSFaceAreaRatio: Currently only volumes of type tetrahedron supported in volume to root-mean-square face area ratio calculation.");
	return NAN;
}

//	Tetrahedron
template <class TAAPosVRT>
number CalculateVolToRMSFaceAreaRatio(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
	number ratio;

	/*
	 * optimal volume to root-mean-square face area ratio of a
	 * regular tetrahedron with edge lengths a:
	 * Q = V/A_rms^(3/2)
	 *
	 * Info: return value is normalized by factor pow(3, 7/4.0) / 2.0 / sqrt(2);
	 * (s. Shewchuk 2002)
	 */

//	Calculate the ratio
	ratio = CalculateTetrahedronVolToRMSFaceAreaRatio(grid, tet, aaPos);

	return ratio;
}

//	Volume
template <class TAAPosVRT>
number CalculateVolToRMSFaceAreaRatio(Grid& grid, Volume* vol, TAAPosVRT& aaPos)
{
	switch (vol->reference_object_id())
	{
		case ROID_TETRAHEDRON:
		{
			return CalculateVolToRMSFaceAreaRatio(grid, static_cast<Tetrahedron*>(vol), aaPos);
		}

		default:
		{
		 	UG_THROW("CalculateVolToRMSFaceAreaRatio: Currently only volumes of type tetrahedron supported in aspect ratio calculation.");
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
	number smallestVolumeVolume = CalculateVolume(smallestVolume, aaPos);
	++volumesBegin;

//	compare all tetrahedrons and find minimal volume
	for(; volumesBegin != volumesEnd; ++volumesBegin)
	{
		Volume* curVolume = *volumesBegin;
		number curVolumeVolume = CalculateVolume(curVolume, aaPos);

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
	number largestVolumeVolume = CalculateVolume(largestVolume, aaPos);
	++volumesBegin;

//	compare all tetrahedrons and find minimal volume
	for(; volumesBegin != volumesEnd; ++volumesBegin)
	{
		Volume* curVolume = *volumesBegin;
		number curVolumeVolume = CalculateVolume(curVolume, aaPos);

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
