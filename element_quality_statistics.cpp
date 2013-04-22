/*
 * element_quality_statistics.cpp
 *
 *  Created on: 17.04.2012
 *      Author: Martin Stepniewski
 */


#include "element_quality_statistics.h"
#include "common/util/table.h"


namespace ug
{



////////////////////////////////////////////////////////////////////////////////////////////
//	CollectAssociatedSides
////////////////////////////////////////////////////////////////////////////////////////////

///	Collects all edges (= 2) which exist in the given face and which share the given vertex.
/**	This method uses Grid::mark **/
UG_API
inline void CollectAssociatedSides(EdgeBase* sidesOut[2], Grid& grid, Face* f, VertexBase* vrt)
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
		VertexBase* vrt = f->vertex(vrtIter);

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
			(s.	"Qualitäts-Metriken und Optimierung von Tetraeder-Netzen",
				 Seminararbeit von Johannes Ahlmann, Universität Karlsruhe)
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
		VertexBase* vrt = f->vertex(vrtIter);

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
			(s.	"Qualitäts-Metriken und Optimierung von Tetraeder-Netzen",
				 Seminararbeit von Johannes Ahlmann, Universität Karlsruhe)
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
//	CalculateSubsetSurfaceArea
number CalculateSubsetSurfaceArea(MultiGrid& mg, int subsetIndex, MGSubsetHandler& sh)
{
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	number subsetSurfaceArea = 0.0;
	for(FaceIterator fIter = sh.begin<Face>(subsetIndex, 0); fIter != sh.end<Face>(subsetIndex, 0); ++fIter)
	{
		Face* f = *fIter;
		subsetSurfaceArea += FaceArea(f, aaPos);
	}

	return subsetSurfaceArea;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateSubsetVolume
number CalculateSubsetVolume(MultiGrid& mg, int subsetIndex, MGSubsetHandler& sh)
{
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	number subsetVolume = 0.0;
	for(VolumeIterator vIter = sh.begin<Volume>(subsetIndex, 0); vIter != sh.end<Volume>(subsetIndex, 0); ++vIter)
	{
		Volume* v = *vIter;
		subsetVolume += CalculateVolume(*v, aaPos);
	}

	return subsetVolume;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	ElementQualityStatistics
////////////////////////////////////////////////////////////////////////////////////////////

//	Wrapper for multigrids
//void ElementQualityStatistics(MultiGrid& mg, int level)
void ElementQualityStatistics(MultiGrid& mg)
{
	if(mg.num_volumes() > 0)
		ElementQualityStatistics3d(mg, mg.get_geometric_objects());
	else
		ElementQualityStatistics2d(mg, mg.get_geometric_objects());
}

//	Wrapper for grids
void ElementQualityStatistics(Grid& grid)
{
	if(grid.num_volumes() > 0)
		ElementQualityStatistics3d(grid, grid.get_geometric_objects());
	else
		ElementQualityStatistics2d(grid, grid.get_geometric_objects());
}

//	Actual procedures
void ElementQualityStatistics2d(Grid& grid, GeometricObjectCollection goc)
{
	//PROFILE_FUNC();
	Grid::VertexAttachmentAccessor<APosition2> aaPos(grid, aPosition2);


//	Numbers
	number n_minEdge = 0.0;
	number n_maxEdge = 0.0;
	number n_minFace = 0.0;
	number n_maxFace = 0.0;
	number n_minFaceAngle = 0.0;
	number n_maxFaceAngle = 0.0;
	number n_minTriAspectRatio = 0.0;
	number n_maxTriAspectRatio = 0.0;


//	Elements
	EdgeBase* minEdge;
	EdgeBase* maxEdge;
	Face* minFace;
	Face* maxFace;
	Face* minAngleFace;
	Face* maxAngleFace;
	Face* minAspectRatioTri;
	Face* maxAspectRatioTri;


//	Basic grid properties on level i
	UG_LOG(endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl);
	UG_LOG("GRID QUALITY STATISTICS" << endl << endl);
	UG_LOG("*** Output info:" << endl);
	UG_LOG("    - The 'aspect ratio' (AR) represents the ratio of minimal height and " << endl <<
		   "      maximal edge length of a triangle or tetrahedron respectively." << endl);
	UG_LOG("    - The MinAngle-Histogram lists the number of minimal element angles in " << endl <<
		   "      different degree ranges." << endl << endl);
	for(uint i = 0; i < goc.num_levels(); ++i)
	{
		//PROFILE_BEGIN(eqs_qualityStatistics2d);
	//	----------
	//	2D section
	//	----------
		minEdge = FindShortestEdge(goc.begin<EdgeBase>(i), goc.end<EdgeBase>(i), aaPos);
		maxEdge = FindLongestEdge(goc.begin<EdgeBase>(i), goc.end<EdgeBase>(i), aaPos);
		minFace = FindSmallestFace(goc.begin<Face>(i), goc.end<Face>(i), aaPos);
		maxFace = FindLargestFace(goc.begin<Face>(i), goc.end<Face>(i), aaPos);
		minAngleFace = FindElementWithSmallestMinAngle(	grid,
														goc.begin<Face>(i),
														goc.end<Face>(i),
														aaPos);
		maxAngleFace = FindElementWithLargestMaxAngle(	grid,
														goc.begin<Face>(i),
														goc.end<Face>(i),
														aaPos);
		if(minEdge != NULL)
			n_minEdge = EdgeLength(minEdge, aaPos);
		if(maxEdge != NULL)
			n_maxEdge = EdgeLength(maxEdge, aaPos);
		if(minFace != NULL)
			n_minFace = FaceArea(minFace, aaPos);
		if(maxFace != NULL)
			n_maxFace = FaceArea(maxFace, aaPos);
		if(minAngleFace != NULL)
			n_minFaceAngle = CalculateMinAngle(grid, minAngleFace, aaPos);
		if(maxAngleFace != NULL)
			n_maxFaceAngle = CalculateMaxAngle(grid, maxAngleFace, aaPos);

	//	Check for triangles
		if(goc.num<Triangle>(i) > 0)
		{
			minAspectRatioTri = FindElementWithSmallestAspectRatio(	grid,
																	goc.begin<Face>(i),
																	goc.end<Face>(i),
																	aaPos);
			maxAspectRatioTri = FindElementWithLargestAspectRatio(	grid,
																	goc.begin<Face>(i),
																	goc.end<Face>(i),
																	aaPos);
			if(minAspectRatioTri != NULL)
				n_minTriAspectRatio = CalculateAspectRatio(grid, minAspectRatioTri, aaPos);
			if(maxAspectRatioTri != NULL)
				n_maxTriAspectRatio = CalculateAspectRatio(grid, maxAspectRatioTri, aaPos);
		}
		//PROFILE_END();


		//PROFILE_BEGIN(eqs_qualityStatisticsOutput);
	//	Table summary
		ug::Table<std::stringstream> table(8, 4);
		table(0, 0) << "Number of volumes"; 	table(0, 1) << goc.num_volumes(i);
		table(1, 0) << "Number of faces"; 		table(1, 1) << goc.num_faces(i);
		table(2, 0) << "Number of vertices";	table(2, 1) << goc.num_vertices(i);

		table(3, 0) << " "; table(3, 1) << " ";
		table(3, 2) << " "; table(3, 3) << " ";

		table(4, 0) << "Shortest edge";	table(4, 1) << n_minEdge;
		table(4, 2) << "Longest edge";	table(4, 3) << n_maxEdge;

		table(5, 0) << "Smallest face angle";	table(5, 1) << n_minFaceAngle;
		table(5, 2) << "Largest face angle";	table(5, 3) << n_maxFaceAngle;

		table(6, 0) << "Smallest face";	table(6, 1) << n_minFace;
		table(6, 2) << "Largest face";	table(6, 3) << n_maxFace;

		if(goc.num<Triangle>(i) > 0)
		{
			table(7, 0) << "Smallest triangle AR";	table(7, 1) << n_minTriAspectRatio;
			table(7, 2) << "Largest triangle AR";	table(7, 3) << n_maxTriAspectRatio;
		}


	//	Output section
		UG_LOG("+++++++++++++++++" << endl);
		UG_LOG(" Grid level " << i << ":" << endl);
		UG_LOG("+++++++++++++++++" << endl << endl);
		UG_LOG(table);


		//PROFILE_BEGIN(eqs_minAngleHistogram);
		MinAngleHistogram(grid, goc.begin<Face>(i), goc.end<Face>(i), aaPos, 10);


		//PROFILE_END();
		UG_LOG(endl);
	}

	UG_LOG(endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl);
}

void ElementQualityStatistics3d(Grid& grid, GeometricObjectCollection goc)
{
	//PROFILE_FUNC();
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);


//	Numbers
	number n_minEdge = 0.0;
	number n_maxEdge = 0.0;
	//number n_minFace;
	//number n_maxFace;
	number n_minFaceAngle = 0.0;
	number n_maxFaceAngle = 0.0;
	number n_minTriAspectRatio = 0.0;
	number n_maxTriAspectRatio = 0.0;

	number n_minVolume = 0.0;
	number n_maxVolume = 0.0;
	number n_minVolAngle = 0.0;
	number n_maxVolAngle = 0.0;
	number n_minVolDihedral = 0.0;
	number n_maxVolDihedral = 0.0;
	number n_minTetAspectRatio = 0.0;
	number n_maxTetAspectRatio = 0.0;


//	Elements
	EdgeBase* minEdge;
	EdgeBase* maxEdge;
	//Face* minAreaFace;
	//Face* maxAreaFace;
	Face* minAngleFace;
	Face* maxAngleFace;
	Face* minAspectRatioTri;
	Face* maxAspectRatioTri;

	Volume* minVolume;
	Volume* maxVolume;
	Volume* minAngleVol;
	Volume* maxAngleVol;
	Volume* minDihedralVol;
	Volume* maxDihedralVol;
	Tetrahedron* minAspectRatioTet;
	Tetrahedron* maxAspectRatioTet;


//	Basic grid properties on level i
	UG_LOG(endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl);
	UG_LOG("GRID QUALITY STATISTICS" << endl << endl);
	UG_LOG("*** Output info:" << endl);
	UG_LOG("    - The 'aspect ratio' (AR) represents the ratio of minimal height and " << endl <<
		   "      maximal edge length of a triangle or tetrahedron respectively." << endl);
	UG_LOG("    - The MinAngle-Histogram lists the number of minimal element angles in " << endl <<
		   "      different degree ranges." << endl << endl);
	for(uint i = 0; i < goc.num_levels(); ++i)
	{
		//PROFILE_BEGIN(eqs_qualityStatistics2d);
	//	----------
	//	2D section
	//	----------
		minEdge = FindShortestEdge(goc.begin<EdgeBase>(i), goc.end<EdgeBase>(i), aaPos);
		maxEdge = FindLongestEdge(goc.begin<EdgeBase>(i), goc.end<EdgeBase>(i), aaPos);
		//minFace =
		//maxFace =
		minAngleFace = FindElementWithSmallestMinAngle(	grid,
														goc.begin<Face>(i),
														goc.end<Face>(i),
														aaPos);
		maxAngleFace = FindElementWithLargestMaxAngle(	grid,
														goc.begin<Face>(i),
														goc.end<Face>(i),
														aaPos);
		if(minEdge != NULL)
			n_minEdge = EdgeLength(minEdge, aaPos);
		if(maxEdge != NULL)
			n_maxEdge = EdgeLength(maxEdge, aaPos);
		if(minAngleFace != NULL)
			n_minFaceAngle = CalculateMinAngle(grid, minAngleFace, aaPos);
		if(maxAngleFace != NULL)
			n_maxFaceAngle = CalculateMaxAngle(grid, maxAngleFace, aaPos);

	//	Check for triangles
		if(goc.num<Triangle>(i) > 0)
		{
			minAspectRatioTri = FindElementWithSmallestAspectRatio(	grid,
																	goc.begin<Face>(i),
																	goc.end<Face>(i),
																	aaPos);
			maxAspectRatioTri = FindElementWithLargestAspectRatio(	grid,
																	goc.begin<Face>(i),
																	goc.end<Face>(i),
																	aaPos);
			if(minAspectRatioTri != NULL)
				n_minTriAspectRatio = CalculateAspectRatio(grid, minAspectRatioTri, aaPos);
			if(maxAspectRatioTri != NULL)
				n_maxTriAspectRatio = CalculateAspectRatio(grid, maxAspectRatioTri, aaPos);
		}
		//PROFILE_END();

	//	----------
	//	3D section
	//	----------
		if(goc.num<Volume>(i) > 0)
		{
			//PROFILE_BEGIN(eqs_qualityStatistics3d);
			minVolume = FindSmallestVolume(	goc.begin<Volume>(i),
											goc.end<Volume>(i),
											aaPos);

			maxVolume = FindLargestVolume(	goc.begin<Volume>(i),
											goc.end<Volume>(i),
											aaPos);

			minAngleVol = FindElementWithSmallestMinAngle(	grid,
															goc.volumes_begin(i),
															goc.volumes_end(i),
															aaPos);
			maxAngleVol = FindElementWithLargestMaxAngle(	grid,
															goc.volumes_begin(i),
															goc.volumes_end(i),
				 											aaPos);
			minDihedralVol = FindVolumeWithSmallestMinDihedral(	grid,
																goc.volumes_begin(i),
																goc.volumes_end(i),
																aaPos);
			maxDihedralVol = FindVolumeWithLargestMaxDihedral(	grid,
																goc.volumes_begin(i),
																goc.volumes_end(i),
																aaPos);
			if(minVolume != NULL)
				n_minVolume = CalculateVolume(*minVolume, aaPos);
			if(maxVolume != NULL)
				n_maxVolume = CalculateVolume(*maxVolume, aaPos);
			if(minAngleVol != NULL)
				n_minVolAngle = CalculateMinAngle(grid, minAngleVol, aaPos);
			if(maxAngleVol != NULL)
				n_maxVolAngle = CalculateMaxAngle(grid, maxAngleVol, aaPos);
			if(minDihedralVol != NULL)
				n_minVolDihedral = CalculateMinDihedral(grid, minDihedralVol, aaPos);
			if(maxDihedralVol != NULL)
				n_maxVolDihedral = CalculateMaxDihedral(grid, maxDihedralVol, aaPos);

		//	Tetrahedron section
			if(goc.num<Tetrahedron>(i) > 0)
			{
				minAspectRatioTet = FindElementWithSmallestAspectRatio(	grid,
																		goc.begin<Tetrahedron>(i),
																		goc.end<Tetrahedron>(i),
																		aaPos);
				maxAspectRatioTet = FindElementWithLargestAspectRatio(	grid,
																		goc.begin<Tetrahedron>(i),
																		goc.end<Tetrahedron>(i),
																		aaPos);
				if(minAspectRatioTet != NULL)
					n_minTetAspectRatio = CalculateAspectRatio(grid, minAspectRatioTet, aaPos);
				if(maxAspectRatioTet != NULL)
					n_maxTetAspectRatio = CalculateAspectRatio(grid, maxAspectRatioTet, aaPos);
			}
		}

		//PROFILE_BEGIN(eqs_qualityStatisticsOutput);
	//	Table summary
		ug::Table<std::stringstream> table(11, 4);
		table(0, 0) << "Number of volumes"; 	table(0, 1) << goc.num_volumes(i);
		table(1, 0) << "Number of faces"; 		table(1, 1) << goc.num_faces(i);
		table(2, 0) << "Number of vertices";	table(2, 1) << goc.num_vertices(i);

		table(3, 0) << " "; table(3, 1) << " ";
		table(3, 2) << " "; table(3, 3) << " ";

		table(4, 0) << "Shortest edge";	table(4, 1) << n_minEdge;
		table(4, 2) << "Longest edge";	table(4, 3) << n_maxEdge;

		table(5, 0) << "Smallest face angle";	table(5, 1) << n_minFaceAngle;
		table(5, 2) << "Largest face angle";	table(5, 3) << n_maxFaceAngle;

		if(goc.num_volumes(i) > 0)
		{
			table(6, 0) << "Smallest volume";		table(6, 1) << n_minVolume;
			table(6, 2) << "Largest volume";		table(6, 3) << n_maxVolume;
			table(7, 0) << "Smallest volume angle";	table(7, 1) << n_minVolAngle;
			table(7, 2) << "Largest volume angle";	table(7, 3) << n_maxVolAngle;
			table(8, 0) << "Smallest volume dihedral";	table(8, 1) << n_minVolDihedral;
			table(8, 2) << "Largest volume dihedral";	table(8, 3) << n_maxVolDihedral;

			if(goc.num<Triangle>(i) > 0)
			{
				table(9, 0) << "Smallest triangle AR"; table(9, 1) << n_minTriAspectRatio;
				table(9, 2) << "Largest triangle AR"; table(9, 3) << n_maxTriAspectRatio;
			}

			if(goc.num<Tetrahedron>(i) > 0)
			{
				table(10, 0) << "Smallest tetrahedron AR";	table(10, 1) << n_minTetAspectRatio;
				table(10, 2) << "Largest tetrahedron AR";	table(10, 3) << n_maxTetAspectRatio;
			}
		}


	//	Output section
		UG_LOG("+++++++++++++++++" << endl);
		UG_LOG(" Grid level " << i << ":" << endl);
		UG_LOG("+++++++++++++++++" << endl << endl);
		UG_LOG(table);

		//PROFILE_BEGIN(eqs_minAngleHistogram);
		if(goc.num_volumes(i) > 0)
		{
			MinAngleHistogram(grid, goc.begin<Volume>(i), goc.end<Volume>(i), aaPos, 10);
		}
		else
		{
			MinAngleHistogram(grid, goc.begin<Face>(i), goc.end<Face>(i), aaPos, 10);
		}
		//PROFILE_END();

		UG_LOG(endl);
	}

	UG_LOG(endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	DistributeNPointsOnASphere
////////////////////////////////////////////////////////////////////////////////////////////
void GetNEvenlyDistributedSphereCoords(vector<vector3>& coords, int N, double radius)
{
/*
 * Check, if N defines a platonic solid
 */

/*
//	Tetrahedron
	if(N == 4)
	{
		const number A = sqrt(2);
		const number B = sqrt(6);
		const number vrts[4][3] = {	{ 0,        0      ,  1},
									{ 2.0/3*A,  0      , -1.0/3},
									{-1.0/3*A,  1.0/3*B, -1.0/3},
									{-1.0/3*A, -1.0/3*B, -1.0/3}};

		vector3 tmpCoords;
		for(int i = 0; i < N; i++)
		{
			tmpCoords = vector3(vrts[i][0], vrts[i][1], vrts[i][2]);
			tmpCoords.x *= radius;
			tmpCoords.y *= radius;
			tmpCoords.z *= radius;
			coords.push_back(tmpCoords);
		}
	}


//	Octahedron
	else if(N == 6)
	{
		const number vrts[6][3] = {	{1,0,0},
									{-1,0,0},
									{0,1,0},
									{0,0,1},
									{0,-1,0},
									{0,0,-1}
								  };

		vector3 tmpCoords;
		for(int i = 0; i < N; i++)
		{
			tmpCoords = vector3(vrts[i][0], vrts[i][1], vrts[i][2]);
			tmpCoords.x *= radius;
			tmpCoords.y *= radius;
			tmpCoords.z *= radius;
			coords.push_back(tmpCoords);
		}
	}


//	Cube
	else if(N == 8)
	{
		const number c = 1.0/sqrt(3);
		const number vrts[8][3] = {	{c,c,c},
									{-c,c,c},
									{c,-c,c},
									{c,c,-c},
									{-c,-c,-c},
									{-c,-c,c},
									{-c,c,-c},
									{c,-c,-c}
								  };

		vector3 tmpCoords;
		for(int i = 0; i < N; i++)
		{
			tmpCoords = vector3(vrts[i][0], vrts[i][1], vrts[i][2]);
			tmpCoords.x *= radius;
			tmpCoords.y *= radius;
			tmpCoords.z *= radius;
			coords.push_back(tmpCoords);
		}
	}


//	Icosahedron
	else if(N == 12)
	{
		const number A = 0.85065080835204;
		const number B = 0.525731112119134;
		const number vrts[12][3] = {	{-B, A, 0},
										{0, B, A},
										{B, A, 0},
										{0, B, -A},
										{-A, 0, B},
										{A, 0, B},
										{A, 0, -B},
										{-A, 0, -B},
										{-B, -A, 0},
										{0, -B, A},
										{B, -A, 0},
										{0, -B, -A}
									};

		vector3 tmpCoords;
		for(int i = 0; i < N; i++)
		{
			tmpCoords = vector3(vrts[i][0], vrts[i][1], vrts[i][2]);
			tmpCoords.x *= radius;
			tmpCoords.y *= radius;
			tmpCoords.z *= radius;
			coords.push_back(tmpCoords);
		}
	}


//	Dodecahedron
	else if(N == 20)
	{
		const number c = 1.0/sqrt(3);
		const number x = (sqrt(5)+1)/(2*sqrt(3));
		const number y = (sqrt(5)-1)/(2*sqrt(3));
		const number vrts[20][3] = {	{c,c,c}, {c,-c,c}, {c,c,-c}, {c,-c,-c},
										{-c,c,c}, {-c,-c,c}, {-c,c,-c}, {-c,-c,-c},
										{x,y,0}, {-x,y,0}, {x,-y,0}, {-x,-y,0},
										{0,x,y}, {0,-x,y}, {0,x,-y}, {0,-x,-y},
										{y,0,x}, {-y,0,x}, {y,0,-x}, {-y,0,-x}};

		vector3 tmpCoords;
		for(int i = 0; i < N; i++)
		{
			tmpCoords = vector3(vrts[i][0], vrts[i][1], vrts[i][2]);
			tmpCoords.x *= radius;
			tmpCoords.y *= radius;
			tmpCoords.z *= radius;
			coords.push_back(tmpCoords);
		}
	}
*/

/*
 * otherwise use heuristic
 */

	//else
	//{
		/***********************
		 * Golden spiral section
		 * (2nd best algorithm)
		 **********************/
		/*
		int n = N;
		double r, phi;
		vector3 tmpCoords;

		for(int i = 0; i < n; i++)
		{
			tmpCoords.y 	= i * (double)2/n - 1 + (double)1/n;
			r				= sqrt(1.0 - tmpCoords.y * tmpCoords.y);
			phi 			= i * M_PI * (3.0 - sqrt(5));
			tmpCoords.x		= r * cos(phi);
			tmpCoords.z 	= r * sin(phi);

			tmpCoords.x *= radius;
			tmpCoords.y *= radius;
			tmpCoords.z *= radius;

			coords.push_back(tmpCoords);
		}
		*/


		/****************************************
		 * "Distributing many points on a sphere"
		 * by E.B. Saff and A.B.J. Kuijlaars
		 ***************************************/

		/*
		int n = N;
		double theta[n+1], phi[n+1], h;
		vector3 tmpCoords;

		for(int i = 1; i <= n; i++)
		{
			h = -1 + 2 * (double)(i-1)/(n-1);
			theta[i] = acos(h);

			if(i == 1 || i == n)
				phi[i] = 0;
			else
				phi[i] = fmod(phi[i-1] + 3.6/sqrt(n*(1.0-h*h)), (2*M_PI));

			tmpCoords.x = sin(theta[i])*cos(phi[i]);
			tmpCoords.y = sin(phi[i])*sin(theta[i]);
			tmpCoords.z = cos(theta[i]);

			tmpCoords.x *= radius;
			tmpCoords.y *= radius;
			tmpCoords.z *= radius;

			coords.push_back(tmpCoords);
		}
		*/


		/****************************************
		 * "Distributing many points on a sphere"
		 * by E.B. Saff and A.B.J. Kuijlaars
		 *
		 * Modification #1
		 ***************************************/

		/*
		int n = N;
		vector3 tmpCoords;
		double p = 0.5;
		double a = 1 - 2 * p / (n-3);
		double b = p * (double)(n+1)/(n-1);
		double theta[n+1], phi[n+1], h[n+1], r[n+1];

		for(int i = 1; i <= n; i++)
		{
			if(i == 1)
			{
				r[1] 		= 0.0;
				theta[1] 	= M_PI;
				phi[1] 		= 0.0;
			}

			else if (i == n)
			{
				theta[n] 	= 0.0;
				phi[n] 		= 0.0;
			}

			else
			{
				double k 	= a * i + b;
				h[i] 		= -1 + 2 * (k-1)/(n-1);
				r[i] 		= sqrt(1 - h[i]*h[i]);
				theta[i] 	= acos(h[i]);
				phi[i] 		= fmod(phi[i-1] + 3.6/sqrt(n)*(double)2/(r[i-1]+r[i]), (2*M_PI));
			}

			tmpCoords.x = sin(theta[i])*cos(phi[i]);
			tmpCoords.y = sin(phi[i])*sin(theta[i]);
			tmpCoords.z = cos(theta[i]);

			tmpCoords.x *= radius;
			tmpCoords.y *= radius;
			tmpCoords.z *= radius;

			coords.push_back(tmpCoords);
		}
		*/


		/****************************************
		 * "Distributing many points on a sphere"
		 * by E.B. Saff and A.B.J. Kuijlaars
		 *
		 * Modification #2
		 * (best algorithm)
		 ***************************************/

		int n = N;
		vector3 tmpCoords;
		double theta[n+1], phi[n+1], h[n+1];

		for(int i = 1; i <= n; i++)
		{
			h[i] 		= -1 + (double)(2*i-1)/n;
			theta[i] 	= acos(h[i]);
			phi[i] 		= sqrt(n*M_PI)*theta[i];

			tmpCoords.x = sin(theta[i])*cos(phi[i]);
			tmpCoords.y = sin(phi[i])*sin(theta[i]);
			tmpCoords.z = cos(theta[i]);

			tmpCoords.x *= radius;
			tmpCoords.y *= radius;
			tmpCoords.z *= radius;

			coords.push_back(tmpCoords);
		}
	//}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	BuildBouton
////////////////////////////////////////////////////////////////////////////////////////////
void BuildBouton(number radius, int numRefinements, int numReleaseSites)
{
//	Initial grid management setup
	Grid grid;
	grid.attach_to_vertices(aPosition);
	grid.attach_to_vertices(aNormal);
	//grid.attach_to_faces(aNormal);
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	Grid::VertexAttachmentAccessor<ANormal> aaNorm(grid, aNormal);
	Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);
	//Grid::FaceAttachmentAccessor<ANormal> aaFaceNorm(grid, aNormal);

	SubsetHandler sh(grid);
	sh.set_default_subset_index(0);

	Selector sel(grid);


//	Generate raw icosphere
	vector3 center(0.0, 0.0, 0.0);
	GenerateIcosphere(grid, center, radius, numRefinements, aPosition);
	sh.set_subset_name("bouton_bnd", 0);


//	Get evenly distributed sphere coordinates
	vector<vector3> coords;
	GetNEvenlyDistributedSphereCoords(coords, numReleaseSites, radius);


//	Testwise creation of evenly distributed vertices
	/*
	VertexBase* vrts[numReleaseSites];
	for(size_t i = 0; i < numReleaseSites; ++i)
	{
		vrts[i] = *grid.create<Vertex>();
		aaPos[vrts[i]] = coords[i];
	}
	*/


/*
 * 	Find corresponding vertices on the icospere and assign
 * 	a first set of evenly distributed vertices for
 * 	the mature release sites.
 */
	number minDist, tmpMinDist;

	sel.clear();
	for(VertexBaseIterator vIter = grid.vertices_begin(); vIter != grid.vertices_end(); ++vIter)
	{
		VertexBase* vrt = *vIter;
		sel.select(vrt);
	}

	if(grid.num<Vertex>() > 0)
	{
		for(size_t i = 0; i < coords.size(); ++i)
		{
			bool gotOne = false;
			VertexBase* tmpVrt;

			for(VertexBaseIterator vIter = grid.vertices_begin(); vIter != grid.vertices_end(); ++vIter)
			{
				VertexBase* vrt = *vIter;
				tmpMinDist = VecDistance(aaPos[vrt], coords[i]);

				if(((!gotOne) || (tmpMinDist < minDist)) && sel.is_selected(vrt))
				{
					minDist = tmpMinDist;
					tmpVrt = vrt;
					gotOne = true;
				}
			}

			sel.deselect(tmpVrt);
			if(i % 2 == 0)
			{
				sh.assign_subset(tmpVrt, 4);
			}
			else
			{
				sh.assign_subset(tmpVrt, 2);
			}
		}
	}

	sh.set_subset_name("T-bars_bnd", 4);
	sh.set_subset_name("immature_AZ", 2);

	for(VertexBaseIterator vIter = sh.begin<VertexBase>(2); vIter != sh.end<VertexBase>(2); ++vIter)
	{
		VertexBase* vrt = *vIter;
		for(Grid::AssociatedFaceIterator fIter = grid.associated_faces_begin(vrt); fIter != grid.associated_faces_end(vrt); ++fIter)
		{
			Face* f = *fIter;
			sh.assign_subset(f, 2);
		}
	}

	for(FaceIterator fIter = sh.begin<Face>(2); fIter != sh.end<Face>(2); ++fIter)
	{
		Face* f = *fIter;
		for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(f); eIter != grid.associated_edges_end(f); ++eIter)
		{
			EdgeBase* e = *eIter;
			sh.assign_subset(e, 2);
		}

		for(size_t i = 0; i < f->num_vertices(); ++i)
		{
			VertexBase* vrt = f->vertex(i);
			sh.assign_subset(vrt, 2);
		}
	}



/*
//	Rotate evenly distributed sphere coordinates around x axes by "a" degrees
 	double a = 45.0;
	for(int i = 0; i < coords.size(); ++i)
	{
		double y, z;
		y = coords[i].y;
		z = coords[i].z;
		coords[i].y = y * cos(a) - z * sin(a);
		coords[i].z = y * sin(a) + z * cos(a);
	}
*/

/*
 * 	Find corresponding vertices on the icospere and assign
 * 	a second set of evenly distributed vertices for
 * 	the immature release sites.
 */
/*
	sel.clear();
	for(VertexBaseIterator vIter = grid.vertices_begin(); vIter != grid.vertices_end(); ++vIter)
	{
		VertexBase* vrt = *vIter;
		sel.select(vrt);
	}

	if(grid.num<Vertex>() > 0)
	{
		for(size_t i = 0; i < coords.size(); ++i)
		{
			bool gotOne = false;
			VertexBase* tmpVrt;

			for(VertexBaseIterator vIter = grid.vertices_begin(); vIter != grid.vertices_end(); ++vIter)
			{
				VertexBase* vrt = *vIter;
				tmpMinDist = VecDistance(aaPos[vrt], coords[i]);

				if(((!gotOne) || (tmpMinDist < minDist)) && sel.is_selected(vrt))
				{
					minDist = tmpMinDist;
					tmpVrt = vrt;
					gotOne = true;
				}
			}

			sel.deselect(tmpVrt);
			sh.assign_subset(tmpVrt, 2);

			*//*
			for(Grid::AssociatedFaceIterator fIter = grid.associated_faces_begin(tmpVrt); fIter != grid.associated_faces_end(tmpVrt); ++fIter)
			{
				Face* f = *fIter;
				sh.assign_subset(f, 2);
			}
			*/
/*
		}
	}
*/

//	Divide evenly distributed vertices into two separate subsets
	/*
	sel.clear();
	int counter = 0;
	for(VertexBaseIterator vIter = sh.begin<VertexBase>(1); vIter != sh.end<VertexBase>(1); ++vIter)
	{
		VertexBase* vrt = *vIter;
		if(counter % 2 == 0)
		{
			sel.select(vrt);
		}

		counter++;
	}

	for(VertexBaseIterator vIter = sel.begin<VertexBase>(); vIter != sel.end<VertexBase>(); ++vIter)
	{
		VertexBase* vrt = *vIter;
		sh.assign_subset(vrt, 2);

		for(Grid::AssociatedFaceIterator fIter = grid.associated_faces_begin(vrt); fIter != grid.associated_faces_end(vrt); ++fIter)
		{
			Face* f = *fIter;
			sh.assign_subset(f, 2);
		}
	}
	*/


//	Create mature release sites
	sel.clear();
	CalculateVertexNormals(grid, aaPos, aaNorm);
	vector<VertexBase*> vrts;
	for(VertexBaseIterator vIter = sh.begin<VertexBase>(4); vIter != sh.end<VertexBase>(4); ++vIter)
	{
		VertexBase* vrt = *vIter;
		vrts.push_back(vrt);
	}

	/*
	 * "mature_AZ"
	 */
	for(size_t i = 0; i < vrts.size(); ++i)
	{
		vector3 negNorm;
		VertexBase* vrt = vrts[i];
		VecScale(negNorm, aaNorm[vrt], -1.0);
		AdaptSurfaceGridToCylinder(sel, grid, vrt, negNorm, 0.15, 0.03);

		for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
		{
			Face* f = *fIter;
			sh.assign_subset(f, 1);
			for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(f); eIter != grid.associated_edges_end(f); ++eIter)
			{
				EdgeBase* e = *eIter;
				sh.assign_subset(e, 1);
			}

			for(size_t i = 0; i < f->num_vertices(); ++i)
			{
				VertexBase* vrt = f->vertex(i);
				sh.assign_subset(vrt, 1);
			}
		}

		Refine(grid, sel);
	}

	sh.set_subset_name("mature_AZ", 1);


	/*
	 * CChCl
	 */
	for(size_t i = 0; i < vrts.size(); ++i)
	{
		vector3 negNorm;
		VertexBase* vrt = vrts[i];
		VecScale(negNorm, aaNorm[vrt], -1.0);
		AdaptSurfaceGridToCylinder(sel, grid, vrt, negNorm, 0.04, 0.02);

		for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
		{
			Face* f = *fIter;
			sh.assign_subset(f, 3);
			for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(f); eIter != grid.associated_edges_end(f); ++eIter)
			{
				EdgeBase* e = *eIter;
				sh.assign_subset(e, 3);
			}

			for(size_t i = 0; i < f->num_vertices(); ++i)
			{
				VertexBase* vrt = f->vertex(i);
				sh.assign_subset(vrt, 3);
			}
		}

		Refine(grid, sel);
	}

	sh.set_subset_name("CChCl", 3);


	/*
	 * T-bars_bnd
	 */
	//number TBarHeight = 0.06;
	//number TBarHeight = 0.02;
	number TBarHeight = 0.015;
	number TableLegRadius = 0.02;
	number TableTopRadius = 5 * TableLegRadius;
	number TableTopHeight = 0.01;

	vector<EdgeBase*> vExtrusionEdges;
	vector<EdgeBase*> vTmpExtrusionEdges;
	vector3 vZero(0.0,0.0,0.0);

	vector3 vExtrDir;
	vector3 vUnitExtrDir;

	Selector tmpSel(grid);

	for(size_t i = 0; i < vrts.size(); ++i)
	{
		tmpSel.clear();
		vExtrusionEdges.clear();
		vTmpExtrusionEdges.clear();

	//	create cylinder surface
		VertexBase* vrt = vrts[i];
		vUnitExtrDir = aaNorm[vrt];
		VecScale(vExtrDir, vUnitExtrDir, -1.0 * TBarHeight);
		AdaptSurfaceGridToCylinder(sel, grid, vrt, vExtrDir, TableLegRadius, 0.01);

	//	assign closure of cylinder surface to subset "5 = T-bars_post"
		for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
		{
			Face* f = *fIter;
			sh.assign_subset(f, 5);

			for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(f); eIter != grid.associated_edges_end(f); ++eIter)
			{
				EdgeBase* e = *eIter;
				sh.assign_subset(e, 5);
			}

			for(size_t i = 0; i < f->num_vertices(); ++i)
			{
				VertexBase* vrt = f->vertex(i);
				sh.assign_subset(vrt, 5);
			}
		}

	//	select T-bars_bottom boundary edges and extrude them and assign extrusion to subset "4 = T-bars_bnd"
		SelectAreaBoundary(tmpSel, sel.begin<Face>(), sel.end<Face>());

		for(EdgeBaseIterator eIter = tmpSel.begin<EdgeBase>(); eIter != tmpSel.end<EdgeBase>(); ++eIter)
		{
			EdgeBase* e = *eIter;
			vExtrusionEdges.push_back(e);
			vTmpExtrusionEdges.push_back(e);

			sh.assign_subset(e, 4);
			sh.assign_subset(e->vertex(0), 4);
			sh.assign_subset(e->vertex(1), 4);
		}

		Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);
		Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);
		Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);
		Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);

	//	assign lower ring of T-bar_post to "3 = CChCl"
		for(size_t i = 0; i < vTmpExtrusionEdges.size(); i++)
		{
			EdgeBase* e = vTmpExtrusionEdges[i];
			sh.assign_subset(e, 3);
			sh.assign_subset(e->vertex(0), 3);
			sh.assign_subset(e->vertex(1), 3);
		}

	//	assign upper ring of T-bar posts to subset "6 = T-bars_bottom"
		vTmpExtrusionEdges.clear();
		for(size_t i = 0; i < vExtrusionEdges.size(); i++)
		{
			EdgeBase* e = vExtrusionEdges[i];
			sh.assign_subset(e, 6);

			vTmpExtrusionEdges.push_back(vExtrusionEdges[i]);
		}

	//	extrude and scale table top horizontally
		tmpSel.clear();
		Extrude(grid, NULL, &vExtrusionEdges, NULL, vZero, EO_CREATE_FACES);

		for(size_t i = 0; i < vExtrusionEdges.size(); i++)
		{
			EdgeBase* e = vExtrusionEdges[i];
			tmpSel.select(e->vertex(0));
			tmpSel.select(e->vertex(1));
			tmpSel.select(e);

			sh.assign_subset(e, 6);
			sh.assign_subset(e->vertex(0), 6);
			sh.assign_subset(e->vertex(1), 6);
		}

	//	seperate post from table bottom (reassign upper ring of T-bar table post to subset "4 = T-bars_bnd")
		for(size_t i = 0; i < vTmpExtrusionEdges.size(); i++)
		{
			EdgeBase* e = vTmpExtrusionEdges[i];
			sh.assign_subset(e, 4);
			sh.assign_subset(e->vertex(0), 4);
			sh.assign_subset(e->vertex(1), 4);
		}

		center = CalculateBarycenter(tmpSel.begin<VertexBase>(), tmpSel.end<VertexBase>(), aaPos);

		for(VertexBaseIterator vIter = tmpSel.begin<VertexBase>(); vIter != tmpSel.end<VertexBase>(); ++vIter)
		{
			VertexBase* vrt = *vIter;
			VecSubtract(vExtrDir, aaPos[vrt], center);
			VecNormalize(vExtrDir, vExtrDir);
			VecScale(vExtrDir, vExtrDir, TableTopRadius);
			VecAdd(aaPos[vrt], aaPos[vrt], vExtrDir);
		}

	//	select edges to extrude vertically and assign table top side edges to subset "7 = T-bars_sides"
		tmpSel.clear();
		vExtrusionEdges.clear();
		SelectAreaBoundary(tmpSel, sh.begin<Face>(6), sh.end<Face>(6));
		for(EdgeBaseIterator eIter = tmpSel.begin<EdgeBase>(); eIter != tmpSel.end<EdgeBase>(); ++eIter)
		{
			EdgeBase* e = *eIter;
			if(NumAssociatedFaces(grid, e) == 1)
			{
				vExtrusionEdges.push_back(e);
				sh.assign_subset(e, 7);
				sh.assign_subset(e->vertex(0), 7);
				sh.assign_subset(e->vertex(1), 7);
			}
		}

	//	extrude table top vertically and assign upper edge ring to subset "8 = T-bars_tabletop"
		VecScale(vExtrDir, vUnitExtrDir, -1.0 * TableTopHeight);
		Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);

		tmpSel.clear();
		for(size_t i = 0; i < vExtrusionEdges.size(); i++)
		{
			EdgeBase* e = vExtrusionEdges[i];
			tmpSel.select(e);
			sh.assign_subset(e, 8);
		}

		sh.set_default_subset_index(8);
		TriangleFill_SweepLine(grid, tmpSel.begin<EdgeBase>(), tmpSel.end<EdgeBase>(), aPosition, aInt, &sh, 8);

		sh.set_default_subset_index(-1);

	}

	sh.set_subset_name("T-bars_post", 5);
	sh.set_subset_name("T-bars_bottom", 6);
	sh.set_subset_name("T-bars_sides", 7);
	sh.set_subset_name("T-bars_tabletop", 8);



	/************************************
	 * Optimization of the triangulation
	 ************************************/


	/*
	 * table bottom
	 */
	tmpSel.clear();
	Triangulate(grid, sh.begin<Quadrilateral>(6), sh.end<Quadrilateral>(6));

	QualityGridGeneration(grid, sh.begin<Face>(6), sh.end<Face>(6), aaPos, 10.0);

	for(FaceIterator fIter = sh.begin<Face>(6); fIter != sh.end<Face>(6); ++fIter)
	{
		Face* f = *fIter;
		tmpSel.select(f);
	}
	sh.set_default_subset_index(6);
	Refine(grid, tmpSel);
	QualityGridGeneration(grid, sh.begin<Face>(6), sh.end<Face>(6), aaPos, 20.0);


	/*
	 * table top
	 */
	tmpSel.clear();
	for(FaceIterator fIter = sh.begin<Face>(8); fIter != sh.end<Face>(8); ++fIter)
	{
		Face* f = *fIter;
		tmpSel.select(f);
	}
	sh.set_default_subset_index(8);
	Refine(grid, tmpSel);
	QualityGridGeneration(grid, sh.begin<Face>(8), sh.end<Face>(8), aaPos, 30.0);


	/*
	 * table sides
	 */
	sel.clear();
	sel.enable_autoselection(true);
	tmpSel.clear();

	ABool aBoolMarked = false;
	grid.attach_to_faces(aBoolMarked);
	Grid::FaceAttachmentAccessor<ABool> aaBoolMarked(grid, aBoolMarked);

	int counter = 1;

	sh.set_default_subset_index(7);
	Triangulate(grid, sh.begin<Quadrilateral>(7), sh.end<Quadrilateral>(7));

	for(FaceIterator fIter = sh.begin<Face>(7); fIter != sh.end<Face>(7); ++fIter)
	{
		Face* f = *fIter;
		sel.select(f);
	}

	for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
	{
		Face* f = *fIter;

		if(aaBoolMarked[f] == false)
		{
			tmpSel.select(f);
			SelectLinkedFlatFaces(tmpSel, 1.0);

			for(FaceIterator fMarkedIter = tmpSel.begin<Face>(); fMarkedIter != tmpSel.end<Face>(); ++fMarkedIter)
			{
				aaBoolMarked[*fMarkedIter] = true;
				sh.assign_subset(*fMarkedIter, counter + 8);
			}

			tmpSel.clear();
			counter++;
		}
	}

	for(int i = 0; i < counter; i++)
		QualityGridGeneration(grid, sh.begin<Face>(i + 8), sh.end<Face>(i + 8), aaPos, 10.0);


	/*
	 * table post
	 */
	Triangulate(grid, sh.begin<Quadrilateral>(4), sh.end<Quadrilateral>(4));


	/*
	 * RR
	 */
	QualityGridGeneration(grid, sh.begin<Face>(1), sh.end<Face>(1), aaPos, 30.0);
	QualityGridGeneration(grid, sh.begin<Face>(0), sh.end<Face>(0), aaPos, 25.0);


	/*
	 * Join redundant t-bar subsets
	 */
	while(sh.num_subsets() > 5)
	{
		sh.join_subsets(4, 4, 5, true);
	}



	/************************************
	 * Volume grid generation
	 ***********************************/
	//Tetrahedralize(grid, sh, 18.0, true, true);

	/*
	for(size_t i = 0; i < sh.num_subsets(); i++)
	{
		sel.clear();

		for(FaceIterator fIter = sh.begin<Face>(i); fIter != sh.end<Face>(i); ++fIter)
		{
			Face* f = *fIter;
			sel.select(f);
			sh.assign_subset(f, i);
		}

		SelectAssociatedGeometricObjects(sel);

		for(VertexBaseIterator vIter = sel.begin<VertexBase>(); vIter != sel.end<VertexBase>(); ++vIter)
		{
			VertexBase* vrt = *vIter;
			sh.assign_subset(vrt, i);
		}

		for(EdgeBaseIterator eIter = sel.begin<EdgeBase>(); eIter != sel.end<EdgeBase>(); ++eIter)
		{
			EdgeBase* e = *eIter;
			sh.assign_subset(e, i);
		}
	}
	*/

	//SeparateSubsetsByLowerDimSubsets<Volume>(grid, sh, true);
	//AdjustSubsetsForSimulation(sh, true);

	AssignSubsetColors(sh);

//	VertexBase* vrt;
//	Grid::edge_traits::secure_container edges;
//	grid.associated_elements(edges, vrt);

	stringstream ss;
	ss << "bouton" << numReleaseSites << ".ugx";
	string outfile = ss.str();

	SaveGridToUGX(grid, sh, outfile.c_str());
}




void SaveSelectionStatesToFile(Grid& mg, Selector& msel, const char* filename)
{
//	create a subset handler which holds different subsets for the different selection states
	//MultiGrid& mg = *msel.multi_grid();
	SubsetHandler sh(mg);


	for(Selector::traits<Volume>::iterator iter = msel.begin<Volume>();
		iter != msel.end<Volume>(); ++iter)
	{
		sh.assign_subset(*iter, msel.get_selection_status(*iter));
	}

	for(Selector::traits<Face>::iterator iter = msel.begin<Face>();
		iter != msel.end<Face>(); ++iter)
	{
		sh.assign_subset(*iter, msel.get_selection_status(*iter));
	}

	for(Selector::traits<EdgeBase>::iterator iter = msel.begin<EdgeBase>();
		iter != msel.end<EdgeBase>(); ++iter)
	{
		sh.assign_subset(*iter, msel.get_selection_status(*iter));
	}

	for(Selector::traits<VertexBase>::iterator iter = msel.begin<VertexBase>();
		iter != msel.end<VertexBase>(); ++iter)
	{
		sh.assign_subset(*iter, msel.get_selection_status(*iter));
	}


	const char* subsetNames[] = {"one", "two"};

	for(int i = 0; i < 2; ++i)
		sh.subset_info(i).name = subsetNames[i];

	AssignSubsetColors(sh);
	EraseEmptySubsets(sh);
	//SaveGridHierarchyTransformed(mg, sh, filename, LG_DISTRIBUTION_Z_OUTPUT_TRANSFORM);

	SaveGridToUGX(mg, sh, filename);
}





}




