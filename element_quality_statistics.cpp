/*
 * element_quality_statistics.cpp
 *
 *  Created on: 17.04.2012
 *      Author: Martin Stepniewski
 */


#include "element_quality_statistics.h"
#include "common/util/table.h"
#include "lib_grid/grid_objects/tetrahedron_rules.h"


namespace ug
{

////////////////////////////////////////////////////////////////////////////////////////////
//	MoveVertexToSmoothTetGridSubdivisionPosition
void MoveVertexToSmoothTetGridSubdivisionPosition(MultiGrid& mg, Vertex* vrt, 	Grid::VertexAttachmentAccessor<APosition>& aaPos,
																				Grid::VertexAttachmentAccessor<APosition>& aaSmoothPos)
{
//	Declare centroid coordinate vector
	typedef typename APosition::ValueType pos_type;
	pos_type p;

//	Declare vertex volume valence
	size_t valence = 0;

//	Collect associated volumes
	std::vector<Volume*> volumes;
	CollectAssociated(volumes, mg, vrt);

//	Iterate over all associated volumes
	for(Grid::AssociatedVolumeIterator vIter = mg.associated_volumes_begin(vrt); vIter != mg.associated_volumes_end(vrt); ++vIter)
	{
		VecSet(p, 0);
		Volume* vol = *vIter;
		++valence;

	//	TETRAHEDRON CASE
		if(vol->reference_object_id() == ROID_TETRAHEDRON)
		{
		//	Iterate over all associated vertices inside tetrahedron

			//
			// Alternative iteration
			//
			/*
			for(Grid::AssociatedEdgeIterator eIter = mg.associated_edges_begin(vol); eIter != mg.associated_edges_end(vol); ++eIter)
			{
				Edge* e = *eIter;
				if(GetConnectedVertex(e, vrt) != NULL)
					VecAdd(p, p, aaPos[GetConnectedVertex(e, vrt)]);
			}
			*/

			for(size_t i = 0; i < vol->num_vertices(); ++i)
			{
				if(i != GetVertexIndex(vol, vrt))
				{
					VecAdd(p, p, aaPos[vol->vertex(i)]);
				}
			}

		//	TODO: refer to subdivision rules object
			number centerWgt 	= -1.0/16;
			number nbrWgt 		= 17.0/48;

			VecScaleAppend(aaSmoothPos[vrt], centerWgt, aaPos[vrt], nbrWgt, p);
		}

	//	OCTAHEDRON CASE
		else if(vol->reference_object_id() == ROID_OCTAHEDRON)
		{
		//	Iterate over all vertices inside octahedron, first associate ones and last the opposing one
			for(Grid::AssociatedEdgeIterator eIter = mg.associated_edges_begin(vol); eIter != mg.associated_edges_end(vol); ++eIter)
			{
				Edge* e = *eIter;
				if(GetConnectedVertex(e, vrt) != NULL)
				{
					VecAdd(p, p, aaPos[GetConnectedVertex(e, vrt)]);
				}
			}

			Vertex* oppVrt = vol->vertex(vol->get_opposing_object(vrt).second);

		//	TODO: refer to subdivision rules object
			number centerWgt 	= 3.0/8;
			number nbrWgt 		= 1.0/12;
			number oppNbrWgt 	= 7.0/24;

			VecScaleAppend(aaSmoothPos[vrt], centerWgt, aaPos[vrt], nbrWgt, p, oppNbrWgt, aaPos[oppVrt]);
		}

	//	UNSUPPORTED VOLUME ELEMENT CASE
		else
		{
			UG_THROW("Volume type not supported for subdivision volumes refinement.");
		}
	}

//	Scale vertex position by the number of associated volume elements
	VecScale(aaSmoothPos[vrt],  aaSmoothPos[vrt], 1.0/valence);
}

////////////////////////////////////////////////////////////////////////////////////////////
//	SubdivisionTetGridSmooth
void SubdivisionTetGridSmooth(MultiGrid& mg, MGSubsetHandler& sh)
{
	typedef typename APosition::ValueType pos_type;

//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);

	APosition aSmoothPosition(0.0);
	mg.attach_to_vertices(aSmoothPosition);
	Grid::VertexAttachmentAccessor<APosition> aaSmoothPos(mg, aSmoothPosition);

//	Load subdivision surfaces rules
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();

//	Check, if volumes are included in input grid
	bool volumesExist = mg.num<Volume>() > 0;
	if(!volumesExist)
		UG_THROW("SubdivisionTetGridSmooth: No volumes included in input grid for smooth TetGrid subdivision refinement.");

//	Initialize new SmoothPosition of vertices with 0.0
	/*
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

		aaSmoothPos[vrt].x() = 0.0;
		aaSmoothPos[vrt].y() = 0.0;
		aaSmoothPos[vrt].z() = 0.0;

		//UG_LOG(aaSmoothPos[vrt].x() << ", " << aaSmoothPos[vrt].y() << ", " << aaSmoothPos[vrt].z() << endl);
	}
	*/

// 	Loop all vertices of top_level
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

	//	Even vertex
		if(mg.get_parent(vrt)->reference_object_id() == ROID_VERTEX)
		{
			Vertex* parentVrt = dynamic_cast<Vertex*>(mg.get_parent(vrt));

		//	Boundary vertex
			if(IsBoundaryVertex3D(mg, vrt))
			{
			//	perform loop subdivision on even surface vertices
			//	first get neighboured vertices
				size_t valence = 0;
				pos_type p;
				VecSet(p, 0);

				for(Grid::AssociatedEdgeIterator iter = mg.associated_edges_begin(parentVrt); iter != mg.associated_edges_end(parentVrt); ++iter)
				{
					if((!volumesExist) || IsBoundaryEdge3D(mg, *iter))
					{
						VecAdd(p, p, aaPos[GetConnectedVertex(*iter, parentVrt)]);
						++valence;
					}
				}

				number centerWgt 	= subdiv.ref_even_inner_center_weight(valence);
				number nbrWgt 		= subdiv.ref_even_inner_nbr_weight(valence);

				VecScaleAdd(aaSmoothPos[vrt], centerWgt, aaPos[parentVrt], nbrWgt, p);
			}

		//	Inner vertex
			else
				MoveVertexToSmoothTetGridSubdivisionPosition(mg, vrt, aaPos, aaSmoothPos);

		}

	//	Odd vertex
		if(mg.get_parent(vrt)->reference_object_id() == ROID_EDGE)
		{
		//	Get parent edge
			Edge* parentEdge = dynamic_cast<Edge*>(mg.get_parent(vrt));

		//	Boundary vertex
			if(IsBoundaryVertex3D(mg, vrt))
			{
			//	apply loop-subdivision on inner elements
			//	get the neighboured triangles
				Face* f[2];
				int numAssociatedBndFaces = 0;

				std::vector<Face*> faces;
				CollectAssociated(faces, mg, parentEdge);

				for(size_t i = 0; i < faces.size(); ++i)
				{
					if(IsBoundaryFace3D(mg, faces[i]))
					{
						if(numAssociatedBndFaces < 2)
						{
							f[numAssociatedBndFaces] = faces[i];
						}
						++numAssociatedBndFaces;
					}
				}

				if(numAssociatedBndFaces == 2)
				{
					if(f[0]->num_vertices() == 3 && f[1]->num_vertices() == 3)
					{
					//	the 4 vertices that are important for the calculation
						Vertex* v[4];
						v[0] = parentEdge->vertex(0); v[1] = parentEdge->vertex(1);
						v[2] = GetConnectedVertex(parentEdge, f[0]);
						v[3] = GetConnectedVertex(parentEdge, f[1]);

						vector4 wghts;

						wghts = subdiv.ref_odd_inner_weights();

						VecScaleAdd(aaSmoothPos[vrt],
									wghts.x(), aaPos[v[0]], wghts.y(), aaPos[v[1]],
									wghts.z(), aaPos[v[2]], wghts.w(), aaPos[v[3]]);
					}
					else
						UG_THROW("Non triangular faces included in grid.");
				}
				else
					UG_THROW("numAssociatedBndFaces != 2.");
			}
			else
				MoveVertexToSmoothTetGridSubdivisionPosition(mg, vrt, aaPos, aaSmoothPos);
		}

	//	Volume vertex
		if(mg.get_parent(vrt)->reference_object_id() == ROID_OCTAHEDRON)
			MoveVertexToSmoothTetGridSubdivisionPosition(mg, vrt, aaPos, aaSmoothPos);
	}

//	Move vertices to their smoothed position
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;
		VecScale(aaPos[vrt], aaSmoothPos[vrt], 1.0);
	}

}

////////////////////////////////////////////////////////////////////////////////////////////
//	SubdivisionTetGridSmoothBasic
void SubdivisionTetGridSmoothBasic(MultiGrid& mg, MGSubsetHandler& sh)
{
//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);

	APosition aTmpPosition;
	mg.attach_to_vertices(aTmpPosition);
	Grid::VertexAttachmentAccessor<APosition> aaTmpPos(mg, aTmpPosition);

//	Select all elements for linear refinement
	//sel.select(mg.begin<Vertex>(0), mg.end<Vertex>(0));
	//sel.select(mg.begin<Face>(0), mg.end<Face>(0));
	//sel.select(mg.begin<Volume>(0), mg.end<Volume>(0));

	//tet_rules::SetRefinementRule(tet_rules::HYBRID_TET_OCT);
	//Refine(mg, sel);

//	Initialize new TmpPosition of vertices with 0
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	//for(VertexIterator vrtIter = mg.begin<Vertex>(i); vrtIter != mg.end<Vertex>(i); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

		aaTmpPos[vrt].x() = 0.0;
		aaTmpPos[vrt].y() = 0.0;
		aaTmpPos[vrt].z() = 0.0;
	}

//	Smooth tetrahedral vertices (see Schaefer et al, "Smooth subdivision of tetrahedral meshes")
	for(VolumeIterator volIter = mg.begin<Tetrahedron>(mg.top_level()); volIter != mg.end<Tetrahedron>(mg.top_level()); ++volIter)
	//for(VolumeIterator volIter = mg.begin<Tetrahedron>(i); volIter != mg.end<Tetrahedron>(i); ++volIter)
	{
		Volume* vol = *volIter;

		VecScaleAppend(	aaTmpPos[vol->vertex(0)],
						-1.0/16, aaPos[vol->vertex(0)],
						17.0/48, aaPos[vol->vertex(1)],
						17.0/48, aaPos[vol->vertex(2)],
						17.0/48, aaPos[vol->vertex(3)]);

		VecScaleAppend(	aaTmpPos[vol->vertex(1)],
						-1.0/16, aaPos[vol->vertex(1)],
						17.0/48, aaPos[vol->vertex(0)],
						17.0/48, aaPos[vol->vertex(2)],
						17.0/48, aaPos[vol->vertex(3)]);

		VecScaleAppend(	aaTmpPos[vol->vertex(2)],
						-1.0/16, aaPos[vol->vertex(2)],
						17.0/48, aaPos[vol->vertex(0)],
						17.0/48, aaPos[vol->vertex(1)],
						17.0/48, aaPos[vol->vertex(3)]);

		VecScaleAppend(	aaTmpPos[vol->vertex(3)],
						-1.0/16, aaPos[vol->vertex(3)],
						17.0/48, aaPos[vol->vertex(0)],
						17.0/48, aaPos[vol->vertex(1)],
						17.0/48, aaPos[vol->vertex(2)]);
	}

//	Smooth octahedral vertices (see Schaefer et al, "Smooth subdivision of tetrahedral meshes")
	for(VolumeIterator volIter = mg.begin<Octahedron>(mg.top_level()); volIter != mg.end<Octahedron>(mg.top_level()); ++volIter)
	//for(VolumeIterator volIter = mg.begin<Octahedron>(i); volIter != mg.end<Octahedron>(i); ++volIter)
	{
		Volume* vol = *volIter;

		Vertex* vrt1 = vol->vertex(0);
		Vertex* vrt2 = vol->vertex(1);
		Vertex* vrt3 = vol->vertex(2);
		Vertex* vrt4 = vol->vertex(3);
		Vertex* vrt5 = vol->vertex(4);
		Vertex* vrt6 = vol->vertex(5);

	//	1
		VecScaleAppend(	aaTmpPos[vrt1],
						1.0/12, aaPos[vrt2],
						1.0/12, aaPos[vrt3],
						1.0/12, aaPos[vrt4],
						1.0/12, aaPos[vrt5]);

		VecScaleAppend(	aaTmpPos[vrt1],
						3.0/8,  aaPos[vrt1],
						7.0/24, aaPos[vrt6]);

	//	2
		VecScaleAppend(	aaTmpPos[vrt2],
						1.0/12, aaPos[vrt1],
						1.0/12, aaPos[vrt3],
						1.0/12, aaPos[vrt5],
						1.0/12, aaPos[vrt6]);

		VecScaleAppend(	aaTmpPos[vrt2],
						3.0/8,  aaPos[vrt2],
						7.0/24, aaPos[vrt4]);

	//	3
		VecScaleAppend(	aaTmpPos[vrt3],
						1.0/12, aaPos[vrt1],
						1.0/12, aaPos[vrt2],
						1.0/12, aaPos[vrt4],
						1.0/12, aaPos[vrt6]);

		VecScaleAppend(	aaTmpPos[vrt3],
						3.0/8,  aaPos[vrt3],
						7.0/24, aaPos[vrt5]);

	//	4
		VecScaleAppend(	aaTmpPos[vrt4],
						1.0/12, aaPos[vrt1],
						1.0/12, aaPos[vrt3],
						1.0/12, aaPos[vrt5],
						1.0/12, aaPos[vrt6]);

		VecScaleAppend(	aaTmpPos[vrt4],
						3.0/8,  aaPos[vrt4],
						7.0/24, aaPos[vrt2]);

	//	5
		VecScaleAppend(	aaTmpPos[vrt5],
						1.0/12, aaPos[vrt1],
						1.0/12, aaPos[vrt2],
						1.0/12, aaPos[vrt4],
						1.0/12, aaPos[vrt6]);

		VecScaleAppend(	aaTmpPos[vrt5],
						3.0/8,  aaPos[vrt5],
						7.0/24, aaPos[vrt3]);

	//	6
		VecScaleAppend(	aaTmpPos[vrt6],
						1.0/12, aaPos[vrt2],
						1.0/12, aaPos[vrt3],
						1.0/12, aaPos[vrt4],
						1.0/12, aaPos[vrt5]);

		VecScaleAppend(	aaTmpPos[vrt6],
						3.0/8,  aaPos[vrt6],
						7.0/24, aaPos[vrt1]);
	}

//	Calculate cell valencies of each vertex (:= of associated volumes)
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	//for(VertexIterator vrtIter = mg.begin<Vertex>(i); vrtIter != mg.end<Vertex>(i); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;
		int num_associated_volumes = 0;

	//	Calculate number of associated volumes
		for(Grid::AssociatedVolumeIterator vIter = mg.associated_volumes_begin(vrt); vIter != mg.associated_volumes_end(vrt); ++vIter)
			num_associated_volumes++;

		//UG_LOG(aaTmpPos[vrt].x() << "; " << aaTmpPos[vrt].y() << "; " << aaTmpPos[vrt].z() << "; " << num_associated_volumes << endl);

		VecScale(aaTmpPos[vrt], aaTmpPos[vrt], 1.0/num_associated_volumes);
	}

//	Move vertices to their smoothed position
	for(VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level()); vrtIter != mg.end<Vertex>(mg.top_level()); ++vrtIter)
	//for(VertexIterator vrtIter = mg.begin<Vertex>(i); vrtIter != mg.end<Vertex>(i); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;
		VecScale(aaPos[vrt], aaTmpPos[vrt], 1.0);
	}

//	Export grid
	//SaveGridToUGX(mg, sh, "test.ugx");
	//SaveGridToUGX(mg, sh, "test.ugx");
}



////////////////////////////////////////////////////////////////////////////////////////////
//	ElementQualityStatistics
////////////////////////////////////////////////////////////////////////////////////////////

//	Wrapper for multigrids
//void ElementQualityStatistics(MultiGrid& mg, int level)
void ElementQualityStatistics(MultiGrid& mg)
{
	if(mg.num_volumes() > 0)
		ElementQualityStatistics3d(mg, mg.get_grid_objects());
	else
		ElementQualityStatistics2d(mg, mg.get_grid_objects());
}

//	Wrapper for grids
void ElementQualityStatistics(Grid& grid)
{
	if(grid.num_volumes() > 0)
		ElementQualityStatistics3d(grid, grid.get_grid_objects());
	else
		ElementQualityStatistics2d(grid, grid.get_grid_objects());
}

//	Actual procedures
void ElementQualityStatistics2d(Grid& grid, GridObjectCollection goc)
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
	Edge* minEdge;
	Edge* maxEdge;
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
		minEdge = FindShortestEdge(goc.begin<Edge>(i), goc.end<Edge>(i), aaPos);
		maxEdge = FindLongestEdge(goc.begin<Edge>(i), goc.end<Edge>(i), aaPos);
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

void ElementQualityStatistics3d(Grid& grid, GridObjectCollection goc)
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
	Edge* minEdge;
	Edge* maxEdge;
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
		minEdge = FindShortestEdge(goc.begin<Edge>(i), goc.end<Edge>(i), aaPos);
		maxEdge = FindLongestEdge(goc.begin<Edge>(i), goc.end<Edge>(i), aaPos);
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




}




