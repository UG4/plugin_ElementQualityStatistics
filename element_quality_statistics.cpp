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


	//	Calculate and output standard deviation for tetrahedral/hexahedral angles
		if(goc.num_volumes(i) > 0)
		{
			number sd_tet = 0.0;
			number sd_hex = 0.0;
			number mean_tet = 0.0;
			number mean_hex = 0.0;
			int numTets = 0;
			int numHex 	= 0;
			number regVolDihedral = 90.0;

			vector<number> vDihedralsOut;

			for(VolumeIterator vIter = goc.begin<Volume>(i); vIter != goc.end<Volume>(i); ++vIter)
			{
				vDihedralsOut.clear();
				Volume* vol = *vIter;

				if(vol->reference_object_id() == ROID_TETRAHEDRON)
				{
					numTets += 1;
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
					numHex += 1;
					regVolDihedral = 90.0;
					CalculateVolumeDihedrals(vDihedralsOut, grid, vol, aaPos);

					for(size_t k = 0; k < vDihedralsOut.size(); ++k)
					{
						sd_hex += (regVolDihedral-vDihedralsOut[k])*(regVolDihedral-vDihedralsOut[k]);
						mean_hex += vDihedralsOut[k];
					}
				}
			}

			if(numTets > 0)
			{
				sd_tet *= (1.0/(6*numTets));
				sd_tet = sqrt(sd_tet);
				mean_tet *= (1.0/(6*numTets));
			}

			if(numHex > 0)
			{
				sd_hex *= (1.0/(12*numHex));
				sd_hex = sqrt(sd_hex);
				mean_hex = (1.0/(12*numHex));
			}

			UG_LOG("Standard deviation of dihedral angles to regular case" << endl);
			UG_LOG("(70.5288° for tetrahedrons, 90° for hexahedrons)" << endl);
			UG_LOG(endl);
			UG_LOG("	Tetrahedrons (" << numTets << "):" << endl);
			if(numTets > 0)
			{
				UG_LOG("		sd   = " << sd_tet << endl);
				UG_LOG("		mean = " << mean_tet << endl);
			}
			UG_LOG(endl);
			UG_LOG("	Hexahedrons (" << numHex << "):" << endl);
			if(numHex > 0)
			{
				UG_LOG("		sd   = " << sd_hex << endl);
				UG_LOG("		mean = " << mean_tet << endl);
			}
			UG_LOG(endl);
		}
	}

	UG_LOG(endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl);
}




}




