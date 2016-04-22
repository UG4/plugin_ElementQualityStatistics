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
//	CreateElementQualityHistogram
void CreateElementQualityHistogram(vector<int>& histOut, const std::vector<number>& vQualities, int numSections)
{
//	Determine min and max quality values and range
	number minVal = numeric_limits<number>::max();
	number maxVal = numeric_limits<number>::min();

	for(size_t i = 0; i < vQualities.size(); ++i)
	{
		minVal = min(minVal, vQualities[i]);
		maxVal = max(maxVal, vQualities[i]);
	}

	number range;
	range = maxVal - minVal;

	if(range <= 0)
	{
		histOut.resize(vQualities.size(), 0);
		return;
	}

	histOut.resize(vQualities.size());

//	Calculate correspondence of the i-th quality value to its adequate histogram section
	for(size_t i = 0; i < vQualities.size(); ++i)
	{
		int section = (int)((number)numSections * (vQualities[i] - minVal) / range);

		if(section < 0)
			section = 0;

		if(section >= numSections)
			section = numSections - 1;

		histOut[i] = section;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	AssignSubsetsByElementQuality3d
void AssignSubsetsByElementQuality3d(MultiGrid& grid, MGSubsetHandler& sh, int numSecs)
{
//	Position attachment management
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	Management for VolumeAttachment aElemID
	AInt aElemID;
	grid.attach_to_volumes(aElemID);
	Grid::VolumeAttachmentAccessor<AInt> aaElemID(grid, aElemID);

//	Init vector for element qualities to be stored
	std::vector<number> vQualities;

	int numElems = 0;

//	Calculate aspect ratio for every element and attach element id
	for(VolumeIterator vIter = grid.begin<Tetrahedron>(); vIter != grid.end<Tetrahedron>(); ++vIter)
	{
		Tetrahedron* tet = static_cast<Tetrahedron*>(*vIter);
		number q = CalculateAspectRatio(grid, tet, aaPos);
		vQualities.push_back(q);

		aaElemID[tet] = numElems;
		numElems++;
	}

//	Create element quality histogram by specified number of sections
	std::vector<int> hist;
	CreateElementQualityHistogram(hist, vQualities, numSecs);

//	Assign elements to subsets by their aspect ration and set subset name by quality
	for(VolumeIterator vIter = grid.begin<Tetrahedron>(); vIter != grid.end<Tetrahedron>(); ++vIter)
	{
		Tetrahedron* tet = static_cast<Tetrahedron*>(*vIter);
		sh.assign_subset(tet, hist[aaElemID[tet]]);

		const char* siName = (std::to_string(vQualities[aaElemID[tet]])).c_str();
		sh.set_subset_name(siName, hist[aaElemID[tet]]);
	}

//	Set color range for the assigned subsets and name subsets by their qualities
	for(int i = 0; i < sh.num_subsets(); ++i)
	{
		number ia = 0;

		if(sh.num_subsets() > 1)
			ia = (number)i/ (number)(sh.num_subsets()-1);

		number r = max<number>(0, -1 + 2 * ia);
		number g = 1. - fabs(2 * (ia - 0.5));
		number b = max<number>(0, 1. - 2 * ia);

		SubsetInfo& si = sh.subset_info(i);

		si.color.x() = r;
		si.color.y() = g;
		si.color.z() = b;
		si.color.w() = 1.f;
	}

//	Detach VolumeAttachment from grid
	grid.detach_from_volumes(aElemID);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	AssignSubsetsByElementQuality
void AssignSubsetsByElementQuality(MultiGrid& mg, MGSubsetHandler& sh, int dim, int numSecs)
{
	if(dim == 2){}
	else if(dim == 3)
		AssignSubsetsByElementQuality3d(mg, sh, numSecs);
	else
		UG_THROW("Only dimensions 2 or 3 supported.");
}


////////////////////////////////////////////////////////////////////////////////////////////
//	AssignSubsetsByElementQuality
void AssignSubsetsByElementQuality(Grid& grid, SubsetHandler& sh, int dim, int numSecs)
{
	if(dim == 2){}
	else if(dim == 3){}
	else
		UG_THROW("Only dimensions 2 or 3 supported.");
}


////////////////////////////////////////////////////////////////////////////////////////////
//	ElementQualityStatistics
////////////////////////////////////////////////////////////////////////////////////////////

//	Wrapper for multigrids
//void ElementQualityStatistics(MultiGrid& mg, int level)
void ElementQualityStatistics(MultiGrid& mg, int dim)
{
	if(dim == 2)
		ElementQualityStatistics2d(mg, mg.get_grid_objects());
	else if(dim == 3)
		ElementQualityStatistics3d(mg, mg.get_grid_objects());
	else
		UG_THROW("Only dimensions 2 or 3 supported.");
}

//	Wrapper for grids
void ElementQualityStatistics(Grid& grid, int dim)
{
	if(dim == 2)
		ElementQualityStatistics2d(grid, grid.get_grid_objects());
	else if(dim == 3)
		ElementQualityStatistics3d(grid, grid.get_grid_objects());
	else
		UG_THROW("Only dimensions 2 or 3 supported.");
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
		MaxAngleHistogram(grid, goc.begin<Face>(i), goc.end<Face>(i), aaPos, 10);
		//PROFILE_END();

		UG_LOG(endl);

		PrintAngleStatistics2d(grid, goc, i, aaPos);
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
	number n_minFace;
	number n_maxFace;
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
	Face* minFace = NULL;
	Face* maxFace = NULL;
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
		if(goc.num_volumes(i) == 0)
		{
			minFace = FindSmallestFace(goc.begin<Face>(i), goc.end<Face>(i), aaPos);
			maxFace = FindLargestFace(goc.begin<Face>(i), goc.end<Face>(i), aaPos);
		}
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
				n_minVolume = CalculateVolume(minVolume, aaPos);
			if(maxVolume != NULL)
				n_maxVolume = CalculateVolume(maxVolume, aaPos);
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

		if(goc.num<Triangle>(i) > 0)
		{
			table(6, 0) << "Smallest triangle AR"; table(6, 1) << n_minTriAspectRatio;
			table(6, 2) << "Largest triangle AR"; table(6, 3) << n_maxTriAspectRatio;
		}

		if(goc.num_volumes(i) == 0)
		{
			table(7, 0) << "Smallest face";	table(7, 1) << n_minFace;
			table(7, 2) << "Largest face";	table(7, 3) << n_maxFace;
		}

		if(goc.num_volumes(i) > 0)
		{
			table(8, 0) << "Smallest volume";		table(8, 1) << n_minVolume;
			table(8, 2) << "Largest volume";		table(8, 3) << n_maxVolume;
			table(9, 0) << "Smallest volume angle";	table(9, 1) << n_minVolAngle;
			table(9, 2) << "Largest volume angle";	table(9, 3) << n_maxVolAngle;
			table(10, 0) << "Smallest volume dihedral";	table(10, 1) << n_minVolDihedral;
			table(10, 2) << "Largest volume dihedral";	table(10, 3) << n_maxVolDihedral;

			if(goc.num<Tetrahedron>(i) > 0)
			{
				table(11, 0) << "Smallest tetrahedron AR";	table(11, 1) << n_minTetAspectRatio;
				table(11, 2) << "Largest tetrahedron AR";	table(11, 3) << n_maxTetAspectRatio;
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
			MaxAngleHistogram(grid, goc.begin<Volume>(i), goc.end<Volume>(i), aaPos, 10);
		}
		else
		{
			MinAngleHistogram(grid, goc.begin<Face>(i), goc.end<Face>(i), aaPos, 10);
			MaxAngleHistogram(grid, goc.begin<Face>(i), goc.end<Face>(i), aaPos, 10);
		}
		//PROFILE_END();

		UG_LOG(endl);

		PrintAngleStatistics3d(grid, goc, i, aaPos);
	}

	UG_LOG(endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl);
}




}




