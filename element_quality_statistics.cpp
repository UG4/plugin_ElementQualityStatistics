/*
 * Copyright (c) 2013-2019:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
 *
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


#include "common/util/table.h"
#include "common/util/stringify.h"
#include "element_quality_statistics.h"
#include "lib_grid/grid_objects/tetrahedron_rules.h"
#include "pcl/pcl_base.h"


namespace ug
{


////////////////////////////////////////////////////////////////////////////////////////////
//	PrintVertexVolumeValence
void PrintVertexVolumeValence(MultiGrid& mg, SubsetHandler& sh, int subsetIndex)
{
	AInt aNumElems;
	mg.attach_to_vertices_dv(aNumElems, 0);
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	Grid::VertexAttachmentAccessor<AInt> aaNumElems(mg, aNumElems);

	for(VolumeIterator vIter = mg.begin<Volume>(mg.top_level()); vIter != mg.end<Volume>(mg.top_level()); ++vIter)
	{
		Volume* vol = *vIter;

		for(size_t i = 0; i < vol->num_vertices(); ++i)
		{
			++aaNumElems[vol->vertex(i)];
		}
	}

	for(VertexIterator vIter = mg.begin<Vertex>(mg.top_level()); vIter != mg.end<Vertex>(mg.top_level()); ++vIter)
	{
		Vertex* vrt = *vIter;

		if(sh.get_subset_index(vrt) == subsetIndex)
			UG_LOG("Vertex " << aaPos[vrt] << ": aaNumElems[vrt] = " << aaNumElems[vrt] << std::endl);
	}
}


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
//	AssignSubsetToElementWithSmallestMinAngle
void AssignSubsetToElementWithSmallestMinAngle(MultiGrid& grid, MGSubsetHandler& sh, int dim, const char* roid, int si)
{
	if(dim == 2)
		AssignSubsetToElementWithSmallestMinAngle2d(grid, sh, roid, si);
	else if(dim == 3)
		AssignSubsetToElementWithSmallestMinAngle3d(grid, sh, roid, si);
	else
		UG_THROW("ERROR in AssignSubsetToElementWithSmallestMinAngle: Only dimensions 2 or 3 supported.");
}


////////////////////////////////////////////////////////////////////////////////////////////
//	AssignSubsetToElementWithSmallestMinAngle2d
void AssignSubsetToElementWithSmallestMinAngle2d(MultiGrid& grid, MGSubsetHandler& sh, const char* roid, int si)
{
	Selector sel(grid);
	GridObjectCollection goc = grid.get_grid_objects();
	Grid::VertexAttachmentAccessor<APosition2> aaPos(grid, aPosition2);
	uint i = goc.num_levels() - 1;

	if(strcmp(roid, "triangle") == 0 || strcmp(roid, "quadrilateral") == 0)
	{
		Face* minAngleElement;
		minAngleElement = FindElementWithSmallestMinAngle(	grid,
															goc.begin<Face>(i),
															goc.end<Face>(i),
															aaPos);

		sel.select(minAngleElement);
		CloseSelection(sel);

		sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), si);
		sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), si);
		sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), si);

		sh.set_subset_name("face_with_smallest_minAngle", si);
	}
	else
		UG_THROW("ERROR in AssignSubsetToElementWithSmallestMinAngle2d: only 'triangle' and 'quadrilateral' supported.");

	AssignSubsetColors(sh);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	AssignSubsetToElementWithSmallestMinAngle3d
void AssignSubsetToElementWithSmallestMinAngle3d(MultiGrid& grid, MGSubsetHandler& sh, const char* roid, int si)
{
	Selector sel(grid);
	GridObjectCollection goc = grid.get_grid_objects();
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	uint i = goc.num_levels() - 1;

	if(strcmp(roid, "triangle") == 0 || strcmp(roid, "quadrilateral") == 0)
	{
		Face* minAngleElement;
		minAngleElement = FindElementWithSmallestMinAngle(	grid,
															goc.begin<Face>(i),
															goc.end<Face>(i),
															aaPos);

		sel.select(minAngleElement);
		CloseSelection(sel);

		sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), si);
		sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), si);
		sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), si);

		sh.set_subset_name("face_with_smallest_minAngle", si);
	}
	else if(strcmp(roid, "tetrahedron") == 0)
	{
		Tetrahedron* minAngleElement;
		minAngleElement = FindElementWithSmallestMinAngle(	grid,
															goc.begin<Tetrahedron>(i),
															goc.end<Tetrahedron>(i),
															aaPos);

		sel.select(minAngleElement);
		CloseSelection(sel);

		sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), si);
		sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), si);
		sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), si);
		sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), si);

		sh.set_subset_name("tet_with_smallest_minAngle", si);
	}
	else
		UG_THROW("ERROR in AssignSubsetToElementWithSmallestMinAngle: only 'triangle', 'quadrilateral' and 'tetrahedron' supported.");

	AssignSubsetColors(sh);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	MeasureTetrahedronWithSmallestMinAngle
void MeasureTetrahedronWithSmallestMinAngle(MultiGrid& grid)
{
	Selector sel(grid);
	MGSubsetHandler sh(grid);
	GridObjectCollection goc = grid.get_grid_objects();
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	uint i = goc.num_levels() - 1;

	Tetrahedron* minAngleElement;
	minAngleElement = FindElementWithSmallestMinAngle(	grid,
														goc.begin<Tetrahedron>(i),
														goc.end<Tetrahedron>(i),
														aaPos);

	sel.select(minAngleElement);
	CloseSelection(sel);

//	consecutive subset assignment of volume sub-elements
	int subsetIndex = 0;

//	Edges
	for(EdgeIterator eIter = sel.begin<Edge>(); eIter != sel.end<Edge>(); ++eIter)
	{
		Edge* e = *eIter;
		sh.assign_subset(e, subsetIndex);
		subsetIndex++;
	}

//	Faces
	for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
	{
		Face* f = *fIter;
		sh.assign_subset(f, subsetIndex);
		subsetIndex++;
	}

//	Volume (assign vertices to same subset as volume element for valence calculation)
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), subsetIndex);
	sh.assign_subset(minAngleElement, subsetIndex);
	sh.set_subset_name("tet_with_smallest_minAngle", subsetIndex);

//	Taking measures
	number measure = 0.0;
	double faceAreaNormSquared_a = 0.0;
	double edgeLengthNormSq_b = 0.0;

	UG_LOG(std::endl);
	UG_LOG("-------------------------------------" << std::endl);
	UG_LOG("ElementWithSmallestMinAngle measures:" << std::endl);

	cout.precision(17);

	for(int i = 0; i <= 5; ++i)
	{
		measure = CalculateVolume(sh.begin<Edge>(i, grid.top_level()), sh.end<Edge>(i, grid.top_level()), aaPos);
		edgeLengthNormSq_b += measure*measure;
		UG_LOG("Edge " << i << " length : " << measure << std::endl);
	}

	for(int i = 6; i <= 9; ++i)
	{
		measure = CalculateVolume(sh.begin<Face>(i, grid.top_level()), sh.end<Face>(i, grid.top_level()), aaPos);
		faceAreaNormSquared_a += measure*measure;
		UG_LOG("Face " << i << " area : " << measure << std::endl);
	}

	measure = CalculateVolume(minAngleElement, aaPos);
	UG_LOG("Volume V = " << measure << std::endl);

	UG_LOG("Area norm squared a = " << faceAreaNormSquared_a << std::endl);
	UG_LOG("Length norm squared b = " << edgeLengthNormSq_b << std::endl);
	UG_LOG("-------------------------------------" << std::endl);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindBoundsForStiffnesMatrixMaxEigenvalue
void FindBoundsForStiffnesMatrixMaxEigenvalue(MultiGrid& mg, MGSubsetHandler& shOut)
{
	Selector sel(mg);

//	Vertex valence defitions
	AInt aNumElems;
	mg.attach_to_vertices_dv(aNumElems, 0);
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	Grid::VertexAttachmentAccessor<AInt> aaNumElems(mg, aNumElems);

//	Bound defining volume
	Volume* bdv = *mg.begin<Volume>(mg.top_level());

//	Measurement values
	number volume = 0.0;
	number faceArea = 0.0;
	number faceAreaNormSquared_a = 0.0;
	number edgeLength = 0.0;
	number edgeLengthNormSq_b = 0.0;
	number upperBound = 0.0;
	number upperBoundTmp = 0.0;
	int maxVertexValence = 0;

//	Calculate vertex valences
	for(VolumeIterator vIter = mg.begin<Volume>(mg.top_level()); vIter != mg.end<Volume>(mg.top_level()); ++vIter)
	{
		Volume* vol = *vIter;

		for(size_t i = 0; i < vol->num_vertices(); ++i)
		{
			++aaNumElems[vol->vertex(i)];
		}
	}

	sel.clear();

//	Determine bound defining volume
	for(VolumeIterator vIter = mg.begin<Volume>(mg.top_level()); vIter != mg.end<Volume>(mg.top_level()); ++vIter)
	{
		Volume* vol = *vIter;

		sel.select(vol);
		CloseSelection(sel);

	//	Volume
		volume = CalculateVolume(vol, aaPos);

	//	Faces
		for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
		{
			Face* f = *fIter;
			faceArea = CalculateVolume(f, aaPos);
			faceAreaNormSquared_a += faceArea*faceArea;
		}

		upperBoundTmp = faceAreaNormSquared_a/9.0/volume;

		if(upperBound < upperBoundTmp)
		{
			upperBound = upperBoundTmp;
			bdv = vol;
		}

		sel.clear();
		faceAreaNormSquared_a = 0.0;
	}

//	Complementary TetrahedronVolToRMSFaceAreaRatio calculation
	Tetrahedron* bdt = static_cast<Tetrahedron*>(bdv);
	number volToRMSFaceAreaRatio = CalculateTetrahedronVolToRMSFaceAreaRatio(mg, bdt, aaPos);

//	Identification of element with smallest TetrahedronVolToRMSFaceAreaRatio for debugging purposes
//	GridObjectCollection goc = mg.get_grid_objects();
//	Tetrahedron* bdt_eqs = FindElementWithSmallestVolToRMSFaceAreaRatio(mg,
//															goc.begin<Tetrahedron>(mg.top_level()),
//															goc.end<Tetrahedron>(mg.top_level()),
//															aaPos);
//	number volToRMSFaceAreaRatio_eqs = CalculateTetrahedronVolToRMSFaceAreaRatio(mg, bdt_eqs, aaPos);

	//sel.select(bdt_eqs);
	sel.select(bdv);
	CloseSelection(sel);

//	Take measurements of bdv

//	Volume
	//volume = CalculateVolume(bdt_eqs, aaPos);
	volume = CalculateVolume(bdv, aaPos);

//	Faces
	for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
	{
		Face* f = *fIter;
		faceArea = CalculateVolume(f, aaPos);
		faceAreaNormSquared_a += faceArea*faceArea;
	}

//	Edges
	for(EdgeIterator eIter = sel.begin<Edge>(); eIter != sel.end<Edge>(); ++eIter)
	{
		Edge* e = *eIter;
		edgeLength = CalculateVolume(e, aaPos);
		edgeLengthNormSq_b += edgeLength*edgeLength;
	}

//	maximal element vertex valence
	for(VertexIterator vrtIter = sel.begin<Vertex>(); vrtIter != sel.end<Vertex>(); ++vrtIter)
	{
		Vertex* vrt = *vrtIter;

		if(maxVertexValence < aaNumElems[vrt])
		{
			maxVertexValence = aaNumElems[vrt];
		}
	}

	cout.precision(17);

	number volToMeanSquareFaceAreaRatio = 4.0/(upperBound*9.0);

	UG_LOG(std::endl);
	UG_LOG("-----------------------------------------------------------------------------" << std::endl);
	UG_LOG("Bound defining tetrahedron measures (for stiffness matrix maximal eigenvalue):" << std::endl);
	UG_LOG(std::endl);
	UG_LOG("volToMeanSquareFaceAreaRatio = " << volToMeanSquareFaceAreaRatio << std::endl);
	UG_LOG("volToRMSFaceAreaRatio        = " << volToRMSFaceAreaRatio << std::endl);
//	UG_LOG("volToRMSFaceAreaRatio_eqs    = " << volToRMSFaceAreaRatio_eqs << std::endl);
	UG_LOG("Lower local bound            = " << upperBound/3.0 << std::endl);
	UG_LOG("Upper local bound            = " << upperBound << std::endl);
	UG_LOG("Lower global bound           = " << upperBound/3.0 << std::endl);
	UG_LOG("Upper global bound           = " << upperBound*maxVertexValence << std::endl);
	UG_LOG("Max valence                  = " << maxVertexValence << std::endl);
	UG_LOG("Volume V                     = " << volume << std::endl);
	UG_LOG("Area norm squared a          = " << faceAreaNormSquared_a << std::endl);
	UG_LOG("Length norm squared b        = " << edgeLengthNormSq_b << std::endl);
	UG_LOG("-----------------------------------------------------------------------------" << std::endl);

	shOut.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), 0);
	shOut.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), 0);
	shOut.assign_subset(sel.begin<Face>(), sel.end<Face>(), 0);
	shOut.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), 0);

	AssignSubsetColors(shOut);
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
		number q = CalculateMinDihedral(grid, tet, aaPos);
		//ALTERNATIVELY:
		//number q = CalculateAspectRatio(grid, tet, aaPos);
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

		const char* siName = (mkstr(vQualities[aaElemID[tet]])).c_str();
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
void ElementQualityStatistics(MultiGrid& mg, int dim, number angleHistStepSize, number aspectRatioHistStepSize, bool bWriteHistograms)
{
	if(dim == 2)
		ElementQualityStatistics2d(mg, mg.get_grid_objects(), angleHistStepSize, aspectRatioHistStepSize, bWriteHistograms);
	else if(dim == 3)
		ElementQualityStatistics3d(mg, mg.get_grid_objects(), angleHistStepSize, aspectRatioHistStepSize, bWriteHistograms);
	else
		UG_THROW("Only dimensions 2 or 3 supported.");
}

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
void ElementQualityStatistics(Grid& grid, int dim, number angleHistStepSize, number aspectRatioHistStepSize, bool bWriteHistograms)
{
	if(dim == 2)
		ElementQualityStatistics2d(grid, grid.get_grid_objects(), angleHistStepSize, aspectRatioHistStepSize, bWriteHistograms);
	else if(dim == 3)
		ElementQualityStatistics3d(grid, grid.get_grid_objects(), angleHistStepSize, aspectRatioHistStepSize, bWriteHistograms);
	else
		UG_THROW("Only dimensions 2 or 3 supported.");
}

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
void ElementQualityStatistics2d(Grid& grid, GridObjectCollection goc, number angleHistStepSize, number aspectRatioHistStepSize, bool bWriteHistograms)
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

	number n_minQuadAspectRatio = 0.0;
	number n_maxQuadAspectRatio = 0.0;


//	Elements
	Edge* minEdge;
	Edge* maxEdge;
	Face* minFace;
	Face* maxFace;
	Face* minAngleFace;
	Face* maxAngleFace;

	Triangle* minAspectRatioTri;
	Triangle* maxAspectRatioTri;

	Quadrilateral* minAspectRatioQuad;
	Quadrilateral* maxAspectRatioQuad;


//	Basic grid properties on level i
	UG_LOG(endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl);
	UG_LOG("GRID QUALITY STATISTICS" << endl << endl);
	UG_LOG("*** Output info:" << endl);
	UG_LOG("    - The 'aspect ratio' (AR) represents the ratio of minimal height and " << endl <<
		   "      maximal edge length of a triangle or tetrahedron respectively." << endl);
	UG_LOG("    - The Min- and MaxAngle-Histogram lists the number of min/max element angles in " << endl <<
		   "      different degree ranges (dihedrals for volumes!)." << endl << endl);
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
																	goc.begin<Triangle>(i),
																	goc.end<Triangle>(i),
																	aaPos);
			maxAspectRatioTri = FindElementWithLargestAspectRatio(	grid,
																	goc.begin<Triangle>(i),
																	goc.end<Triangle>(i),
																	aaPos);
			if(minAspectRatioTri != NULL)
				n_minTriAspectRatio = CalculateAspectRatio(grid, minAspectRatioTri, aaPos);
			if(maxAspectRatioTri != NULL)
				n_maxTriAspectRatio = CalculateAspectRatio(grid, maxAspectRatioTri, aaPos);
		}

	//	Check for quadrilaterals
		if (goc.num<Quadrilateral>(i) > 0)
		{
			minAspectRatioQuad = FindElementWithSmallestAspectRatio(grid,
																	goc.begin<Quadrilateral>(i),
																	goc.end<Quadrilateral>(i),
																	aaPos);
			maxAspectRatioQuad = FindElementWithLargestAspectRatio(grid,
																   goc.begin<Quadrilateral>(i),
																   goc.end<Quadrilateral>(i),
																   aaPos);
			if(minAspectRatioQuad != NULL)
				n_minQuadAspectRatio = CalculateAspectRatio(grid, minAspectRatioQuad, aaPos);
			if(maxAspectRatioQuad != NULL)
				n_maxQuadAspectRatio = CalculateAspectRatio(grid, maxAspectRatioQuad, aaPos);
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

		if(goc.num<Quadrilateral>(i) > 0)
		{
			table(8, 0) << "Smallest quad AR";	table(8, 1) << n_minQuadAspectRatio;
			table(8, 2) << "Largest quad AR";	table(8, 3) << n_maxQuadAspectRatio;
		}

	//	Output section
		UG_LOG("+++++++++++++++++" << endl);
		UG_LOG(" Grid level " << i << ":" << endl);
		UG_LOG("+++++++++++++++++" << endl << endl);
		UG_LOG(table);


		//PROFILE_BEGIN(eqs_minAngleHistogram);
	//	TODO:
		//MinAngleHistogram(grid, goc.begin<Face>(i), goc.end<Face>(i), aaPos, angleHistStepSize, i);
		//MaxAngleHistogram(grid, goc.begin<Face>(i), goc.end<Face>(i), aaPos, angleHistStepSize, i);
		//PROFILE_END();

		UG_LOG(endl);
		PrintAngleStatistics2d(grid, goc, i, aaPos);
	}

	UG_LOG(endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl);
}

void ElementQualityStatistics3d(Grid& grid, GridObjectCollection goc, number angleHistStepSize, number aspectRatioHistStepSize, bool bWriteHistograms)
{
	//PROFILE_FUNC();
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	DistributedGridManager* dgm = grid.distributed_grid_manager();

//	Output
	vector<number> volMinAngles;
	vector<number> volMaxAngles;
	vector<number> volAspectRatios;
	vector<number> volToRMSFaceAreaRatios;
	vector<number> faceMinAngles;
	vector<number> faceMaxAngles;

	//ug::Table<std::stringstream> volMinAngleTable(numRanges, 2);
	ug::Table<std::stringstream> volMinAngleTable;
	ug::Table<std::stringstream> volMaxAngleTable;
	ug::Table<std::stringstream> volAspectRatioTable;
	ug::Table<std::stringstream> volToRMSFaceAreaRatioTable;

//	Numbers
	size_t numVertices = 0;
	size_t numEdges = 0;
	size_t numFaces = 0;
	size_t numVolumes = 0;

	number n_minEdge = std::numeric_limits<double>::max();
	number n_maxEdge = std::numeric_limits<double>::min();
	number n_minFace = std::numeric_limits<double>::max();
	number n_maxFace = std::numeric_limits<double>::min();
	number n_minFaceAngle = std::numeric_limits<double>::max();
	number n_maxFaceAngle = std::numeric_limits<double>::min();

	number n_minTriAspectRatio = std::numeric_limits<double>::max();
	number n_maxTriAspectRatio = std::numeric_limits<double>::min();

	number n_minQuadAspectRatio = std::numeric_limits<double>::max();
	number n_maxQuadAspectRatio = std::numeric_limits<double>::min();

	number n_minVolume = std::numeric_limits<double>::max();
	number n_maxVolume = std::numeric_limits<double>::min();
	number n_minVolAngle = std::numeric_limits<double>::max();
	number n_maxVolAngle = std::numeric_limits<double>::min();

	number n_minTetAspectRatio = std::numeric_limits<double>::max();
	number n_maxTetAspectRatio = std::numeric_limits<double>::min();
	number n_minTetVolToRMSFaceAreaRatio = std::numeric_limits<double>::max();
	number n_maxTetVolToRMSFaceAreaRatio = std::numeric_limits<double>::min();

	number n_minHexAspectRatio = std::numeric_limits<double>::max();
	number n_maxHexAspectRatio = std::numeric_limits<double>::min();
	number n_minHexVolToRMSFaceAreaRatio = std::numeric_limits<double>::max();
	number n_maxHexVolToRMSFaceAreaRatio = std::numeric_limits<double>::min();


//	Elements
	Edge* minEdge;
	Edge* maxEdge;
	Face* minFace = NULL;
	Face* maxFace = NULL;
	Face* minAngleFace;
	Face* maxAngleFace;

	Triangle* minAngleTri;
	Triangle* maxAngleTri;
	Triangle* minAspectRatioTri;
	Triangle* maxAspectRatioTri;

	Quadrilateral* minAngleQuad;
	Quadrilateral* maxAngleQuad;
	Quadrilateral* minAspectRatioQuad;
	Quadrilateral* maxAspectRatioQuad;

	Volume* minVolume;
	Volume* maxVolume;
	Volume* minAngleVol;
	Volume* maxAngleVol;

	Tetrahedron* minAspectRatioTet;
	Tetrahedron* maxAspectRatioTet;
	Tetrahedron* minVolToRMSFaceAreaRatioTet;
	Tetrahedron* maxVolToRMSFaceAreaRatioTet;

	Hexahedron* minAspectRatioHex;
	Hexahedron* maxAspectRatioHex;


//	Basic grid properties on level i
	UG_LOG(endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl);
	UG_LOG("GRID QUALITY STATISTICS" << endl << endl);
	UG_LOG("*** Output info:" << endl);
	UG_LOG("    - The 'aspect ratio' (AR) represents the ratio of minimal height and " << endl <<
		   "      maximal edge length of a triangle or tetrahedron respectively." << endl);
	UG_LOG("    - The Min- and MaxAngle-Histogram lists the number of min/max element angles in " << endl <<
		   "      different degree ranges (dihedrals for volumes!)." << endl << endl);
	for(uint i = 0; i < goc.num_levels(); ++i)
	{
		//PROFILE_BEGIN(eqs_qualityStatistics2d);

	//	--------------------
	//	Number of elements
	//	--------------------
		numVertices = 0;
		numEdges = 0;
		numFaces = 0;
		numVolumes = 0;

		for(VertexIterator vrtIter = goc.begin<Vertex>(i); vrtIter != goc.end<Vertex>(i); ++vrtIter)
		{
			Vertex* vrt = *vrtIter;

			#ifdef UG_PARALLEL
			//	ghosts (vertical masters) as well as horizontal slaves (low dimensional elements only) have to be ignored,
			//	since they have a copy on another process and
			//	since we already consider that copy...
				if(dgm->is_ghost(vrt) || dgm->contains_status(vrt, ES_H_SLAVE))
					continue;
			#endif

			numVertices++;
		}

		for(EdgeIterator eIter = goc.begin<Edge>(i); eIter != goc.end<Edge>(i); ++eIter)
		{
			Edge* e = *eIter;

			#ifdef UG_PARALLEL
			//	ghosts (vertical masters) as well as horizontal slaves (low dimensional elements only) have to be ignored,
			//	since they have a copy on another process and
			//	since we already consider that copy...
				if(dgm->is_ghost(e) || dgm->contains_status(e, ES_H_SLAVE))
					continue;
			#endif

			numEdges++;
		}

		for(FaceIterator fIter = goc.begin<Face>(i); fIter != goc.end<Face>(i); ++fIter)
		{
			Face* f = *fIter;

			#ifdef UG_PARALLEL
			//	ghosts (vertical masters) as well as horizontal slaves (low dimensional elements only) have to be ignored,
			//	since they have a copy on another process and
			//	since we already consider that copy...
				if(dgm->is_ghost(f) || dgm->contains_status(f, ES_H_SLAVE))
					continue;
			#endif

			numFaces++;
		}

		for(VolumeIterator vIter = goc.begin<Volume>(i); vIter != goc.end<Volume>(i); ++vIter)
		{
			Volume* v = *vIter;

			#ifdef UG_PARALLEL
			//	ghosts (vertical masters) have to be ignored,
			//	since they have a copy on another process and
			//	since we already consider that copy...
				if(dgm->is_ghost(v))
					continue;
			#endif

			numVolumes++;
		}

		#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1){
			//	sum the numbers of all involved processes. Since we ignored ghosts,
			//	each process contributes the numbers of a unique part of the grid.
				pcl::ProcessCommunicator pc;
				numVertices = pc.allreduce(numVertices, PCL_RO_SUM);
				numEdges = pc.allreduce(numEdges, PCL_RO_SUM);
				numFaces = pc.allreduce(numFaces, PCL_RO_SUM);
				numVolumes = pc.allreduce(numVolumes, PCL_RO_SUM);
			}
		#endif

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
			minAspectRatioTri = FindElementWithSmallestAspectRatio(grid,
																	goc.begin<Triangle>(i),
																	goc.end<Triangle>(i),
																	aaPos);
			maxAspectRatioTri = FindElementWithLargestAspectRatio(grid,
																	goc.begin<Triangle>(i),
																	goc.end<Triangle>(i),
																	aaPos);
			if(minAspectRatioTri != NULL)
				n_minTriAspectRatio = CalculateAspectRatio(grid, minAspectRatioTri, aaPos);
			if(maxAspectRatioTri != NULL)
				n_maxTriAspectRatio = CalculateAspectRatio(grid, maxAspectRatioTri, aaPos);
		}

	//	Check for quadrilaterals
		if (goc.num<Quadrilateral>(i) > 0)
		{
			minAspectRatioQuad = FindElementWithSmallestAspectRatio(grid,
																	goc.begin<Quadrilateral>(i),
																	goc.end<Quadrilateral>(i),
																	aaPos);
			maxAspectRatioQuad = FindElementWithLargestAspectRatio(	grid,
																	goc.begin<Quadrilateral>(i),
																	goc.end<Quadrilateral>(i),
																	aaPos);
			if(minAspectRatioQuad != NULL)
				n_minQuadAspectRatio = CalculateAspectRatio(grid, minAspectRatioQuad, aaPos);
			if(maxAspectRatioQuad != NULL)
				n_maxQuadAspectRatio = CalculateAspectRatio(grid, maxAspectRatioQuad, aaPos);
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

			if(minVolume != NULL)
				n_minVolume = CalculateVolume(minVolume, aaPos);
			if(maxVolume != NULL)
				n_maxVolume = CalculateVolume(maxVolume, aaPos);
			if(minAngleVol != NULL)
				n_minVolAngle = CalculateMinAngle(grid, minAngleVol, aaPos);
			if(maxAngleVol != NULL)
				n_maxVolAngle = CalculateMaxAngle(grid, maxAngleVol, aaPos);

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

				minVolToRMSFaceAreaRatioTet = FindElementWithSmallestVolToRMSFaceAreaRatio(grid,
																		goc.begin<Tetrahedron>(i),
																		goc.end<Tetrahedron>(i),
																		aaPos);
				maxVolToRMSFaceAreaRatioTet = FindElementWithLargestVolToRMSFaceAreaRatio(grid,
																		goc.begin<Tetrahedron>(i),
																		goc.end<Tetrahedron>(i),
																		aaPos);



				if(minAspectRatioTet != NULL)
					n_minTetAspectRatio = CalculateAspectRatio(grid, minAspectRatioTet, aaPos);
				if(maxAspectRatioTet != NULL)
					n_maxTetAspectRatio = CalculateAspectRatio(grid, maxAspectRatioTet, aaPos);
				if(minVolToRMSFaceAreaRatioTet != NULL)
					n_minTetVolToRMSFaceAreaRatio = CalculateVolToRMSFaceAreaRatio(grid, minVolToRMSFaceAreaRatioTet, aaPos);
				if(maxVolToRMSFaceAreaRatioTet != NULL)
					n_maxTetVolToRMSFaceAreaRatio = CalculateVolToRMSFaceAreaRatio(grid, maxVolToRMSFaceAreaRatioTet, aaPos);
			}

		//	Hexahedron section
			if (goc.num<Hexahedron>(i) > 0) {
				minAspectRatioHex = FindElementWithSmallestAspectRatio(	grid,
																		goc.begin<Hexahedron>(i),
																		goc.end<Hexahedron>(i),
																		aaPos);
				maxAspectRatioHex = FindElementWithLargestAspectRatio(	grid,
																		goc.begin<Hexahedron>(i),
																		goc.end<Hexahedron>(i),
																		aaPos);

				if(minAspectRatioHex != NULL)
					n_minHexAspectRatio = CalculateAspectRatio(grid, minAspectRatioHex, aaPos);
				if(maxAspectRatioHex != NULL)
					n_maxHexAspectRatio = CalculateAspectRatio(grid, maxAspectRatioHex, aaPos);
			}
		}

		#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1){
			//	sum the numbers of all involved processes. Since we ignored ghosts,
			//	each process contributes the numbers of a unique part of the grid.
				pcl::ProcessCommunicator pc;
				n_minEdge = pc.allreduce(n_minEdge, PCL_RO_MIN);
				n_maxEdge = pc.allreduce(n_maxEdge, PCL_RO_MAX);
				n_minFace = pc.allreduce(n_minFace, PCL_RO_MIN);
				n_maxFace = pc.allreduce(n_maxFace, PCL_RO_MAX);
				n_minFaceAngle = pc.allreduce(n_minFaceAngle, PCL_RO_MIN);
				n_maxFaceAngle = pc.allreduce(n_maxFaceAngle, PCL_RO_MAX);
				n_minTriAspectRatio = pc.allreduce(n_minTriAspectRatio, PCL_RO_MIN);
				n_maxTriAspectRatio = pc.allreduce(n_maxTriAspectRatio, PCL_RO_MAX);
//				n_minQuadAngle = pc.allreduce(n_minQuadAngle, PCL_RO_MIN);
//				n_maxQuadAngle = pc.allreduce(n_maxQuadAngle, PCL_RO_MAX);
				n_minQuadAspectRatio = pc.allreduce(n_minQuadAspectRatio, PCL_RO_MIN);
				n_maxQuadAspectRatio = pc.allreduce(n_maxQuadAspectRatio, PCL_RO_MAX);
				n_minVolume = pc.allreduce(n_minVolume, PCL_RO_MIN);
				n_maxVolume = pc.allreduce(n_maxVolume, PCL_RO_MAX);
				n_minVolAngle = pc.allreduce(n_minVolAngle, PCL_RO_MIN);
				n_maxVolAngle = pc.allreduce(n_maxVolAngle, PCL_RO_MAX);
				n_minTetAspectRatio = pc.allreduce(n_minTetAspectRatio, PCL_RO_MIN);
				n_maxTetAspectRatio = pc.allreduce(n_maxTetAspectRatio, PCL_RO_MAX);
				n_minTetVolToRMSFaceAreaRatio = pc.allreduce(n_minTetVolToRMSFaceAreaRatio, PCL_RO_MIN);
				n_maxTetVolToRMSFaceAreaRatio = pc.allreduce(n_maxTetVolToRMSFaceAreaRatio, PCL_RO_MAX);
			}
		#endif

		//PROFILE_BEGIN(eqs_qualityStatisticsOutput);
	//	Table summary
		ug::Table<std::stringstream> table(11, 4);
		table(0, 0) << "Number of volumes"; 	table(0, 1) << numVolumes;
		table(1, 0) << "Number of faces"; 		table(1, 1) << numFaces;
		table(2, 0) << "Number of vertices";	table(2, 1) << numVertices;

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

		if (goc.num<Quadrilateral>(i) > 0)
		{
			table(7, 0) << "Smallest quadrilateral AR"; table(7, 1) << n_minQuadAspectRatio;
			table(7, 2) << "Largest quadrilateral AR"; table(7, 3) << n_maxQuadAspectRatio;
		}

		table(8, 0) << "Smallest face";	table(8, 1) << n_minFace;
		table(8, 2) << "Largest face";	table(8, 3) << n_maxFace;

		if(goc.num_volumes(i) > 0)
		{
			table(9, 0) << "Smallest volume";		table(9, 1) << n_minVolume;
			table(9, 2) << "Largest volume";		table(9, 3) << n_maxVolume;
			table(10, 0) << "Smallest volume dihedral";	table(10, 1) << n_minVolAngle;
			table(10, 2) << "Largest volume dihedral";	table(10, 3) << n_maxVolAngle;

			if(goc.num<Tetrahedron>(i) > 0)
			{
				table(11, 0) << "Smallest tet AR";	table(11, 1) << n_minTetAspectRatio;
				table(11, 2) << "Largest tet AR";	table(11, 3) << n_maxTetAspectRatio;
				table(12, 0) << "Smallest tet Vol/FaceAreaRatio";	table(12, 1) << n_minTetVolToRMSFaceAreaRatio;
				table(12, 2) << "Largest tet Vol/FaceAreaRatio";	table(12, 3) << n_maxTetVolToRMSFaceAreaRatio;
			}

			if (goc.num<Hexahedron>(i) > 0) {
				table(13, 0) << "Smallest hex AR";	table(13, 1) << n_minHexAspectRatio;
				table(13, 2) << "Largest hex AR";	table(13, 3) << n_maxHexAspectRatio;
				table(14, 0) << "Smallest hex Vol/FaceAreaRatio";	table(14, 1) << n_minHexVolToRMSFaceAreaRatio;
				table(14, 2) << "Largest hex Vol/FaceAreaRatio";	table(14, 3) << n_maxHexVolToRMSFaceAreaRatio;
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
			volMinAngles.clear();
			volMaxAngles.clear();
			volAspectRatios.clear();
			volToRMSFaceAreaRatios.clear();

			CollectMinAngles(grid, goc.begin<Volume>(i), goc.end<Volume>(i), aaPos, volMinAngles);
			CollectMaxAngles(grid, goc.begin<Volume>(i), goc.end<Volume>(i), aaPos, volMaxAngles);
			CollectAspectRatios(grid, goc.begin<Volume>(i), goc.end<Volume>(i), aaPos, volAspectRatios);
			CollectVolToRMSFaceAreaRatios(grid, goc.begin<Volume>(i), goc.end<Volume>(i), aaPos, volToRMSFaceAreaRatios);
		}
		else
		{
		//	TODO:
//			CollectMinAngles(grid, goc.begin<Face>(i), goc.end<Face>(i), aaPos, faceMinAngles);
//			CollectMaxAngles(grid, goc.begin<Face>(i), goc.end<Face>(i), aaPos, faceMaxAngles);
//			CollectAspectRatios(grid, goc.begin<Face>(i), goc.end<Face>(i), aaPos, faceAspectRatios);
		}
		//PROFILE_END();

		UG_LOG(endl << "(*) MinAngle-Histogram for '" << "3d' elements");
		UG_LOG(endl);
		volMinAngleTable.clear();
		PrintAngleHistogram(volMinAngles, angleHistStepSize, volMinAngleTable);

		UG_LOG(endl << "(*) MaxAngle-Histogram for '" << "3d' elements");
		UG_LOG(endl);
		volMaxAngleTable.clear();
		PrintAngleHistogram(volMaxAngles, angleHistStepSize, volMaxAngleTable);

		UG_LOG(endl << "(*) AspectRatio-Histogram for '" << "3d' elements");
		UG_LOG(endl);
		volAspectRatioTable.clear();
		PrintAspectRatioHistogram(volAspectRatios, aspectRatioHistStepSize, volAspectRatioTable);

		UG_LOG(endl << "(*) VolToRMSFaceAreaRatioHistogram-Histogram for '" << "3d' elements");
		UG_LOG(endl);
		volToRMSFaceAreaRatioTable.clear();
		PrintAspectRatioHistogram(volToRMSFaceAreaRatios, aspectRatioHistStepSize, volToRMSFaceAreaRatioTable);

	//	----------------------------------------
	//	Histogram table file output section
	//	----------------------------------------
		if(bWriteHistograms)
		{
			int procRank = 0;
			#ifdef UG_PARALLEL
				procRank = pcl::ProcRank();
			#endif
			if(procRank == 0)
			{
				ofstream ofstr;
				std::stringstream ss;
				ss << "volMinAngles_lvl_" << i << ".csv";
				ofstr.open(ss.str().c_str());
				ofstr << volMinAngleTable.to_csv(";");
				ofstr.close();

				ofstr.clear();
				ss.str("");
				ss << "volMaxAngles_lvl_" << i << ".csv";
				ofstr.open(ss.str().c_str());
				ofstr << volMaxAngleTable.to_csv(";");
				ofstr.close();

				ofstr.clear();
				ss.str("");
				ss << "volAspectRatios_lvl_" << i << ".csv";
				ofstr.open(ss.str().c_str());
				ofstr << volAspectRatioTable.to_csv(";");
				ofstr.close();

				ofstr.clear();
				ss.str("");
				ss << "volToRMSFaceAreaRatios_lvl_" << i << ".csv";
				ofstr.open(ss.str().c_str());
				ofstr << volToRMSFaceAreaRatioTable.to_csv(";");
				ofstr.close();
			}
		}

		UG_LOG(endl);
		PrintAngleStatistics3d(grid, goc, i, aaPos);
	}

	UG_LOG(endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	PrintHistograms
////////////////////////////////////////////////////////////////////////////////////////////
void PrintAngleHistogram(vector<number>& locAngles, number stepSize, ug::Table<std::stringstream>& outTable)
{
	vector<number> angles;

	#ifdef UG_PARALLEL
	//	propagate to all processes
		pcl::ProcessCommunicator pc;
		pc.allgatherv(angles, locAngles);
	#else
		angles = locAngles;
	#endif

//	Sort the calculated angles in an ascending way
	sort (angles.begin(), angles.end());

//	Evaluate the minimal and maximal degree rounding to 10
	int minDeg = round(number(angles.front()) / 10.0) * 10;
	int maxDeg = round(number(angles.back()) / 10.0) * 10;

//	Expand minDeg and maxDeg by plus minus 10 degrees or at least to 0 or 180 degrees
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
	for(uint i = 0; i < angles.size(); ++i)
	{
		number angle = angles[i];
		for (uint range = 0; range < numRanges; range++)
		{
			if (angle < minDeg + (range+1)*stepSize)
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
	ug::Table<std::stringstream> angleTable(numRows, 6);

//	First third
	uint i = 0;
	for(; i < numRows; ++i)
	{
		angleTable(i, 0) << minDeg + i*stepSize << " - " << minDeg + (i+1)*stepSize << " deg : ";
		angleTable(i, 1) << counter[i];
	}

//	Second third
//	Check, if second third of table is needed
	if(i < counter.size())
	{
		for(; i < 2*numRows; ++i)
		{
			angleTable(i-numRows, 2) << minDeg + i*stepSize << " - " << minDeg + (i+1)*stepSize << " deg : ";
			angleTable(i-numRows, 3) << counter[i];
		}
	}

//	Third third
	if(i < counter.size())
//	Check, if third third of table is needed
	{
		for(; i < numRanges; ++i)
		{
			angleTable(i-2*numRows, 4) << minDeg + i*stepSize << " - " << minDeg + (i+1)*stepSize << " deg : ";
			angleTable(i-2*numRows, 5) << counter[i];
		}
	}

//	Output table
	UG_LOG(endl << angleTable);


//	----------------------------------------
//	Histogram table file output section
//	----------------------------------------
	numRanges = floor(180.0/stepSize);
	counter.clear();
	counter.resize(numRanges, 0.0);
	int numElems = 0;

//	Count the elements in their corresponding minAngle range
	for(uint i = 0; i < angles.size(); ++i)
	{
		number angle = angles[i];
		for (uint range = 0; range < numRanges; range++)
		{
			if (angle < (range+1)*stepSize)
			{
				++counter[range];
				++numElems;
				break;
			}
		}
	}

	outTable.add_rows(numRanges);
	outTable.add_cols(2);

	for(uint i = 0; i < numRanges; ++i)
	{
		outTable(i, 0) << i*stepSize << " - " << (i+1)*stepSize;
		outTable(i, 1) << 100.0/numElems*counter[i];
	}
}


void PrintAspectRatioHistogram(vector<number>& locAspectRatios, number stepSize, ug::Table<std::stringstream>& outTable)
{
	vector<number> aspectRatios;

	#ifdef UG_PARALLEL
	//	propagate to all processes
		pcl::ProcessCommunicator pc;
		pc.allgatherv(aspectRatios, locAspectRatios);
	#else
		aspectRatios = locAspectRatios;
	#endif

//	Sort the calculated aspectRatios in an ascending way
	sort(aspectRatios.begin(), aspectRatios.end());

//	Evaluate the minimal and maximal aspectRatio rounding to 0.01
	number minAspectRatio = round(number(aspectRatios.front()) * 10.0) / 10.0;
	number maxAspectRatio = round(number(aspectRatios.back()) * 10.0) / 10.0;

//	Expand minAspectRatio and maxAspectRatio by plus minus 0.1 or at least to 0 or 1.0
	if((minAspectRatio-0.1) > 0)
		minAspectRatio = minAspectRatio - 0.1;
	else
		minAspectRatio = 0;

	if((maxAspectRatio+0.1) < 1.0)
		maxAspectRatio = maxAspectRatio + 0.1;
	else
		maxAspectRatio = 1.0;

//	Evaluate the number of ranges in respect to the specified step size
	uint numRanges = round((maxAspectRatio-minAspectRatio) / stepSize);
	vector<uint> counter(numRanges, 0);

//	Count the elements in their corresponding aspectRatio range
	for(uint i = 0; i < aspectRatios.size(); ++i)
	{
		number aspectRatio = aspectRatios[i];
		for (uint range = 0; range < numRanges; range++)
		{
			if (aspectRatio < minAspectRatio + (range+1)*stepSize)
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
	ug::Table<std::stringstream> aspectRatioTable(numRows, 6);

//	First third
	uint i = 0;
	for(; i < numRows; ++i)
	{
		aspectRatioTable(i, 0) << minAspectRatio + i*stepSize << " - " << minAspectRatio + (i+1)*stepSize << " : ";
		aspectRatioTable(i, 1) << counter[i];
	}

//	Second third
//	Check, if second third of table is needed
	if(i < counter.size())
	{
		for(; i < 2*numRows; ++i)
		{
			aspectRatioTable(i-numRows, 2) << minAspectRatio + i*stepSize << " - " << minAspectRatio + (i+1)*stepSize << " : ";
			aspectRatioTable(i-numRows, 3) << counter[i];
		}
	}

//	Third third
	if(i < counter.size())
//	Check, if third third of table is needed
	{
		for(; i < numRanges; ++i)
		{
			aspectRatioTable(i-2*numRows, 4) << minAspectRatio + i*stepSize << " - " << minAspectRatio + (i+1)*stepSize << " : ";
			aspectRatioTable(i-2*numRows, 5) << counter[i];
		}
	}

//	Output table
	UG_LOG(endl << aspectRatioTable);


//	----------------------------------------
//	Histogram table file output section
//	----------------------------------------
	numRanges = floor(1.0/stepSize);
	counter.clear();
	counter.resize(numRanges, 0.0);
	int numElems = 0;

//	Count the elements in their corresponding aspectRatio range
	for(uint i = 0; i < aspectRatios.size(); ++i)
	{
		number aspectRatio = aspectRatios[i];
		for (uint range = 0; range < numRanges; range++)
		{
			if (aspectRatio < (range+1)*stepSize)
			{
				++counter[range];
				++numElems;
				break;
			}
		}
	}

	outTable.add_rows(numRanges);
	outTable.add_cols(2);

	for(uint i = 0; i < numRanges; ++i)
	{
		outTable(i, 0) << i*stepSize << " - " << (i+1)*stepSize;
		outTable(i, 1) << 100.0/numElems*counter[i];
	}
}
}
