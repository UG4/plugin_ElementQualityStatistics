/*
 * element_quality_statistics.cpp
 *
 *  Created on: 17.04.2012
 *      Author: Martin Stepniewski
 */


#ifndef __ELEMENT_QUALITY_STATISTICS_H__
#define __ELEMENT_QUALITY_STATISTICS_H__

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
#include "elem_stat_util.h"
#include "lib_grid/algorithms/element_angles.h"
#include "lib_grid/algorithms/element_aspect_ratios.h"



using namespace std;


namespace ug {


////////////////////////////////////////////////////////////////////////////////////////////
//	CollectMinAngles
template <class TIterator, class TAAPosVRT>
void CollectMinAngles(Grid& grid, TIterator elementsBegin,
					  TIterator elementsEnd,
					  TAAPosVRT& aaPos,
					  vector<number>& minAngles)
{
	//PROFILE_FUNC();
	DistributedGridManager* dgm = grid.distributed_grid_manager();

//	if elementsBegin equals elementsEnd, then the list is empty and we can
//	immediately return NULL
	if(elementsBegin == elementsEnd)
	{
		UG_LOG("ERROR in CollectMinAngles: elementsBegin == elementsEnd." << std::endl);
		return;
	}

//	Calculate the minAngle of every element
	for(TIterator iter = elementsBegin; iter != elementsEnd; ++iter)
	{
		#ifdef UG_PARALLEL
		//	ghosts (vertical masters) as well as horizontal slaves (low dimensional elements only) have to be ignored,
		//	since they have a copy on another process and
		//	since we already consider that copy...
			if(dgm->is_ghost(*iter) || dgm->contains_status(*iter, ES_H_SLAVE))
				continue;
		#endif

		number curMinAngle = CalculateMinAngle(grid, *iter, aaPos);
		minAngles.push_back(curMinAngle);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CollectMaxAngles
template <class TIterator, class TAAPosVRT>
void CollectMaxAngles(Grid& grid, TIterator elementsBegin,
					  TIterator elementsEnd,
					  TAAPosVRT& aaPos,
					  vector<number>& maxAngles)
{
	//PROFILE_FUNC();
	DistributedGridManager* dgm = grid.distributed_grid_manager();

//	if elementsBegin equals elementsEnd, then the list is empty and we can
//	immediately return NULL
	if(elementsBegin == elementsEnd)
	{
		UG_LOG("ERROR in CollectMaxAngles: elementsBegin == elementsEnd." << std::endl);
		return;
	}

//	Calculate the minAngle of every element
	for(TIterator iter = elementsBegin; iter != elementsEnd; ++iter)
	{
		#ifdef UG_PARALLEL
		//	ghosts (vertical masters) as well as horizontal slaves (low dimensional elements only) have to be ignored,
		//	since they have a copy on another process and
		//	since we already consider that copy...
			if(dgm->is_ghost(*iter) || dgm->contains_status(*iter, ES_H_SLAVE))
				continue;
		#endif

		number curMaxAngle = CalculateMaxAngle(grid, *iter, aaPos);
		maxAngles.push_back(curMaxAngle);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	AspectRatioHistogram
template <class TIterator, class TAAPosVRT>
void CollectAspectRatios(Grid& grid, TIterator elementsBegin,
						 TIterator elementsEnd,
						 TAAPosVRT& aaPos,
						 vector<number>& aspectRatios)
{
	//PROFILE_FUNC();
	DistributedGridManager* dgm = grid.distributed_grid_manager();

//	if elementsBegin equals elementsEnd, then the list is empty and we can
//	immediately return NULL
	if(elementsBegin == elementsEnd)
	{
		UG_LOG("ERROR in CollectAspectRatios: elementsBegin == elementsEnd." << std::endl);
		return;
	}

//	Calculate the aspectRatio of every element
	for(TIterator iter = elementsBegin; iter != elementsEnd; ++iter)
	{
		#ifdef UG_PARALLEL
		//	ghosts (vertical masters) as well as horizontal slaves (low dimensional elements only) have to be ignored,
		//	since they have a copy on another process and
		//	since we already consider that copy...
			if(dgm->is_ghost(*iter) || dgm->contains_status(*iter, ES_H_SLAVE))
				continue;
		#endif

		number curAspectRatio = CalculateAspectRatio(grid, *iter, aaPos);
		aspectRatios.push_back(curAspectRatio);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	VolToRMSFaceAreaRatioHistogram
template <class TIterator, class TAAPosVRT>
void CollectVolToRMSFaceAreaRatios(Grid& grid, TIterator elementsBegin,
						  	  	   TIterator elementsEnd,
								   TAAPosVRT& aaPos,
								   vector<number>& ratios)
{
	//PROFILE_FUNC();
	DistributedGridManager* dgm = grid.distributed_grid_manager();

//	if elementsBegin equals elementsEnd, then the list is empty and we can
//	immediately return NULL
	if(elementsBegin == elementsEnd)
	{
		UG_LOG("ERROR in CollectVolToRMSFaceAreaRatios: elementsBegin == elementsEnd." << std::endl);
		return;
	}

//	Calculate the ratio of every element
	for(TIterator iter = elementsBegin; iter != elementsEnd; ++iter)
	{
		#ifdef UG_PARALLEL
		//	ghosts (vertical masters) as well as horizontal slaves (low dimensional elements only) have to be ignored,
		//	since they have a copy on another process and
		//	since we already consider that copy...
			if(dgm->is_ghost(*iter) || dgm->contains_status(*iter, ES_H_SLAVE))
				continue;
		#endif

		number curRatio = CalculateVolToRMSFaceAreaRatio(grid, *iter, aaPos);
		ratios.push_back(curRatio);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	PrintAngleStatistics2d
template <class TAAPosVRT>
void PrintAngleStatistics2d(Grid& grid, GridObjectCollection& goc, int level, TAAPosVRT& aaPos)
{
	size_t numTriangles = 0;
	size_t numQuadrilaterals = 0;

	number sd_tri = 0.0;
	number sd_quad = 0.0;
	number mean_tri = 0.0;
	number mean_quad = 0.0;
	number regAngle = 90.0;

	DistributedGridManager* dgm = grid.distributed_grid_manager();
	int i = level;

//	Calculate and output standard deviation for triangular/quadrilateral angles
	if(goc.num<Triangle>(i) > 0 || goc.num<Quadrilateral>(i) > 0)
	{
		vector<number> vAnglesOut;

		for(FaceIterator fIter = goc.begin<Face>(i); fIter != goc.end<Face>(i); ++fIter)
		{
			vAnglesOut.clear();
			Face* f = *fIter;

			#ifdef UG_PARALLEL
			//	ghosts (vertical masters) as well as horizontal slaves (low dimensional elements only) have to be ignored,
			//	since they have a copy on another process and
			//	since we already consider that copy...
				if(dgm->is_ghost(f) || dgm->contains_status(f, ES_H_SLAVE))
					continue;
			#endif

			if(f->reference_object_id() == ROID_TRIANGLE)
			{
				regAngle = 60.0;
				CalculateAngles(vAnglesOut, grid, f, aaPos);
				numTriangles++;

				for(size_t k = 0; k < vAnglesOut.size(); ++k)
				{
					sd_tri += (regAngle-vAnglesOut[k])*(regAngle-vAnglesOut[k]);
					mean_tri += vAnglesOut[k];
				}
			}

			if(f->reference_object_id() == ROID_QUADRILATERAL)
			{
				regAngle = 90.0;
				CalculateAngles(vAnglesOut, grid, f, aaPos);
				numQuadrilaterals++;

				for(size_t k = 0; k < vAnglesOut.size(); ++k)
				{
					sd_quad += (regAngle-vAnglesOut[k])*(regAngle-vAnglesOut[k]);
					mean_quad += vAnglesOut[k];
				}
			}
		}
	}

	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
		//	sum the numbers of all involved processes. Since we ignored ghosts,
		//	each process contributes the numbers of a unique part of the grid.
			pcl::ProcessCommunicator pc;
			numTriangles = pc.allreduce(numTriangles, PCL_RO_SUM);
			numQuadrilaterals = pc.allreduce(numQuadrilaterals, PCL_RO_SUM);
			sd_tri = pc.allreduce(sd_tri, PCL_RO_SUM);
			mean_tri = pc.allreduce(mean_tri, PCL_RO_SUM);
			sd_quad = pc.allreduce(sd_quad, PCL_RO_SUM);
			mean_quad = pc.allreduce(mean_quad, PCL_RO_SUM);
		}
	#endif

//	Calculate and output standard deviation for triangular/quadrilateral angles
	if(goc.num<Triangle>(i) > 0 || goc.num<Quadrilateral>(i) > 0)
	{
		if(goc.num<Triangle>(i) > 0)
		{
			sd_tri *= (1.0/(3*numTriangles));
			sd_tri = sqrt(sd_tri);
			mean_tri *= (1.0/(3*numTriangles));
		}

		if(goc.num<Quadrilateral>(i) > 0)
		{
			sd_quad *= (1.0/(4*numQuadrilaterals));
			sd_quad = sqrt(sd_quad);
			mean_quad *= (1.0/(4*numQuadrilaterals));
		}

		UG_LOG("(*) Standard deviation of face angles to regular case" << endl);
		UG_LOG("	(60° for triangles, 90° for quadrilaterals)" << endl);
		UG_LOG(endl);
		UG_LOG("	Triangles (" << numTriangles << "):" << endl);
		if(goc.num<Triangle>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_tri << endl);
			UG_LOG("		mean = " << mean_tri << endl);
		}
		UG_LOG(endl);
		UG_LOG("	Quadrilaterals (" << numQuadrilaterals << "):" << endl);
		if(goc.num<Quadrilateral>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_quad << endl);
			UG_LOG("		mean = " << mean_quad << endl);
		}
		UG_LOG(endl);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	PrintAngleStatistics3d
template <class TAAPosVRT>
void PrintAngleStatistics3d(Grid& grid, GridObjectCollection& goc, int level, TAAPosVRT& aaPos)
{
	size_t numTriangles = 0;
	size_t numQuadrilaterals = 0;
	size_t numTetrahedrons = 0;
	size_t numHexahedrons = 0;
	size_t numOctahedrons = 0;

	number sd_tri = 0.0;
	number sd_quad = 0.0;
	number mean_tri = 0.0;
	number mean_quad = 0.0;
	number regAngle = 90.0;

	number sd_tet = 0.0;
	number sd_hex = 0.0;
	number sd_oct = 0.0;
	number mean_tet = 0.0;
	number mean_hex = 0.0;
	number mean_oct = 0.0;
	number regVolDihedral = 90.0;

	DistributedGridManager* dgm = grid.distributed_grid_manager();
	int i = level;

//	Calculate and output standard deviation for triangular/quadrilateral angles
	if(goc.num<Triangle>(i) > 0 || goc.num<Quadrilateral>(i) > 0)
	{
		vector<number> vAnglesOut;

		for(FaceIterator fIter = goc.begin<Face>(i); fIter != goc.end<Face>(i); ++fIter)
		{
			vAnglesOut.clear();
			Face* f = *fIter;

			#ifdef UG_PARALLEL
			//	ghosts (vertical masters) as well as horizontal slaves (low dimensional elements only) have to be ignored,
			//	since they have a copy on another process and
			//	since we already consider that copy...
				if(dgm->is_ghost(f) || dgm->contains_status(f, ES_H_SLAVE))
					continue;
			#endif

			if(f->reference_object_id() == ROID_TRIANGLE)
			{
				regAngle = 60.0;
				CalculateAngles(vAnglesOut, grid, f, aaPos);
				numTriangles++;

				for(size_t k = 0; k < vAnglesOut.size(); ++k)
				{
					sd_tri += (regAngle-vAnglesOut[k])*(regAngle-vAnglesOut[k]);
					mean_tri += vAnglesOut[k];
				}
			}

			if(f->reference_object_id() == ROID_QUADRILATERAL)
			{
				regAngle = 90.0;
				CalculateAngles(vAnglesOut, grid, f, aaPos);
				numQuadrilaterals++;

				for(size_t k = 0; k < vAnglesOut.size(); ++k)
				{
					sd_quad += (regAngle-vAnglesOut[k])*(regAngle-vAnglesOut[k]);
					mean_quad += vAnglesOut[k];
				}
			}
		}
	}

	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
		//	sum the numbers of all involved processes. Since we ignored ghosts,
		//	each process contributes the numbers of a unique part of the grid.
			pcl::ProcessCommunicator pc;
			numTriangles = pc.allreduce(numTriangles, PCL_RO_SUM);
			numQuadrilaterals = pc.allreduce(numQuadrilaterals, PCL_RO_SUM);
			sd_tri = pc.allreduce(sd_tri, PCL_RO_SUM);
			mean_tri = pc.allreduce(mean_tri, PCL_RO_SUM);
			sd_quad = pc.allreduce(sd_quad, PCL_RO_SUM);
			mean_quad = pc.allreduce(mean_quad, PCL_RO_SUM);
		}
	#endif

//	Calculate and output standard deviation for triangular/quadrilateral angles
	if(goc.num<Triangle>(i) > 0 || goc.num<Quadrilateral>(i) > 0)
	{
		if(goc.num<Triangle>(i) > 0)
		{
			sd_tri *= (1.0/(3*numTriangles));
			sd_tri = sqrt(sd_tri);
			mean_tri *= (1.0/(3*numTriangles));
		}

		if(goc.num<Quadrilateral>(i) > 0)
		{
			sd_quad *= (1.0/(4*numQuadrilaterals));
			sd_quad = sqrt(sd_quad);
			mean_quad *= (1.0/(4*numQuadrilaterals));
		}

		UG_LOG("(*) Standard deviation of face angles to regular case" << endl);
		UG_LOG("	(60° for triangles, 90° for quadrilaterals)" << endl);
		UG_LOG(endl);
		UG_LOG("	Triangles (" << numTriangles << "):" << endl);
		if(goc.num<Triangle>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_tri << endl);
			UG_LOG("		mean = " << mean_tri << endl);
		}
		UG_LOG(endl);
		UG_LOG("	Quadrilaterals (" << numQuadrilaterals << "):" << endl);
		if(goc.num<Quadrilateral>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_quad << endl);
			UG_LOG("		mean = " << mean_quad << endl);
		}
		UG_LOG(endl);
	}

//	Calculate and output standard deviation for tetrahedral/hexahedral angles
	if(goc.num<Tetrahedron>(i) > 0 || goc.num<Hexahedron>(i) > 0 || goc.num<Octahedron>(i) > 0)
	{
		vector<number> vDihedralsOut;

		for(VolumeIterator vIter = goc.begin<Volume>(i); vIter != goc.end<Volume>(i); ++vIter)
		{
			vDihedralsOut.clear();
			Volume* vol = *vIter;

			#ifdef UG_PARALLEL
			//	ghosts (vertical masters) have to be ignored,
			//	since they have a copy on another process and
			//	since we already consider that copy...
				if(dgm->is_ghost(vol))
					continue;
			#endif

			if(vol->reference_object_id() == ROID_TETRAHEDRON)
			{
				regVolDihedral = 70.52877937;
				CalculateAngles(vDihedralsOut, grid, vol, aaPos);
				numTetrahedrons++;

				for(size_t k = 0; k < vDihedralsOut.size(); ++k)
				{
					sd_tet += (regVolDihedral-vDihedralsOut[k])*(regVolDihedral-vDihedralsOut[k]);
					mean_tet += vDihedralsOut[k];
				}
			}

			if(vol->reference_object_id() == ROID_HEXAHEDRON)
			{
				regVolDihedral = 90.0;
				CalculateAngles(vDihedralsOut, grid, vol, aaPos);
				numHexahedrons++;

				for(size_t k = 0; k < vDihedralsOut.size(); ++k)
				{
					sd_hex += (regVolDihedral-vDihedralsOut[k])*(regVolDihedral-vDihedralsOut[k]);
					mean_hex += vDihedralsOut[k];
				}
			}

			if(vol->reference_object_id() == ROID_OCTAHEDRON)
			{
				regVolDihedral = 109.4712206;
				CalculateAngles(vDihedralsOut, grid, vol, aaPos);
				numOctahedrons++;

				for(size_t k = 0; k < vDihedralsOut.size(); ++k)
				{
					sd_oct += (regVolDihedral-vDihedralsOut[k])*(regVolDihedral-vDihedralsOut[k]);
					mean_oct += vDihedralsOut[k];
				}
			}
		}
	}

	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
		//	sum the numbers of all involved processes. Since we ignored ghosts,
		//	each process contributes the numbers of a unique part of the grid.
			pcl::ProcessCommunicator pc;
			numTetrahedrons = pc.allreduce(numTetrahedrons, PCL_RO_SUM);
			numHexahedrons = pc.allreduce(numHexahedrons, PCL_RO_SUM);
			numOctahedrons = pc.allreduce(numOctahedrons, PCL_RO_SUM);
			sd_tet = pc.allreduce(sd_tet, PCL_RO_SUM);
			mean_tet = pc.allreduce(mean_tet, PCL_RO_SUM);
			sd_hex = pc.allreduce(sd_hex, PCL_RO_SUM);
			mean_hex = pc.allreduce(mean_hex, PCL_RO_SUM);
			sd_oct = pc.allreduce(sd_oct, PCL_RO_SUM);
			mean_oct = pc.allreduce(mean_oct, PCL_RO_SUM);
		}
	#endif

//	Calculate and output standard deviation for tetrahedral/hexahedral angles
	if(goc.num<Tetrahedron>(i) > 0 || goc.num<Hexahedron>(i) > 0 || goc.num<Octahedron>(i) > 0)
	{
		if(goc.num<Tetrahedron>(i) > 0)
		{
			sd_tet *= (1.0/(6*numTetrahedrons));
			sd_tet = sqrt(sd_tet);
			mean_tet *= (1.0/(6*numTetrahedrons));
		}

		if(goc.num<Hexahedron>(i) > 0)
		{
			sd_hex *= (1.0/(12*numHexahedrons));
			sd_hex = sqrt(sd_hex);
			mean_hex *= (1.0/(12*numHexahedrons));
		}

		if(goc.num<Octahedron>(i) > 0)
		{
			sd_oct *= (1.0/(12*numOctahedrons));
			sd_oct = sqrt(sd_oct);
			mean_oct *= (1.0/(12*numOctahedrons));
		}

		UG_LOG("(*) Standard deviation of dihedral angles to regular case" << endl);
		UG_LOG("	(70.5288° for tetrahedrons, 90° for hexahedrons, 109.471° for Octahedrons)" << endl);
		UG_LOG(endl);
		UG_LOG("	Tetrahedrons (" << numTetrahedrons << "):" << endl);
		if(goc.num<Tetrahedron>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_tet << endl);
			UG_LOG("		mean = " << mean_tet << endl);
		}
		UG_LOG(endl);
		UG_LOG("	Hexahedrons (" << numHexahedrons << "):" << endl);
		if(goc.num<Hexahedron>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_hex << endl);
			UG_LOG("		mean = " << mean_hex << endl);
		}
		UG_LOG(endl);
		UG_LOG("	Octahedrons (" << numOctahedrons << "):" << endl);
		if(goc.num<Octahedron>(i) > 0)
		{
			UG_LOG("		sd   = " << sd_oct << endl);
			UG_LOG("		mean = " << mean_oct << endl);
		}
		UG_LOG(endl);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	PrintAngleHistograms
void PrintAngleHistogram(vector<number>& locAngles, number stepSize, ug::Table<std::stringstream>& outTable);
void PrintAspectRatioHistogram(vector<number>& locAspectRatios, number stepSize, ug::Table<std::stringstream>& outTable);


////////////////////////////////////////////////////////////////////////////////////////////
//	PrintVertexVolumeValence
void PrintVertexVolumeValence(MultiGrid& mg, SubsetHandler& sh, int subsetIndex);


////////////////////////////////////////////////////////////////////////////////////////////
//	AssignSubsetToElementWithSmallestMinAngle
void AssignSubsetToElementWithSmallestMinAngle(MultiGrid& grid, MGSubsetHandler& sh, int dim, const char* roid);
void AssignSubsetToElementWithSmallestMinAngle2d(MultiGrid& grid, MGSubsetHandler& sh, const char* roid);
void AssignSubsetToElementWithSmallestMinAngle3d(MultiGrid& grid, MGSubsetHandler& sh, const char* roid);


////////////////////////////////////////////////////////////////////////////////////////////
//	AssignSubsetsByElementQuality3d
void AssignSubsetsByElementQuality3d(MultiGrid& mg, MGSubsetHandler& sh, int numSecs);


////////////////////////////////////////////////////////////////////////////////////////////
//	AssignSubsetsByElementQuality
void AssignSubsetsByElementQuality(MultiGrid& mg, MGSubsetHandler& sh, int dim, int numSecs);


////////////////////////////////////////////////////////////////////////////////////////////
//	AssignSubsetsByElementQuality
void AssignSubsetsByElementQuality(Grid& grid, SubsetHandler& sh, int dim, int numSecs);


////////////////////////////////////////////////////////////////////////////////////////////
//	ElementQualityStatistics
////////////////////////////////////////////////////////////////////////////////////////////

//	Wrapper
void ElementQualityStatistics(MultiGrid& mg, int dim, number angleHistStepSize, number aspectRatioHistStepSize, bool bWriteHistograms);
void ElementQualityStatistics(MultiGrid& mg, int dim);
void ElementQualityStatistics(Grid& grid, int dim, number angleHistStepSize, number aspectRatioHistStepSize, bool bWriteHistograms);
void ElementQualityStatistics(Grid& grid, int dim);

//	Actual procedures
void ElementQualityStatistics2d(Grid& grid, GridObjectCollection goc, number angleHistStepSize = 10.0, number aspectRatioHistStepSize = 0.1, bool bWriteHistograms = true);
void ElementQualityStatistics3d(Grid& grid, GridObjectCollection goc, number angleHistStepSize = 10.0, number aspectRatioHistStepSize = 0.1, bool bWriteHistograms = true);


}	 
#endif  //__ELEMENT_QUALITY_STATISTICS_H__

