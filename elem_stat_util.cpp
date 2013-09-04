/*
 * eqs_util.cpp
 *
 *  Created on: 03.09.2013
 *      Author: Martin Stepniewski
 */


#include "elem_stat_util.h"


namespace ug
{


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateSubsetSurfaceArea
number CalculateSubsetSurfaceArea(MultiGrid& mg, int subsetIndex, MGSubsetHandler& sh)
{
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	number subsetSurfaceArea = 0.0;

	DistributedGridManager* dgm = mg.distributed_grid_manager();
	for(FaceIterator fIter = sh.begin<Face>(subsetIndex, 0); fIter != sh.end<Face>(subsetIndex, 0); ++fIter)
	{
		Face* f = *fIter;
		#ifdef UG_PARALLEL
			if(dgm->is_ghost(f))
				continue;
		#endif
		subsetSurfaceArea += FaceArea(f, aaPos);
	}

	#ifdef UG_PARALLEL
		if(pcl::GetNumProcesses() > 1){
		//	sum the volumes of all involved processes. Since we ignored ghosts,
		//	each process contributes the volume of a unique part of the grid.
			pcl::ProcessCommunicator pc;
			subsetSurfaceArea = pc.allreduce(subsetSurfaceArea, PCL_RO_SUM);
		}
	#endif

	return subsetSurfaceArea;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateSubsetVolume
number CalculateSubsetVolume(MultiGrid& mg, int subsetIndex, MGSubsetHandler& sh)
{
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	number subsetVolume = 0.0;

	DistributedGridManager* dgm = mg.distributed_grid_manager();//NULL if not parallel
	for(VolumeIterator vIter = sh.begin<Volume>(subsetIndex, 0); vIter != sh.end<Volume>(subsetIndex, 0); ++vIter)
	{
		Volume* v = *vIter;
		#ifdef UG_PARALLEL
		//	ghosts have to be ignored, since they have a copy on another process and
		//	since we alread consider that copy...
			if(dgm->is_ghost(v))
				continue;
		#endif
		subsetVolume += CalculateVolume(*v, aaPos);
	}

	#ifdef UG_PARALLEL
		if(pcl::GetNumProcesses() > 1){
		//	sum the volumes of all involved processes. Since we ignored ghosts,
		//	each process contributes the volume of a unique part of the grid.
			pcl::ProcessCommunicator pc;
			subsetVolume = pc.allreduce(subsetVolume, PCL_RO_SUM);
		}
	#endif

	return subsetVolume;
}







}




