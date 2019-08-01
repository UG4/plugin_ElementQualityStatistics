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
		//	ghosts (vertical slaves) as well as horizontal slaves (low dimensional elements only) have to be ignored,
		//	since they have a copy on another process and
		//	since we already consider that copy...
			if(dgm->is_ghost(f) || dgm->contains_status(f, ES_H_SLAVE))
				continue;
		#endif
		subsetSurfaceArea += FaceArea(f, aaPos);
	}

	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
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
		//	since we already consider that copy...
			if(dgm->is_ghost(v))
				continue;
		#endif
		subsetVolume += CalculateVolume(v, aaPos);
	}

	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
		//	sum the volumes of all involved processes. Since we ignored ghosts,
		//	each process contributes the volume of a unique part of the grid.
			pcl::ProcessCommunicator pc;
			subsetVolume = pc.allreduce(subsetVolume, PCL_RO_SUM);
		}
	#endif

	return subsetVolume;
}







}




