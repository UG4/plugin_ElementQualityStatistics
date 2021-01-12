/*
 * Copyright (c) 2013-2021:  G-CSC, Goethe University Frankfurt
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

#include "registry/registry.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"
#include "element_quality_statistics.h"
#include "elem_stat_util.h"

#include <string>

using namespace std;

extern "C" UG_API void
InitUGPlugin_ElementQualityStatistics(ug::bridge::Registry* reg, string parentGroup)
{
	string grp(parentGroup); grp.append("ElementQualityStatistics/");

//	Register ElementQualityStatistics wrapper functions
	reg->add_function(	"ElementQualityStatistics",
						(void (*)(ug::Grid&, int, number, number, bool)) (&ug::ElementQualityStatistics),
						grp, "", "grid#dim#angleHiststepSize,aspectRatioHiststepSize", "Prints element quality statistics for a grid object");
	reg->add_function(	"ElementQualityStatistics",
						(void (*)(ug::Grid&, int)) (&ug::ElementQualityStatistics),
						grp, "", "grid#dim", "Prints element quality statistics for a grid object");
	reg->add_function(	"ElementQualityStatistics",
						(void (*)(ug::MultiGrid&, int, number, number, bool)) (&ug::ElementQualityStatistics),
						grp, "", "mg#dim#stepSize", "Prints element quality statistics for a multigrid object");
	reg->add_function(	"ElementQualityStatistics",
						(void (*)(ug::MultiGrid&, int)) (&ug::ElementQualityStatistics),
						grp, "", "mg#dim", "Prints element quality statistics for a multigrid object");

//	Register CalculateSubsetSurfaceArea
	reg->add_function(	"get_subset_surface_area", &ug::CalculateSubsetSurfaceArea,
						grp, "Subset surface area", "mg#subsetIndex#sh", "Returns subset surface area.");
	reg->add_function(	"get_subset_volume", &ug::CalculateSubsetVolume,
						grp, "Subset volume", "mg#subsetIndex#sh", "Returns subset volume.");

//	Register AssignSubsetsByElementQuality
	reg->add_function(	"AssignSubsetsByElementQuality",
						(void (*)(ug::Grid&, ug::SubsetHandler&, int, int)) (&ug::AssignSubsetsByElementQuality),
						grp, "", "grid#dim", "");
	reg->add_function(	"AssignSubsetsByElementQuality",
						(void (*)(ug::MultiGrid&, ug::MGSubsetHandler&, int, int)) (&ug::AssignSubsetsByElementQuality),
						grp, "", "mg#sh", "");
	reg->add_function(	"AssignSubsetToElementWithSmallestMinAngle",
						(void (*)(ug::MultiGrid&, ug::MGSubsetHandler&, int, const char*, int)) (&ug::AssignSubsetToElementWithSmallestMinAngle),
						grp, "", "mg#sh#roid", "");

	reg->add_function(	"MeasureTetrahedronWithSmallestMinAngle",
						(void (*)(ug::MultiGrid&)) (&ug::MeasureTetrahedronWithSmallestMinAngle),
						grp, "", "mg", "");

	reg->add_function(	"FindBoundsForStiffnesMatrixMaxEigenvalue",
						(void (*)(ug::MultiGrid&, ug::MGSubsetHandler&)) (&ug::FindBoundsForStiffnesMatrixMaxEigenvalue),
						grp, "", "mg#sh", "");

	reg->add_function(	"PrintVertexVolumeValence",
						(void (*)(ug::MultiGrid&, ug::MGSubsetHandler&, int)) (&ug::PrintVertexVolumeValence),
						grp, "", "mg#sh#si", "");
}
