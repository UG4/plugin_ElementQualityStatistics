/*
 * plugin_main.cpp
 *
 *  Created on: 17.04.2012
 *      Author: Martin Stepniewski
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
						(void (*)(ug::Grid&)) (&ug::ElementQualityStatistics),
						grp, "", "grid", "Prints element quality statistics for a grid object");
	reg->add_function(	"ElementQualityStatistics",
						(void (*)(ug::MultiGrid&/*, int level*/)) (&ug::ElementQualityStatistics),
						grp, "", "mg", "Prints element quality statistics for a multigrid object");

//	Register CalculateSubsetSurfaceArea
	reg->add_function(	"get_subset_surface_area", &ug::CalculateSubsetSurfaceArea,
						grp, "Subset surface area", "mg#subsetIndex#sh", "Returns subset surface area.");
	reg->add_function(	"get_subset_volume", &ug::CalculateSubsetVolume,
						grp, "Subset volume", "mg#subsetIndex#sh", "Returns subset volume.");


//	Register RefineTetVolumeSmoothly
	reg->add_function(	"SubdivisionTetGridSmooth", &ug::SubdivisionTetGridSmooth,
			grp, "SubdivisionTetGridSmooth", "mg", "SubdivisionTetGridSmooth");

	reg->add_function(	"SubdivisionTetGridSmoothBasic", &ug::SubdivisionTetGridSmoothBasic,
			grp, "SubdivisionTetGridSmoothBasic", "mg", "SubdivisionTetGridSmoothBasic");

}
