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
						(void (*)(ug::Grid&, int)) (&ug::ElementQualityStatistics),
						grp, "", "grid#dim", "Prints element quality statistics for a grid object");
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
						(void (*)(ug::Grid&, ug::SubsetHandler& sh, int, int)) (&ug::AssignSubsetsByElementQuality),
						grp, "", "grid#dim", "");
	reg->add_function(	"AssignSubsetsByElementQuality",
						(void (*)(ug::MultiGrid&, ug::MGSubsetHandler& sh, int, int)) (&ug::AssignSubsetsByElementQuality),
						grp, "", "mg#dim", "");
}
