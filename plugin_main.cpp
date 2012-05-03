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

#include <string>

using namespace std;

extern "C" UG_API void InitUGPlugin(ug::bridge::Registry* reg, string parentGroup)
{
	string grp(parentGroup); grp.append("ElementQualityStatistics/");

//	Register ElementQualityStatistics wrapper functions
	reg->add_function(	"ElementQualityStatistics",
						(void (*)(ug::Grid&)) (&ug::ElementQualityStatistics), grp);
	reg->add_function(	"ElementQualityStatistics",
						(void (*)(ug::MultiGrid&, int level)) (&ug::ElementQualityStatistics), grp);


}
