// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#include "registry/registry.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"
#include "element_quality_statistics.h"

#include <string>

using namespace std;

void SayHello()
{
	UG_LOG("HELLO From Grid Statistics!\n");
}

extern "C" UG_API void InitUGPlugin(ug::bridge::Registry* reg, string parentGroup)
{
	string grp(parentGroup); grp.append("ElementQualityStatistics/");

	reg->add_function("SayHello", &SayHello, grp);
	reg->add_function("element_quality_statistics", &ug::ElementQualityStatistics, grp);
}
