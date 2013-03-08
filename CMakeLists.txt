################################################################################
#
################################################################################

cmake_minimum_required(VERSION 2.6)

project(ELEMENT_QUALITY_STATISTICS)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


# sources
set(SOURCES	plugin_main.cpp
			element_quality_statistics.cpp
			volume_calculation.cpp)



if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/element_quality_statistics ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(ElementQualityStatistics SHARED ${SOURCES})
	target_link_libraries (ElementQualityStatistics ug4)
endif(buildEmbeddedPlugins)
