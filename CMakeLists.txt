################################################################################
#
################################################################################

cmake_minimum_required(VERSION 2.6)

project(ELEMENT_QUALITY_STATISTICS)

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")


# sources
set(SOURCES	plugin_main.cpp
			element_quality_statistics.cpp
			volume_calculation.cpp)


# create a shared library from the sources and link it against ug4.
add_library(ElementQualityStatistics SHARED ${SOURCES})

if(STATIC)
	target_link_libraries (ElementQualityStatistics ug4_s)
else(STATIC)
	target_link_libraries (ElementQualityStatistics ug4)
endif(STATIC)