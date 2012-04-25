/*
 * volume_calculation.h
 *
 *  Created on: 11.01.2012
 *      Author: marscher
 */

#ifndef VOLUME_CALCULATION_H_
#define VOLUME_CALCULATION_H_

#include "lib_grid/lib_grid.h"

namespace ug {

// wrappers for high level volume classes
number CalculateVolume(const Tetrahedron&,
		Grid::VertexAttachmentAccessor<APosition>&);
number CalculateVolume(const Prism&,
		Grid::VertexAttachmentAccessor<APosition>&);
number CalculateVolume(const Pyramid&,
		Grid::VertexAttachmentAccessor<APosition>&);
number CalculateVolume(const Hexahedron&,
		Grid::VertexAttachmentAccessor<APosition>&);
number CalculateVolume(const Volume&,
		Grid::VertexAttachmentAccessor<APosition>&);
number CalculateVolume(geometry_traits<Volume>::iterator,
		geometry_traits<Volume>::iterator,
		Grid::VertexAttachmentAccessor<APosition>&);

// tetrahedron volume calculation from msteppnie used
// pyramid
number CalculateVolumePyramid(const vector3& v0, const vector3& v1,
		const vector3& v2, const vector3& v3, const vector3& v4);

// prism
number CalculateVolumePrism(const vector3& v0, const vector3& v1,
		const vector3& v2, const vector3& v3, const vector3& v4,
		const vector3& v5);

// hexahedron
number CalculateVolumeHexahedron(const vector3& v0, const vector3& v1,
		const vector3& v2, const vector3& v3, const vector3& v4,
		const vector3& v5, const vector3& v6, const vector3& v8);

} // end namespace ug

#endif /* VOLUME_CALCULATION_H_ */
