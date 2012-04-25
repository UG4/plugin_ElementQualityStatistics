/*
 * volume_calculation.c
 *
 *  Created on: 11.01.2012
 *      Author: marscher
 *
 *  uses algorithm for hexahedron volume calculation from
 *  J. Grandy. October 30, 1997. Efficient Computation of Volume of. Hexahedral Cells.
 */

#include "volume_calculation.h"
#include "lib_grid/lib_grid.h"

#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

namespace ug {

number CalculateVolume(const Volume& vol,
		Grid::VertexAttachmentAccessor<APosition>& aaPos) {
	switch (vol.reference_object_id()) {
	case ROID_TETRAHEDRON:
		return CalculateVolume(static_cast<Tetrahedron>(vol), aaPos);
	case ROID_PRISM:
		return CalculateVolume(static_cast<Prism>(vol), aaPos);
	case ROID_PYRAMID:
		return CalculateVolume(static_cast<Pyramid>(vol), aaPos);
	case ROID_HEXAHEDRON:
		return CalculateVolume(static_cast<Hexahedron>(vol), aaPos);
	default:
		UG_ASSERT(false, "dont know how to calculate given volume.");
	}

	return NAN;
}

number CalculateVolume(const Tetrahedron& tet,
		Grid::VertexAttachmentAccessor<APosition>& aaPos) {
	vector3& a = aaPos[tet.vertex(0)];
	vector3& b = aaPos[tet.vertex(1)];
	vector3& c = aaPos[tet.vertex(2)];
	vector3& d = aaPos[tet.vertex(3)];

	return CalculateTetrahedronVolume(a, b, c, d);
}

// prism
number CalculateVolumePrism(const vector3& a, const vector3& b,
		const vector3& c, const vector3& d, const vector3& e,
		const vector3& f) {
	number result = 0;
	vector3 center;
	vector3 arr[] = { a, b, c, d, e, f };
	CalculateCenter(center, arr, 6);
//	UG_LOG("center prism: " << center << endl)

	result += CalculateTetrahedronVolume(a, b, c, center);
//	UG_LOG("t1: " << result << endl)
	result += CalculateTetrahedronVolume(d, e, f, center);
//	UG_LOG("t1+t2: " << result << endl)

	result +=CalculateVolumePyramid(a, b, e, d, center);
	result +=CalculateVolumePyramid(b, c, f, e, center);
	result +=CalculateVolumePyramid(c, a, d, f, center);

	return result;
}

number CalculateVolume(const Prism& prism,
		Grid::VertexAttachmentAccessor<APosition>& aaPos) {
	vector3& a = aaPos[prism.vertex(0)];
	vector3& b = aaPos[prism.vertex(1)];
	vector3& c = aaPos[prism.vertex(2)];
	vector3& d = aaPos[prism.vertex(3)];
	vector3& e = aaPos[prism.vertex(4)];
	vector3& f = aaPos[prism.vertex(5)];

	return CalculateVolumePrism(a, b, c, d, e, f);
}

// pyramid
number CalculateVolumePyramid(const vector3& a, const vector3& b,
		const vector3& c, const vector3& d, const vector3& e) {
	number result = 0;
//	UG_LOG("a: " << a << " b: " << b << " c: " << c << " d: " << d << " e: " << e << endl)
	//fixme check for a set of volumes that this condition is met!
	// check a,b,c,d are in same plane
	// fixme why does this check only work in x,y plane?!
//	vector3 r, n, x;
//	VecCross(n, a, b);
//	VecNormalize(n, n);
//
//	VecSubtract(r, a, b);
//	VecSubtract(x, c, r);
//	number dot = VecDot(x, n);
//
//	UG_LOG("dot pyra: " << dot<< endl)
//	UG_ASSERT(dot < SMALL,
//			"pyramid volume calculation needs all base are points in one plane!");

	vector3 center, h_, h, c1, c2, da, ba, cb, cd;
	VecSubtract(da, d, a);
	VecSubtract(ba, b, a);
	VecSubtract(cb, c, b);
	VecSubtract(cd, c, d);
	VecCross(c1, da, ba);
	VecCross(c2, cb, cd);
	number A = 0.5 * VecLength(c1) + 0.5 * VecLength(c2);
//	UG_LOG("A pyra: " << A <<endl)

	vector3 arr[] = { a, b, c, d };
	CalculateCenter(center, arr, 4);

	number height = DistancePointToPlane(e, center, c1);

//	VecSubtract(h_, e, center);
//	VecAdd(h, h_, center);
//	VecSubtract(h, e, center);
//	 VecLength(h);
//	UG_LOG("pyra h: " << height << endl)

	result = 1.0 / 3.0 * A * height;
//	UG_LOG("pyra vol: " << result << endl)

	return result;
}

number CalculateVolume(const Pyramid& pyramid,
		Grid::VertexAttachmentAccessor<APosition>& aaPos) {
	vector3& a = aaPos[pyramid.vertex(0)];
	vector3& b = aaPos[pyramid.vertex(1)];
	vector3& c = aaPos[pyramid.vertex(2)];
	vector3& d = aaPos[pyramid.vertex(3)];
	vector3& top = aaPos[pyramid.vertex(4)];

	return CalculateVolumePyramid(a, b, c, d, top);
}

/**
 * Algorithm (14) from J.Grandy "Efficient Computation of Volume of Hexahedral Cells"
 * which allows volume calculation of irregular hexahedrons
 */
number CalculateVolume(const Hexahedron& hexahedron,
		Grid::VertexAttachmentAccessor<APosition>& aaPos) {
	vector3& a = aaPos[hexahedron.vertex(0)];
	vector3& b = aaPos[hexahedron.vertex(1)];
	vector3& c = aaPos[hexahedron.vertex(2)];
	vector3& d = aaPos[hexahedron.vertex(3)];
	vector3& e = aaPos[hexahedron.vertex(4)];
	vector3& f = aaPos[hexahedron.vertex(5)];
	vector3& g = aaPos[hexahedron.vertex(6)];
	vector3& h = aaPos[hexahedron.vertex(7)];
	return CalculateVolumeHexahedron(a, b, c, d, e, f, g, h);
//	number result = 0;
//
//	vector3 v[] =
//			{ aaPos[hexa.vertex(0)], aaPos[hexa.vertex(1)],
//					aaPos[hexa.vertex(2)], aaPos[hexa.vertex(3)],
//					aaPos[hexa.vertex(4)], aaPos[hexa.vertex(5)],
//					aaPos[hexa.vertex(6)], aaPos[hexa.vertex(7)] };
//
//	// matrices to calculate determinant
//	matrix33 m1, m2, m3;
//
//	// determine long diagonal
//	vector<number> diagonalLength;
//	vector<number>::iterator iter;
//
//	// diagonals
//	vector3 v6v0, v7v1, v4v2, v3v5;
//	// longest diagonal
//	vector3* ld = NULL;
//	pair<int, int> LD;
//	// triangulation vectors
//	pair<int, int> T1, T2, T3;
//	// edge
//	pair<int, int> E;
//
//	VecSubtract(v6v0, v[6], v[0]);
//	VecSubtract(v7v1, v[7], v[1]);
//	VecSubtract(v4v2, v[4], v[2]);
//	VecSubtract(v3v5, v[3], v[5]);
//
//	diagonalLength.push_back(VecLength(v6v0));
//	diagonalLength.push_back(VecLength(v7v1));
//	diagonalLength.push_back(VecLength(v4v2));
//	diagonalLength.push_back(VecLength(v3v5));
//	// determine longest diagonal
//	iter = max_element(diagonalLength.begin(), diagonalLength.end());
//	uint pos = distance(diagonalLength.begin(), iter);
//	switch (pos) {
//	case 0: {
//		LD = make_pair(6, 0);
//		ld = &v6v0;
//		T1 = make_pair(2, 5);
//		T2 = make_pair(7, 5);
//		T3 = make_pair(2, 7);
//		E = make_pair(1, 0);
//		break;
//	}
//	case 1: {
//		LD = make_pair(7, 1);
//		ld = &v7v1;
//		T1 = make_pair(4, 3);
////		T2 = make_pair();
//		E = make_pair(1, 0);
//		break;
//	}
//	case 2: {
//		LD = make_pair(4, 2);
//		ld = &v4v2;
//		T1 = make_pair(5, 0);
//		E = make_pair(2, 1);
//		break;
//	}
//	case 3: {
//		LD = make_pair(3, 5);
//		ld = &v3v5;
//		T1 = make_pair(4, 1);
//		E = make_pair(3, 0);
//	}
//	}
//
//	// fill determinants
//	m1.assign(*ld, 0);
//	m2.assign(*ld, 0);
//	m3.assign(*ld, 0);
//
//	// shuffle t and e
////	vector3 t, e;
////	VecSubtract(t, v[T1.first], v[T1.second]);
////	m1.assign(t, 1);
////	m1.assign(e, 2);
//
//	result = 1. / 6 * (Determinant(m1) + Determinant(m2) + Determinant(m3));
//	return result;
}

number CalculateVolumeHexahedron(const vector3& a, const vector3& b,
		const vector3& c, const vector3& d, const vector3& e, const vector3& f,
		const vector3& g, const vector3& h) {
	number result = 0;
	vector3 arr[] = { a, b, c, d, e, f, g, h };
	vector3 center;
	CalculateCenter(center, arr, 8);

	// top and bottom
	result += CalculateVolumePyramid(a, b, c, d, center);
	result += CalculateVolumePyramid(e, f, g, h, center);

	// sides
	result += CalculateVolumePyramid(a, b, f, e, center);
	result += CalculateVolumePyramid(b, c, g, f, center);

	result += CalculateVolumePyramid(c, d, g, h, center);
	result += CalculateVolumePyramid(a, d, h, e, center);

	return result;
}

number CalculateVolume(geometry_traits<Volume>::iterator begin,
		geometry_traits<Volume>::iterator end,
		Grid::VertexAttachmentAccessor<APosition>& aaPos) {
	number result = 0;
	geometry_traits<Volume>::iterator iter = begin;
	while (iter != end) {
		Volume* vol = *iter;
		result += CalculateVolume(*vol, aaPos);
		iter++;
	}
	return result;
}

} // end of namespace ug
