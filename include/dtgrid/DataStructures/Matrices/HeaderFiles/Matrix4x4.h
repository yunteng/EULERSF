/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _matrix4x4_h
#define _matrix4x4_h

#include <iostream>
#include <string>
#include <sstream>


namespace Matrix
{

    template <typename Real> class Vector3;
    template <typename Real> class Vector4;

    template <typename Real>
    class Matrix4x4
    {
    public:
	Matrix4x4();
	Matrix4x4(std::string s);
	Matrix4x4(const Real m[4][4]);

	void load(std::istream& is);
	void save(std::ostream& os) const;


	Matrix4x4 transpose() const;
	// Input: i (row), j (column)
	Real &operator()(unsigned int i, unsigned int j);
	Real *operator[](unsigned int i) { return m[i]; }
	Vector3<Real> operator*(const Vector3<Real>& vec3) const;
	Vector4<Real> operator*(const Vector4<Real>& vec4) const;
	Matrix4x4 operator*(const Matrix4x4& m2) const;

	Vector3<Real> getTranslation() const;
	void setTranslation(const Vector3<Real> t);


	// Static methods

	static Matrix4x4 scale(Real scale);
	static Matrix4x4 scale(Real sx, Real sy, Real sz);
	// Angles in radians
	static Matrix4x4 rotationXYZ(Real rx, Real ry, Real rz);
	static Matrix4x4 translation(Real tx, Real ty, Real tz);
	static Matrix4x4 identity();

    protected:
	Real m[4][4];

	static Matrix4x4 identityMatrix;
    };

}


namespace std
{
    template <typename Real>
	ostream& operator<<(ostream& os, const Matrix::Matrix4x4<Real>& m)
    {
	m.save(os);
	return os;
    }

    template <typename Real>
	istream& operator>>(istream& is, Matrix::Matrix4x4<Real>& m)
    {
	m.load(is);
	return is;
    }
}



#include "Matrix4x4_Impl.h"

#endif
