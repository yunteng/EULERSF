/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _matrix4x4_impl_h
#define _matrix4x4_impl_h

#include <cmath>
#include "../../Vectors/HeaderFiles/Vector3.h"
#include "../../Vectors/HeaderFiles/Vector4.h"



using namespace std;


namespace Matrix
{

    template <typename Real> Matrix4x4<Real> Matrix4x4<Real>::identityMatrix = Matrix4x4<Real>(); 



    template <typename Real>
	Matrix4x4<Real>::Matrix4x4()
    {
	m[0][0] = 1; m[0][1] = 0; m[0][2] = 0; m[0][3] = 0;
	m[1][0] = 0; m[1][1] = 1; m[1][2] = 0; m[1][3] = 0;
	m[2][0] = 0; m[2][1] = 0; m[2][2] = 1; m[2][3] = 0;
	m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 1;
    }


    template <typename Real>
	Matrix4x4<Real>::Matrix4x4(const Real m[4][4])
    {
	this->m[0][0] = m[0][0]; this->m[0][1] = m[0][1]; this->m[0][2] = m[0][2]; this->m[0][3] = m[0][3];
	this->m[1][0] = m[1][0]; this->m[1][1] = m[1][1]; this->m[1][2] = m[1][2]; this->m[1][3] = m[1][3];
	this->m[2][0] = m[2][0]; this->m[2][1] = m[2][1]; this->m[2][2] = m[2][2]; this->m[2][3] = m[2][3];
	this->m[3][0] = m[3][0]; this->m[3][1] = m[3][1]; this->m[3][2] = m[3][2]; this->m[3][3] = m[3][3];
    }


    template <typename Real>
	Matrix4x4<Real>::Matrix4x4(std::string s)
    {
	std::stringstream ss;
	ss << s;
	load(ss);
    }


    template <typename Real>
	void Matrix4x4<Real>::load(std::istream& is)
    {
	is >> m[0][0]; is >> m[0][1]; is >> m[0][2]; is >> m[0][3];   // row 1
	is >> m[1][0]; is >> m[1][1]; is >> m[1][2]; is >> m[1][3];   // row 2
	is >> m[2][0]; is >> m[2][1]; is >> m[2][2]; is >> m[2][3];   // row 3
	is >> m[3][0]; is >> m[3][1]; is >> m[3][2]; is >> m[3][3];   // row 4
    }


    template <typename Real>
	void Matrix4x4<Real>::save(std::ostream& os) const
    {
	os << m[0][0]; os << " "; os << m[0][1]; os << " "; os << m[0][2]; os << " "; os << m[0][3]; os << " ";   // row 1
	os << m[1][0]; os << " "; os << m[1][1]; os << " "; os << m[1][2]; os << " "; os << m[1][3]; os << " ";   // row 2
	os << m[2][0]; os << " "; os << m[2][1]; os << " "; os << m[2][2]; os << " "; os << m[2][3]; os << " ";   // row 3
	os << m[3][0]; os << " "; os << m[3][1]; os << " "; os << m[3][2]; os << " "; os << m[3][3]; os << " ";   // row 4
    }


    template <typename Real>
	Vector3<Real> Matrix4x4<Real>::getTranslation() const
    {
	Vector3<Real> v3;
	v3[0] = m[0][3];
	v3[1] = m[1][3];
	v3[2] = m[2][3];
	return v3;
    }


    template <typename Real>
	void Matrix4x4<Real>::setTranslation(const Vector3<Real> t)
    {
	m[0][3] = t[0];
	m[1][3] = t[1];
	m[2][3] = t[2];
    }

	
    template <typename Real>
	Matrix4x4<Real> Matrix4x4<Real>::identity()
    {
	return identityMatrix;
    }


    template <typename Real>
	Vector3<Real> Matrix4x4<Real>::operator*(const Vector3<Real>& vec3) const
    {
	Vector4<Real> vec4 = Vector4<Real>(vec3, Real(1.0));

	return (*this) * vec4;
    }


    template <typename Real>
	Vector4<Real> Matrix4x4<Real>::operator*(const Vector4<Real>& vec4) const
    {
	Vector4<Real> res;

	res[0] = vec4[0] * m[0][0] + vec4[1] * m[0][1] + vec4[2] * m[0][2] + vec4[3] * m[0][3];
	res[1] = vec4[0] * m[1][0] + vec4[1] * m[1][1] + vec4[2] * m[1][2] + vec4[3] * m[1][3];
	res[2] = vec4[0] * m[2][0] + vec4[1] * m[2][1] + vec4[2] * m[2][2] + vec4[3] * m[2][3];
	res[3] = vec4[0] * m[3][0] + vec4[1] * m[3][1] + vec4[2] * m[3][2] + vec4[3] * m[3][3];

	return res;
    }


    template <typename Real>
	Matrix4x4<Real> Matrix4x4<Real>::operator*(const Matrix4x4& m2) const
    {
	Matrix4x4 res;
	unsigned int i, j, k;

	for (i=0; i<4; i++)
	{
	    for (j=0; j<4; j++)
	    {
		res.m[i][j] = Real(0.0);
		for (k=0; k<4; k++)
		{
		    res.m[i][j] += m[i][k] * m2.m[k][j];
		}
	    }
	}

	return res;
    }


    // Input: i (row), j (column)
    template <typename Real>
	Real &Matrix4x4<Real>::operator()(unsigned int i, unsigned int j)
    {
	return m[i][j];
    }


    template <typename Real>
	Matrix4x4<Real> Matrix4x4<Real>::transpose() const
    {
	int i, j;
	Matrix4x4 a = Matrix4x4();  // the transpose

	for (i=0; i<4; i++)
	{
	    for (j=0; j<4; j++)
	    {
		a.m[i][j] = m[j][i];
	    }
	}

	return a;
    }




    // Static methods

    template <typename Real>
	Matrix4x4<Real> Matrix4x4<Real>::scale(Real sx, Real sy, Real sz)
    {
	Real m[4][4];

	m[0][0] = sx; m[0][1] = 0; m[0][2] = 0; m[0][3] = 0;
	m[1][0] = 0; m[1][1] = sy; m[1][2] = 0; m[1][3] = 0;
	m[2][0] = 0; m[2][1] = 0; m[2][2] = sz; m[2][3] = 0;
	m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 1;

	return Matrix4x4(m);
    }


    template <typename Real>
	Matrix4x4<Real> Matrix4x4<Real>::scale(Real scale)
    {
	Real m[4][4];

	m[0][0] = scale; m[0][1] = 0; m[0][2] = 0; m[0][3] = 0;
	m[1][0] = 0; m[1][1] = scale; m[1][2] = 0; m[1][3] = 0;
	m[2][0] = 0; m[2][1] = 0; m[2][2] = scale; m[2][3] = 0;
	m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 1;

	return Matrix4x4(m);
    }


    template <typename Real>
	Matrix4x4<Real> Matrix4x4<Real>::rotationXYZ(Real rx, Real ry, Real rz)
    {
	Real m[4][4];

	m[0][0] = cos(ry) * cos(rz); 
	m[0][1] = cos(ry) * sin(rz); 
	m[0][2] = -sin(ry); 
	m[0][3] = 0;
	m[1][0] = sin(rx) * sin(ry) * cos(rz) - cos(rx) * sin(rz); 
	m[1][1] = sin(rx) * sin(ry) * sin(rz) + cos(rx) * cos(rz); 
	m[1][2] = cos(ry) * sin(rx); 
	m[1][3] = 0;
	m[2][0] = cos(rx) * sin(ry) * cos(rz) + sin(rx) * sin(rz);
	m[2][1] = cos(rx) * sin(ry) * sin(rz) - sin(rx) * cos(rz); 
	m[2][2] = cos(ry) * cos(rx); 
	m[2][3] = 0;
	m[3][0] = 0; 
	m[3][1] = 0; 
	m[3][2] = 0; 
	m[3][3] = 1;

	return Matrix4x4(m);
    }

    template <typename Real>
	Matrix4x4<Real> Matrix4x4<Real>::translation(Real tx, Real ty, Real tz)
    {
	Real m[4][4];

	m[0][0] = 1; m[0][1] = 0; m[0][2] = 0; m[0][3] = tx;
	m[1][0] = 0; m[1][1] = 1; m[1][2] = 0; m[1][3] = ty;
	m[2][0] = 0; m[2][1] = 0; m[2][2] = 1; m[2][3] = tz;
	m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 1;

	return Matrix4x4(m);
    }


}
#endif
