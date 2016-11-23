/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _vector3_h
#define _vector3_h

#include "Vector4.h"
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>


namespace Matrix
{

    template <typename Real> class Vector2;

    template <typename Real>
    class Vector3
    {
    public:

	typedef Real ComponentType;


    public:
	Vector3() { v[0] = v[1] = v[2] = Real(0.0); }
	Vector3(const Real values[3]) {v[0] = values[0]; v[1] = values[1]; v[2] = values[2]; }
	Vector3(Real v0) { v[0] = v0; v[1] = v0; v[2] = v0; }
	template <typename Real2> Vector3(Real2 v0) { v[0] = static_cast<Real>(v0); v[1] = static_cast<Real>(v0); v[2] = static_cast<Real>(v0); }
	Vector3(Real v0, Real v1, Real v2) { v[0] = v0; v[1] = v1; v[2] = v2; }
	Vector3(const Vector4<Real>& vec4) { v[0]=vec4[0]/vec4[3]; v[1]=vec4[1]/vec4[3]; v[2]=vec4[2]/vec4[3]; }
	Vector3(const Vector2<Real>& vec2, Real v2) { v[0]=vec2[0]; v[1]=vec2[1]; v[2] = v2; }
	template <typename Real2> Vector3(const Vector3<Real2>& vec3) { *this = vec3; }
	Vector3(std::string str) { std::stringstream s; s << str; load(s); } 
	~Vector3() { }

	void load(std::istream& is) { is >> v[0] >> v[1] >> v[2]; }
	void save(std::ostream& os) const { os << v[0] << " " << v[1] << " " << v[2]; }

	Real& operator[](unsigned int i) { return v[i]; }
	Real operator[](unsigned int i) const { return v[i]; }
	Real& operator()(unsigned int i) { return v[i]; }
	Real operator()(unsigned int i) const { return v[i]; }

	Vector3 operator-(const Vector3& vec3) const { return Vector3(v[0]-vec3.v[0],v[1]-vec3.v[1],v[2]-vec3.v[2]); }
	Vector3 operator+(const Vector3& vec3) const { return Vector3(v[0]+vec3.v[0],v[1]+vec3.v[1],v[2]+vec3.v[2]); }

	Vector3& operator+=(const Vector3& vec3) {
		v[0] += vec3[0]; 
		v[1] += vec3[1]; 
		v[2] += vec3[2];
		return *this; 
	}

	template<typename Scalar> Vector3& operator/=(Scalar s) {
		v[0] /= s; 
		v[1] /= s; 
		v[2] /= s;
		return *this; 
	}

	Real length() const { return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ); }
	Real lengthSquared() const { return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; }
	Real operator*(const Vector3<Real>& vec3) const { return v[0]*vec3[0] + v[1]*vec3[1] + v[2]*vec3[2]; }
	Vector3& normalize() { Real len = length(); if (len!=0) { v[0] /= len; v[1] /= len; v[2] /= len; } return *this; }

	template <typename Real2> Vector3 operator=(const Vector3<Real2>& vec3) { v[0]=(Real)vec3[0]; v[1]=(Real)vec3[1]; v[2]=(Real)vec3[2]; return *this; }
	template <typename Real2> bool operator==(const Vector3<Real2>& vec3) const { return v[0]==(Real)vec3[0] && v[1]==(Real)vec3[1] && v[2]==(Real)vec3[2]; }

	template<typename Scalar> Vector3<Real> operator*(Scalar s) const { return Vector3<Real>((Real)(s*v[0]), (Real)(s*v[1]), (Real)(s*v[2])); }
	template<typename Scalar> Vector3<Real> operator/(Scalar s) const { return Vector3<Real>((Real)(v[0]/s), (Real)(v[1]/s), (Real)(v[2]/s)); }
	Vector3 entryMult(const Vector3& vec3) const { return Vector3(v[0]*vec3.v[0],v[1]*vec3.v[1],v[2]*vec3.v[2]); }

	bool operator<( const Vector3<Real>& vec3 ) const //lexicographic
	{
		if ( v[0]<vec3[0] ) 
			return true;
		else if ( v[0]>vec3[0] )
			return false;
		else if ( v[1]<vec3[1] )
			return true;
		else if ( v[1]>vec3[1] )
			return false;
		else if ( v[2]<vec3[2] )
			return true;
		else
			return false;
	}

	Real *getArrayPtr() { return v; }

	// static methods

	static Vector3 floor(const Vector3& vec3) { return Vector3(::floor(vec3[0]),::floor(vec3[1]),::floor(vec3[2])); }
	static Vector3 ceil(const Vector3& vec3) { return Vector3(::ceil(vec3[0]),::ceil(vec3[1]),::ceil(vec3[2])); }
	static Vector3 round(const Vector3& vec3) { return Vector3(::round(vec3[0]),::round(vec3[1]),::round(vec3[2]));}

	// conversion operators

	Real *conversionToRealPtr() { return v; }

    protected:
	Real v[3];


	// friends
	template<typename Real2> friend Vector3<Real2> cross(const Vector3<Real2>& v1, const Vector3<Real2>& v2);
    };


    template <typename Real2>
	Vector3<Real2> cross(const Vector3<Real2>& v1, const Vector3<Real2>& v2)
    {
	Vector3<Real2> res;

	res[0] = v1[1] * v2[2] - v1[2] * v2[1];
	res[1] = v1[2] * v2[0] - v1[0] * v2[2];
	res[2] = v1[0] * v2[1] - v1[1] * v2[0];

	return res;
    }


};


namespace std
{
    template <typename Real1, typename Real2>
	Matrix::Vector3<Real2> operator*(Real1 r, const Matrix::Vector3<Real2>& v)
    {
	Matrix::Vector3<Real2> res;

	res[0] = v[0]*r;
	res[1] = v[1]*r;
	res[2] = v[2]*r;

	return res;
    }


    template <typename Real>
	ostream& operator<<(ostream& os, const Matrix::Vector3<Real>& m)
    {
	m.save(os);
	return os;
    }


    template <typename Real>
	istream& operator>>(istream& is, Matrix::Vector3<Real>& m)
    {
	m.load(is);
	return is;
    }

}


#endif
