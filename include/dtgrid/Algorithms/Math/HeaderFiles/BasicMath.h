/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _math_basicmath_h
#define _math_basicmath_h

#define _USE_MATH_DEFINES   // M_PI etc in math.h
#include <stdlib.h>
#include <math.h>


namespace Math
{

   template<typename Type>
   struct MathToolsBinSearchDummyTransform
   {
       Type operator()(const Type& t) const { return t; }
   };

   template<typename Array, typename Type>
	int binarySearch(const Array& a, int start, int end, Type t)
    {
	MathToolsBinSearchDummyTransform<Type> transform;
	return binarySearch(a, start, end, t, transform);
    }

    /*! \brief Logarithmic search in the array, a,  delimited by the indices [start;end). 
    *  If element t exists in a, the index of t in a is returned. 
    *  If t does not exist in a, and t lies within the range of the elements in a,
    *  the index of the smallest element that is larger
    *  than t is returned.
    *  Otherwise the index of the closest boundary (start or end) is returned.
    */
    template<typename Array, typename Transform, typename Type>
	int binarySearch(const Array& a, int start, int end, Type t, const Transform& transform)
    {
	int in1, in2, it;	

	if ( start == end )
	{
	    return end;
	}
	if ( transform(a[end-1]) < t )
	{
	    return end;
	}
	if ( transform(a[start]) > t )
	{
	    return start;
	}


	in1 = 0;
	in2 = end-start-1;

	// binary search
	while (in2 > in1)
	{
	    it = (in1+in2)>>1;

	    if (t < transform(a[it+start]))
	    {
		in2 = it-1;
	    }
	    else if (t > transform(a[it+start]))
	    {
		in1 = it+1;
	    }
	    else
	    {
		in1 = it;
		break;
	    }
	}

	if (t > transform(a[in1+start]))
	{
	    // this is safe as we know that t is within the limits of the array
	    in1++;
	}


	return in1+start;
    }


    inline unsigned int getSmallestPowerOfTwoLargerThanOrEqualTo(unsigned int x)
    {
	unsigned int xOrig = x;
	unsigned int xLog = 0;
	unsigned int y;

	x = (x >> 1);
	while(x > 0) { x = (x >> 1); xLog++; }
	y = (1<<xLog);
	if ( y < xOrig ) { y = (y<<1); }
	return y;
    }



    template<typename Scalar>
	inline Scalar computeOrientation(Scalar v1x, Scalar v1y, Scalar v2x, Scalar v2y)
    {
	Scalar value = v1x * v2y - v1y * v2x;

	return (Scalar)(sign3(value));
    }


    template<typename Scalar>
	inline int sign3(Scalar s)
    {
	return ( s<0 ? -1 : ( s>0 ? 1 : 0 ) );
    }



    template<typename Scalar>
	inline Scalar sign(Scalar s)
    {
	return ((s)<=0 ? Scalar(-1) : Scalar(1));
    }


    template<typename Scalar>
	inline Scalar fabs2(Scalar a, int *sign)
    {
	if (a<0)
	{
	    *sign = -1;
	    return -a;
	}
	else
	{
	    *sign = 1;
	    return a;
	}
    }

    template<typename Scalar>
	inline Scalar pow2(Scalar s)
    {
	return s*s;
    }

    template<typename Scalar>
	inline Scalar pow3(Scalar s)
    {
	return s*s*s;
    }

    template<typename Scalar>
	inline Scalar log2(Scalar s)
    {
	return log10(s) * ( 1.0 / 0.30102999566398 );  //   log10(s) / log10(2)
    }


    template<typename Scalar>
    inline double dist2D(Scalar p1[2], Scalar p2[2])
    {
	return sqrt( (double)( pow2(p2[0]-p1[0]) + pow2(p2[1]-p1[1]) ) );
    }


    template<typename Scalar>
    inline double dist2DSquared(Scalar p1[2], Scalar p2[2])
    {
	return (double)(pow2(p2[0]-p1[0]) + pow2(p2[1]-p1[1]));
    }


    template<typename Scalar> 
	inline bool solveQuadraticEquation(Scalar A, Scalar B, Scalar C, Scalar *alpha, Scalar *beta)
    {
	Scalar D;

	D = pow2(B) - 4*A*C;

	if ( D >= 0 )
	{
	    *alpha = (-B + sqrt(D))/(2*A);
	    *beta = (-B - sqrt(D))/(2*A);
	    return true;
	}
	else
	{
	    return false;
	}
    }


    inline bool isEven(int i)
    {
	return i%2 == 0; 
    }


    inline double rand0to1()
    {
	return (double)rand()/RAND_MAX;
    }


    template<typename Scalar>
	inline Scalar length(Scalar x, Scalar y, Scalar z)
    {
	return sqrt( pow2(x) + pow2(y) + pow2(z) );
    }


    template<typename Scalar>
	inline Scalar length2(Scalar x, Scalar y, Scalar z)
    {
	return pow2(x) + pow2(y) + pow2(z);
    }


    template<typename Scalar>
	inline void normalize(Scalar *x, Scalar *y, Scalar *z)
    {
	Scalar len = length(*x, *y, *z);

	*x = *x / len;
	*y = *y / len;
	*z = *z / len;
    }

    template<typename Scalar>
	inline void normalize(Scalar n[3])
    {
	Scalar len = length(n[0], n[1], n[2]);

	n[0] = n[0] / len;
	n[1] = n[1] / len;
	n[2] = n[2] / len;
    }

    template<typename Scalar>
	inline Scalar fround(Scalar s)
    {
	return floor(s+Scalar(0.5));
    }

    template<typename Scalar>
	inline void crossVec(Scalar x1, Scalar y1, Scalar z1, Scalar x2, Scalar y2, Scalar z2, Scalar *xn, Scalar *yn, Scalar *zn)
    {
	*xn = y1 * z2 - z1 * y2;
	*yn = z1 * x2 - x1 * z2;
	*zn = x1 * y2 - y1 * x2;
    }


}

#endif
