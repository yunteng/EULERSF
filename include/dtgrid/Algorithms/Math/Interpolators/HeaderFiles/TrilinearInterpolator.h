/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _math_trilinearinterpolator_h
#define _math_trilinearinterpolator_h

#include "DataStructures/Vectors/HeaderFiles/Vector3.h"

namespace Math
{

    template<typename Data, typename Real>
    class TrilinearInterpolator
    {
    public:

	// Input:
	//
	// 'values' array, ordering:
	//     X:         X + 1:
	// (and Y to the right, Z up) 
	//
	//  1 -- 3       5 -- 7
	//  |    |       |    |
	//  0 -- 2       4 -- 6
	//
	// unitOffset:
	// offset in "unit cube coordinates" of point to be interpolated
	//
	// incorporate more tests here!! should not be negative. take fabs and max of value and zero
	static Data interp(Data values[8], Matrix::Vector3<Real> unitOffset)
	{
	    return 
		values[0] * (Real(1.0)-unitOffset[0]) * (Real(1.0)-unitOffset[1]) * (Real(1.0)-unitOffset[2]) +
		values[4] * unitOffset[0] * (Real(1.0)-unitOffset[1]) * (Real(1.0)-unitOffset[2]) +
		values[2] * (Real(1.0)-unitOffset[0]) * unitOffset[1] * (Real(1.0)-unitOffset[2]) +
		values[1] * (Real(1.0)-unitOffset[0]) * (Real(1.0)-unitOffset[1]) * unitOffset[2] +
		values[5] * unitOffset[0] * (Real(1.0)-unitOffset[1]) * unitOffset[2] +
		values[3] * (Real(1.0)-unitOffset[0]) * unitOffset[1] * unitOffset[2] +
		values[6] * unitOffset[0] * unitOffset[1] * (Real(1.0)-unitOffset[2]) +
		values[7] * unitOffset[0] * unitOffset[1] * unitOffset[2];
	}



	// Input:
	//
	// 'values' array, ordering:
	//     X:         X + 1:
	// (and Y to the right, Z up) 
	//
	//  1 -- 3       5 -- 7
	//  |    |       |    |
	//  0 -- 2       4 -- 6
	//
	// unitOffset:
	// offset in "unit cube coordinates" of point to be interpolated
	//
	// incorporate more tests here!! should not be negative. take fabs and max of value and zero
	static Data interp(Data values[8], Matrix::Vector3<Real> unitOffset, bool mask[8])
	{
	    Data finalVal;
	    Real wSum = 0;
	    Real w[8];

	    if (mask[0]) { w[0] = (Real(1.0)-unitOffset[0]) * (Real(1.0)-unitOffset[1]) * (Real(1.0)-unitOffset[2]); wSum+=w[0]; } else { w[0]=0; } 
	    if (mask[1]) { w[1] = (Real(1.0)-unitOffset[0]) * (Real(1.0)-unitOffset[1]) * unitOffset[2]; wSum+=w[1]; } else { w[1]=0; } 
	    if (mask[2]) { w[2] = (Real(1.0)-unitOffset[0]) * unitOffset[1] * (Real(1.0)-unitOffset[2]); wSum+=w[2]; } else { w[2]=0; } 
	    if (mask[3]) { w[3] = (Real(1.0)-unitOffset[0]) * unitOffset[1] * unitOffset[2]; wSum+=w[3]; } else { w[3]=0; } 
	    if (mask[4]) { w[4] = unitOffset[0] * (Real(1.0)-unitOffset[1]) * (Real(1.0)-unitOffset[2]); wSum+=w[4]; } else { w[4]=0; } 
	    if (mask[5]) { w[5] = unitOffset[0] * (Real(1.0)-unitOffset[1]) * unitOffset[2]; wSum+=w[5]; } else { w[5]=0; } 
	    if (mask[6]) { w[6] = unitOffset[0] * unitOffset[1] * (Real(1.0)-unitOffset[2]); wSum+=w[6]; } else { w[6]=0; } 
	    if (mask[7]) { w[7] = unitOffset[0] * unitOffset[1] * unitOffset[2]; wSum+=w[7]; } else { w[7]=0; } 

	    if (wSum==0)
		return static_cast<Data>(0);
	    
	    finalVal = values[0] * w[0];
	    finalVal += values[1] * w[1];
	    finalVal += values[2] * w[2];
	    finalVal += values[3] * w[3];
	    finalVal += values[4] * w[4];
	    finalVal += values[5] * w[5];
	    finalVal += values[6] * w[6];
	    finalVal += values[7] * w[7];

	    finalVal /= wSum;

	    return finalVal;
	}


    };

}

#endif
