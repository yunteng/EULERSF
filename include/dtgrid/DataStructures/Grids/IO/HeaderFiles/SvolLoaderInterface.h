/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_svolloaderinterface_h_
#define _grids_svolloaderinterface_h_

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <Core/Exception/HeaderFiles/DefaultException.h>

namespace Grids
{

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    class SvolLoaderInterface
    {
    public:

	class InitParameters
	{
	public:
	    UIntType numIndexEndMarkers;
	    UIntType numAAEndMarkers;
	    UIntType numValueEndMarkers;
	    std::string fileName;
	};

	class LoadOutput
	{
	public:

	    UIntType *va1D;         // indices that point to the start of each (x)-column in yIndex
	    UIntType numVa1D;    
	    IndexType *xIndex;              // start and end x indices of the connected components contained in the projection of the narrow band onto the x axis
	    UIntType numXIndex; 
	    UIntType lastXIndex;
	    UIntType *aa1D;         // points into va1D. gives the start (and in this case also end) index of each connected component in the x direction

	    UIntType *va2D;         // indices that point to the start of each (x,y)-column in zIndex
	    UIntType numVa2D;
	    IndexType *yIndex;              // start and end y indices of the connected components in the y direction contained in the projection of the narrow band onto the (x,y) plane
	    UIntType numYIndex;
	    UIntType lastYIndex;
	    UIntType *aa2D;         // points into va2D. gives the start (and in this case also end) index of each connected component in the y direction

	    IndexType *zIndex;              // start and end z index of each connected component in the z direction
	    UIntType numZIndex;
	    UIntType lastZIndex;
	    UIntType *aa3D;         // points into va3D. gives the start (and in this case also end) index of each connected component in the z direction

	    UIntType numVa3D;
	    DataType *va3D;				   // values of grid points in narrow band (will always correspond to values[0])

	    // prefs

	    RealType beta, gamma, dx;
	    DataType insideConstant, outsideConstant;

	    IndexType bbox[3][2];                      // (min,max) for x, y and z in grid coordinates
	    RealType translation[3];
	    RealType rotation[3];

	};


    public:

	virtual ~SvolLoaderInterface() { }

	virtual void init(const typename SvolLoaderInterface<IndexType, RealType, DataType, UIntType>::InitParameters& initParameters) = 0;

	virtual void load(LoadOutput *output) = 0;

    protected:

	inline void reportError(long code, const std::string& file, int line)
	{
	    if (code == -1L)
	    {
		Core::throwDefaultException(strerror(errno), file, line);
	    }
	}

    };
}


#endif
