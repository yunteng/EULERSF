/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_svolsaverinterface_h_
#define _grids_svolsaverinterface_h_

#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <Core/Exception/HeaderFiles/DefaultException.h>
#include <errno.h>

namespace Grids
{

    /*
     * The SvolSaverInterface allows pushing the various components of a dtgrid data structure in random order.
     * ostream is not supported, since it does not allow for seek operations using 64 bit which is required for large svols.
     *
     */
    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    class SvolSaverInterface
    {
    public:

	static const std::string fileNameExtension;

    public:

	class InitParameters
	{
	public:
	    UIntType numXIndex, numYIndex, numZIndex;
	    UIntType numVa1D, numVa2D, numVa3D;
	    RealType dx, beta, gamma, insideConstant, outsideConstant;
	    RealType translation[3];
	    RealType rotation[3];
	    std::string fileName;
	    IndexType bbox[3][2];

	    void printASCII(std::ostream& os)
	    {
		os << "numXIndex = " << numXIndex << std::endl;
		os << "numYIndex = " << numYIndex << std::endl;
		os << "numZIndex = " << numZIndex << std::endl;
		os << "numVa1D = " << numVa1D << std::endl;
		os << "numVa2D = " << numVa2D << std::endl;
		os << "numVa3D = " << numVa3D << std::endl;
		os << "dx = " << dx << std::endl;
		os << "beta = " << beta << std::endl;
		os << "gamma = " << gamma << std::endl;
		os << "insideConstant = " << insideConstant << std::endl;
		os << "outsideConstant = " << outsideConstant << std::endl;
		os << "fileName = " << fileName << std::endl;
		os << "bbox = " << bbox[0][0] << " " << bbox[0][1] << " " << bbox[1][0] << " " << bbox[1][1] << " " << bbox[2][0] << " " << bbox[2][1] << std::endl;
	    }
	};

    public:

	virtual ~SvolSaverInterface() { }

	virtual void init(const typename SvolSaverInterface<IndexType, RealType, DataType, UIntType>::InitParameters& initParameters) = 0;
	virtual void finalize() = 0;
	virtual void flushAll() = 0;

	virtual void pushVa3D(DataType& va3D) = 0;
	virtual void pushVa2D(UIntType va2D) = 0;
	virtual void pushVa1D(UIntType va1D) = 0;
	virtual void pushAA1D(UIntType aa1D) = 0;
	virtual void pushAA2D(UIntType aa2D) = 0;
	virtual void pushAA3D(UIntType aa3D) = 0;
	virtual void pushXIndex(IndexType x) = 0;
	virtual void pushYIndex(IndexType y) = 0;
	virtual void pushZIndex(IndexType z) = 0;

	virtual IndexType& xIndexBack() = 0;
	virtual IndexType& yIndexBack() = 0;
	virtual IndexType& zIndexBack() = 0;
	virtual UIntType& aa1DBack() = 0;
	virtual UIntType& aa2DBack() = 0;
	virtual UIntType& aa3DBack() = 0;
	virtual UIntType va1DSize() = 0;
	virtual UIntType va2DSize() = 0;
	virtual UIntType va3DSize() = 0;
	virtual UIntType xIndexSize() = 0;
	virtual UIntType yIndexSize() = 0;
	virtual UIntType zIndexSize() = 0;
	virtual DataType& va3DBack() = 0;

	virtual bool pushTwoAA() = 0;


	// debug
	virtual UIntType getNumZIndexFromFile() { return 0; }

	virtual const std::string& getFileNameExtension() { return fileNameExtension; }


    protected:

	inline void reportError(long code, const std::string& file, int line)
	{
	    if (code == -1L)
	    {
#ifdef WIN32		
		Core::throwDefaultException(_sys_errlist[errno], file, line);
#else
		Core::throwDefaultException(sys_errlist[errno], file, line);
#endif		
	    }
	}

    };


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    const std::string SvolSaverInterface<IndexType, RealType, DataType, UIntType>::fileNameExtension = "svol";

}



#endif
