/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_svolloaderdttop1_h_
#define _grids_svolloaderdttop1_h_

#include <string>
#include <stdlib.h>
#include <Core/Exception/HeaderFiles/DefaultException.h>
#include "SvolLoaderInterface.h"
#include "SvolSaver_DTTop1.0.h"

namespace Grids
{

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    class SvolLoaderDTTop1 : public SvolLoaderInterface<IndexType, RealType, DataType, UIntType>
    {
    public:

	const static std::string idString;
	const static std::string header;


    public:

	void init(const typename SvolLoaderInterface<IndexType, RealType, DataType, UIntType>::InitParameters& initParameters);

	void load(typename SvolLoaderInterface<IndexType, RealType, DataType, UIntType>::LoadOutput *output);


    protected:

	inline bool isBigEndian();
	inline bool isLittleEndian();
	template<typename MyReal>
	inline void byteSwap ( MyReal& real );

	template<typename MyType1, typename MyType2>
	void loadArray(std::istream& is, MyType1 *a, unsigned int numEntries, bool swapEndian);
	template<typename MyIntType>
	void loadIntegerArray(std::istream& is, MyIntType *a, unsigned int numEntries, unsigned int intSize, bool swapEndian);
	template<typename MyRealType>
	void loadRealArray(std::istream& is, MyRealType *a, unsigned int numEntries, unsigned int realSize, bool swapEndian);

	void buildAAData(UIntType *aa, IndexType *index, UIntType numIndex);

    protected:

	UIntType numIndexEndMarkers;
	UIntType numAAEndMarkers;
	UIntType numValueEndMarkers;
	std::string fileName;
    };


    //////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION
    //////////////////////////////////////////////////////////////////////////


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    const std::string SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::idString = "DTG_v1.0  ";

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    const std::string SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::header = "DTG v1.0  ";


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    bool SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::isBigEndian() 
    {
	short int word = 0x001;
	char *byte = reinterpret_cast<char *> ( &word );
	return byte[0] == 0;
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    bool SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::isLittleEndian() 
    {
	short int word = 0x001;
	char *byte = reinterpret_cast<char *> ( &word );
	return byte[0] != 0;
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    template<typename MyReal>
    void SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::byteSwap ( MyReal& real ) 
    {
	unsigned char *ptr = reinterpret_cast<unsigned char *> ( &real );
	register int i = 0;
	register int j = sizeof ( MyReal ) - 1;
	unsigned char tmp;
	while ( i < j ) 
	{
	    tmp = ptr[i];
	    ptr[i] = ptr[j];
	    ptr[j] = tmp;
	    i++; j--;
	}
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    template<typename MyType1, typename MyType2>
    void SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::loadArray(std::istream& is, MyType1 *a, unsigned int numEntries, bool swapEndian)
    {
	unsigned int i;
	MyType2 tmp;

	if (swapEndian)
	{
	    for (i=0; i<numEntries; i++)
	    {
		is.read((char *)&tmp, sizeof(MyType2));
		byteSwap(tmp);
		a[i] = (MyType1)tmp;
	    }
	}
	else
	{
	    for (i=0; i<numEntries; i++)
	    {
		is.read((char *)&tmp, sizeof(MyType2));
		a[i] = (MyType1)tmp;
	    }
	}
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    template<typename MyIntType>
    void SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::loadIntegerArray(std::istream& is, MyIntType *a, unsigned int numEntries, unsigned int intSize, bool swapEndian)
    {
	if (intSize == sizeof(MyIntType) && !swapEndian)
	{
	    is.read((char *)a, sizeof(MyIntType) * numEntries);
	}
	else
	{
	    switch (intSize)
	    {
	    case sizeof(char):
		loadArray<MyIntType, char>(is, a, numEntries, swapEndian);
		break;
	    case sizeof(short):
		loadArray<MyIntType, short>(is, a, numEntries, swapEndian);
		break;
	    case sizeof(int):
		loadArray<MyIntType, int>(is, a, numEntries, swapEndian);
		break;
	    case sizeof(long long):
		loadArray<MyIntType, long long>(is, a, numEntries, swapEndian);
		break;
	    }
	}
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    template<typename MyRealType>
    void SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::loadRealArray(std::istream& is, MyRealType *a, unsigned int numEntries, unsigned int realSize, bool swapEndian)
    {
	if (realSize == sizeof(MyRealType) && !swapEndian)
	{
	    is.read((char *)a, sizeof(MyRealType) * numEntries);
	}
	else
	{
	    switch (realSize)
	    {
	    case sizeof(float):
		loadArray<MyRealType, float>(is, a, numEntries, swapEndian);
		break;
	    case sizeof(double):
		loadArray<MyRealType, double>(is, a, numEntries, swapEndian);
		break;
	    }
	}
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    void SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::buildAAData(UIntType *aa, IndexType *index, UIntType numIndex)
    {
	UIntType count = 0;

	for (UIntType i=0, j=0; i<numIndex; i+=2, j++)
	{
	    aa[j] = count;
	    count += static_cast<UIntType>(index[i+1]-index[i]+1);
	}
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    void SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::init(const typename SvolLoaderInterface<IndexType, RealType, DataType, UIntType>::InitParameters& initParameters)
    {
	numIndexEndMarkers = initParameters.numIndexEndMarkers;
	numAAEndMarkers = initParameters.numAAEndMarkers;
	numValueEndMarkers = initParameters.numValueEndMarkers;
	fileName = initParameters.fileName;
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    void SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::load(typename SvolLoaderInterface<IndexType, RealType, DataType, UIntType>::LoadOutput *output)
    {
	std::ifstream is(fileName.c_str(), ios_base::binary);

	if (!is)
	{
	    Core::throwDefaultException("Unable to open file "+fileName+" for loading!", __FILE__, __LINE__);
	}


	// load topology
	{
	    std::string line;
	    char header[10];
	    unsigned int indexSize, realSize, intSize;
	    bool swapEndian;

	    is.read(header, 10);

	    // From this point we assume that the file is in the correct format

	    getline(is, line);  // ignore: IndexTypeTextDescription
	    is >> line;         // IndexTypeSizeInBytes
	    is >> indexSize;
	    getline(is, line);  // eat the newline

	    getline(is, line);  // ignore: IntTypeTextDescription
	    is >> line;         // IntTypeSizeInBytes
	    is >> intSize;
	    getline(is, line);  // eat the newline

	    getline(is, line);  // ignore: RealTypeTextDescription
	    is >> line;         // RealTypeSizeInBytes
	    is >> realSize;
	    getline(is, line);  // eat the newline

	    is >> line;         // Endian
	    is >> line;
	    if ( (line == "LittleEndian" && isBigEndian()) || (line == "BigEndian" && isLittleEndian()) )
	    {
		swapEndian = true;
	    }
	    else
	    {
		swapEndian = false;
	    }
	    getline(is, line);  // eat newline


	    if (indexSize != sizeof(IndexType))
	    {
		cerr << "Warning: Loading file, but be aware that Index size is " << sizeof(IndexType) << ", Index size in file is " << indexSize << std::endl;
	    }

	    if (intSize != sizeof(unsigned int))

	    {
		cerr << "Warning: Loading file, but be aware that Int size is " << sizeof(unsigned int) << ", Int size in file is " << intSize << std::endl;
	    }

	    if (realSize != sizeof(RealType))
	    {
		cerr << "Warning: Loading file, but be aware that Real size is " << sizeof(RealType) << ", Real size in file is " << realSize << std::endl;
	    }


	    loadIntegerArray(is, (IndexType *)output->bbox, 6, indexSize, swapEndian);
	    loadRealArray(is, (RealType *)output->translation, 3, realSize, swapEndian);
	    loadRealArray(is, (RealType *)output->rotation, 3, realSize, swapEndian);


	    // 1D
	    loadIntegerArray(is, &output->numXIndex, 1, intSize, swapEndian);
	    output->xIndex = new IndexType[output->numXIndex+numIndexEndMarkers];
	    loadIntegerArray(is, output->xIndex, output->numXIndex, indexSize, swapEndian);
	    output->aa1D = new unsigned int[(output->numXIndex>>1)+numAAEndMarkers];
	    if (SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::saveAAData)
	    {
		loadIntegerArray(is, output->aa1D, output->numXIndex>>1, intSize, swapEndian);
	    }
	    else
	    {
		buildAAData(output->aa1D, output->xIndex, output->numXIndex);
	    }
	    loadIntegerArray(is, &output->numVa1D, 1, intSize, swapEndian);
	    output->va1D = new unsigned int[output->numVa1D+numValueEndMarkers];
	    loadIntegerArray(is, output->va1D, output->numVa1D, intSize, swapEndian);

	    // 2D
	    loadIntegerArray(is, &output->numYIndex, 1, intSize, swapEndian);
	    output->yIndex = new IndexType[output->numYIndex+numIndexEndMarkers];
	    loadIntegerArray(is, output->yIndex, output->numYIndex, indexSize, swapEndian);
	    output->aa2D = new unsigned int[(output->numYIndex>>1)+numAAEndMarkers];
	    if (SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::saveAAData)
	    {
		loadIntegerArray(is, output->aa2D, output->numYIndex>>1, intSize, swapEndian);
	    }
	    else
	    {
		buildAAData(output->aa2D, output->yIndex, output->numYIndex);
	    }
	    loadIntegerArray(is, &output->numVa2D, 1, intSize, swapEndian);
	    output->va2D = new unsigned int[output->numVa2D+numValueEndMarkers];
	    loadIntegerArray(is, output->va2D, output->numVa2D, intSize, swapEndian);

	    // 3D
	    loadIntegerArray(is, &output->numZIndex, 1, intSize, swapEndian);
	    output->zIndex = new IndexType[output->numZIndex+numIndexEndMarkers];
	    loadIntegerArray(is, output->zIndex, output->numZIndex, indexSize, swapEndian);
	    output->aa3D = new unsigned int[(output->numZIndex>>1)+numAAEndMarkers];
	    if (SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::saveAAData)
	    {
		loadIntegerArray(is, output->aa3D, output->numZIndex>>1, intSize, swapEndian);
	    }
	    else
	    {
		buildAAData(output->aa3D, output->zIndex, output->numZIndex);
	    }
	    loadIntegerArray(is, &output->numVa3D, 1, intSize, swapEndian);
	    output->va3D = new DataType[output->numVa3D+numValueEndMarkers];
	}


	// load values as a separate scalar field
	{
	    std::string line;
	    unsigned int dataSize, tmpInt, j;
	    bool swapEndian;
	    unsigned int numberOfScalarFields;

	    // Load the scalar fields header
	    getline(is, line);
	    if (line != "Auxiliary Fields v1.0")
	    {
		Core::throwDefaultException("Unknown auxiliary fields format: "+line, __FILE__, __LINE__);
	    }

	    // from this point on we assume the file format is correct

	    is >> line;    // NumberOfScalarFields
	    is >> numberOfScalarFields;
	    getline(is, line);  // eat the newline

	    getline(is, line);  // ignore "Header DTGridTopology Auxiliary Field v1.0"
	    getline(is, line);  // ignore DataTextDescription
	    if (line != "DataTextDescription SignedDistanceField")
	    {
		Core::throwDefaultException("The first scalar field should be 'DataTextDescription SignedDistanceField' but is: "+line, __FILE__, __LINE__);
	    }
	    is >> line;         // Endian
	    is >> line;
	    if ( (line == "LittleEndian" && isBigEndian()) || (line == "BigEndian" && isLittleEndian()) )
	    {
		swapEndian = true;
	    }
	    else
	    {
		swapEndian = false;
	    }
	    getline(is, line);  // eat newline
	    getline(is, line);  // ignore DataTypeTextDescription
	    is >> line;         // DataTypeSizeInBytes
	    is >> dataSize;
	    getline(is, line);  // eat newline
	    getline(is, line);  // ignore NumberOfDataItems, should be the same as the number of values specified by the topology
	    output->dx = output->gamma = output->beta = output->insideConstant = output->outsideConstant = 0;  // default values
	    is >> line;  // NumberOfConstants
	    is >> tmpInt;
	    for (j=0; j<tmpInt; j++)
	    {
		is >> line;  // Constant
		is >> line;  // Name
		if (line == "dx")
		{
		    is >> output->dx;
		}
		else if (line == "gamma")
		{
		    is >> output->gamma;
		}
		else if (line == "beta")
		{
		    is >> output->beta;
		}
		else if (line == "insideConstant")
		{
		    is >> output->insideConstant;
		}
		else if (line == "outsideConstant")
		{
		    is >> output->outsideConstant;
		}

		getline(is, line);  // eat newline
	    }

	    // load numerical data
	    this->loadRealArray(is, output->va3D, output->numVa3D, dataSize, swapEndian);
	}

	output->lastXIndex = output->numXIndex-1;
	output->lastYIndex = output->numYIndex-1;
	output->lastZIndex = output->numZIndex-1;
    }



}
#endif
