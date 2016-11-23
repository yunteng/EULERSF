/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_svolsaver_dttop10_h_
#define _grids_svolsaver_dttop10_h_

#include <string>
#include <fcntl.h>
#include "SvolSaverInterface.h"
#include <Core/Exception/HeaderFiles/DefaultException.h>

#ifndef WIN32

#ifndef _LARGEFILE64_SOURCE
  #define _LARGEFILE64_SOURCE
#endif

  #include <sys/types.h>
  #include <unistd.h>
  #define _open open
  #define _close close
  #define _read read
  #define _write write

#ifdef __APPLE__
  #define O_LARGEFILE 0
  #define _lseeki64 lseek
#else
  #define _lseeki64 lseek64
#endif

#else
  #include <io.h>
  #include <sys/stat.h>
#endif


namespace Grids
{
    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    class SvolSaverDTTop1 : public SvolSaverInterface<IndexType, RealType, DataType, UIntType>
    {
    public:

	const static std::string idString;
	const static std::string header;
	const static bool saveAAData = false;

    public:

	SvolSaverDTTop1();
	~SvolSaverDTTop1();

	void init(const typename SvolSaverInterface<IndexType, RealType, DataType, UIntType>::InitParameters& initParameters);
	void finalize();
	void flushAll();

	void pushVa3D(DataType& va3D);
	void pushVa2D(UIntType va2D);
	void pushVa1D(UIntType va1D);
	void pushAA1D(UIntType aa1D);
	void pushAA2D(UIntType aa2D);
	void pushAA3D(UIntType aa3D);
	void pushXIndex(IndexType x);
	void pushYIndex(IndexType y);
	void pushZIndex(IndexType z);

	IndexType& xIndexBack() { return bufferXIndex[bufferCountXIndex-1]; }
	IndexType& yIndexBack() { return bufferYIndex[bufferCountYIndex-1]; }
	IndexType& zIndexBack() { return bufferZIndex[bufferCountZIndex-1]; }
	UIntType& aa1DBack() { return bufferAA1D[(bufferCountXIndex>>1)-1]; }
	UIntType& aa2DBack() { return bufferAA2D[(bufferCountYIndex>>1)-1]; }
	UIntType& aa3DBack() { return bufferAA3D[(bufferCountZIndex>>1)-1]; }
	UIntType va1DSize() { return va1DSizeVar; }
	UIntType va2DSize() { return va2DSizeVar; }
	UIntType va3DSize() { return va3DSizeVar; }
	UIntType xIndexSize() { return xIndexSizeVar; }
	UIntType yIndexSize() { return yIndexSizeVar; }
	UIntType zIndexSize() { return zIndexSizeVar; }
	DataType& va3DBack() { return bufferVa3D[bufferCountVa3D-1]; }

	bool pushTwoAA() { return false; }

	template<typename ValType>
	void flush(unsigned int bufferCount, long long offset, ValType *buffer);
	template<typename ValType>
	void push(ValType& val, unsigned int *bufferCount, unsigned int bufferSize, long long *offset, ValType *buffer);


    protected:

	inline bool isBigEndian();
	inline bool isLittleEndian();
	template<typename MyReal>
	inline void byteSwap ( MyReal& real );


    protected:
	bool finalized;
	int file;
	unsigned int bufferCountVa1D, bufferCountVa2D, bufferCountVa3D, bufferCountAA1D, bufferCountAA2D, bufferCountAA3D, bufferCountXIndex, bufferCountYIndex, bufferCountZIndex;
	unsigned int bufferSizeVa1D, bufferSizeVa2D, bufferSizeVa3D, bufferSizeAA1D, bufferSizeAA2D, bufferSizeAA3D, bufferSizeXIndex, bufferSizeYIndex, bufferSizeZIndex;
	UIntType *bufferVa1D, *bufferVa2D, *bufferAA1D, *bufferAA2D, *bufferAA3D;
	UIntType xIndexSizeVar, yIndexSizeVar, zIndexSizeVar, va1DSizeVar, va2DSizeVar, va3DSizeVar;
	DataType *bufferVa3D;
	IndexType *bufferXIndex, *bufferYIndex, *bufferZIndex;
	// current offsets into file
	long long offsetVa1D, offsetVa2D, offsetVa3D, offsetAA1D, offsetAA2D, offsetAA3D, offsetXIndex, offsetYIndex, offsetZIndex;
    };



    //////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION
    //////////////////////////////////////////////////////////////////////////


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    const std::string SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::idString = "DTG_v1.0  ";
    
    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    const std::string SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::header = "DTG v1.0  ";


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::SvolSaverDTTop1()
    {
	bufferSizeVa1D = 1024*1024 * sizeof(UIntType);
	bufferSizeVa2D = 1024*1024 * sizeof(UIntType);
	bufferSizeVa3D = 1024*1024 * sizeof(RealType);
	bufferSizeAA1D = 1024*1024 * sizeof(UIntType);
	bufferSizeAA2D = 1024*1024 * sizeof(UIntType);
	bufferSizeAA3D = 1024*1024 * sizeof(UIntType);
	bufferSizeXIndex = 1024*1024 * sizeof(IndexType);
	bufferSizeYIndex = 1024*1024 * sizeof(IndexType);
	bufferSizeZIndex = 1024*1024 * sizeof(IndexType);
	bufferCountVa1D = 0;
	bufferCountVa2D = 0;
	bufferCountVa3D = 0;
	bufferCountAA1D = 0;
	bufferCountAA2D = 0;
	bufferCountAA3D = 0;
	bufferCountXIndex = 0;
	bufferCountYIndex = 0;
	bufferCountZIndex = 0;
	bufferVa1D = new UIntType[(size_t)bufferSizeVa1D];
	bufferVa2D = new UIntType[(size_t)bufferSizeVa2D];
	bufferVa3D = new DataType[(size_t)bufferSizeVa3D];

	if (saveAAData)
	{
	    bufferAA1D = new UIntType[(size_t)bufferSizeAA1D];
	    bufferAA2D = new UIntType[(size_t)bufferSizeAA2D];
	    bufferAA3D = new UIntType[(size_t)bufferSizeAA3D];
	}
	else
	{
	    bufferAA1D = NULL;
	    bufferAA2D = NULL;
	    bufferAA3D = NULL;
	}

	bufferXIndex = new IndexType[(size_t)bufferSizeXIndex];
	bufferYIndex = new IndexType[(size_t)bufferSizeYIndex];
	bufferZIndex = new IndexType[(size_t)bufferSizeZIndex];
	xIndexSizeVar = 0;
	yIndexSizeVar = 0;
	zIndexSizeVar = 0;
	va1DSizeVar = 0;
	va2DSizeVar = 0;
	va3DSizeVar = 0;
	finalized = false;
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::~SvolSaverDTTop1()
    {
	if (!finalized)
	{
	    finalize();
	}
	delete[] bufferVa1D;
	delete[] bufferVa2D;
	delete[] bufferVa3D;
	delete[] bufferAA1D;
	delete[] bufferAA2D;
	delete[] bufferAA3D;
	delete[] bufferXIndex;
	delete[] bufferYIndex;
	delete[] bufferZIndex;
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    bool SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::isBigEndian() 
    {
	short int word = 0x001;
	char *byte = reinterpret_cast<char *> ( &word );
	return byte[0] == 0;
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    bool SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::isLittleEndian() 
    {
	short int word = 0x001;
	char *byte = reinterpret_cast<char *> ( &word );
	return byte[0] != 0;
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    template<typename MyReal>
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::byteSwap ( MyReal& real ) 
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
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::finalize()
    {
	flush(bufferCountVa1D, offsetVa1D, bufferVa1D);
	flush(bufferCountVa2D, offsetVa2D, bufferVa2D);
	flush(bufferCountVa3D, offsetVa3D, bufferVa3D);
	if (saveAAData)
	{
	    flush(bufferCountAA1D, offsetAA1D, bufferAA1D);
	    flush(bufferCountAA2D, offsetAA2D, bufferAA2D);
	    flush(bufferCountAA3D, offsetAA3D, bufferAA3D);
	}
	flush(bufferCountXIndex, offsetXIndex, bufferXIndex);
	flush(bufferCountYIndex, offsetYIndex, bufferYIndex);
	flush(bufferCountZIndex, offsetZIndex, bufferZIndex);
	_close(file);
	finalized = true;
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::flushAll()
    {
	flush(bufferCountVa1D, offsetVa1D, bufferVa1D);
	flush(bufferCountVa2D, offsetVa2D, bufferVa2D);
	flush(bufferCountVa3D, offsetVa3D, bufferVa3D);
	if (saveAAData)
	{
	    flush(bufferCountAA1D, offsetAA1D, bufferAA1D);
	    flush(bufferCountAA2D, offsetAA2D, bufferAA2D);
	    flush(bufferCountAA3D, offsetAA3D, bufferAA3D);
	}
	flush(bufferCountXIndex, offsetXIndex, bufferXIndex);
	flush(bufferCountYIndex, offsetYIndex, bufferYIndex);
	flush(bufferCountZIndex, offsetZIndex, bufferZIndex);
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::init(const typename SvolSaverInterface<IndexType, RealType, DataType, UIntType>::InitParameters& initParameters)
    {
	// Alternative: _sopen_s
#ifdef WIN32	
	file = _open(initParameters.fileName.c_str(), _O_WRONLY | _O_BINARY | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE );
#else
	file = _open(initParameters.fileName.c_str(), O_LARGEFILE | O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );
#endif	

	if (file == -1)
	{
	    Core::throwDefaultException(std::string("SvolSaverDTTop1: Unable to open file ") + initParameters.fileName + std::string(" for writing."), __FILE__, __LINE__ );
	}

	_write(file, (const void *)header.c_str(), 10);

	_write(file, (const void *)"IndexTypeTextDescription ", 25);
	if ( sizeof(IndexType)==sizeof(char) )
	{
	    _write(file, (const void *)"char\n", 5);
	}
	else if ( sizeof(IndexType)==sizeof(short) )
	{
	    _write(file, (const void *)"short\n", 6);
	}
	else if ( sizeof(IndexType)==sizeof(int) )
	{
	    _write(file, (const void *)"int\n", 4);
	}
	else if ( sizeof(IndexType)==sizeof(long long) )
	{
	    _write(file, (const void *)"long long\n", 10);
	}
	_write(file, (const void *)"IndexTypeSizeInBytes ", 21);
	if ( sizeof(IndexType)==sizeof(char) )
	{
	    std::stringstream ss;
	    ss << sizeof(char) << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}
	else if ( sizeof(IndexType)==sizeof(int) )
	{
	    std::stringstream ss;
	    ss << sizeof(int) << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}
	else if ( sizeof(IndexType)==sizeof(short) )
	{
	    std::stringstream ss;
	    ss << sizeof(short) << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}
	else if ( sizeof(IndexType)==sizeof(long long) )
	{
	    std::stringstream ss;
	    ss << sizeof(long long) << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}

	_write(file, (const void *)"IntTypeTextDescription ", 23);
	if ( sizeof(UIntType)==sizeof(char) )
	{
	    _write(file, (const void *)"char\n", 5);
	}
	else if ( sizeof(UIntType)==sizeof(short) )
	{
	    _write(file, (const void *)"short\n", 6);
	}
	else if ( sizeof(UIntType)==sizeof(int) )
	{
	    _write(file, (const void *)"int\n", 4);
	}
	else if ( sizeof(UIntType)==sizeof(long long) )
	{
	    _write(file, (const void *)"long long\n", 10);
	}
	_write(file, (const void *)"IntTypeSizeInBytes ", 19);
	if ( sizeof(UIntType)==sizeof(char) )
	{
	    std::stringstream ss;
	    ss << sizeof(char) << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}
	else if ( sizeof(UIntType)==sizeof(int) )
	{
	    std::stringstream ss;
	    ss << sizeof(int) << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}
	else if ( sizeof(UIntType)==sizeof(short) )
	{
	    std::stringstream ss;
	    ss << sizeof(short) << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}
	else if ( sizeof(UIntType)==sizeof(long long) )
	{
	    std::stringstream ss;
	    ss << sizeof(long long) << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}

	_write(file, (const void *)"RealTypeTextDescription ", 24);
	if ( sizeof(RealType)==sizeof(float) )
	{
	    _write(file, (const void *)"float\n", 6);
	}
	else if ( sizeof(RealType)==sizeof(double) )
	{
	    _write(file, (const void *)"double\n", 7);
	}
	else
	{
	    _write(file, (const void *)"custom\n", 7);
	}
	_write(file, (const void *)"RealTypeSizeInBytes ", 20);
	{
	    std::stringstream ss;
	    ss << sizeof(DataType) << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}

	if (isLittleEndian())
	{
	    _write(file, (const void *)"Endian LittleEndian\n", 20);
	}
	else
	{
	    _write(file, (const void *)"Endian BigEndian\n", 17);
	}

	
	_write(file, (const char *)initParameters.bbox, sizeof(IndexType)*6);
	_write(file, (const char *)initParameters.translation, sizeof(RealType)*3);
	_write(file, (const char *)initParameters.rotation, sizeof(RealType)*3);

	_write(file, (const char *)&initParameters.numXIndex, sizeof(UIntType));
	//offsetXIndex = _telli64(file);
	offsetXIndex = _lseeki64(file, 0, SEEK_CUR);
	offsetAA1D = offsetXIndex + initParameters.numXIndex * sizeof(IndexType);
	if (saveAAData)
	{
	    offsetVa1D = offsetAA1D + (initParameters.numXIndex>>1) * sizeof(UIntType);
	}
	else
	{
	    offsetVa1D = offsetAA1D;
	}
	_lseeki64(file, offsetVa1D, SEEK_SET);
	_write(file, (const void *)&initParameters.numVa1D, sizeof(UIntType));
	offsetVa1D += sizeof(UIntType);

	offsetYIndex = offsetVa1D + initParameters.numVa1D * sizeof(UIntType);
	_lseeki64(file, offsetYIndex, SEEK_SET);
	_write(file, (const char *)&initParameters.numYIndex, sizeof(UIntType));
	offsetYIndex += sizeof(UIntType);
	offsetAA2D = offsetYIndex + initParameters.numYIndex * sizeof(IndexType);
	if (saveAAData)
	{
	    offsetVa2D = offsetAA2D + (initParameters.numYIndex>>1) * sizeof(UIntType);
	}
	else
	{
	    offsetVa2D = offsetAA2D;
	}
	_lseeki64(file, offsetVa2D, SEEK_SET);
	_write(file, (const void *)&initParameters.numVa2D, sizeof(UIntType));
	offsetVa2D += sizeof(UIntType);

	offsetZIndex = offsetVa2D + initParameters.numVa2D * sizeof(UIntType);
	_lseeki64(file, offsetZIndex, SEEK_SET);
	_write(file, (const char *)&initParameters.numZIndex, sizeof(UIntType));
	offsetZIndex += sizeof(UIntType);
	offsetAA3D = offsetZIndex + initParameters.numZIndex * sizeof(IndexType);
	if (saveAAData)
	{
	    offsetVa3D = offsetAA3D + (initParameters.numZIndex>>1) * sizeof(UIntType);
	}
	else
	{
	    offsetVa3D = offsetAA3D;
	}
	_lseeki64(file, offsetVa3D, SEEK_SET);
	_write(file, (const void *)&initParameters.numVa3D, sizeof(UIntType));

	_write(file, (const void *)"Auxiliary Fields v1.0\n", 22);
	_write(file, (const void *)"NumberOfScalarFields 1\n", 23);
	_write(file, (const void *)"Header DTGridTopology Auxiliary Field v1.0\n", 43);
	_write(file, (const void *)"DataTextDescription SignedDistanceField\n", 40);
	if (isLittleEndian())
	{
	    _write(file, (const void *)"Endian LittleEndian\n", 20);
	}
	else
	{
	    _write(file, (const void *)"Endian BigEndian\n", 17);
	}
	_write(file, (const void *)"DataTypeTextDescription ", 24);
	if ( sizeof(DataType)==sizeof(float) )
	{
	    _write(file, (const void *)"float\n", 6);
	}
	else if ( sizeof(DataType)==sizeof(double) )
	{
	    _write(file, (const void *)"double\n", 7);
	}
	else
	{
	    _write(file, (const void *)"custom\n", 7);
	}
	_write(file, (const void *)"DataTypeSizeInBytes ", 20);
	{
	    std::stringstream ss;
	    ss << sizeof(DataType) << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}

	_write(file, (const void *)"NumberOfDataItems ", 18);
	{
	    std::stringstream ss;
	    ss << initParameters.numVa3D << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}

	{
	    std::stringstream ss;
	    ss << "NumberOfConstants 5" << std::endl;
	    ss << "Constant dx " << initParameters.dx << std::endl;
	    ss << "Constant gamma " << initParameters.gamma << std::endl;
	    ss << "Constant beta " << initParameters.beta << std::endl;
	    ss << "Constant insideConstant " << initParameters.insideConstant << std::endl;
	    ss << "Constant outsideConstant " << initParameters.outsideConstant << std::endl;
	    _write(file, (const void *)ss.str().c_str(), (unsigned int)ss.str().length());
	}

	//offsetVa3D = _telli64(file);
	offsetVa3D = _lseeki64(file, 0, SEEK_CUR);
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    template<typename ValType>
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::flush(unsigned int bufferCount, long long offset, ValType *buffer)
    {
	_lseeki64(file, offset, SEEK_SET);
	_write(file, (const void *)buffer, (unsigned int)(bufferCount*sizeof(ValType)));
    }



    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    template<typename ValType>
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::push(ValType& val, unsigned int *bufferCount, unsigned int bufferSize, long long *offset, ValType *buffer)
    {
	if (*bufferCount == bufferSize)
	{
	    _lseeki64(file, *offset, SEEK_SET);
	    _write(file, (const void *)buffer, (unsigned int)(bufferSize*sizeof(ValType)));
	    *bufferCount = 0;
	    *offset += bufferSize*sizeof(ValType);
	}

	buffer[*bufferCount] = val;
	(*bufferCount)++;
    }


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::pushVa3D(DataType& va3D)
    {
	push(va3D, &bufferCountVa3D, bufferSizeVa3D, &offsetVa3D, bufferVa3D);
	va3DSizeVar++;
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::pushVa2D(UIntType va2D)
    {
	push(va2D, &bufferCountVa2D, bufferSizeVa2D, &offsetVa2D, bufferVa2D);
	va2DSizeVar++;
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::pushVa1D(UIntType va1D)
    {
	push(va1D, &bufferCountVa1D, bufferSizeVa1D, &offsetVa1D, bufferVa1D);
	va1DSizeVar++;
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::pushAA1D(UIntType aa1D)
    {
	if (saveAAData)
	    push(aa1D, &bufferCountAA1D, bufferSizeAA1D, &offsetAA1D, bufferAA1D);
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::pushAA2D(UIntType aa2D)
    {
	if (saveAAData)
	    push(aa2D, &bufferCountAA2D, bufferSizeAA2D, &offsetAA2D, bufferAA2D);
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::pushAA3D(UIntType aa3D)
    {
	if (saveAAData)
	    push(aa3D, &bufferCountAA3D, bufferSizeAA3D, &offsetAA3D, bufferAA3D);
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::pushXIndex(IndexType x)
    {
	push(x, &bufferCountXIndex, bufferSizeXIndex, &offsetXIndex, bufferXIndex);
	xIndexSizeVar++;
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::pushYIndex(IndexType y)
    {
	push(y, &bufferCountYIndex, bufferSizeYIndex, &offsetYIndex, bufferYIndex);
	yIndexSizeVar++;
    }

    template<typename IndexType, typename RealType, typename DataType, typename UIntType>							  
    void SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::pushZIndex(IndexType z)
    {
	push(z, &bufferCountZIndex, bufferSizeZIndex, &offsetZIndex, bufferZIndex);
	zIndexSizeVar++;
    }

}

#endif
