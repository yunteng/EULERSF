/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _arrays_array3d_h
#define _arrays_array3d_h

#include <stdlib.h>
#include <fstream>
#include <iostream>

#include <Core/Exception/HeaderFiles/DefaultException.h>


namespace Arrays
{

    template<class Data> 
    class Array3D
    {
    public:

	class Iterator
	{
	public:
	    Iterator(Array3D<Data> *parent=NULL)
		: parent(parent)
	    {
		i = 0;
		x = 0;
		y = 0;
		z = 0;
	    }

	    int getI()
	    {
		return x;
	    }

	    int getJ()
	    {
		return y;
	    }

	    int getK()
	    {
		return z;
	    }

	    Data getValue()
	    {
		return parent->a[i];
	    }

	    bool hasNext()
	    {
		return i < parent->getStorageSize();
	    }

	    Iterator next()
	    {
		i++;

		if (z < parent->dim[2]-1)
		{
		    z++;
		}
		else
		{
		    z = 0;
		    if (y < parent->dim[1]-1)
		    {
			y++;
		    }
		    else
		    {
			y = 0;
			x++;
		    }
		}

		return *this;
	    }

	protected:

	    unsigned int i;
	    int x, y, z;
	    Array3D<Data> *parent;
	};



    protected:
	Data *a;
	unsigned int d1md2;
	unsigned int dim[3];
	unsigned int storageSize;

    protected:

	void copy(const Array3D& array)
	{
	    if ( storageSize != array.storageSize )
	    {
		if (a != NULL)
		{
		    delete[] a;
		}
		a = new Data[array.storageSize];
	    }
	    dim[0] = array.dim[0];
	    dim[1] = array.dim[1];
	    dim[2] = array.dim[2];
	    d1md2 = array.d1md2;
	    storageSize = array.storageSize;
	    memcpy(a, array.a, storageSize*sizeof(Data));
	}

    public:

	template<typename Index>
	Array3D(const Index dim[3], bool initialize, Data initVal)
	{
	    unsigned int i;

	    this->dim[0] = dim[0];
	    this->dim[1] = dim[1];
	    this->dim[2] = dim[2];
	    d1md2= dim[1] * dim[2];
	    storageSize = dim[0] * dim[1] * dim[2];
	    a = new Data[storageSize];

	    if (initialize)
	    {
		for (i=0; i<storageSize; i++)
		{
		    a[i] = initVal;
		}
	    }
	}

	Array3D(const Array3D& array)
	{
	    storageSize = 0;
	    a = NULL;
	    copy(array);
	}

	Array3D()
	{
	    storageSize = 0;
	    a = NULL;
	    dim[0] = dim[1] = dim[2] = 0;
	}

	~Array3D(void)
	{
	    delete[] a;
	}

	void eraseStorage()
	{
	    dim[0] = dim[1] = dim[2] = 0;
	    d1md2 = 0;
	    storageSize = 0;
	    delete[] a;
	    a = NULL;
	}

	Data& operator()(unsigned int i, unsigned int j, unsigned int k)
	{
	    return a[i*d1md2 + j*dim[2] + k];		
	}

	Data operator()(unsigned int i, unsigned int j, unsigned int k) const
	{
	    return a[i*d1md2 + j*dim[2] + k];		
	}

	template<typename UIntType>
	Data& operator()(UIntType i, UIntType j, UIntType k)
	{
	    return a[i*d1md2 + j*dim[2] + k];		
	}

	template<typename UIntType>
	Data operator()(UIntType i, UIntType j, UIntType k) const
	{
	    return a[i*d1md2 + j*dim[2] + k];		
	}

	template<typename UIntType>
	unsigned int computeLinearIndex(UIntType i, UIntType j, UIntType k) const
	{
	    return i*d1md2 + j*dim[2] + k;
	}

	unsigned int computeLinearIndex(unsigned int i, unsigned int j, unsigned int k) const
	{
	    return i*d1md2 + j*dim[2] + k;
	}

	void computeIndex(unsigned int i, unsigned int *x, unsigned int *y, unsigned int *z) const
	{
	    *x = i / d1md2;
	    *y = (i % d1md2) / dim[2];
	    *z = (i % d1md2) % dim[2];
	}

	// as above, but instead uses a linear index
	Data& operator()(unsigned int i)
	{
	    return a[i];		
	}

	Data operator()(unsigned int i) const
	{
	    return a[i];		
	}

	void operator=(const Array3D& array)
	{
	    copy(array);
	}

	bool operator==(const Array3D& array)
	{
	    if (storageSize != array.storageSize)
		return false;

	    for (unsigned int i=0; i<storageSize; i++)
	    {
		if (a[i] != array.a[i])
		    return false;
	    }

	    return true;
	}

	template<typename Index>
	void size(Index dim[3]) const
	{
	    dim[0] = this->dim[0];
	    dim[1] = this->dim[1];
	    dim[2] = this->dim[2];
	}

	const unsigned int *getSize() const
	{
	    return dim;
	}

	unsigned int getStorageSize() const
	{
	    return storageSize;
	}

	void clear()
	{
	    unsigned int i;

	    for (i=0; i<storageSize; i++)
	    {
		a[i] = 0;
	    }
	}


	/*! Only clear a rectangular sub-region of the array.
	*  The end-points are non-inclusive. */
	template<typename Index2>
	void clear(Index2 start[3], Index2 end[3])
	{
	    Index2 x, y, z;

	    for (x=start[0]; x<end[0]; x++)
	    {
		for (y=start[1]; y<end[1]; y++)
		{
		    for (z=start[2]; z<end[2]; z++)
		    {
			(*this)(x,y,z) = 0;
		    }
		}
	    }
	}

	void set(Data v)
	{
	    unsigned int i;

	    for (i=0; i<storageSize; i++)
	    {
		a[i] = v;
	    }
	}

	unsigned int getMemUsed()
	{
	    return storageSize * sizeof(Data) + sizeof(this);
	}

	unsigned int getMemUsage()
	{
	    return getMemUsed();
	}

	void save(const std::string& fileName) const
	{
	    std::ofstream os(fileName.c_str(), ios_base::binary);
	    if (!os) { Core::throwDefaultException(std::string("Unable to open file "+fileName, __FILE__, __LINE__); }
	    save(os);
	    os.close();
	}

	void save(std::ostream& os) const
	{
	    os.write((const char *)dim, 3*sizeof(unsigned int));
	    os.write((const char *)&storageSize, sizeof(unsigned int));
	    if (storageSize!=0 && a!=NULL)
	    {
		os.write((const char *)a, sizeof(Data)*storageSize);
	    }
	}

	void load(std::istream& is)
	{
	    is.read((char *)dim, 3*sizeof(unsigned int));
	    is.read((char *)&storageSize, sizeof(unsigned int));
	    d1md2= dim[1] * dim[2];
	    delete[] a;
	    if (storageSize!=0)
	    {
		a = new Data[storageSize];
		is.read((char *)a, sizeof(Data)*storageSize);
	    }
	    else
	    {
		a = NULL;
	    }
	}

	Iterator iterator()
	{
	    return Iterator(this);
	}


	friend class Iterator;
    };

    namespace std
    {
	template <class Data>
	ostream& operator<<(ostream& os, const Array3D<Data>& array3D)
	{
	    array3D.save(os);
	    return os;
	}

	template <class Data>
	istream& operator>>(istream& is, Array3D<Data>& array3D)
	{
	    array3D.load(is);
	    return is;
	}
    }

}

#endif
