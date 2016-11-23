/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _arrays_array1d_h
#define _arrays_array1d_h

#include <iostream>
#include <fstream>

namespace Arrays
{

    /**
    * This array is based on page tables.
    * Iterating through the array is fast.
    * Pushing a new element is fast.
    * Only lookups are slower since they require an extra page table lookup.
    * Template parameters:
    * - Note that pageSize and pageExp must agree
    */
    template<class Type, unsigned int pageSize=4096, unsigned int pageExp=12, unsigned int pageMask=((1<<12)-1)>
    class Array1D
    {
    public:
	typedef Type *TypePtr;

	// The iterator supports two different interfaces
	class Iterator
	{
	public:
	    typedef int difference_type;
	    typedef Type* pointer;
	    typedef Type& reference;
	    typedef Type value_type;
	    typedef struct std::forward_iterator_tag iterator_category;
	public:
	    Iterator()
	    {
		this->parent = NULL;
	    }

	    Iterator(Array1D *parent, bool begin=true)
	    {
		this->parent = parent;

		if (begin)
		{
		    // begin iterator
		    if (parent->mySize == 0)
		    {
			currentPtr = parentCurrentPtr = NULL;
		    }
		    else
		    {
			parentCurrentPtr = parent->currentPtr;
			internalPageCount = 0;
			count = 0;
			pageCount = 0;
			currentPtr = parent->pages[0];
		    }
		}
		else
		{
		    // end iterator
		    if ( parent->mySize == 0 || (parent->mySize % pageSize) )
		    {
			pageCount = parent->mySize / pageSize;
		    }
		    else
		    {
			pageCount = parent->mySize / pageSize - 1;
		    }
		    count = parent->mySize;
		    currentPtr = parent->currentPtr;
		    parentCurrentPtr = parent->currentPtr;
		    internalPageCount = parent->mySize % pageSize;
		}
	    }

	    Iterator(Array1D *parent, Type *parentCurrentPtr, Type *currentPtr, unsigned int internalPageCount, unsigned int count, unsigned int pageCount)
		: parent(parent), currentPtr(currentPtr), parentCurrentPtr(parentCurrentPtr), internalPageCount(internalPageCount), count(count), pageCount(pageCount)
	    {
	    }


	    void setValue(Type& t)
	    {
		*currentPtr = t;
	    }


	    Type& operator*()
	    {
		return *currentPtr;
	    }

	    Type *getPtr()
	    {
		return currentPtr;
	    }

	    // Interface I

	    bool hasNext() const
	    {
		return currentPtr != parentCurrentPtr;    		
	    }

	    void next()
	    {
		count++;
		internalPageCount++;
		if (internalPageCount == pageSize)
		{
		    // parentCurrentPtr points to the element just beyond the end
		    // so here we check if currentPtr points to the end element
		    if (currentPtr != parentCurrentPtr-1)
		    {
			internalPageCount = 0;
			pageCount++;
			currentPtr = parent->pages[pageCount];
		    }
		    else
		    {
			// here we increment past the end
			currentPtr++;
		    }
		}
		else
		{
		    currentPtr++;
		}	    
	    }


	    /*! We assume that this method is not used when the iterator points to the first element and also that it is not used when it points past the end */
	    void prev()
	    {
		count--;
		if ( internalPageCount == 0 )
		{
		    internalPageCount = pageSize-1;
		    pageCount--;
		    currentPtr = parent->pages[pageCount] + (pageSize-1);
		}
		else
		{
		    internalPageCount--;
		    currentPtr--;
		}	    
	    }

	    ////////////////////


	    // Interface II

	    Iterator& operator++()
	    {
		next();
		return *this;
	    }


	    Iterator& operator++(int)
	    {
		next();
		return *this;
	    }


	    Iterator& operator--()
	    {
		prev();
		return *this;
	    }


	    Iterator& operator--(int)
	    {
		prev();
		return *this;
	    }

	    bool operator==(const Iterator& oi) const
	    {
		return oi.currentPtr == currentPtr;
	    }


	    bool operator!=(const Iterator& oi) const
	    {
		return oi.currentPtr != currentPtr;
	    }


	    Iterator operator+(int i) const
	    {
		unsigned int oCount = count + i;
		unsigned int oInternalPageCount = oCount & pageMask;
		unsigned int oPageCount = oCount >> pageExp;
		Type *oCurrentPtr = parent->pages[oPageCount] + oInternalPageCount;

		return Iterator(parent, parentCurrentPtr, oCurrentPtr, oInternalPageCount, oCount, oPageCount);
	    }
	    ////////////////////

	protected:
	    Array1D *parent;
	    Type *currentPtr;
	    Type *parentCurrentPtr;
	    unsigned int internalPageCount;
	    unsigned int count;
	    unsigned int pageCount;
	};

    public:
	Array1D(unsigned int initialSize=0)
	{
	    mySize = 0;
	    myAllocatedSize = 0;
	    numPages = 0;
	    pageArraySize = 0;
	    pages = NULL;
	    currentPtr = NULL;
	    currentEndPtr = NULL;
	    if (initialSize != 0)
	    {
		allocatePages( (initialSize >> pageExp) + ((initialSize & pageMask) == 0 ? 0 : 1) );
		currentPtr = pages[0];
		currentEndPtr = currentPtr + pageSize;
	    }
	}

	~Array1D()
	{
	    clear();
	}


	void init(std::string fileName, unsigned int numElements)
	{
	    unsigned int numIterations = (numElements / pageSize + ( (numElements % pageSize) == 0 ? 0 : 1 ) );
	    unsigned int i;
	    reserve(numElements);
	    std::ifstream is;
	    is.open(fileName.c_str(), std::ios_base::binary);
	    for (i=0; i<numIterations; i++)
	    {
		is.read((char *)pages[i], pageSize * sizeof(Type));
	    }
	    is.close();
	}


	/* This method sets currentPtr correctly */
	void reserve(unsigned int newSize)
	{
	    resize(newSize);
	    if ( mySize < newSize )
	    {
		mySize = newSize;
		currentPtr = &pages[mySize >> pageExp][mySize & pageMask];
		currentEndPtr = currentPtr + pageSize - ( mySize % pageSize );
	    }
	}

	void resize(unsigned int newAllocatedSize)
	{
	    // only resize if current size is smaller
	    if ( newAllocatedSize > myAllocatedSize )
	    {
		unsigned int diff = newAllocatedSize - myAllocatedSize;
		unsigned int numNewPages = diff / pageSize + ( (diff % pageSize)==0?0:1 );
		allocatePages(numNewPages);
	    }
	}

	void clear()
	{
	    unsigned int i;

	    for (i=0; i<numPages; i++)
	    {
		delete[] pages[i];
	    }
	    delete[] pages;
	    pageArraySize = 0;
	    numPages = 0;
	    pages = NULL;
	    currentPtr = NULL;
	    currentEndPtr = NULL;
	    mySize = 0;
	    myAllocatedSize = 0;
	}


	long long getMemUsed()
	{
	    return getMemUsage();
	}

	long long getMemUsage()
	{
	    long long memUsage = sizeof(this);
	    memUsage += numPages * pageSize * sizeof(Type);
	    memUsage += pageArraySize * sizeof(TypePtr);
	    return memUsage;
	}


	/*! this method assumes that the array is non-empty */
	Type& back()
	{
	    return (*this)(mySize-1);
	}


	void pop_back()
	{
	    mySize--;
	    currentPtr--;
	    if ( currentPtr == &pages[mySize >> pageExp][0])
	    {
		if ( mySize == 0 )
		{
		    currentPtr = NULL;
		    currentEndPtr = NULL;
		}
		else
		{
		    currentPtr = pages[(mySize-1) >> pageExp] + pageSize;
		    currentEndPtr = currentPtr;
		}
	    }
	    // currentPtr is needed for iteration!
	}


	void push_back(const Type &t)
	{
	    if ( currentPtr == currentEndPtr )
	    {
		if (mySize == myAllocatedSize)
		{
		    allocatePages(1);
		    currentPtr = pages[numPages-1];
		    currentEndPtr = currentPtr + pageSize;
		}
		else
		{
		    currentPtr = &pages[mySize >> pageExp][mySize & pageMask];
		    currentEndPtr = currentPtr + pageSize;
		}
	    }
	    // currentPtr is needed for iteration!
	    *currentPtr = t;
	    mySize++;
	    currentPtr++;
	}

	/*! \brief Allocate 'i' number of consecutive data elements. 
	*  It is assumed that this is always possible (e.g. for Octrees 'i' is always 8 and the page size is a multiple of 8) */ 
	Type *allocate(unsigned int i)
	{
	    Type *res;

	    if (mySize == myAllocatedSize)
	    {
		allocatePages(1);
	    }

	    res = currentPtr;
	    currentPtr+=i;
	    mySize+=i;
	    return res;    
	}

	Iterator iterator()
	{
	    return Iterator(this);
	}

	Iterator begin()
	{
	    return Iterator(this, true);	
	}

	Iterator end()
	{
	    return Iterator(this, false);
	}

	Type& operator[](int i) const
	{
	    return pages[i >> pageExp][i & pageMask];
	}

	Type& operator()(int i) const
	{
	    return pages[i >> pageExp][i & pageMask];
	}

	unsigned int size() const
	{
	    return mySize;
	}

	bool empty() const
	{
	    return mySize==0;
	}

	void serialize(std::ostream& os) const
	{
	    unsigned int writtenSize = 0;
	    unsigned int numWrittenPages = 0;
	    unsigned int writeSize;

	    os.write((const char *)&mySize, sizeof(unsigned int));

	    while (writtenSize < mySize)
	    {
		writeSize = std::min((mySize - writtenSize), pageSize) * sizeof(Type);
		os.write((const char *)(pages[numWrittenPages]), writeSize);
		writtenSize += pageSize;
		numWrittenPages++;
	    }
	}


	void deserialize(std::istream& is)
	{
	    unsigned int i;
	    unsigned int readSize;
	    is.read((char *)&mySize, sizeof(unsigned int));
	    unsigned int numNewPages = mySize / pageSize;
	    numPages = 0;
	    pageArraySize = 0;
	    if ( mySize % pageSize != 0)
	    {
		numNewPages++;
	    }
	    allocatePages(numNewPages);

	    for (i=0; i<numPages; i++)
	    {
		readSize = std::min(mySize-i*pageSize, pageSize) * sizeof(Type);
		is.read((char *)(pages[i]), readSize);
	    }

	    currentPtr = &pages[mySize >> pageExp][mySize & pageMask];
	    currentEndPtr = currentPtr + pageSize - ( mySize % pageSize ); 
	}


    protected:

	void allocatePages(unsigned int num)
	{
	    unsigned int i;
	    TypePtr *newPages;

	    if (numPages+num > pageArraySize)
	    {
		// reallocate page array
		// this is done in powers of two in order to have amortized constant time allocation
		if (pageArraySize == 0)
		{
		    //pageArraySize = 1;
		    pageArraySize = num;
		}
		else
		{
		    while (pageArraySize < numPages+num)
		    {
			pageArraySize = (pageArraySize<<1);
		    }
		}

		newPages = new TypePtr[pageArraySize];
		for (i=0; i<numPages; i++)
		{
		    newPages[i] = pages[i];
		}
		delete[] pages;
		pages = newPages;
	    }

	    for (i=0; i<num; i++)
	    {
		pages[numPages++] = new Type[pageSize];
		myAllocatedSize += pageSize;
	    }
	}

    protected:
	unsigned int mySize;
	unsigned int myAllocatedSize;

	// page based storage
	unsigned int numPages;
	unsigned int pageArraySize;
	TypePtr *pages;
	Type *currentPtr;
	Type *currentEndPtr;


	friend class Iterator;
    };


}

#endif
