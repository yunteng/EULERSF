/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _containers_indexvaluecontainer_h
#define _containers_indexvaluecontainer_h

#include <DataStructures/Arrays/HeaderFiles/Array1D.h>
#include <iostream>
#include <fstream>


namespace Containers
{

    template<typename Real, typename Index>
    class IndexValueContainer
    {
    public:
	struct IndexValue
	{
	    Index a, b;
	    Real v;
	};

	enum LoadMode { LM_AT_ONCE, LM_SEQUENTIALLY };

	// the iterator can only be acquired if the IndexValueContainer is sorted
	class Iterator
	{
	public:
	    Iterator();
	    Iterator(IndexValueContainer *parent, bool clear);
	    Iterator(const Iterator& oi);
	    ~Iterator();

	    Real getValue() const;
	    void setValue(Real v);
	    void getIndex(Index *x, Index *y, Index *z) const;
	    Index getI() const;
	    Index getJ() const;
	    Index getK() const;
	    void next();
	    bool hasNext() const;
	    void commit() { }

	    Iterator& operator=(const Iterator& oi);


	protected:
	    unsigned int i;
	    IndexValueContainer *parent;
	    IndexValue *current;
	    IndexValue privateCurrent;
	    typename Arrays::Array1D<IndexValue>::Iterator iter;
	    unsigned int dim0;
	    bool isSorted;
	    bool clear;
	    bool parentLoaded;

	    std::ifstream *privateInStream;
	    unsigned int privateStreamSize;
	    unsigned int privateStreamCount;
	};

    public:
	IndexValueContainer(const std::string& fileName, LoadMode loadMode);
	IndexValueContainer();
	template<typename DimType>
	IndexValueContainer(DimType dim[3]);
	~IndexValueContainer();

	void init(const std::string& fileName, LoadMode loadMode);
	void clear();

	void add(IndexValueContainer *oivc, bool clearOther);
	void add(Index x, Index y, Index z);
	void add(Index x, Index y, Index z, Real value);

	Iterator iterator(bool clear=false);

	void sort();
	void size(unsigned int dim[3]) const;
	unsigned int getNumElements() const { return numElements; }


	void serialize(const std::string& fileName) const;
	void serialize(std::ostream& os) const;
	void deserialize(std::istream& is);

	unsigned int getMemUsed();


    protected:
	unsigned int dim[3];
	Arrays::Array1D< IndexValue > *iva;
	bool isSorted;
	unsigned int numElements;
	bool loaded;
	std::string fileName;

	// friends
	friend class Iterator;
    };

}

#include "IndexValueContainerImpl.h"


#endif
