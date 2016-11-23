/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _containers_indexvaluecontainerimpl_h
#define _containers_indexvaluecontainerimpl_h


namespace Containers
{

    template<typename Real, typename Index>
    IndexValueContainer<Real, Index>::IndexValueContainer(const std::string& fileName, LoadMode loadMode)
    {
	init(fileName, loadMode);
    }


    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::init(const std::string& fileName, LoadMode loadMode)
    {
	iva = NULL;
	this->fileName = fileName;
	if ( loadMode == LM_AT_ONCE )
	{
	    std::ifstream *privateInStream;
	    privateInStream = new std::ifstream(fileName.c_str(), std::ios_base::binary);
	    deserialize(*privateInStream);
	    loaded = true;
	}
	else if (loadMode == LM_SEQUENTIALLY )
	{
	    std::ifstream *privateInStream;
	    privateInStream = new std::ifstream(fileName.c_str(), std::ios_base::binary);
	    privateInStream->read((char *)dim, 3 * sizeof(unsigned int));
	    privateInStream->read((char *)&isSorted, sizeof(bool));
	    privateInStream->read((char *)&numElements, sizeof(unsigned int));
	    loaded = false;
	    privateInStream->close();
	    delete privateInStream;
	}
	else
	{
	    IndexValueContainer();
	}
    }


    template<typename Real, typename Index>
    IndexValueContainer<Real, Index>::IndexValueContainer()
    {
	this->dim[0] = 0;
	this->dim[1] = 0;
	this->dim[2] = 0;
	iva = NULL;
	numElements = 0;
	isSorted = false;
	loaded = true;
    }


    template<typename Real, typename Index>
    template<typename DimType>
    IndexValueContainer<Real, Index>::IndexValueContainer(DimType dim[3])
    {
	this->dim[0] = dim[0];
	this->dim[1] = dim[1];
	this->dim[2] = dim[2];
	//iva = new vector< IndexValue >[ dim[2] ];
	iva = new Arrays::Array1D< IndexValue >[ dim[2] ];
	numElements = 0;
	isSorted = false;
	loaded = true;
    }

    template<typename Real, typename Index>
    IndexValueContainer<Real, Index>::~IndexValueContainer()
    {
	delete[] iva;
    }


    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::clear()
    {
	unsigned int i, size;

	if (isSorted)
	{
	    size = dim[0];
	}
	else
	{
	    size = dim[2];
	}


	for (i=0; i<size; i++)
	{
	    iva[i].clear();
	}

	isSorted = false;
	numElements = 0;
    }


    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::size(unsigned int dim[3]) const
    {
	dim[0] = this->dim[0];
	dim[1] = this->dim[1];
	dim[2] = this->dim[2];
    }


    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::add(Index x, Index y, Index z)
    {
	IndexValue iv;
	iv.a = x;
	iv.b = y;
	iva[z].push_back(iv);
	numElements++;
    }


    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::add(Index x, Index y, Index z, Real value)
    {
	IndexValue iv;
	iv.a = x;
	iv.b = y;
	iv.v = value;
	iva[z].push_back(iv);
	numElements++;
    }


    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::add(IndexValueContainer *oivc, bool clearOther)
    {
	unsigned int i, j, size;

	for (i=0; i<this->dim0; i++)
	{
	    Arrays::Array1D< typename IndexValueContainer<Real, Index>::IndexValue >& ova = oivc->iva[i];
	    typename Arrays::Array1D< typename IndexValueContainer<Real, Index>::IndexValue >::Iterator iter = ova.iterator();
	    while (iter.hasNext())
	    {
		iva[i].push_back( *iter );
		iter.next();
	    }


	    if (clearOther)
	    {
		ova.clear();
	    }
	}
    }



    template<typename Real, typename Index>
    unsigned int IndexValueContainer<Real, Index>::getMemUsed()
    {
	if ( isSorted )
	{
	    return sizeof(this) + dim[0] * sizeof(Arrays::Array1D<IndexValue>) + numElements * sizeof(IndexValue);
	}
	else
	{
	    return sizeof(this) + dim[2] * sizeof(Arrays::Array1D<IndexValue>) + numElements * sizeof(IndexValue);
	}
    }



    template<typename Real, typename Index>
    typename IndexValueContainer<Real, Index>::Iterator IndexValueContainer<Real, Index>::iterator(bool clear)
    {
	return Iterator(this, clear);
    }


    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::sort()
    {
	// do two stable bucket sorts
	// during insertion, a bucket sort wrt Z was performed
	unsigned int i, t;

	// do bucket sort in Y
	Arrays::Array1D<IndexValue> *ivay = new Arrays::Array1D<IndexValue>[dim[1]];

	for (i=0; i<dim[2]; i++)
	{
	    Arrays::Array1D< typename IndexValueContainer<Real, Index>::IndexValue >& iv = iva[i];
	    typename Arrays::Array1D< typename IndexValueContainer<Real, Index>::IndexValue >::Iterator iter = iv.iterator();
	    while (iter.hasNext())
	    {
		IndexValue& current = *iter;	
		t = current.b;
		current.b = i;

		ivay[t].push_back( current );
		iter.next();
	    }

	    iv.clear();
	}
	delete[] iva;
	iva = ivay;


	// do bucket sort in X
	Arrays::Array1D<IndexValue> *ivax = new Arrays::Array1D<IndexValue>[dim[0]];

	for (i=0; i<dim[1]; i++)
	{
	    Arrays::Array1D< IndexValue >& iv = iva[i];
	    typename Arrays::Array1D< IndexValue >::Iterator iter = iv.iterator();
	    while (iter.hasNext())
	    {
		IndexValue& current = *iter;	
		t = current.a;
		current.a = i;

		ivax[t].push_back( current );
		iter.next();
	    }

	    iv.clear();
	}
	delete[] iva;
	iva = ivax;


	isSorted = true;
    }



    // IMPLEMENTATION OF CLASS ITERATOR

    template<typename Real, typename Index>
    IndexValueContainer<Real, Index>::Iterator::Iterator()
	: parent(NULL)
    {
	clear = false;
	this->dim0 = 0;
	i = 0;
	//j = 0;
	parentLoaded = true;
	privateInStream = NULL;
    }

    template<typename Real, typename Index>
    IndexValueContainer<Real, Index>::Iterator::~Iterator()
    {
	delete privateInStream;
    }

    template<typename Real, typename Index>
    IndexValueContainer<Real, Index>::Iterator::Iterator(IndexValueContainer *parent, bool clear)
	: parent(parent), clear(clear)
    {
	parentLoaded = parent->loaded;
	isSorted = parent->isSorted;

	if (isSorted)
	{
	    dim0 = (*parent).dim[0];
	}
	else
	{
	    dim0 = (*parent).dim[2];
	}
	i = 0;

	if (parentLoaded)
	{
	    privateInStream = NULL;
	    while ( i<dim0 && (*parent).iva[i].size()==0 )
	    {
		i++;
	    }
	    if ( i<dim0 )
	    {
		iter = parent->iva[i].iterator();
		current = iter.getPtr();
	    }
	}
	else
	{
	    current = &privateCurrent;
	    privateInStream = new std::ifstream(parent->fileName.c_str(), std::ios_base::binary);
	    privateInStream->seekg(3 * sizeof(unsigned int) + sizeof(bool) + sizeof(unsigned int));
	    privateStreamCount = 0;
	    privateInStream->read((char *)&privateStreamSize, sizeof(unsigned int));
	    // find the first coordinate for which there is a non-empty array
	    while ( i<dim0 && privateStreamSize==0 )
	    {
		if (i != dim0-1)
		{
		    privateInStream->read((char *)&privateStreamSize, sizeof(unsigned int));
		}
		i++;
	    }
	    if ( i<dim0 )
	    {
		privateInStream->read((char *)&privateCurrent, sizeof(IndexValue));
		privateStreamCount++;
	    }
	}
    }


    template<typename Real, typename Index>
    IndexValueContainer<Real, Index>::Iterator::Iterator(const Iterator& oi)
    {
	i = oi.i;
	parent = oi.parent;
	dim0 = oi.dim0;
	isSorted = oi.isSorted;
	clear = oi.clear;
	parentLoaded = oi.parentLoaded;

	if (oi.privateInStream != NULL)
	{
	    privateCurrent = oi.privateCurrent;
	    current = &privateCurrent;
	    privateInStream = new std::ifstream(parent->fileName.c_str(), std::ios_base::binary);
	    privateInStream->seekg(oi.privateInStream->tellg());
	    privateStreamSize = oi.privateStreamSize;
	    privateStreamCount = oi.privateStreamCount;
	}
	else
	{
	    iter = oi.iter;
	    current = iter.getPtr();    
	    privateInStream = NULL;
	}
    }


    template<typename Real, typename Index>
    typename IndexValueContainer<Real, Index>::Iterator& IndexValueContainer<Real, Index>::Iterator::operator=(const Iterator& oi)
    {
	i = oi.i;
	parent = oi.parent;
	dim0 = oi.dim0;
	isSorted = oi.isSorted;
	clear = oi.clear;
	parentLoaded = oi.parentLoaded;

	delete privateInStream;
	if (oi.privateInStream != NULL)
	{
	    privateCurrent = oi.privateCurrent;
	    current = &privateCurrent;
	    privateInStream = new std::ifstream(parent->fileName.c_str(), std::ios_base::binary);
	    privateInStream->seekg(oi.privateInStream->tellg());
	    privateStreamSize = oi.privateStreamSize;
	    privateStreamCount = oi.privateStreamCount;
	}
	else
	{
	    iter = oi.iter;
	    current = iter.getPtr();    
	    privateInStream = NULL;
	}

	return *this;
    }


    template<typename Real, typename Index>
    Real IndexValueContainer<Real, Index>::Iterator::getValue() const
    {
	return current->v;
    }

    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::Iterator::setValue(Real v)
    {
	current->v = v;
    }

    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::Iterator::getIndex(Index *x, Index *y, Index *z) const
    {
	if (isSorted)
	{
	    *x = i;
	    *y = current->a;
	    *z = current->b;
	}
	else
	{
	    *z = i;
	    *x = current->a;
	    *y = current->b;
	}
    }


    template<typename Real, typename Index>
    Index IndexValueContainer<Real, Index>::Iterator::getI() const
    {
	if (isSorted)
	{
	    return i;
	}
	else
	{
	    return current->a;
	}
    }


    template<typename Real, typename Index>
    Index IndexValueContainer<Real, Index>::Iterator::getJ() const
    {
	if (isSorted)
	{
	    return current->a;
	}
	else
	{
	    return current->b;
	}
    }


    template<typename Real, typename Index>
    Index IndexValueContainer<Real, Index>::Iterator::getK() const
    {
	if (isSorted)
	{
	    return current->b;
	}
	else
	{
	    return i;
	}
    }


    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::Iterator::next()
    {
	if (parentLoaded)
	{
	    iter.next();
	    if ( iter.hasNext() )
	    {
		current = iter.getPtr();
	    }
	    else
	    {
		if (clear)
		{
		    parent->iva[i].clear();
		}
		i++;
		while ( i<dim0 && (*parent).iva[i].size()==0 )
		{
		    i++;
		}
		if ( i<dim0 )
		{
		    iter = parent->iva[i].iterator();
		    current = iter.getPtr();
		}
	    }
	}
	else
	{
	    if (privateStreamCount < privateStreamSize)
	    {
		// continue reading for this coordinate
		// the actual indexvalue is read into privateCurrent
		privateInStream->read((char *)&privateCurrent, sizeof(IndexValue));	    
		privateStreamCount++;
	    }
	    else
	    {
		privateStreamCount = 0;
		privateInStream->read((char *)&privateStreamSize, sizeof(unsigned int));
		i++;
		// find next coordinate for which there is a non-empty array
		while ( i<dim0 && privateStreamSize==0 )
		{
		    if (i != dim0-1)
		    {
			privateInStream->read((char *)&privateStreamSize, sizeof(unsigned int));
		    }
		    i++;
		}
		if ( i<dim0 )
		{
		    privateInStream->read((char *)&privateCurrent, sizeof(IndexValue));
		    privateStreamCount++;
		}
	    }
	}
    }

    template<typename Real, typename Index>
    bool IndexValueContainer<Real, Index>::Iterator::hasNext() const
    {
	return i!=dim0;
    }



    // SERIALIZATION / DESERIALIZATION




    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::serialize(const std::string& fileName) const
    {
	std::ofstream os(fileName.c_str(), std::ios_base::binary);
	serialize(os);
    }

    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::serialize(std::ostream& os) const
    {
	unsigned int i, e;
	os.write((const char *)dim, 3 * sizeof(unsigned int));
	os.write((const char *)&isSorted, sizeof(bool));
	os.write((const char *)&numElements, sizeof(unsigned int));
	if (isSorted)
	{
	    e = dim[0];
	}
	else
	{
	    e = dim[2];
	}

	for (i=0; i<e; i++)
	{
	    iva[i].serialize(os);
	}
    }


    template<typename Real, typename Index>
    void IndexValueContainer<Real, Index>::deserialize(std::istream& is)
    {
	unsigned int i, e;
	is.read((char *)dim, 3 * sizeof(unsigned int));
	is.read((char *)&isSorted, sizeof(bool));
	is.read((char *)&numElements, sizeof(unsigned int));
	if (isSorted)
	{
	    e = dim[0];
	}
	else
	{
	    e = dim[2];
	}

	delete[] iva;
	iva = new Arrays::Array1D< IndexValue >[ e ];
	for (i=0; i<e; i++)
	{
	    iva[i].deserialize(is);
	}    
    }

}

#endif
