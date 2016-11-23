/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_svolloaderfactory_h_
#define _grids_svolloaderfactory_h_

#include <iostream>
#include <fstream>
#include <Core/Exception/HeaderFiles/DefaultException.h>
#include "SvolLoaderInterface.h"
#include "SvolLoader_DTTop1.0.h"


namespace Grids
{

    class SvolLoaderFactory
    {
    public:

	template <typename IndexType, typename RealType, typename DataType, typename UIntType>
	static SvolLoaderInterface<IndexType, RealType, DataType, UIntType> *getInstance(const std::string& fileName, std::string *header)
	{
	    std::ifstream is(fileName.c_str());

	    if (!is)
	    {
		Core::throwDefaultException("Unable to open file "+fileName+" for loading!", __FILE__, __LINE__);
	    }

	    // Load header: Assumption: Header is always 10 bytes long!
	    char headerLocal[11];
	    is.read(headerLocal, 10);
	    headerLocal[10] = 0;

	    *header = string(headerLocal);

	    if (*header == SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>::header)
	    {
		return new SvolLoaderDTTop1<IndexType, RealType, DataType, UIntType>();
	    }

	    return NULL;
	}

    };


}



#endif
