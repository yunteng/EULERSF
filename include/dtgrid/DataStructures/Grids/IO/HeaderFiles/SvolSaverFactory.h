/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_svolsaverfactory_h_
#define _grids_svolsaverfactory_h_

#include "SvolSaverInterface.h"
#include "SvolSaver_DTTop1.0.h"


namespace Grids
{

    class SvolSaverFactory
    {
    public:

	template <typename IndexType, typename RealType, typename DataType, typename UIntType>
	    static SvolSaverInterface<IndexType, RealType, DataType, UIntType> *getInstance(const std::string& idString)
	{
	    if (idString == SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>::idString)
	    {
		return new SvolSaverDTTop1<IndexType, RealType, DataType, UIntType>();
	    }

	    return NULL;
	}

    };


}



#endif
