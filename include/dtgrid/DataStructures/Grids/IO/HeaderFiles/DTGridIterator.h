/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_dtgriditerator_h_
#define _grids_dtgriditerator_h_


namespace Grids
{
    template<typename IndexType, typename RealType, typename DataType, typename UIntType, typename IndexTypeIterator, typename DataTypeIterator, typename UIntTypeIterator>
    class DTGridIterator
    {
    public:

	class Input
	{
	public:

	    /* indices that point to the start of each (x)-column in yIndex */
	    UIntTypeIterator va1DBegin, va1DEnd;
	    /* start and end x indices of the connected components contained in the projection of the narrow band onto the x axis */
	    IndexTypeIterator xIndexBegin, xIndexEnd;              

	    /* indices that point to the start of each (x,y)-column in zIndex */
	    UIntType va2DBegin, va2DEnd;
	    /* start and end y indices of the connected components in the y direction contained in the projection of the narrow band onto the (x,y) plane */
	    IndexType yIndexBegin, yIndexEnd;              

	    /* start and end z index of each connected component in the z direction */
	    IndexTypeIterator zIndexBegin, zIndexEnd;              

	    /* values of grid points in narrow band (will always correspond to values[0]) */
	    DataTypeIterator va3DBegin, va3DEnd;				   
	};


    public:

	DTGridIterator(const Input& input);




    };

}

#endif
