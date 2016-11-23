/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_h
#define _grids_h


namespace Grids
{

    enum IterationOrder { IteratorOrder_XYZ, IteratorOrder_ZYX, IterationOrder_YX, IterationOrder_XY };

    enum StencilFormat { SF_NONE, SF_FIRSTORDER, SF_FIRSTORDER_CURVATURE, SF_WENO, SF_WENO_CURVATURE, SF_VOXEL, SF_VOXEL_GRAD, SF_BOX, SF_FIRSTORDER_NEIGHBORS };

    enum TubeType { NO_TUBE=0x01, ZERO_CROSSING_TUBE=0x02, BETA_TUBE=0x04, GAMMA_TUBE=0x08, ENTIRE_TUBE=0x10, INSIDE_TUBE=0x20, OUTSIDE_TUBE=0x40, ENTIRE_VOLUME=0x80, INSIDE_VOLUME=0x100, OUTSIDE_VOLUME=0x200};

    enum IteratorType { IT_TUBE, IT_VOLUME };


    struct StencilTraits
    {
	template<StencilFormat stencilFormat>
	static int getStencilLength()
	{
	    return ( (stencilFormat==SF_BOX ? 27 : (stencilFormat==SF_VOXEL_GRAD ? 32 : (stencilFormat==SF_VOXEL ? 8 : (stencilFormat==SF_FIRSTORDER ? 7 : ( (stencilFormat==SF_FIRSTORDER_CURVATURE || stencilFormat==SF_WENO) ? 19 : ((stencilFormat==SF_WENO_CURVATURE) ? 31 : ( stencilFormat==SF_FIRSTORDER_NEIGHBORS ? 25 : 1 ) ) ) ) ) ) ) );
	}

	/** The maximum number of _Iterator instances in the Iterator stencil. See the specification of the StencilFormat below. */
	static const int maxNumberOfStencilIterators = 32;
    };


}


#endif
