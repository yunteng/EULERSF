/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _levelset_enrighttestvelocityfield3d_h
#define _levelset_enrighttestvelocityfield3d_h

#include "VelocityField3D.h"
#include <Algorithms/Math/HeaderFiles/BasicMath.h>


namespace LevelSet
{

    /**
    *  Designed for use in a unit computational domain: [0;1]^3
    */
    template<class Real, class Index>
    class EnrightTestVelocityField3D : public VelocityField3D<Real, Index>
    {
    public:
	EnrightTestVelocityField3D(Real dx);

	EnrightTestVelocityField3D(Index dim[3]);

	~EnrightTestVelocityField3D(void) { }

	void advect(Real dt);

	// input:   p : position
	// output:  v : Deformation at p
	void operator()(const Index p[3], Real v[3]) const;

	void operator()(Index i, Index j, Index k, Real v[3]) const;

    protected:
	unsigned int dim[3];
	Real dxr, dyr, dzr, c1;
	Real t;
	static const Real T;	
    };

}


#include "EnrightTestVelocityField3DImpl.h"


#endif
