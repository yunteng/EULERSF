/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _levelset_velocityfield3d_h
#define _levelset_velocityfield3d_h

#include <Core/Exception/HeaderFiles/DefaultException.h>

namespace LevelSet
{

    template<class Real, class Index>
    class VelocityField3D
    {
    public:

	virtual ~VelocityField3D(void) { }

	virtual void advect(Real dt) = 0;

	// input:   p : position
	// output:  v : velocity at p
	virtual void operator()(const Index p[3], Real v[3]) const = 0;

	virtual void operator()(Index i, Index j, Index k, Real v[3]) const = 0;

	virtual void setTime(Real t) { this->t = t; }

	virtual bool maxAbsConstantInNarrowBand() { return false; }

	virtual Real getMaxAbs() const { Core::throwDefaultException("Method not implemented!", __FILE__, __LINE__); return -1; }

    protected:

	Real t;

    };

}

#endif
