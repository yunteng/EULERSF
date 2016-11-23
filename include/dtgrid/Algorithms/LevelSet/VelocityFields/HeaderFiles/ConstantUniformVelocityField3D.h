/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _levelset_ConstantUniformVelocityfield3d_h
#define _levelset_ConstantUniformVelocityfield3d_h

#include "VelocityField3D.h"

namespace LevelSet
{

    template<class Real, class Index>
    class ConstantUniformVelocityField3D : public VelocityField3D<Real, Index>
    {
    public:
	ConstantUniformVelocityField3D(Real v[3]);

	~ConstantUniformVelocityField3D(void) { }

	void advect(Real dt);

	// input:   p : position
	// output:  v : ConstantUniformVelocity at p
	void operator()(const Index p[3], Real v[3]) const;

	void operator()(Index i, Index j, Index k, Real v[3]) const;

	Real getMaxAbs() const;

	bool maxAbsConstantInNarrowBand() { return true; }

    protected:
	Real v[3];
	Real maxAbs;
    };

}

#include "ConstantUniformVelocityField3DImpl.h"

#endif
