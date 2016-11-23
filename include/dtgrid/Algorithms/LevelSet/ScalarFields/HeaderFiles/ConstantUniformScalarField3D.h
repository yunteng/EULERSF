/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _levelset_constantuniformscalarfield3d_h
#define _levelset_constantuniformscalarfield3d_h

#include "ScalarField3D.h"


namespace LevelSet
{

    template<class Real, class Index, class Grid, class Iterator>
    class ConstantUniformScalarField3D : public ScalarField3D<Real, Index, Grid, Iterator>
    {
    public:
	ConstantUniformScalarField3D();
	ConstantUniformScalarField3D(Real value);
	~ConstantUniformScalarField3D(void);

	void init(Grid *phi);

	inline void propagate(Real dt, Grid *phi);

	inline Real computeHyperbolicTerm(const Iterator& iter);

	inline Real computeParabolicTerm(const Iterator& iter, Real kappa);

	inline bool isHyperBolic() const { return true; }

	inline bool isParabolic() const { return false; }

    protected:
	Real maxAbs;
	Real value;
    };


}

#include "ConstantUniformScalarField3DImpl.h"

#endif
