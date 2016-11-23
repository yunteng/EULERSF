/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _levelset_meancurvatureflowscalarfield3d_h
#define _levelset_meancurvatureflowscalarfield3d_h

#include "ScalarField3D.h"

namespace LevelSet
{

    template<class Real, class Index, class Grid, class Iterator>
    class MeanCurvatureFlowScalarField3D : public ScalarField3D<Real, Index, Grid, Iterator>
    {
    public:
	MeanCurvatureFlowScalarField3D();
	~MeanCurvatureFlowScalarField3D(void);

	void init(Grid *phi);

	inline void propagate(Real dt, Grid *phi);

	inline Real computeHyperbolicTerm(const Iterator& iter);

	inline Real computeParabolicTerm(const Iterator& iter, Real kappa);

	inline bool isHyperBolic() const { return false; }

	inline bool isParabolic() const { return true; }

    };


}

#include "MeanCurvatureFlowScalarField3D_Impl.h"

#endif
