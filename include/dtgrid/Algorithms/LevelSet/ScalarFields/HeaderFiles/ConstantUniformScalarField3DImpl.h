/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#include <math.h>


namespace LevelSet
{

    template<class Real, class Index, class Grid, class Iterator>
    ConstantUniformScalarField3D<Real, Index, Grid, Iterator>::ConstantUniformScalarField3D()
    {
	this->value = 0;
	maxAbs = 0;
    }


    template<class Real, class Index, class Grid, class Iterator>
    ConstantUniformScalarField3D<Real, Index, Grid, Iterator>::ConstantUniformScalarField3D(Real value)
    {
	this->value = value;
	maxAbs = fabs(value);
    }

    template<class Real, class Index, class Grid, class Iterator>
    ConstantUniformScalarField3D<Real, Index, Grid, Iterator>::~ConstantUniformScalarField3D(void)
    {
    }

    template<class Real, class Index, class Grid, class Iterator>
    void ConstantUniformScalarField3D<Real, Index, Grid, Iterator>::init(Grid *phi)
    {
    }

    template<class Real, class Index, class Grid, class Iterator>
    void ConstantUniformScalarField3D<Real, Index, Grid, Iterator>::propagate(Real dt, Grid *phi)
    {
    }

  
    template<class Real, class Index, class Grid, class Iterator>
    Real ConstantUniformScalarField3D<Real, Index, Grid, Iterator>::computeHyperbolicTerm(const Iterator& iter)
    {
	return value;
    }

    template<class Real, class Index, class Grid, class Iterator>
    Real ConstantUniformScalarField3D<Real, Index, Grid, Iterator>::computeParabolicTerm(const Iterator& iter, Real kappa)
    {
	return 0;
    }

}
