/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

namespace LevelSet
{

    template<class Real, class Index, class Grid, class Iterator>
    MeanCurvatureFlowScalarField3D<Real, Index, Grid, Iterator>::MeanCurvatureFlowScalarField3D()
    {
    }

    template<class Real, class Index, class Grid, class Iterator>
    MeanCurvatureFlowScalarField3D<Real, Index, Grid, Iterator>::~MeanCurvatureFlowScalarField3D(void)
    {
    }

    template<class Real, class Index, class Grid, class Iterator>
    void MeanCurvatureFlowScalarField3D<Real, Index, Grid, Iterator>::init(Grid *phi)
    {
    }


    template<class Real, class Index, class Grid, class Iterator>
    void MeanCurvatureFlowScalarField3D<Real, Index, Grid, Iterator>::propagate(Real dt, Grid *phi)
    {
    }


    template<class Real, class Index, class Grid, class Iterator>
    Real MeanCurvatureFlowScalarField3D<Real, Index, Grid, Iterator>::computeHyperbolicTerm(const Iterator& iter)
    {
	return 0;
    }

    template<class Real, class Index, class Grid, class Iterator>
    Real MeanCurvatureFlowScalarField3D<Real, Index, Grid, Iterator>::computeParabolicTerm(const Iterator& iter, Real kappa)
    {
	return 1;
    }

}
