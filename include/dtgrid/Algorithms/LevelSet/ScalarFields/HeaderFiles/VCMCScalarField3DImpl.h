/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#define _USE_MATH_DEFINES   // M_PI etc in math.h
#include <math.h>
#include <DataStructures/Grids/HeaderFiles/Grids.h>
#include <Algorithms/Math/HeaderFiles/BasicMath.h>

namespace LevelSet
{

    template<class Real, class Index, class Grid, class Iterator>
    VCMCScalarField3D<Real, Index, Grid, Iterator>::VCMCScalarField3D(bool useSpecifiedAMC, Real AMC)
    {
	this->useSpecifiedAMC = useSpecifiedAMC;
	this->averageMeanCurvature = AMC;
	this->absAMC = fabs(AMC);
    }

    template<class Real, class Index, class Grid, class Iterator>
    VCMCScalarField3D<Real, Index, Grid, Iterator>::~VCMCScalarField3D(void)
    {
    }

    template<class Real, class Index, class Grid, class Iterator>
    void VCMCScalarField3D<Real, Index, Grid, Iterator>::init(Grid *phi)
    {
	if (!useSpecifiedAMC)
	{

	    // compute average mean curvature
	    // epsilon = 1.5*dx according to osher, fedkiw page 15
#ifdef WIN32
	    typename Grid::template StencilTubeIterator<Grids::BETA_TUBE> iter = phi->beginStencilTubeIterator<Grids::BETA_TUBE, Grids::SF_FIRSTORDER_CURVATURE, false>();
	    typename Grid::template StencilTubeIterator<Grids::BETA_TUBE> iend = phi->endStencilTubeIterator<Grids::BETA_TUBE, Grids::SF_FIRSTORDER_CURVATURE, false>();
#else
	    typename Grid::template StencilTubeIterator<Grids::BETA_TUBE> iter = phi->template beginStencilTubeIterator<Grids::BETA_TUBE, Grids::SF_FIRSTORDER_CURVATURE, false>();
	    typename Grid::template StencilTubeIterator<Grids::BETA_TUBE> iend = phi->template endStencilTubeIterator<Grids::BETA_TUBE, Grids::SF_FIRSTORDER_CURVATURE, false>();
#endif
	    dx = phi->getDx();
	    Real epsilon = static_cast<Real>(1.5) * dx;
	    Real volElement = Math::pow3(dx);
	    Real deltaVal;
	    Real c1 = static_cast<Real>(1.0) / (static_cast<Real>(2.0) * epsilon);
	    Real c2 = static_cast<Real>(M_PI) / epsilon;
	    Real H;
	    Real gradLen;
	    Real area;
	    Real cellArea;
	    Real phiVal;

	    // compute surface integral
	    averageMeanCurvature = 0;
	    area = 0;

	    while (iter != iend)
	    {
		phiVal = iter.getValue();
		if (fabs(phiVal) <= epsilon)
		{
		    deltaVal = c1 + c1 * cos( phiVal * c2 );
		    H = iter.template meanCurvature<Grids::SF_FIRSTORDER_CURVATURE>();
		    gradLen = iter.gradientLength();
		    cellArea = deltaVal * gradLen * volElement;
		    area += cellArea;
		    averageMeanCurvature +=  H * cellArea; 
		}

#ifdef WIN32
		iter.operator++<Grids::SF_FIRSTORDER_CURVATURE, false>();
#else
		iter.template operator++<Grids::SF_FIRSTORDER_CURVATURE, false>();
#endif
	    }
	    iter.commit();

	    averageMeanCurvature /= area;

	    absAMC = fabs(averageMeanCurvature);
	}
    }

    template<class Real, class Index, class Grid, class Iterator>
    void VCMCScalarField3D<Real, Index, Grid, Iterator>::propagate(Real dt, Grid *phi)
    {
    }


    template<class Real, class Index, class Grid, class Iterator>
    Real VCMCScalarField3D<Real, Index, Grid, Iterator>::computeHyperbolicTerm(const Iterator& iter)
    {
	return averageMeanCurvature;
    }


    template<class Real, class Index, class Grid, class Iterator>
    Real VCMCScalarField3D<Real, Index, Grid, Iterator>::computeParabolicTerm(const Iterator& iter, Real kappa)
    {
	return 1;
    }

}
