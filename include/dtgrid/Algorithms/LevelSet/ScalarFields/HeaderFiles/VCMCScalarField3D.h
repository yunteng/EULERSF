/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _levelset_VCMCscalarfield3d_h
#define _levelset_VCMCscalarfield3d_h

#include "ScalarField3D.h"


namespace LevelSet
{

    template<class Real, class Index, class Grid, class Iterator>
    class VCMCScalarField3D : public ScalarField3D<Real, Index, Grid, Iterator>
    {
    public:
	VCMCScalarField3D(bool useSpecifiedAMC=false, Real AMC=0);
	~VCMCScalarField3D(void);

	void init(Grid *phi);

	inline void propagate(Real dt, Grid *phi);

	inline Real computeHyperbolicTerm(const Iterator& iter);

	inline Real computeParabolicTerm(const Iterator& iter, Real kappa);

	inline bool isHyperBolic() const { return true; }

	inline bool isParabolic() const { return true; }

	inline Real getAverageMeanCurvature() { return averageMeanCurvature; }

	inline void setAverageMeanCurvature(Real averageMeanCurvature) { this->averageMeanCurvature = averageMeanCurvature; absAMC = fabs(averageMeanCurvature); }

    protected:
	bool useSpecifiedAMC;
	Real averageMeanCurvature;
	Real absAMC;
	Real dx;
    };

}

#include "VCMCScalarField3DImpl.h"


#endif
