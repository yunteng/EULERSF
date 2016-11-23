/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _levelset_scalarfield3d_h
#define _levelset_scalarfield3d_h


namespace LevelSet
{


    template<class Real, class Index, class Grid, class Iterator>
    class ScalarField3D
    {
    public:
	virtual ~ScalarField3D(void) { }

	virtual void init(Grid *phi) = 0;

	virtual void propagate(Real dt, Grid *phi) = 0;

	virtual Real computeHyperbolicTerm(const Iterator& iter) = 0;

	/* In the equation dPhi/dt = b * kappa * norm(grad(Phi)), where kappa is the curvature, the method computeParabolicTerm returns 'b' */
	virtual Real computeParabolicTerm(const Iterator& iter, Real kappa) = 0;

	virtual bool isHyperBolic() const = 0;

	virtual bool isParabolic() const = 0;

	virtual void setTime(Real t) { this->t = t; }

    protected:
	Real t;
    };


}

#endif
