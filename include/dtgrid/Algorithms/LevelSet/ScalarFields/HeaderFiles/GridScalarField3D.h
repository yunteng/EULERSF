/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _levelset_gridscalarfield3d_h
#define _levelset_gridscalarfield3d_h

#include "ScalarField3D.h"

namespace LevelSet
{

    template<class Real, class Index, class Grid, class Iterator>
    class GridScalarField3D : public ScalarField3D<Real, Index, Grid, Iterator>
    {
    public:
	GridScalarField3D();
	/** GridScalar3D owns 'g' */
	GridScalarField3D(Grid *g);
	~GridScalarField3D(void);

	void init(Grid *phi);

	inline void propagate(Real dt, Grid *phi);

	inline Real computeHyperbolicTerm(const Iterator& otherIter);

	inline Real computeParabolicTerm(const Iterator& iter, Real kappa);

	inline bool isHyperBolic() const { return true; }

	inline bool isParabolic() const { return false; }

    protected:
	Grid *g;
	typename Grid::TubeIterator iter, iend;
	Real minusGamma, plusGamma, value;
    };


}

#include "GridScalarField3D_Impl.h"

#endif
