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
    GridScalarField3D<Real, Index, Grid, Iterator>::GridScalarField3D()
    {
	this->g = NULL;
    }

    template<class Real, class Index, class Grid, class Iterator>
    GridScalarField3D<Real, Index, Grid, Iterator>::GridScalarField3D(Grid *g)
    {
	this->g = g;
	iter = g->beginTubeIterator();
	iend = g->endTubeIterator();
	value = *iter;
	minusGamma = -g->getGamma();
	plusGamma = g->getGamma();
    }

    template<class Real, class Index, class Grid, class Iterator>
    GridScalarField3D<Real, Index, Grid, Iterator>::~GridScalarField3D(void)
    {
	iter.commit();
	delete g;
    }

    template<class Real, class Index, class Grid, class Iterator>
    void GridScalarField3D<Real, Index, Grid, Iterator>::init(Grid *phi)
    {
    }


    template<class Real, class Index, class Grid, class Iterator>
    void GridScalarField3D<Real, Index, Grid, Iterator>::propagate(Real dt, Grid *phi)
    {
	iter.commit();
	iter = g->beginTubeIterator();
	iend = g->endTubeIterator();
	value = *iter;
    }


    template<class Real, class Index, class Grid, class Iterator>
    Real GridScalarField3D<Real, Index, Grid, Iterator>::computeHyperbolicTerm(const Iterator& otherIter)
    {
	Index i, j, k;

	if ( iter == iend )
	{
	    return minusGamma;
	}


	otherIter.getIndex(&i, &j, &k);


	if ( iter.getI() == i && iter.getJ() == j && iter.getK() == k )
	{
	    return -value;
	}
	else
	{
	    // we are lexicographically smaller or larger
	    while ( iter!=iend && iter.getI() < i )
	    {
		iter++;
	    }
	    if ( iter!=iend && iter.getI() == i )
	    {
		while ( iter!=iend && iter.getI() == i && iter.getJ() < j )
		{
		    iter++;
		}
		if ( iter!=iend && iter.getI() == i && iter.getJ() == j )
		{
		    while ( iter!=iend && iter.getI()==i && iter.getJ()==j && iter.getK() < k )
		    {
			iter++;
		    }
		    if ( iter!=iend && iter.getI() == i && iter.getJ() == j && iter.getK() == k )
		    {
			value = *iter;
			return -value;
		    }
		    else
		    {
			value = *iter;
			return ( value < 0 ? plusGamma : minusGamma );
		    }
		}
		else
		{
		    value = *iter;
		    return ( value < 0 ? plusGamma : minusGamma );
		}
	    }
	    else
	    {
		value = *iter;
		return ( value < 0 ? plusGamma : minusGamma );
	    }
	}

    }

    template<class Real, class Index, class Grid, class Iterator>
    Real GridScalarField3D<Real, Index, Grid, Iterator>::computeParabolicTerm(const Iterator& iter, Real kappa)
    {
	return 0;
    }

}
