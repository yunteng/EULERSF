/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

namespace LevelSet
{

    template<class Real, class Index>
    const Real EnrightTestVelocityField3D<Real, Index>::T = 3;

    template<class Real, class Index>
    EnrightTestVelocityField3D<Real, Index>::EnrightTestVelocityField3D(Real dx)
    {
	t = 0;
	dxr = 1.0 / dx;
	dyr = 1.0 / dx;
	dzr = 1.0 / dx;
	c1 =  (Real)M_PI / T;
    }

    template<class Real, class Index>
    EnrightTestVelocityField3D<Real, Index>::EnrightTestVelocityField3D(Index dim[3])
    {
	this->dim[0] = dim[0];
	this->dim[1] = dim[1];
	this->dim[2] = dim[2];
	t = 0;
	dxr = static_cast<Real>(1.0) / dim[0];
	dyr = static_cast<Real>(1.0) / dim[1];
	dzr = static_cast<Real>(1.0) / dim[2];
	c1 =  (Real)M_PI / T;
    }

    template<class Real, class Index>
    void EnrightTestVelocityField3D<Real, Index>::advect(Real dt)
    {
	t += dt;
    }

    template<class Real, class Index>
    void EnrightTestVelocityField3D<Real, Index>::operator()(const Index p[3], Real v[3]) const
    {
	(*this)(p[0], p[1], p[2], v);
    }

    template<class Real, class Index>
    void EnrightTestVelocityField3D<Real, Index>::operator()(Index x, Index y, Index z, Real v[3]) const
    {
	Real xf = x * dxr;
	Real yf = y * dyr;
	Real zf = z * dzr;

	if (t < T)
	{
	    Real tr = cos(t*c1);
	    Real a = sin(2*(Real)M_PI*yf);
	    Real b = -sin(2*(Real)M_PI*xf);
	    Real c = sin(2*(Real)M_PI*zf);
	    v[0] = tr * ( 2 * Math::pow2(sin((Real)M_PI*xf)) * a * c );
	    v[1] = tr * ( b * Math::pow2(sin((Real)M_PI*yf)) * c );
	    v[2] = tr * ( b * a * Math::pow2(sin((Real)M_PI*zf)) );
	}
	else
	{
	    v[0] = 0;
	    v[1] = 0;
	    v[2] = 0;
	}
    }

}
