/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
namespace LevelSet
{

    template<class Real, class Index>
    ConstantUniformVelocityField3D<Real, Index>::ConstantUniformVelocityField3D(Real v[3])
    {
	this->v[0] = v[0];
	this->v[1] = v[1];
	this->v[2] = v[2];
	this->maxAbs = length(v[0], v[1], v[2]);
    }

    template<class Real, class Index>
    void ConstantUniformVelocityField3D<Real, Index>::advect(Real dt)
    {
    }

    template<class Real, class Index>
    void ConstantUniformVelocityField3D<Real, Index>::operator()(const Index p[3], Real v[3]) const
    {
	(*this)(p[0], p[1], p[2], v);
    }

    template<class Real, class Index>
    void ConstantUniformVelocityField3D<Real, Index>::operator()(Index x, Index y, Index z, Real v[3]) const
    {
	v[0] = this->v[0];
	v[1] = this->v[1];
	v[2] = this->v[2];
    }

    template<class Real, class Index>
    Real ConstantUniformVelocityField3D<Real, Index>::getMaxAbs() const
    {
	return this->maxAbs;
    }

}
