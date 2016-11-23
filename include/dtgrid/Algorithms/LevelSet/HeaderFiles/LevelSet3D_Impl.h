/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/



namespace LevelSet
{


    template<typename Grid>  
    LevelSet3D<Grid>::LevelSet3D(Grid *phi, InitParams initParams, bool ownsPhi)
    {
	setDx(phi->getDx());
	setBeta(phi->getBeta());
	setGamma(phi->getGamma());
	setReinitMaxIter(initParams.numReinitIterations);
	setTDFormat(initParams.tdFormat);
	setSDFormat(initParams.sdFormat);
	setUseVanillaPeng(initParams.useVanillaPeng);
	setReinitCFLNumber(initParams.reinitCFLNumber);
	setPropagateCFLNumber(initParams.propagateCFLNumber);
	setVerbose(initParams.verbose);
	setMaxAbsDeterminationMethod(initParams.maxAbsDeterminationMethod);
	this->ownsPhi = ownsPhi;
	this->phi = phi;
	this->t = 0;
    }



    template<typename Grid>
    LevelSet3D<Grid>::LevelSet3D(LevelSet3D *ls)
    {
	setDx(ls->getDx());
	setBeta(ls->getBeta());
	setGamma(ls->getGamma());
	setReinitMaxIter(ls->getReinitMaxIter());
	setDoFreeBuffer2(ls->getDoFreeBuffer2());
	setTDFormat(ls->getTDFormat());
	setSDFormat(ls->getSDFormat());
	setUseVanillaPeng(ls->getUseVanillaPeng());
	setReinitCFLNumber(ls->getReinitCFLNumber());
	setPropagateCFLNumber(ls->getPropagateCFLNumber());
	setVerbose(ls->getVerbose());
	setMaxAbsDeterminationMethod(ls->getMaxAbsDeterminationMethod());
	setPhi(ls->getPhi);
	this->ownsPhi = false;
	this->t = 0;
    }




    template<typename Grid>  
    LevelSet3D<Grid>::~LevelSet3D()
    {
	if (ownsPhi) { delete phi; }
    }



    template<typename Grid>  
    typename LevelSet3D<Grid>::Real LevelSet3D<Grid>::averageMeanCurvature()
    {
	// compute average mean curvature
	// epsilon = 1.5*dx according to osher, fedkiw page 15
	typename Grid::template StencilTubeIterator<Grids::BETA_TUBE> iter = phi->template beginStencilTubeIterator<Grids::BETA_TUBE, Grids::SF_FIRSTORDER_CURVATURE, false>(0,-1);
	typename Grid::template StencilTubeIterator<Grids::BETA_TUBE> iend = phi->template endStencilTubeIterator<Grids::BETA_TUBE, Grids::SF_FIRSTORDER_CURVATURE, false>(0);
	Real epsilon = 1.5 * this->dx;
	Real volElement = Math::pow3(this->dx);
	Real deltaVal;
	Real c1 = 1.0 / (2.0 * epsilon);
	Real c2 = (Real)M_PI / epsilon;
	Real H;
	Real gradLen;
	Real area;
	Real cellArea;
	Real phiVal;
	Real averageMeanCurvature;

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

	averageMeanCurvature /= area;

	return averageMeanCurvature;
    }


    template<typename Grid>  
    typename LevelSet3D<Grid>::Real LevelSet3D<Grid>::area()
    {
	typename Grid::template StencilTubeIterator<Grids::BETA_TUBE> iter = phi->template beginStencilTubeIterator<Grids::BETA_TUBE, Grids::SF_FIRSTORDER, false>(0,-1);
	typename Grid::template StencilTubeIterator<Grids::BETA_TUBE> iend = phi->template endStencilTubeIterator<Grids::BETA_TUBE, Grids::SF_FIRSTORDER, false>(0);
	const Real epsilon = 1.5 * this->dx;
	const Real volElement = Math::pow3(this->dx);
	Real deltaVal;
	const Real c1 = 1.0 / (2.0 * epsilon);
	const Real c2 = 3.1415926535897931 / epsilon;
	Real gradLen;
	Real surfArea;
	Real cellArea;
	Real phiVal;

	// compute surface area, See fedkiw and Osher book.
	surfArea = 0;
	while (iter != iend)
	{
	    phiVal = iter.getValue();
	    if (fabs(phiVal) <= epsilon)
	    {
		deltaVal = c1 + c1 * cos( phiVal * c2 );
		gradLen = iter.gradientLength();
		cellArea = deltaVal * gradLen * volElement;
		surfArea += cellArea;
	    }

#ifdef WIN32
	    iter.operator++<Grids::SF_FIRSTORDER, false>();
#else
	    iter.template operator++<Grids::SF_FIRSTORDER, false>();
#endif
	}

	return surfArea;
    }


    template<typename Grid>  
    typename LevelSet3D<Grid>::Real LevelSet3D<Grid>::computeC(Real phi) const
    {
	Real phi_fabs = fabs(phi);

	if ( phi_fabs < this->beta )
	{
	    return 1;
	}
	else if ( phi_fabs < this->gamma )
	{
	    return Math::pow2(phi_fabs - this->gamma)*(2*phi_fabs + this->gamma - 3 * this->beta)/Math::pow3(this->gamma - this->beta); 
	}
	else
	{
	    return 0;
	}
    }



    template<typename Grid>  
    template<Grids::StencilFormat stencilFormat, class MyScalarField>
    void LevelSet3D<Grid>::propagate(MyScalarField *sf, const SimulationParameters& simParameters, Real maxTime)
    {
	if (simParameters.useSimBBox)
	{
	    simBBox[0][0] = simParameters.simBBox[0][0];
	    simBBox[0][1] = simParameters.simBBox[0][1];
	    simBBox[1][0] = simParameters.simBBox[1][0];
	    simBBox[1][1] = simParameters.simBBox[1][1];
	    simBBox[2][0] = simParameters.simBBox[2][0];
	    simBBox[2][1] = simParameters.simBBox[2][1];
	    bboxSmoothTransitionWidth[0] = simParameters.bboxSmoothTransitionWidth[0];
	    bboxSmoothTransitionWidth[1] = simParameters.bboxSmoothTransitionWidth[1];
	    propagate<stencilFormat, MyScalarField, true>(sf, maxTime);
	}
	else
	{
	    propagate<stencilFormat, MyScalarField, false>(sf, maxTime);
	}
    }




    template<typename Grid>  
    template<Grids::StencilFormat stencilFormat, class MyScalarField, bool useSimBBox>
    void LevelSet3D<Grid>::propagate(MyScalarField *sf, Real maxTime)
    {
	if (sdFormat == SD_WENO && ( stencilFormat != Grids::SF_WENO && stencilFormat != Grids::SF_WENO_CURVATURE ) )
	{	    
	    Core::throwDefaultException("sdFormat (SD_WENO) does not agree with stencilFormat!", __FILE__, __LINE__);
	}

	switch(LevelSet3D<Grid>::tdFormat)
	{
	case LevelSet3D<Grid>::TD_EULER:
	    {
		Real dt;
		Real dtMax = max((Real)0, maxTime- this->t);
		propagate<stencilFormat, MyScalarField, false, useSimBBox>(sf, 0, 0, dtMax, false, 0, &dt);
		sf->propagate(dt, phi);
		this->t += dt;
		this->dt = dt;
	    }
	    break;

	case LevelSet3D<Grid>::TD_TVD_RUNGE_KRUTTA:
	    {

		// TVD RK METHOD THAT KEEPS SCALARFIELD CONSTANT INBETWEEN GRID POINTS

		Real dt;
		Real dtMax = max((Real)0, maxTime- this->t);
		typename Grid::template StencilTubeIterator<Grids::GAMMA_TUBE> iter, iend;

		// first, step to time 't+dt'
		propagate<stencilFormat, MyScalarField, true, useSimBBox>(sf, 0, 1, dtMax, false, 0, &dt);

		// 0: t+dt
		// 1: t

		// second, step to time 't+2*dt'
		propagate<stencilFormat, MyScalarField, false, useSimBBox>(sf, 0, 0, dtMax, false, 0, &dt);

		// 0: t+2*dt
		// 1: t

		// third, averaging step to produce solution to time 't+0.5*dt'
		iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, false>(0,-1);
		iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, false>(0);
		while (iter != iend)
		{
		    Real phiVal = static_cast<Real>(0.75) * iter.getValue(1) + static_cast<Real>(0.25) * iter.getValue();
		    iter.setValue(phiVal);
#ifdef WIN32		
		    iter.operator++<Grids::SF_NONE, false>();
#else
		    iter.template operator++<Grids::SF_NONE, false>();
#endif
		}

		// 0: t+0.5*dt
		// 1: t

		// fourth, step to time 't+1.5*dt'
		propagate<stencilFormat, MyScalarField, false, useSimBBox>(sf, 0, 0, dtMax, false, 0, &dt);

		// 0: t+1.5*dt
		// 1: t

		// third, averaging step to produce solution to time 't'
		iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, false>(0);
		iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, false>(0);
		while (iter != iend)
		{
		    Real phiVal = static_cast<Real>((1.0/3.0)) * iter.getValue(1) + static_cast<Real>((2.0/3.0)) * iter.getValue();
		    iter.setValue(phiVal);
#ifdef WIN32
		    iter.operator++<Grids::SF_NONE, false>();
#else
		    iter.template operator++<Grids::SF_NONE, false>();
#endif
		}

		// 0: t+dt
		// 1: t

		sf->propagate(dt, phi);

		this->t += dt;
		this->dt = dt;
	    }
	    break;
	}

    }



    template<typename Grid>  
    template<Grids::StencilFormat stencilFormat, class MyScalarField, bool copyElements, bool useSimBBox>
    void LevelSet3D<Grid>::propagate(MyScalarField *sf, unsigned int b1, unsigned int b2, Real maxDt, bool enforceDt, Real enforcedDt, Real *dtCFL, Real *maxSV)
    {
	Real dt, v, t1, t2;
	Real phi_ijk, Dx_m, Dx_p, Dy_m, Dy_p, Dz_m, Dz_p, phi_x_2, phi_y_2, phi_z_2;

	if (sf->isParabolic())
	{

	    if (stencilFormat != Grids::SF_FIRSTORDER_CURVATURE && stencilFormat != Grids::SF_WENO_CURVATURE)
	    {
		Core::throwDefaultException("stencilFormat must be either SF_FIRSTORDER_CURVATURE or SF_WENO_CURVATURE to support parabolic propagation!", __FILE__, __LINE__);
	    }


	    if ( !sf->isHyperBolic() )
	    {
		// The ScalarField is only parabolic
		Real kappa, b;

		typename Grid::template StencilTubeIterator<Grids::GAMMA_TUBE> iter, iend;
		iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, stencilFormat, false>(b1);
		iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, stencilFormat, false>(b1);

		this->maxAbs = 0;

		reserveAuxiliaryBuffer();

		while (iter != iend)
		{
		    kappa = iter.template meanCurvature<stencilFormat>();
		    b = sf->computeParabolicTerm(iter, kappa);

		    // b > 0
		    if ( b > this->maxAbs)
		    {
			this->maxAbs = b;
		    }

		    storeAuxiliaryValue<useSimBBox>( b * kappa * sqrt( Math::pow2(iter.dc2x()) + Math::pow2(iter.dc2y()) + Math::pow2(iter.dc2z()) ), iter.getI(), iter.getJ(), iter.getK(), iter.getArrayIndex() );
#ifdef WIN32
		    iter.operator++<stencilFormat, false>();
#else
		    iter.template operator++<stencilFormat, false>();
#endif
		}

		// Osher, Fedkiw, eq 4.7, pg. 44
		dt = this->propagateCFLNumber * Math::pow2(this->dx) / (3 * 2 * this->maxAbs);


		if (enforceDt)
		{
		    dt = enforcedDt;
		}
		dt = min(dt, maxDt);
		if (dtCFL != NULL)
		{
		    *dtCFL = dt;
		}
		if (maxSV != NULL)
		{
		    *maxSV = this->maxAbs;
		}


		iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, copyElements>(b1, b2);
		iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, copyElements>(b1);

		if (this->useVanillaPeng)
		{
		    while (iter != iend)
		    {
			phi_ijk = iter.getValue(); 
			iter.setValue(phi_ijk + dt * computeC(phi_ijk) * getAuxiliaryValue(iter.getArrayIndex()), b2);
#ifdef WIN32
			iter.operator++<Grids::SF_NONE, copyElements>();
#else
			iter.template operator++<Grids::SF_NONE, copyElements>();
#endif
		    }
		}
		else
		{
		    while (iter != iend)
		    {
			phi_ijk = iter.getValue(); 
			iter.setValue(phi_ijk + dt * getAuxiliaryValue(iter.getArrayIndex()), b2);
#ifdef WIN32
			iter.operator++<Grids::SF_NONE, copyElements>();
#else
			iter.template operator++<Grids::SF_NONE, copyElements>();
#endif
		    }
		}

		freeAuxiliaryBuffer();

		if (b1 != b2)
		{
		    phi->swapBuffers(b1, b2);
		}

	    }
	    else
	    {
		// The Scalar Field is both parabolic and hyperbolic
		// The Scalar Field is a(x,t) - b(x,t)*kappa, where b >= 0
		typename Grid::template StencilTubeIterator<Grids::GAMMA_TUBE> iter, iend;
		Real hyperbolicTerm, parabolicTerm;
		Real maxAbsA, maxAbsB;
		Real kappa, a, b;

		iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, stencilFormat, false>(b1);
		iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, stencilFormat, false>(b1);

		maxAbsA = maxAbsB = 0;

		reserveAuxiliaryBuffer();

		while (iter != iend)
		{
		    kappa = iter.template meanCurvature<stencilFormat>();
		    a = sf->computeHyperbolicTerm(iter);
		    b = sf->computeParabolicTerm(iter, kappa);

		    // discretize hyperbolic term using the godunov scheme
		    switch (stencilFormat)
		    {
		    case Grids::SF_FIRSTORDER:
		    case Grids::SF_FIRSTORDER_CURVATURE:
			Dx_m = iter.dm1x();
			Dx_p = iter.dp1x();
			Dy_m = iter.dm1y();
			Dy_p = iter.dp1y();
			Dz_m = iter.dm1z();
			Dz_p = iter.dp1z();
			break;
		    case Grids::SF_WENO:
		    case Grids::SF_WENO_CURVATURE:
			Dx_m = iter.dm5x();
			Dx_p = iter.dp5x();
			Dy_m = iter.dm5y();
			Dy_p = iter.dp5y();
			Dz_m = iter.dm5z();
			Dz_p = iter.dp5z();
			break;
		    }

		    if ( a > 0 )
		    {
			t1 = max(Dx_m, (Real)0);
			t1 = t1*t1;
			t2 = min(Dx_p, (Real)0);
			t2 = t2*t2;
			phi_x_2 = max(t1,t2);

			t1 = max(Dy_m, (Real)0);
			t1 = t1*t1;
			t2 = min(Dy_p, (Real)0);
			t2 = t2*t2;
			phi_y_2 = max(t1,t2);

			t1 = max(Dz_m, (Real)0);
			t1 = t1*t1;
			t2 = min(Dz_p, (Real)0);
			t2 = t2*t2;
			phi_z_2 = max(t1,t2);

			if ( a > maxAbsA )
			{
			    maxAbsA = a;
			}
		    }
		    else
		    {
			t1 = min(Dx_m, (Real)0);
			t1 = t1*t1;
			t2 = max(Dx_p, (Real)0);
			t2 = t2*t2;
			phi_x_2 = max(t1,t2);

			t1 = min(Dy_m, (Real)0);
			t1 = t1*t1;
			t2 = max(Dy_p, (Real)0);
			t2 = t2*t2;
			phi_y_2 = max(t1,t2);

			t1 = min(Dz_m, (Real)0);
			t1 = t1*t1;
			t2 = max(Dz_p, (Real)0);
			t2 = t2*t2;
			phi_z_2 = max(t1,t2);

			if (-a > maxAbsA)
			{
			    maxAbsA = -a;
			}
		    }

		    hyperbolicTerm = a * sqrt(phi_x_2 + phi_y_2 + phi_z_2);


		    // discretize parabolic term using central differencing
		    parabolicTerm = b * kappa * sqrt( Math::pow2(iter.dc2x()) + Math::pow2(iter.dc2y()) + Math::pow2(iter.dc2z()) );

		    if (b > maxAbsB )
		    {
			maxAbsB = b;
		    }


		    // set value

		    storeAuxiliaryValue<useSimBBox>( hyperbolicTerm - parabolicTerm, iter.getI(), iter.getJ(), iter.getK(), iter.getArrayIndex() );
#ifdef WIN32
		    iter.operator++<stencilFormat, false>();
#else
		    iter.template operator++<stencilFormat, false>();
#endif
		}


		if (maxAbsA == 0)
		{
		    // Osher, Fedkiw eq. 4.7
		    dt = this->propagateCFLNumber * Math::pow2(this->dx) / ( 3 * 2 * maxAbsB);
		}
		else if (maxAbsB == 0)
		{
		    if (maxAbsDeterminationMethod == MAXABS_NORM1)
		    {
			// Osher, Fedkiw eq. 3.10
			dt = this->propagateCFLNumber * this->dx / ( 3 * maxAbsA );
		    }
		    else if (maxAbsDeterminationMethod == MAXABS_NORM2)
		    {
			// Osher, Fedkiw eq. 3.11
			dt = this->propagateCFLNumber * this->dx / maxAbsA ;
		    }
		}
		else
		{
		    if (maxAbsDeterminationMethod == MAXABS_NORM1)
		    {
			// Osher, Fedkiw eq. 4.12
			dt = this->propagateCFLNumber * ( Math::pow2(this->dx) / ( 3 * maxAbsA * this->dx + 3 * 2 * maxAbsB ) );
		    }
		    else
		    {
			// Osher, Fedkiw eq. 4.7 combined with eq. 3.11
			dt = this->propagateCFLNumber * ( Math::pow2(this->dx) / ( maxAbsA * this->dx + 3 * 2 * maxAbsB ) );
		    }
		}


		if (enforceDt)
		{
		    dt = enforcedDt;
		}
		dt = min(dt, maxDt);
		if (dtCFL != NULL)
		{
		    *dtCFL = dt;
		}

		if (maxSV != NULL)
		{
		    *maxSV = maxAbsA + maxAbsB;
		}


		iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, copyElements>(b1, b2);
		iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, copyElements>(b1);

		if (this->useVanillaPeng)
		{
		    while (iter != iend)
		    {
			phi_ijk = iter.getValue();
			iter.setValue( phi_ijk - computeC(phi_ijk) * dt * getAuxiliaryValue(iter.getArrayIndex()), b2 );
#ifdef WIN32
			iter.operator++<Grids::SF_NONE, copyElements>();
#else
			iter.template operator++<Grids::SF_NONE, copyElements>();
#endif
		    }
		}
		else
		{
		    while (iter != iend)
		    {
			phi_ijk = iter.getValue();
			iter.setValue( phi_ijk - dt * getAuxiliaryValue(iter.getArrayIndex()), b2 );
#ifdef WIN32
			iter.operator++<Grids::SF_NONE, copyElements>();
#else
			iter.template operator++<Grids::SF_NONE, copyElements>();
#endif
		    }
		}

		freeAuxiliaryBuffer();

		if (b1 != b2)
		{
		    phi->swapBuffers(b1, b2);
		}

	    }

	}
	else
	{
	    typename Grid::template StencilTubeIterator<Grids::GAMMA_TUBE> iter, iend;

	    iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, stencilFormat, false>(b1);
	    iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, stencilFormat, false>(b1);

	    reserveAuxiliaryBuffer();

	    this->maxAbs = 0;

	    while (iter != iend)
	    {
		v = sf->computeHyperbolicTerm(iter);

		switch (stencilFormat)
		{
		case Grids::SF_FIRSTORDER:
		case Grids::SF_FIRSTORDER_CURVATURE:
		    Dx_m = iter.dm1x();
		    Dx_p = iter.dp1x();
		    Dy_m = iter.dm1y();
		    Dy_p = iter.dp1y();
		    Dz_m = iter.dm1z();
		    Dz_p = iter.dp1z();
		    break;
		case Grids::SF_WENO:
		case Grids::SF_WENO_CURVATURE:
		    Dx_m = iter.dm5x();
		    Dx_p = iter.dp5x();
		    Dy_m = iter.dm5y();
		    Dy_p = iter.dp5y();
		    Dz_m = iter.dm5z();
		    Dz_p = iter.dp5z();
		    break;
		}

		if ( v > 0 )
		{
		    t1 = max(Dx_m, (Real)0);
		    t1 = t1*t1;
		    t2 = min(Dx_p, (Real)0);
		    t2 = t2*t2;
		    phi_x_2 = max(t1,t2);

		    t1 = max(Dy_m, (Real)0);
		    t1 = t1*t1;
		    t2 = min(Dy_p, (Real)0);
		    t2 = t2*t2;
		    phi_y_2 = max(t1,t2);

		    t1 = max(Dz_m, (Real)0);
		    t1 = t1*t1;
		    t2 = min(Dz_p, (Real)0);
		    t2 = t2*t2;
		    phi_z_2 = max(t1,t2);

		    if ( v > this->maxAbs )
		    {
			LevelSet3D<Grid>::maxAbs = v;
		    }
		}
		else
		{
		    t1 = min(Dx_m, (Real)0);
		    t1 = t1*t1;
		    t2 = max(Dx_p, (Real)0);
		    t2 = t2*t2;
		    phi_x_2 = max(t1,t2);

		    t1 = min(Dy_m, (Real)0);
		    t1 = t1*t1;
		    t2 = max(Dy_p, (Real)0);
		    t2 = t2*t2;
		    phi_y_2 = max(t1,t2);

		    t1 = min(Dz_m, (Real)0);
		    t1 = t1*t1;
		    t2 = max(Dz_p, (Real)0);
		    t2 = t2*t2;
		    phi_z_2 = max(t1,t2);

		    Real vabs = fabs(v);
		    if (vabs > this->maxAbs)
		    {
			this->maxAbs = vabs;
		    }
		}

		storeAuxiliaryValue<useSimBBox>( v * sqrt(phi_x_2 + phi_y_2 + phi_z_2), iter.getI(), iter.getJ(), iter.getK(), iter.getArrayIndex() );

#ifdef WIN32
		iter.operator++<stencilFormat, false>();
#else
		iter.template operator++<stencilFormat, false>();
#endif
	    }


	    if (maxAbsDeterminationMethod == MAXABS_NORM1)
	    {
		// Osher, Fedkiw eq. 3.10
		dt = this->propagateCFLNumber * this->dx / ( 3 * this->maxAbs );
	    }
	    else if (maxAbsDeterminationMethod == MAXABS_NORM2)
	    {
		// Osher, Fedkiw eq. 3.11
		dt = this->propagateCFLNumber * this->dx / ( this->maxAbs );
	    }


	    if (enforceDt)
	    {
		dt = enforcedDt;
	    }
	    dt = min(dt, maxDt);
	    if (dtCFL != NULL)
	    {
		*dtCFL = dt;
	    }

	    if (maxSV != NULL)
	    {
		*maxSV = this->maxAbs;
	    }


	    iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, copyElements>(b1, b2);
	    iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, copyElements>(b1);

	    if (this->useVanillaPeng)
	    {
		while (iter != iend)
		{
		    phi_ijk = iter.getValue();
		    iter.setValue( phi_ijk - computeC(phi_ijk) * dt * getAuxiliaryValue(iter.getArrayIndex()), b2 );
#ifdef WIN32
		    iter.operator++<Grids::SF_NONE, copyElements>();
#else
		    iter.template operator++<Grids::SF_NONE, copyElements>();
#endif
		}
	    }
	    else
	    {
		while (iter != iend)
		{
		    phi_ijk = iter.getValue();

		    iter.setValue( phi_ijk - dt * getAuxiliaryValue(iter.getArrayIndex()), b2 );
#ifdef WIN32
		    iter.operator++<Grids::SF_NONE, copyElements>();
#else
		    iter.template operator++<Grids::SF_NONE, copyElements>();
#endif
		}
	    }

	    freeAuxiliaryBuffer();

	    if (b1 != b2)
	    {
		phi->swapBuffers(b1, b2);
	    }
	}
    }





    template<typename Grid>  
    template<Grids::StencilFormat stencilFormat, class MyVelocityField>
    void LevelSet3D<Grid>::advect(MyVelocityField *vf, const SimulationParameters& simParameters, Real maxTime)
    {
	if (simParameters.useSimBBox)
	{
	    simBBox[0][0] = simParameters.simBBox[0][0];
	    simBBox[0][1] = simParameters.simBBox[0][1];
	    simBBox[1][0] = simParameters.simBBox[1][0];
	    simBBox[1][1] = simParameters.simBBox[1][1];
	    simBBox[2][0] = simParameters.simBBox[2][0];
	    simBBox[2][1] = simParameters.simBBox[2][1];
	    bboxSmoothTransitionWidth[0] = simParameters.bboxSmoothTransitionWidth[0];
	    bboxSmoothTransitionWidth[1] = simParameters.bboxSmoothTransitionWidth[1];
	    advect<stencilFormat, MyVelocityField, true>(vf, maxTime);
	}
	else
	{
	    advect<stencilFormat, MyVelocityField, false>(vf, maxTime);
	}
    }



    template<typename Grid>  
    template<Grids::StencilFormat stencilFormat, class MyVelocityField, bool useSimBBox>
    void LevelSet3D<Grid>::advect(MyVelocityField *vf, Real maxTime)
    {

	if (sdFormat == SD_WENO && ( stencilFormat != Grids::SF_WENO && stencilFormat != Grids::SF_WENO_CURVATURE ) )
	{
	    
	    Core::throwDefaultException("sdFormat (SD_WENO) does not agree with stencilFormat!", __FILE__, __LINE__);
	}


	switch(LevelSet3D<Grid>::tdFormat)
	{
	case LevelSet3D<Grid>::TD_EULER:
	    {
		Real dt;
		Real dtMax = max((Real)0, maxTime- this->t);
		advect<stencilFormat, MyVelocityField, false, useSimBBox>(vf, 0, 0, dtMax, false, 0, &dt);
		vf->advect(dt);
		this->t += dt;
		this->dt = dt;
	    }
	    break;

	case LevelSet3D<Grid>::TD_TVD_RUNGE_KRUTTA:
	    {

		// TVD RK METHOD THAT KEEPS VELOCITYFIELD CONSTANT INBETWEEN GRID POINTS

		Real dt;
		Real dtMax = max((Real)0, maxTime- this->t);
		typename Grid::template StencilTubeIterator<Grids::GAMMA_TUBE> iter, iend;

		// first, step to time 't+dt'
		advect<stencilFormat, MyVelocityField, true, useSimBBox>(vf, 0, 1, dtMax, false, 0, &dt);

		// 0: t+dt
		// 1: t

		// second, step to time 't+2*dt'
		advect<stencilFormat, MyVelocityField, false, useSimBBox>(vf, 0, 0, dtMax, false, 0, &dt);

		// 0: t+2*dt
		// 1: t

		// third, averaging step to produce solution to time 't+0.5*dt'
		iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, false>(0);
		iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, false>(0);
		while (iter != iend)
		{
		    Real phiVal = static_cast<Real>(0.75) * iter.getValue(1) + static_cast<Real>(0.25) * iter.getValue();
		    iter.setValue(phiVal);
#ifdef WIN32
		    iter.operator++<Grids::SF_NONE, false>();
#else
		    iter.template operator++<Grids::SF_NONE, false>();
#endif
		}

		// 0: t+0.5*dt
		// 1: t

		// fourth, step to time 't+1.5*dt'
		advect<stencilFormat, MyVelocityField, false, useSimBBox>(vf, 0, 0, dtMax, false, 0, &dt);

		// 0: t+0.5*dt
		// 1: t

		// third, averaging step to produce solution to time 't'
		// OPT: Just iterate over the flat value arrays
		iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, false>(0);
		iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, false>(0);
		while (iter != iend)
		{
		    Real phiVal = static_cast<Real>((1.0/3.0)) * iter.getValue(1) + static_cast<Real>((2.0/3.0)) * iter.getValue();
		    iter.setValue(phiVal);
#ifdef WIN32
		    iter.operator++<Grids::SF_NONE, false>();
#else
		    iter.template operator++<Grids::SF_NONE, false>();
#endif
		}

		// 0: t+dt
		// 1: t

		vf->advect(dt);

		this->t += dt;
		this->dt = dt;
	    }
	    break;
	}

    }




    template<typename Grid>
    template<Grids::StencilFormat stencilFormat, class MyVelocityField, bool copyElements, bool useSimBBox>
    void LevelSet3D<Grid>::advect(MyVelocityField *vf, unsigned int b1, unsigned int b2, Real maxDt, bool enforceDt, Real enforcedDt, Real *dtCFL)
    {
	Real dt;

	Index i, j, k;
	Real v[3];
	Real v_t_phi_x, v_t_phi_y, v_t_phi_z, phi_ijk;
	this->maxAbs = 0;
	Real vSq;
	typename Grid::template StencilTubeIterator<Grids::GAMMA_TUBE> iter, iend;

	iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, stencilFormat, false>(b1);
	iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, stencilFormat, false>(b1);

	reserveAuxiliaryBuffer();

	while (iter != iend)
	{
	    phi_ijk = iter.getValue();
	    iter.getIndex(&i, &j, &k);
	    (*vf)(i, j, k, v);	    

	    // compute differences in x direction
	    if ( v[0] < 0 )
	    {
		switch(stencilFormat)
		{
		case Grids::SF_FIRSTORDER:
		case Grids::SF_FIRSTORDER_CURVATURE:
		    v_t_phi_x = v[0] * iter.dp1x();
		    break;
		case Grids::SF_WENO:
		case Grids::SF_WENO_CURVATURE:
		    v_t_phi_x = v[0] * iter.dp5x();
		    break;
		}
	    }
	    else if ( v[0] > 0)
	    {
		switch(stencilFormat)
		{
		case Grids::SF_FIRSTORDER:
		case Grids::SF_FIRSTORDER_CURVATURE:
		    v_t_phi_x = v[0] * iter.dm1x();
		    break;
		case Grids::SF_WENO:
		case Grids::SF_WENO_CURVATURE:
		    v_t_phi_x = v[0] * iter.dm5x();
		    break;
		}
	    }
	    else
	    {
		v_t_phi_x = 0;
	    }

	    // compute differences in y direction
	    if ( v[1] < 0 )
	    {
		switch(stencilFormat)
		{
		case Grids::SF_FIRSTORDER:
		case Grids::SF_FIRSTORDER_CURVATURE:
		    v_t_phi_y = v[1] * iter.dp1y();
		    break;
		case Grids::SF_WENO:
		case Grids::SF_WENO_CURVATURE:
		    v_t_phi_y = v[1] * iter.dp5y();
		    break;
		}
	    }
	    else if ( v[1] > 0)
	    {
		switch(stencilFormat)
		{
		case Grids::SF_FIRSTORDER:
		case Grids::SF_FIRSTORDER_CURVATURE:
		    v_t_phi_y = v[1] * iter.dm1y();
		    break;
		case Grids::SF_WENO:
		case Grids::SF_WENO_CURVATURE:
		    v_t_phi_y = v[1] * iter.dm5y();
		    break;
		}
	    }
	    else
	    {
		v_t_phi_y = 0;
	    }

	    // compute differences in z direction
	    if ( v[2] < 0 )
	    {
		switch(stencilFormat)
		{
		case Grids::SF_FIRSTORDER:
		case Grids::SF_FIRSTORDER_CURVATURE:
		    v_t_phi_z = v[2] * iter.dp1z();
		    break;
		case Grids::SF_WENO:
		case Grids::SF_WENO_CURVATURE:
		    v_t_phi_z = v[2] * iter.dp5z();
		    break;
		}
	    }
	    else if ( v[2] > 0)
	    {
		switch(stencilFormat)
		{
		case Grids::SF_FIRSTORDER:
		case Grids::SF_FIRSTORDER_CURVATURE:
		    v_t_phi_z = v[2] * iter.dm1z();
		    break;
		case Grids::SF_WENO:
		case Grids::SF_WENO_CURVATURE:
		    v_t_phi_z = v[2] * iter.dm5z();
		    break;
		}
	    }
	    else
	    {
		v_t_phi_z = 0;
	    }



	    storeAuxiliaryValue<useSimBBox>( (v_t_phi_x + v_t_phi_y + v_t_phi_z), iter.getI(), iter.getJ(), iter.getK(), iter.getArrayIndex() );


	    // determine max velocity for CFL condition
	    switch (maxAbsDeterminationMethod)
	    {
	    case MAXABS_NORM1:
		{
		    if ( (vSq = fabs(v[0]) + fabs(v[1]) + fabs(v[2])) > this->maxAbs )
		    {
			this->maxAbs = vSq;
		    }
		}
		break;
	    case MAXABS_NORM2:
		{
		    vSq = Math::pow2(v[0]) + Math::pow2(v[1]) + Math::pow2(v[2]);
		    if (vSq > this->maxAbs)
		    {
			this->maxAbs = vSq;
		    }
		}
		break;
	    }



#ifdef WIN32
	    iter.operator++<stencilFormat, false>();
#else
	    iter.template operator++<stencilFormat, false>();
#endif
	}

	// compute dt based on CFL condition
	// we choose a conservative CFL number of alpha = 0.9
	// osher, fedkiw eq. 3.11

	if (maxAbsDeterminationMethod == MAXABS_NORM2)
	{
	    this->maxAbs = sqrt(this->maxAbs);
	}

	if (this->maxAbs == 0)
	{
	    dt = 0;
	}
	else
	{
	    dt = this->propagateCFLNumber * this->dx / this->maxAbs;
	}

	if (enforceDt)
	{
	    dt = enforcedDt;
	}
	dt = min(dt, maxDt);
	if (dtCFL != NULL)
	{
	    *dtCFL = dt;
	}

	if (verbose)
	    std::cout << "time = " << this->t << ", maxAbs = " << this->maxAbs << std::endl;

	iter = phi->template beginStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, copyElements>(b1, b2);
	iend = phi->template endStencilTubeIterator<Grids::GAMMA_TUBE, Grids::SF_NONE, copyElements>(b1);

	if (this->useVanillaPeng)
	{
	    while (iter != iend)
	    {
		phi_ijk = iter.getValue();
		iter.setValue( phi_ijk - dt * computeC(phi_ijk) * getAuxiliaryValue(iter.getArrayIndex()), b2 );
#ifdef WIN32
		iter.operator++<Grids::SF_NONE, copyElements>();
#else
		iter.template operator++<Grids::SF_NONE, copyElements>();
#endif
	    }
	}
	else
	{
	    while (iter != iend)
	    {
		phi_ijk = iter.getValue();
		iter.setValue( phi_ijk - dt * getAuxiliaryValue(iter.getArrayIndex()), b2 );
#ifdef WIN32
		iter.operator++<Grids::SF_NONE, copyElements>();
#else
		iter.template operator++<Grids::SF_NONE, copyElements>();
#endif
	    }
	}

	freeAuxiliaryBuffer();

	if (b1 != b2)
	{
	    phi->swapBuffers(b1, b2);
	}
    }



    template<typename Grid>  template<Grids::StencilFormat stencilFormat>
    void LevelSet3D<Grid>::reinitialize(const SimulationParameters& simParameters, unsigned int reinitMaxIter)
    {
	if (simParameters.useSimBBox)
	{
	    simBBox[0][0] = simParameters.simBBox[0][0];
	    simBBox[0][1] = simParameters.simBBox[0][1];
	    simBBox[1][0] = simParameters.simBBox[1][0];
	    simBBox[1][1] = simParameters.simBBox[1][1];
	    simBBox[2][0] = simParameters.simBBox[2][0];
	    simBBox[2][1] = simParameters.simBBox[2][1];
	    bboxSmoothTransitionWidth[0] = simParameters.bboxSmoothTransitionWidth[0];
	    bboxSmoothTransitionWidth[1] = simParameters.bboxSmoothTransitionWidth[1];
	    reinitialize<stencilFormat, true>(reinitMaxIter);
	}
	else
	{
	    reinitialize<stencilFormat, false>(reinitMaxIter);
	}
    }

	
    template<typename Grid>  template<Grids::StencilFormat stencilFormat, bool useSimBBox>
    void LevelSet3D<Grid>::reinitialize(unsigned int reinitMaxIter)
    {    
	if (sdFormat == SD_WENO && ( stencilFormat != Grids::SF_WENO && stencilFormat != Grids::SF_WENO_CURVATURE ) )
	{	    
	    Core::throwDefaultException("sdFormat (SD_WENO) does not agree with stencilFormat!", __FILE__, __LINE__);
	}


	unsigned int i;

	if (reinitMaxIter == -1)
	{
	    reinitMaxIter = this->reinitMaxIter;
	}

	reserveAuxiliaryBuffer();


	switch(LevelSet3D<Grid>::tdFormat)
	{
	case LevelSet3D<Grid>::TD_EULER:
	    {
		for (i=0; i<reinitMaxIter; i++)
		{
		    reinitialize<stencilFormat, false, false, useSimBBox>();
		}
	    }
	    break;

	case LevelSet3D<Grid>::TD_TVD_RUNGE_KRUTTA:
	    {
		Real phiVal;
		Real *sgn = new Real[phi->getNumValues()];
		typename Grid::template StencilTubeIterator<Grids::ENTIRE_TUBE> iter, iend;

		for (i=0; i<reinitMaxIter; i++)
		{
		    reinitialize<stencilFormat, true, true, useSimBBox>(0, 1, sgn);
		    // 0: t+dt
		    // 1: t

		    reinitialize<stencilFormat, true, false, useSimBBox>(0, 0, sgn);
		    // 0: t+2*dt
		    // 1: t

		    iter = phi->template beginStencilTubeIterator<Grids::ENTIRE_TUBE, Grids::SF_NONE, false>(0,-1);
		    iend = phi->template endStencilTubeIterator<Grids::ENTIRE_TUBE, Grids::SF_NONE, false>(0);
		    while (iter != iend)
		    {
			phiVal = static_cast<Real>(0.75) * iter.getValue(1) + static_cast<Real>(0.25) * iter.getValue();
			iter.setValue(phiVal);
#ifdef WIN32
			iter.operator++<Grids::SF_NONE, false>();
#else
			iter.template operator++<Grids::SF_NONE, false>();
#endif
		    }
		    // 0: t+0.5*dt
		    // 1: t

		    reinitialize<stencilFormat, true, false, useSimBBox>(0, 0, sgn);
		    // 0: t+0.5*dt
		    // 1: t

		    iter = phi->template beginStencilTubeIterator<Grids::ENTIRE_TUBE, Grids::SF_NONE, false>(0);
		    iend = phi->template endStencilTubeIterator<Grids::ENTIRE_TUBE, Grids::SF_NONE, false>(0);
		    while (iter != iend)
		    {
			phiVal = static_cast<Real>((1.0/3.0)) * iter.getValue(1) + static_cast<Real>((2.0/3.0)) * iter.getValue();
			iter.setValue(phiVal);
#ifdef WIN32
			iter.operator++<Grids::SF_NONE, false>();
#else
			iter.template operator++<Grids::SF_NONE, false>();
#endif
		    }
		    // 0: t+dt
		    // 1: t

		}

		delete[] sgn;

	    }
	    break;
	}


	freeAuxiliaryBuffer();
    }









    template<typename Grid>  
    template<Grids::StencilFormat stencilFormat, bool useSgnBuffer, bool initSgnBuffer, bool useSimBBox>
    void LevelSet3D<Grid>::reinitialize(unsigned int b1, unsigned int b2, Real *sgn)
    {    
	Real dt;
	Real dx2 = this->dx * this->dx;
	Real Dx_m, Dx_p, Dy_m, Dy_p, Dz_m, Dz_p, s, t1, t2, phi_x_2, phi_y_2, phi_z_2, phi_ijk, dPhiLen2;
	UInt arrayIndex;

	typename Grid::template StencilTubeIterator<Grids::ENTIRE_TUBE> iter, iend;

	iter = phi->template beginStencilTubeIterator<Grids::ENTIRE_TUBE, stencilFormat, false>(b1);
	iend = phi->template endStencilTubeIterator<Grids::ENTIRE_TUBE, stencilFormat, false>(b1);

	arrayIndex = 0;

	while (iter != iend)
	{
	    if (stencilFormat == Grids::SF_FIRSTORDER || stencilFormat == Grids::SF_FIRSTORDER_CURVATURE)
	    {
		Dx_m = iter.dm1x();
		Dx_p = iter.dp1x();
		Dy_m = iter.dm1y();
		Dy_p = iter.dp1y();
		Dz_m = iter.dm1z();
		Dz_p = iter.dp1z();
	    }

	    if (stencilFormat == Grids::SF_WENO || stencilFormat == Grids::SF_WENO_CURVATURE)
	    {
		Dx_m = iter.dm5x();
		Dx_p = iter.dp5x();
		Dy_m = iter.dm5y();
		Dy_p = iter.dp5y();
		Dz_m = iter.dm5z();
		Dz_p = iter.dp5z();
	    }


	    phi_ijk = iter.getValue();


	    if ( phi_ijk > 0 )
	    {
		t1 = max(Dx_m, (Real)0);
		t1 = t1*t1;
		t2 = min(Dx_p, (Real)0);
		t2 = t2*t2;
		phi_x_2 = max(t1,t2);

		t1 = max(Dy_m, (Real)0);
		t1 = t1*t1;
		t2 = min(Dy_p, (Real)0);
		t2 = t2*t2;
		phi_y_2 = max(t1,t2);

		t1 = max(Dz_m, (Real)0);
		t1 = t1*t1;
		t2 = min(Dz_p, (Real)0);
		t2 = t2*t2;
		phi_z_2 = max(t1,t2);
	    }
	    else
	    {
		t1 = min(Dx_m, (Real)0);
		t1 = t1*t1;
		t2 = max(Dx_p, (Real)0);
		t2 = t2*t2;
		phi_x_2 = max(t1,t2);

		t1 = min(Dy_m, (Real)0);
		t1 = t1*t1;
		t2 = max(Dy_p, (Real)0);
		t2 = t2*t2;
		phi_y_2 = max(t1,t2);

		t1 = min(Dz_m, (Real)0);
		t1 = t1*t1;
		t2 = max(Dz_p, (Real)0);
		t2 = t2*t2;
		phi_z_2 = max(t1,t2);
	    }

	    dPhiLen2 = phi_x_2 + phi_y_2 + phi_z_2; 

	    if (useSgnBuffer)
	    {
		if (initSgnBuffer)
		{
		    s = sgn[arrayIndex] = phi_ijk / ( sqrt( Math::pow2(phi_ijk) + dPhiLen2 * dx2) );		
		}
		else
		{
		    s = sgn[arrayIndex];
		}
	    }
	    else
	    {
		s = phi_ijk / ( sqrt( Math::pow2(phi_ijk) + dPhiLen2 * dx2) );
	    }

	    storeAuxiliaryValue<useSimBBox>( s * sqrt(dPhiLen2) - s, iter.getI(), iter.getJ(), iter.getK(), iter.getArrayIndex() );

#ifdef WIN32
	    iter.operator++<stencilFormat, false>();
#else
	    iter.template operator++<stencilFormat, false>();
#endif
	    arrayIndex++;
	}

	// Osher, Fedkiw eq. 3.11
	dt = this->reinitCFLNumber * this->dx;

	iter = phi->template beginStencilTubeIterator<Grids::ENTIRE_TUBE, Grids::SF_NONE, false>(b1, b2);
	iend = phi->template endStencilTubeIterator<Grids::ENTIRE_TUBE, Grids::SF_NONE, false>(b1);

	while (iter != iend)
	{
	    phi_ijk = iter.getValue();

	    iter.setValue( phi_ijk - dt * getAuxiliaryValue(iter.getArrayIndex()), b2 );
#ifdef WIN32
	    iter.operator++<Grids::SF_NONE, false>();
#else
	    iter.template operator++<Grids::SF_NONE, false>();
#endif
	}


	if (b1 != b2)
	{
	    phi->swapBuffers(b1, b2);
	}

    }





    template<typename Grid>  
    void LevelSet3D<Grid>::rebuild()
    {
	phi->rebuildNarrowBand();
	phi->freeBuffer(1);
    }



    template<typename Grid>  
    void LevelSet3D<Grid>::reserveAuxiliaryBuffer()
    {
	auxiliaryBuffer = new Real[phi->getNumValues()];
    }

    template<typename Grid>  
    void LevelSet3D<Grid>::freeAuxiliaryBuffer()
    {
	delete[] auxiliaryBuffer;
    }

    template<typename Grid>  template<bool useSimBBox>
    void LevelSet3D<Grid>::storeAuxiliaryValue(Real value, Index x, Index y, Index z, UInt i)
    {
	if (useSimBBox)
	{
	    Index dist;
	    dist = std::min(static_cast<Index>(abs(x-simBBox[0][0])), static_cast<Index>(abs(x-simBBox[0][1])));
	    dist = std::min(static_cast<Index>(abs(y-simBBox[1][0])), dist);
	    dist = std::min(static_cast<Index>(abs(y-simBBox[1][1])), dist);
	    dist = std::min(static_cast<Index>(abs(z-simBBox[2][0])), dist);
	    dist = std::min(static_cast<Index>(abs(z-simBBox[2][1])), dist);

	    if (dist < bboxSmoothTransitionWidth[0])
	    {
		value = 0;
	    }
	    else if (dist < bboxSmoothTransitionWidth[1])
	    {
		value *= static_cast<Real>(dist-bboxSmoothTransitionWidth[0])/(bboxSmoothTransitionWidth[1]-bboxSmoothTransitionWidth[0]);
	    }
	}
	auxiliaryBuffer[i] = value;
    }

    template<typename Grid>  
    typename LevelSet3D<Grid>::Real LevelSet3D<Grid>::getAuxiliaryValue(UInt i)
    {
	return auxiliaryBuffer[i];
    }


    template<typename Grid>  
    template<Grids::StencilFormat advectStencilFormat, Grids::StencilFormat reinitializeStencilFormat, class MyVelocityField>
    void LevelSet3D<Grid>::advectDt(MyVelocityField *vf, const SimulationParameters& simParameters)
    {
	unsigned int i;

	// we need to rebuild before we start the level set advection
	if (simParameters.numIterations > 0)
	{
	    if (verbose) { std::cout << " Rebuilding Initial Level Set..."; std::cout.flush(); }
	    rebuild();
	    std::cout << " done!" << std::endl;
	}

	if (simParameters.useSimBBox)
	{
	    simBBox[0][0] = simParameters.simBBox[0][0];
	    simBBox[0][1] = simParameters.simBBox[0][1];
	    simBBox[1][0] = simParameters.simBBox[1][0];
	    simBBox[1][1] = simParameters.simBBox[1][1];
	    simBBox[2][0] = simParameters.simBBox[2][0];
	    simBBox[2][1] = simParameters.simBBox[2][1];
	    bboxSmoothTransitionWidth[0] = simParameters.bboxSmoothTransitionWidth[0];
	    bboxSmoothTransitionWidth[1] = simParameters.bboxSmoothTransitionWidth[1];
	}

	for (i=0; i<simParameters.numIterations && getTime()<simParameters.maxTime && phi->getNumValues()>0; i++)
	{
	    if (verbose)
	    {
		std::cout << "Iteration " << i << std::endl;
		std::cout << " NumValues = " << phi->getNumValues() << ", Time = " << getTime() << std::endl;
	    }

	    // Time-substepping
	    Real nextTime = getTime() + simParameters.dt;
	    while ( getTime() < nextTime && getTime()<simParameters.maxTime && phi->getNumValues()>0)
	    {
		// ADVECT
		if (verbose) { std::cout << " Advecting..."; std::cout.flush(); }
		if (simParameters.useSimBBox)
		    advect<advectStencilFormat, MyVelocityField, true>(vf, nextTime);
		else
		    advect<advectStencilFormat, MyVelocityField, false>(vf, nextTime);
		if (verbose) { std::cout << "done!" << std::endl; }	

		// REINITIALIZE
		if (verbose) { std::cout << " Reinitializing..."; std::cout.flush(); }
		if (simParameters.useSimBBox)
		    reinitialize<reinitializeStencilFormat, true>();
		else
		    reinitialize<reinitializeStencilFormat, false>();
		if (verbose) { std::cout << "done!" << std::endl; }	

		// REBUILD
		if (verbose) { std::cout << " Rebuilding..."; std::cout.flush(); }
		rebuild();
		if (verbose) { std::cout << "done!" << std::endl; }	
	    }

	    if (simParameters.saveIntermediateSvols && (i % simParameters.saveIntermediateSvolInterval == 0))
	    {
		// Save the Result Svol (Sparse Volume)
		std::stringstream ss;
		ss << simParameters.outputBaseFileName << "_iter_" << i << "_time_" << getTime() << ".svol";
		if (verbose) { std::cout << "Saving the intermediate svol in file " << ss.str() << " ..."; }
		phi->saveSparseVolume(ss.str());
		if (verbose) { std::cout << " done!" << std::endl; }
	    }
	}

    }




    template<typename Grid>  
    template<Grids::StencilFormat propagateStencilFormat, Grids::StencilFormat reinitializeStencilFormat, class MyScalarField>
    void LevelSet3D<Grid>::propagateDt(MyScalarField *sf, const SimulationParameters& simParameters)
    {
	unsigned int i;

	// we need to rebuild before we start the level set propagation
	if (simParameters.numIterations > 0)
	{
	    if (verbose) { std::cout << " Rebuilding Initial Level Set..."; std::cout.flush(); }
	    rebuild();
	    std::cout << " done!" << std::endl;
	}

	if (simParameters.useSimBBox)
	{
	    simBBox[0][0] = simParameters.simBBox[0][0];
	    simBBox[0][1] = simParameters.simBBox[0][1];
	    simBBox[1][0] = simParameters.simBBox[1][0];
	    simBBox[1][1] = simParameters.simBBox[1][1];
	    simBBox[2][0] = simParameters.simBBox[2][0];
	    simBBox[2][1] = simParameters.simBBox[2][1];
	    bboxSmoothTransitionWidth[0] = simParameters.bboxSmoothTransitionWidth[0];
	    bboxSmoothTransitionWidth[1] = simParameters.bboxSmoothTransitionWidth[1];
	}

	for (i=0; i<simParameters.numIterations && getTime()<simParameters.maxTime && phi->getNumValues()>0; i++)
	{
	    if (verbose)
	    {
		std::cout << "Iteration " << i << std::endl;
		std::cout << " NumValues = " << phi->getNumValues() << ", Time = " << getTime() << std::endl;
	    }

	    // Time-substepping
	    Real nextTime = getTime() + simParameters.dt;
	    while ( getTime() < nextTime && getTime()<simParameters.maxTime && phi->getNumValues()>0 )
	    {
		// PROPAGATE
		if (verbose) { std::cout << " Time = " << getTime() << std::endl; }
		if (verbose) { std::cout << " Propagating..."; std::cout.flush(); }
		if (simParameters.useSimBBox)
		    propagate<propagateStencilFormat, MyScalarField, true>(sf, nextTime);
		else
		    propagate<propagateStencilFormat, MyScalarField, false>(sf, nextTime);
		if (verbose) { std::cout << "done!" << std::endl; }	

		// REINITIALIZE
		if (verbose) { std::cout << " Reinitializing..."; std::cout.flush(); }
		if (simParameters.useSimBBox)
		    reinitialize<reinitializeStencilFormat, true>();
		else
		    reinitialize<reinitializeStencilFormat, false>();
		if (verbose) { std::cout << "done!" << std::endl; }	

		// REBUILD
		if (verbose) { std::cout << " Rebuilding..."; std::cout.flush(); }
		rebuild();
		if (verbose) { std::cout << "done!" << std::endl; }	
	    }

	    if (simParameters.saveIntermediateSvols && (i % simParameters.saveIntermediateSvolInterval == 0))
	    {
		// Save the Result Svol (Sparse Volume)
		std::stringstream ss;
		ss << simParameters.outputBaseFileName << "_iter_" << i << "_time_" << getTime() << ".svol";
		if (verbose) { std::cout << "Saving the intermediate svol in file " << ss.str() << " ..."; }
		phi->saveSparseVolume(ss.str());
		if (verbose) { std::cout << " done!" << std::endl; }
	    }
	}

    }





    template<typename Grid>  
    LevelSet3D<Grid>::SimulationParameters::SimulationParameters()
    {
	numIterations = 0;
	saveIntermediateSvols = false;
	saveIntermediateSvolInterval = 1;
	outputBaseFileName = "outputBaseFileName";
	dt = 1;
	maxTime = std::numeric_limits<Real>::max();
	useSimBBox = false;
    }


    template<typename Grid>  
    LevelSet3D<Grid>::InitParams::InitParams()
    {
	maxAbsDeterminationMethod = MAXABS_NORM1;
	useVanillaPeng = false;
	verbose = true;
    }

}


