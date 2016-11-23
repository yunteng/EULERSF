/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _levelset_levelset3d_h
#define _levelset_levelset3d_h

#include <sstream>
#include <string>

#include <Core/Exception/HeaderFiles/DefaultException.h>
#include <DataStructures/Grids/HeaderFiles/Grids.h>
#include <Algorithms/LevelSet/VelocityFields/HeaderFiles/VelocityField3D.h>
#include <Algorithms/LevelSet/ScalarFields/HeaderFiles/ScalarField3D.h>
#include <Algorithms/Math/HeaderFiles/BasicMath.h>


namespace LevelSet
{

    /**
     *   The LevelSet3D class implements a basic finite difference based level set solver including:
     *	 - Advection in an externally generated velocity field
     *   - Motion in the normal direction
     *   - Motion involving mean curvature
     *   - First order, second order and HJ WENO spatial derivative approximations
     *   - Forward Euler and TVD Runge-Kutta temporal derivative approximations
     */
    template<class Grid> 
    class LevelSet3D
    {
    public:

	typedef typename Grid::Real Real;
	typedef typename Grid::Index Index;
	typedef typename Grid::UInt UInt;

	// Temporal Discretization
	enum TDFormat 
	{
	    TD_EULER, TD_TVD_RUNGE_KRUTTA
	};

	// Spatial Discretization
	enum SDFormat 
	{
	    SD_FIRSTORDER, SD_WENO
	};

	enum MaxAbsDeterminationMethod
	{
	    /** Compute maximal velocity in CFL condition using L1 norm, Osher and Fedkiw eq. (3.10) 
	      * This estimate is more conservative than MAXABS_NORM2. */
	    MAXABS_NORM1,
	    /** Compute maximal velocity in CFL condition using L2 norm, Osher and Fedkiw eq. (3.11) */
	    MAXABS_NORM2
	};


	struct SimulationParameters
	{
	    SimulationParameters();

	    unsigned int numIterations;
	    bool saveIntermediateSvols;
	    unsigned int saveIntermediateSvolInterval;
	    std::string outputBaseFileName;
	    Real dt;
	    Real maxTime;
	    /** If useSimBBox is true, the simulation only takes place within the bbox specified
	     *  and the change to the level set goes to zero at the boundary of the bbox, within a
	     *  range of width 'bboxSmoothTransitionWidth'. */
	    bool useSimBBox;
	    Index simBBox[3][2];
	    /** voxels a distance less than bboxSmoothTransitionWidth[0] away from the border are not modified.
	     *  voxels a distance of more than bboxSmoothTransitionWidth[1] away from the border are fully modified. */
	    Index bboxSmoothTransitionWidth[2];
	};




	struct InitParams
	{
	    InitParams();

	    unsigned int numReinitIterations;
	    bool useVanillaPeng;
	    Real reinitCFLNumber;
	    Real propagateCFLNumber;
	    MaxAbsDeterminationMethod maxAbsDeterminationMethod;
	    bool verbose;
	    TDFormat tdFormat;
	    SDFormat sdFormat;
	};


    public:
	LevelSet3D(Grid *phi, InitParams initParams, bool ownsPhi = true);
	LevelSet3D(LevelSet3D *ls);

	~LevelSet3D();


	//////////////////////////////////////
	// Level Set Operations
	//////////////////////////////////////


	// ADVECT

	/** This method only performs advection */
	template<Grids::StencilFormat stencilFormat, class MyVelocityField>
	void advect(MyVelocityField *vf, const SimulationParameters& simParameters, Real maxTime=std::numeric_limits<Real>::max());

	/** This method performs advection, reinitialization and narrow band rebuild */
	template<Grids::StencilFormat advectStencilFormat, Grids::StencilFormat reinitializeStencilFormat, class MyVelocityField>
	void advectDt(MyVelocityField *vf, const SimulationParameters& simParameters);


	// PROPAGATE IN NORMAL DIRECTION    

	/** This method only performs propagation */
	template<Grids::StencilFormat stencilFormat, class MyScalarField>
	void propagate(MyScalarField *sf, const SimulationParameters& simParameters, Real maxTime=std::numeric_limits<Real>::max());

	/** This method performs propagation, reinitialization and narrow band rebuild */
	template<Grids::StencilFormat propagateStencilFormat, Grids::StencilFormat reinitializeStencilFormat, class MyScalarField>
	void propagateDt(MyScalarField *sf, const SimulationParameters& simParameters);


	// REINITIALIZATION

	template<Grids::StencilFormat stencilFormat>
	void reinitialize(const SimulationParameters& simParameters, unsigned int iter=-1);

	
	// NARROW BAND REBUILD

	void rebuild();



	// MISC

	Real area();

	Real averageMeanCurvature();



	// Access/Get/Set Operations

	Grid *getPhi() { return phi; }
	void setPhi(Grid *phi) { this->phi = phi; }
	bool getOwnsPhi() { return ownsPhi; }
	void setOwnsPhi(bool ownsPhi) { this->ownsPhi = ownsPhi; }
	void setTDFormat(TDFormat tdFormat) { this->tdFormat = tdFormat; }
	TDFormat getTDFormat() { return tdFormat; }
	void setSDFormat(SDFormat sdFormat) { this->sdFormat = sdFormat; }
	SDFormat getSDFormat() { return sdFormat; }
	void setTime(Real t) { this->t = t; }
	Real getTime() const { return t; }
	Real getMaxAbs() const { return maxAbs; }
	Real getDx() { return dx; }
	Real getBeta() { return beta; }
	Real getGamma() { return gamma; }
	void setReinitMaxIter(unsigned int reinitMaxIter) { this->reinitMaxIter = reinitMaxIter; }
	unsigned int getReinitMaxIter() { return reinitMaxIter; }
	void setUseVanillaPeng(bool useVanillaPeng) { this->useVanillaPeng = useVanillaPeng; }
	bool getUseVanillaPeng() { return useVanillaPeng; }
	void setReinitCFLNumber(Real reinitCFLNumber) { this->reinitCFLNumber = reinitCFLNumber; }
	Real getReinitCFLNumber() { return reinitCFLNumber; }
	void setPropagateCFLNumber(Real propagateCFLNumber) { this->propagateCFLNumber = propagateCFLNumber; }
	Real getPropagateCFLNumber() { return propagateCFLNumber; }
	void setMaxAbsDeterminationMethod(MaxAbsDeterminationMethod maxAbsDeterminationMethod) { this->maxAbsDeterminationMethod = maxAbsDeterminationMethod; }
	MaxAbsDeterminationMethod getMaxAbsDeterminationMethod() { return maxAbsDeterminationMethod; }
	void setVerbose(bool verbose) { this->verbose = verbose; }
	bool getVerbose() { return verbose; }


    protected:

	void setBeta(Real beta) { this->beta = beta; }
	void setGamma(Real gamma) { this->gamma = gamma; }
	void setDx(Real dx) { this->dx = dx; this->dx2 = dx*dx; }
	inline Real computeC(Real phi) const;

	/** This method only performs advection */
	template<Grids::StencilFormat stencilFormat, class MyVelocityField, bool useSimBBox>
	void advect(MyVelocityField *vf, Real maxTime=std::numeric_limits<Real>::max());

	/** This method only performs propagation */
	template<Grids::StencilFormat stencilFormat, class MyScalarField, bool useSimBBox>
	void propagate(MyScalarField *sf, Real maxTime=std::numeric_limits<Real>::max());

	template<Grids::StencilFormat stencilFormat, bool useSgnBuffer, bool initSgnBuffer, bool useSimBBox>
	void reinitialize(unsigned int b1=0, unsigned int b2=0, Real *sgn=NULL);

	template<Grids::StencilFormat stencilFormat, class MyVelocityField, bool copyElements, bool useSimBBox>
	void advect(MyVelocityField *vf, unsigned int b1, unsigned int b2, Real maxDt, bool enforceDt=false, Real enforcedDt=0, Real *dtCFL=NULL);

	template<Grids::StencilFormat stencilFormat, class MyScalarField, bool copyElements, bool useSimBBox>
	void propagate(MyScalarField *sf, unsigned int b1, unsigned int b2, Real maxDt, bool enforceDt=false, Real enforcedDt=0, Real *dtCFL=NULL, Real *maxSV=NULL);

	template<Grids::StencilFormat stencilFormat, bool useSimBBox>
	void reinitialize(unsigned int iter=-1);

	void reserveAuxiliaryBuffer();
	void freeAuxiliaryBuffer();
	template<bool useSimBBox>
	inline void storeAuxiliaryValue(Real value, Index x, Index y, Index z, UInt i);
	inline Real getAuxiliaryValue(UInt i);

    protected:
	Grid *phi;
	Real dt;
	bool ownsPhi;
	Real maxAbs;
	Real dx, beta, gamma, dx2;
	Real reinitCFLNumber, propagateCFLNumber;
	unsigned int reinitMaxIter;
	bool useVanillaPeng;
	TDFormat tdFormat;
	SDFormat sdFormat;
	Real t;
	MaxAbsDeterminationMethod maxAbsDeterminationMethod; 
	bool verbose;
	Real *auxiliaryBuffer;
	Index simBBox[3][2];
	Index bboxSmoothTransitionWidth[2];
    };


}

#include "LevelSet3D_Impl.h"


#endif

