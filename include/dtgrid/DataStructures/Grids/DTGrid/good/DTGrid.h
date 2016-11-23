/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_dtgrid_h
#define _grids_dtgrid_h


#include <vector>      
#include <iostream>
#include <string>
#include <limits>

#include <DataStructures/Matrices/HeaderFiles/Matrix4x4.h>
#include <Algorithms/Math/Interpolators/HeaderFiles/TrilinearInterpolator.h>
#include <Algorithms/Math/HeaderFiles/BasicMath.h>
#include <Core/Exception/HeaderFiles/DefaultException.h>

#include <DataStructures/Grids/HeaderFiles/Grids.h>
#include <DataStructures/Grids/IO/HeaderFiles/SvolSaverFactory.h>
#include <DataStructures/Grids/IO/HeaderFiles/SvolLoaderFactory.h>
#include <DataStructures/Grids/IO/HeaderFiles/DTGridConstructor.h>
#include <DataStructures/Grids/IO/HeaderFiles/SvolSaver_DTTop1.0.h>
#include "CCUnion.h"


#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif


namespace Grids
{

    /**
	USE_BINARY_SEARCH :        Use binary search for random access
	USE_LINEAR_SEARCH :        Use linear search for random access. In our experience, linear search is a bit faster for typical computer graphics models. */ 
    enum DTGridSearchType { USE_BINARY_SEARCH, USE_LINEAR_SEARCH };


    template<typename Data, typename Index, typename Real, typename UInt>
    struct DTGridTraitsDefault
    {
	typedef Data DataType;
	typedef Index IndexType;
	typedef Real RealType;
	typedef Data *DataPtr;
	typedef UInt UIntType;

	/** The width of dilation used when rebuilding the narrow band (corresponds to the max number of voxels the interface can move between rebuilds) */
	static const UInt rebuildDilationWidth = 1; 
	/** If the region includes an extra band of voxels around the gamma tube (added by a narrowBandRebuild() operation, see below), some optimizations can be applied for iterators that rely on the DTGrid storing a signed distance field, and includeSafeBand should be defined. Note that when using the rebuild method supplied with the DT-Grid, the narrow band will always include a safeband! Furthermore note that when using FMM-based reinitialization, includeSafeBand must be disabled because FMM will not add a safeband. */
	static const bool includeSafeBand = true;
	/** If this is defined, the calls to increment operations in the updateStencil() method are inlined directly in the updateStencilMethod. Improves performance slightly. */
	static const bool inlineUpdateStencilCalls = true;
	/** The maximum number of grid points allowed in a dilation operation. Used to allocate static arrays. */
	static const int maxDilationWidth = 10;
	static const DTGridSearchType randomAccessType = USE_LINEAR_SEARCH;
	/** If this is false, linear array indices in off-center stencil iterators are not maintained (and hence cannot be used). Speeds up iteration */
	static const bool maintainArrayIndexInStencil = false;
	/** The width in grid cells of the zero crossing band */
	static const int zeroCrossingWidth = 1;
	/** If this boolean is true, open level sets are generally supported.
	 *  This mostly incurs a slight slowdown and is therefore disable per default.
	 *  Open level sets are always supported for iteration with a stencil iterator. */
	static const bool openLevelSetsSupport = false;
    };









    /*! 
    *
    *  \brief The DT-Grid is an infinite "Dynamic Tubular Grid" that stores values of arbitrary type.
    *
    *  \section Introduction
    *
    *  Dynamic Tubular Grid means:
    *
    *  <UL>
    *  <LI><B>Grid</B>: The DT-Grid stores arbitrary data defined on an (in principle) infinite uniform grid indexed
    *           by indices (i, j, k).
    *
    *  <LI><B>Tubular</B>: In many cases (e.g. level set simulations and fluid simulations), data is only defined in certain
    *           regions, or tubes, of the volume, and the data outside of these tubes is either undefined
    *           or constant. 
    *	    In the case of narrow band level sets, the values of grid points inside (outside of the narrow band)
    *           are clamped to the width of the narrow band on the inside of the zero-crossing (<i>insideConstant</i>), 
    *           and the values of grid points outside (outside of the narrow band) are clamped 
    *           to the width of the narrow band on the outside of the zero-crossing (<i>outSideConstant</i>).
    *
    *  <LI><B>Dynamic</B>: The DT-Grid supports the tubes to grow and/or shrink dynamically by supporting
    *           a dilation and a rebuildNarrowBand method. The dilation algorithm dilates the tubes
    *           by a certain user-specified width and is generic. The rebuildNarrowBand algorithm is tailored
    *           for level set simulations, but a similar custom algorithm can be constructed for other types
    *           of data.
    *  </UL>
    *
    *  For level set simulations, the DT-Grid is inspired by the paper on narrow band level sets:
    *  "A PDE-based fast local level set method" by Peng et al., Journal of Computational Physics 1999.
    *  Hence it builds an abstraction on top of the data structure itself by partitioning the data into
    *  a set of concentric tubes centered about the zero-crossing. The inner-most tube is the zero-crossing.
    *  The next tube is the beta-tube, then the gamma-tube and finally the entire tube (which may be identical to the gamma-tube). Hence ZERO_CROSSING_WIDTH < beta < gamma.
    *  A tube of width 'w' consists of all the grid points with absolute values less than 'w'.
    *  If the DT-Grid contains all grid points contained in a (gamma+dx)-tube, the DT-Grid is said to include a safe-band. In that case
    *  certain optimizations are possible, and the "#define includeSafeBand" should not be commented out.
    *
    *  \section Template-Arguments
    *
    *  The DT-Grid has the following template arguments:
    *  <OL>
    *  <LI> Data  : The type of data being stored. 
    *  <LI> Index : The indices used to index into the grid. Since the DT-Grid is unbounded (the grid is only limited by the integral
    *               limits of the Index data type) it is most convenient to allow Index to be negative, hence usually int or short is used.
    *               If all grid points can fit within the index-range supported by the short data-type, shorts are preferred over ints since
    *               using shorts will decrease the storage required.
    *  <LI> Real  : Should be float or double. Used for matrix calculations. Also, when the DT-Grid is used for scalar fields and level set
    *               calculations, it is assumed that Data and Real are identical.
    *  </OL>
    *  
    *
    *  \section Construction
    *
    *  Constructing a DT-Grid can be done using one of the supplied constructors or using the method push().
    *  When using the push() method, the number of elements (including compressed indices) need not be known
    *  in advance. In that case, construction proceeds as:
    \verbatim

      grid->clear();
      grid->beginSafePush();
      ...
      push(x, y, z, value);  // called for each grid point that is inserted into the grid.
      ...
      grid->endSafePush();
    \endverbatim
    *
    *  Note that the grid points must be push'ed into the grid in lexicographic (x,y,z) order! E.g. (0,1,2) must be pushed before
    *  (1,1,2). 
    *  If a set of grid points are not in the lexicographic (x,y,z) order and needs to be inserted into a DT-Grid they 
    *  must first be sorted. This can be done using three linear time stable bucket sorts. First sort in z, then in y and finally in x. 
    *
    *  \section Buffers Internal Buffers
    *
    *  DT-Grid can have up to MAX_NUM_BUFFERS (numbered 0 - (MAX_NUM_BUFFERS-1)) internal value-buffers. 
    *  These can be used during advection of e.g. a level set or velocity field.
    *  By default no storage is allocated. One \ref DTGrid_Constructor_Array "constructor" however allocates buffer 0.
    *  The internal value-buffers are manipulated with the methods: 
    *	void allocateBuffer();
    *	void freeBuffer();
    *	void swapBuffers();
    *
    *  \section Access Iterators, Sequential Access and Random Access 
    *
    *  The DT-Grid supports numerous iterators (see below) that provide access to the grid.
    *  Briefly described these iterators are:
    *  <OL>
    *  <LI> TubeIterator : This iterator sequentially visits every grid point in the narrow band (the "entire tube").
    *  <LI> VolumeIterator : This iterator sequentially visits every grid point that is either inside or contained in the narrow band.
    *  <LI> StencilTubeIterator : This iterator sequentially visits every grid point in a given tube (specified by a template argument)
    *       and also iterates a stencil (specified by a template argument) of grid points over the grid hence allowing for constant time
    *       access on average to grid points in the stencil and hence derivatives etc.
    *  <LI> ZeroCrossingIterator, BetaTubeIterator, GammaTubeIterator, EntireTubeIterator : These iterators are derived from the 
    *       StencilTubeIterator by specifying the template argument identifying the tube being iterated over.
    *  </OL>
    *
    *  If the DT-Grid is not accessed sequentially, it supports logarithmic (in the number of connected components along
    *  the compressed directions, see JSC paper for details) random access and logarithmic and constant time neighbor access methods.
    *  The neighbor search methods are based on the idea of a Locator object (see JSC paper for details).
    *  Access to grid points outside of the narrow band or tubes stored by the DT-Grid will return a clamped constant value, identified
    *  by either <i>insideConstant</i> (when the grid point is inside) or <i>outsideConstant</i> (when the grid point is outside).
    *
    *  By using the defines 'USE_BINARY_SEARCH', 'USE_OLD_BINARY_SEARCH' and 'USE_LINEAR_SEARCH', the underlying method for searching
    *  can be changed. In our experience linear search is the fastest for models typically used in computer graphics. Linear search
    *  is the default. If the number of connected components is very large, binary search may be faster. Experimentation is encouraged.
    *  Note that only one of these symbols should be defined at once (otherwise there is a built-in preference).
    *
    *
    *  \section DecoupledArrays Decoupled Arrays
    *
    *  In some situations it may be convenient, e.g. for the fast marching method, to store extra information at each grid point.
    *  One way to do this is to let the DT-Grid store a composite data type, but another, more flexible approach, is to store the
    *  extra information in separate arrays docoupled from the main data-structure. 
    *  Allocating a decoupled array can be done as: <i>MyDataType *myArray = new MyDataType[myDTGrid.getNumValues()];</i>. To access a grid 
    *  point in this array when iterating through the grid can be done as: <i>myArray[myIter.getArrayIndex()]</i>, and can be done as:
    *  <i>myArray[myLocator.iv3D]</i>, when myLocator is a Locator obtained using the getLocator() method or any of the neighbor search methods.
    *
    *
    *  \section Transformations Transformations and Bounding Box
    *
    *  The DT-Grid has an associated bounding box and transformation. <B>It is important to note that the DT-Grid itself is not
    *  limited by any boundaries. It can expand and contract dynamically and infinitely (only limited by the integral limits of
    *  the Index type specified as a template argument). As a result, level set simulations on a DT-Grid are "out-of-the-box".</B> The bounding box is useful for collision detection and CSG operations.
    *  The transformation consists of a rotation, a translation and a scale. This transformation can relate separate DT-Grids in a
    *  common world space and enable CSG operations etc between differently oriented, scaled and translated DT-Grids, e.g. level sets.
    *  Note that the grid-spacing, dx, is identical to the scale. This also means that when the scale changes, all the distance values in the
    *  level set are re-scaled.
    *
    *  \section Assumptions General Assumptions
    *
    *   <OL>
    *   <LI> For level set simulations, Data is usually float or double. In general it can be anything, e.g. a vector for fluid simulations.
    *   <LI> Index should be signed, ie. char, short or int.
    *     It is important that Index is signed as this allows the grid to expand dynamically without 
    *     limits other than those enforced by the numerical limits of the Index used.
    *   <LI> dx == dy == dz
    *   <LI> Note that the only operations that can be called on any iterator, unless it is different from the end()-iterator, are
    *     operator==() and operator!=()
    *   <LI> It is assumed that the entire defined region that is stored constitutes a closed manifold.
    *   </OL>
    *
    *************************************************************************************************/
    template<class Traits>
    class DTGrid
    {
    public:

	typedef typename Traits::DataType Data;
	typedef typename Traits::IndexType Index;
	typedef typename Traits::RealType Real;
	typedef typename Traits::DataType *DataPtr;
	typedef typename Traits::UIntType UInt;


	/*! \brief InitParams is used to initialize the DT-Grid and is passed to the constructor */
	struct InitParams
	{
	    // used only for level set simulations
	    /*! \brief dx is the spacing between grid points. The DT-Grid assumes dx==dy==dz. */
	    Real dx;
	    /*! \brief beta is the width of the beta-tube. dx < beta < gamma. Only used when the DT-Grid stores a scalar field. */
	    Real beta; 
	    /*! \brief gamma is the width of the gamma-tube. dx < beta < gamma. Only used when the DT-Grid stores a scalar field. */
	    Real gamma;
	    // constants for undefined values, used in general
	    /*! \brief The constant value of data outside of the narrow band and inside. Used in general. */
	    Data insideConstant;
	    /*! \brief The constant value of data outside of the narrow band and outside. Used in general. */
	    Data outsideConstant;  

	    InitParams(Real dx=Real(), Real beta=Real(), Real gamma=Real(), 
		Data insideConstant=Data(), Data outsideConstant=Data()) 
		: dx(dx), beta(beta), gamma(gamma), insideConstant(insideConstant), outsideConstant(outsideConstant) { };

	    virtual ~InitParams() { }
	};



	/** The maximum number of internal value buffers in the DT-Grid. These buffers can be used for advection.  */
	static const int maxNumBuffers = 2; 
	/* Number of End-Markers in index-arrays. */
	static const UInt numIndexEndMarkers = 2;
	/** Number of End-Markers in value-arrays. */
	static const UInt numValueEndMarkers = 1;
	/** Number of End-Markers in aa-arrays. */
	static const UInt numAAEndMarkers = 1;




	/*! \brief The Locator contains structural information about a grid point position. 
	 *
	 * It can be used in neighbor search. See the JSC paper 
	 *  A Locator can only point to an existing element (not an element outside the narrow band) 
	 *  Given x, y, z and ic1D, ic2D, ic3D, the indices iv1D, iv2D and iv3D are redundant. 
	 *  In particular:
	 *  <UL>
	 *        <LI>          iv1D = aa1D[ic1D>>1] + x - xIndex[ic1D]
	 *        <LI>          iv2D = aa2D[ic2D>>1] + y - yIndex[ic2D]
	 *        <LI>          iv3D = aa3D[ic3D>>1] + z - zIndex[ic3D]
	 *  </UL>
	 *  Hence iv1D, iv2D, iv3D may be optimized away, if storage is a problem.
	 */
	struct Locator
	{
	    /*! \brief Index into the va1D array of DT-Grid */
	    UInt iv1D;
	    /*! \brief Index into the va2D array of DT-Grid */
	    UInt iv2D;
	    /*! \brief Index into the va3D array of DT-Grid */
	    UInt iv3D;
	    /*! \brief Index into the xIndex array of DT-Grid */
	    UInt ic1D;
	    /*! \brief Index into the yIndex array of DT-Grid */
	    UInt ic2D;
	    /*! \brief Index into the zIndex array of DT-Grid */
	    UInt ic3D;
	    /*! \brief X-index of this grid point */
	    Index x;
	    /*! \brief Y-index of this grid point */
	    Index y;
	    /*! \brief Z-index of this grid point */
	    Index z;
	};


	/*!
	 *  \brief _Iterator sequentially visits all grid points in the data structure.
	 *
	 *  The only thing that _Iterator assumes about the data type stored (Data), is
	 *  that <i>insideConstant</i> and <i>outsideConstant</i> are specified.
	 *  <i>insideConstant</i> corresponds to the constant value in undefined connected components that
	 *  are inside. <i>outsideConstant</i> corresponds to the constant value in undefined connected components
	 *  that are outside.
	 *  _Iterator is protected and only used internally in the data-structure.
	 *  Sequentially visiting all grid points in the data structure can be done from outside
	 *  of the data structure using the public iterators. See below.
	 *
	 *  _Iterator contains a number of methods optimized for fast iteration of a stencil of _Iterator
	 *  instances. These methods are incrementFast, incrementFastZ, incrementUntil.
	 *
	 */
	template <IteratorType it>
	class _Iterator
	{
	public:
	    _Iterator(const DTGrid *parent=NULL, Data *va3D=NULL, UInt numVa3D=0, bool begin=false);

	    _Iterator(const DTGrid *parent, Data *va3D, UInt numVa3D, const Locator& loc, bool valid, bool inCorrectXYColumn);
	    
	    void commit() { }
		
	    inline bool operator==(const _Iterator& iter) const;
	    inline bool operator!=(const _Iterator& iter) const;
	    inline _Iterator &operator++();
	    inline _Iterator &operator++(int) { return operator++(); }
	    inline _Iterator &next() { return operator++(); }
	    inline bool hasNext() const;
	    inline Data operator*() const;
	    inline Data getValue() const { return **this; }
	    inline Data& getValue() { return value; }
	    // setValue sets the value directly in the data structure and in this _Iterator.
	    inline void setValue(Data value) { va3D[iv3D] = this->value = value; } 
	    // setCachedValue only sets the value directly in this _Iterator.
	    inline void setCachedValue(Data value) { this->value = value; } 
	    inline void getIndex(Index *x, Index *y, Index *z) const;
	    // Get linear "array-index" of the current grid point.
	    // This index can be used into additional value arrays decoupled
	    // from the data structure itself.
	    inline UInt getArrayIndex() const { return iv3D; }
	    inline UInt& getArrayIndex() { return iv3D; }
	    inline Index getI() const { return x; }
	    inline Index getJ() const { return y; }
	    inline Index getK() const { return z; }
	    // these methods are used to speed up increment when _Iterator corresponds to a non-center stencil
	    // grid point in Iterator. See the updateStencil method in Iterator.
	    inline void incrementFast();
	    inline void incrementFastZ(const Data& outsideValue, Index zn);
	    inline void incrementUntil(const Data& outsideValue, Index xn, Index yn, Index zn);
	    inline void incrementArrayIndexAndZ() { z++; iv3D++; }
	    inline UInt getIv1D() const { return iv1D; }
	    inline void retrieveValue() { value = va3D[iv3D]; }
	    inline bool isValid() const { return valid; }
            inline void getLocator(Locator& loc) const;

	    //protected:
	public:
	    const DTGrid *parent;
	    Data *va3D;   // cannot be defined as "Data *const va3D" since this will disallow operator=()
	    Index *zIndex;
	    Data value;
	    UInt numVa3D;
	    Index x, y, z;
	    Index zn;
	    UInt iv1D, iv2D, iv3D;
	    UInt ic1D, ic2D, ic3D;
	    Data insideConstant, outsideConstant;
	    bool inCorrectPColumn;  // true, if this _Iterator is in the correct p-column
	    bool valid;  // true, if this _Iterator is at the correct grid point (relative to the center grid point of the stencil) 

	    template<TubeType tt> friend class Iterator;
	    template<TubeType tt> friend class StencilTubeIterator;
	};


	/** 
	 *  The TubeTopologyIterator is like the TubeIterator except that it does not assume a valid 
	 *  array of values exists. Only the topology is considered */
	class TubeTopologyIterator
	{
	public:
	    TubeTopologyIterator(const DTGrid *parent=NULL, bool begin=false);

	    inline bool operator==(const TubeTopologyIterator& iter) const;
	    inline bool operator!=(const TubeTopologyIterator& iter) const;
	    inline TubeTopologyIterator &operator++();
	    inline TubeTopologyIterator &operator++(int) { return operator++(); }
	    inline bool hasNext() const;
	    inline void getIndex(Index *x, Index *y, Index *z) const;
	    /** Get linear "array-index" of the current grid point.
	     * This index can be used into additional value arrays decoupled
	     * from the data structure itself. */
	    inline UInt getArrayIndex() const { return iv3D; }
	    inline UInt& getArrayIndex() { return iv3D; }
	    inline Index getX() const { return x; }
	    inline Index getY() const { return y; }
	    inline Index getZ() const { return z; }
	    /* these methods are used to speed up increment when _Iterator corresponds to a non-center stencil
	     * grid point in Iterator. See the updateStencil method in the StencilIterator object. */
	    inline void incrementFast();
	    inline void incrementFastZ(Index zn);
	    inline void incrementUntil(Index xn, Index yn, Index zn);
	    inline void incrementArrayIndexAndZ() { z++; iv3D++; }

	    inline bool isValid() const { return valid; }
	    inline bool isInCorrectPColumn() const { return inCorrectPColumn; }

	    inline UInt getIv1D() const { return iv1D; }
	    inline UInt getVa1D() const { return parent->va1D[iv1D-1]; }

	protected:
	    const DTGrid *parent;
	    Index *zIndex;
	    UInt numVa3D;
	    Index x, y, z;
	    Index zn;
	    UInt iv1D, iv2D, iv3D;
	    UInt ic1D, ic2D, ic3D;
	    bool inCorrectPColumn;  // true, if this _Iterator is in the correct p-column
	    bool valid;  // true, if this _Iterator is at the correct grid point (relative to the center grid point of the stencil) 
	};


	/*! \brief The CSGType is used in the CSG methods to indicate the desired CSG operation */ 
	enum CSGType { CSG_UNION, CSG_INTERSECTION, CSG_DIFFERENCE };



	/*! \brief TubeIterator sequentially visits each grid point in the DT-Grid (does not include a stencil) */
	typedef _Iterator<IT_TUBE> TubeIterator;
	/*! \brief VolumeIterator sequentially visits all grid points that are either contained in the narrow band or are inside (does not include a stencil) */
	typedef _Iterator<IT_VOLUME> VolumeIterator;


	/*!
	 *
	 *
	 *  \section Introduction
	 *
	 *  \brief StencilTubeIterator sequentially iterates a stencil over all grid points in the tube indicated by template-argument <i>TubeType</i>.
	 *
	 *  StencilTubeIterator is specifically tailored for the situation where the DT-Grid stores
	 *  a scalar field. If you need to store and iterate over other types of data, e.g. velocity fields, you
	 *  can use either the TubeIterator or VolumeIterator classes defined in the DT-Grid.
	 *  You can also define a custom iterator similar to the StencilTubeIterator by wrapping one or several _Iterator
	 *  instances, as _Iterator does not make assumptions on the type of data stored.
         *
	 *  StencilTubeIterator iterates a stencil of _Iterator instances efficiently over a certain tube which makes it possible 
	 *  to gain constant time access to all elements within the stencil on average.
	 *
	 *  Template argument <i>TubeType</i> specifies which tube you wish to visit sequentially. 
	 *  In a signed distance field, a tube of width 'w' is concentric about the zero-crossing and
	 *  consists of all the elements with absolute value less than 'w'.
	 *  There are several possibilities for the 'TubeType':
	 *
	 *  <OL>
	 *   <LI> 'ZERO_CROSSING_TUBE' : Returns all elements close to the zero crossing sequentially.
	 *      (can be specified using the '#define ZERO_CROSSING_WIDTH').
	 *      All elements with absolute value less than ZERO_CROSSING_WIDTH are returned.
	 *   <LI> 'BETA_TUBE' : Returns all elements with absolute value less than 'beta' sequentially.
	 *   <LI> 'GAMMA_TUBE' : Returns all elements with absolute value less than 'gamma' sequentially.
	 *   <LI> 'ENTIRE_TUBE' : Returns all elements sequentially.
         *  </OL>  
	 *
	 *  Likewise it should be easy to add your own tube definitions if needed.
	 *  There is nothing in the data structure itself that relies on or makes use of these tubes.
	 *  Only the iterator imposes this abstraction on the data structure. The abstraction is convenient
	 *  for narrow band level set implementations, see e.g. "A PDE-based fast local level set method"
	 *  by Peng et al., Journal of Computational Physics.
	 *
	 *  To reduce the dimensionality of the number of StencilTubeIterator instances required and increase the
	 *  flexibility, some of the methods are templated as opposed to templating StencilTubeIterator itself further.
	 *  These template arguments are:
	 *
	 *  <OL>
	 *  <LI> <B>StencilFormat</B> : StencilFormat specifies which type of stencil the StencilTubeIterator supports, and
	 *     defines internally the mapping from _Iterator instances (the 'si' array of StencilTubeIterator) to
	 *     grid points of the stencil. The possible StencilFormats and the associated mapping is given below.
	 *     Note that the methods d2xy(), d2xz(), d2yz() and meanCurvature() cannot be called
	 *     unless the StencilFormat is one of the following:
	 *     SF_FIRSTORDER_CURVATURE, SF_ENO_CURVATURE, SF_WENO_CURVATURE.
	 *     The StencilFormat must be used consistently for all operations on the same StencilTubeIterator instance.
	 *     This means that once the StencilTubeIterator has been created using a specific StencilFormat, e.g.
	 *     SF_FIRSTORDER, you cannot simply call meanCurvature() with template argument SF_FIRSTORDER_CURVATURE.
	 *     For efficiency reasons, the StencilTubeIterator does not automatically check and enforce this.
	 *     Also for efficiency reasons you should always use the StencilFormat with the least number
	 *     of stencil grid points possible, as this is the fastest.
	 *     
	 *  <LI> <B>copyElements</B> : When doing level set advection with several buffers inside the DTGrid it is
	 *     sometimes convenient to let the StencilTubeIterator copy elements, that are not touched by the advection directly,
	 *     to the second buffer. If this is not desired, copyElements should be set to false.
	 *     As with StencilFormat above, 'copyElements' must also be used consistently.
	 *  </OL>
	 *
	 *  Note that unfortunately default template arguments are not allowed on template methods.
	 *
	 *  In order to instantiate a StencilTubeIterator you must call one of the following methods defined on
	 *  DTGrid itself:  beginTubeIterator(), endTubeIterator(), beginZeroCrossing(), endZeroCrossing() and so on. 
	 *
	 *  Note that the only operations allowed on an operator returned by any of the endXXX() methods,
	 *  and that the only methods that one is allowed to call on a StencilTubeIterator unless it is not equal
	 *  to the endXXX() StencilTubeIterator is 'operator==' and 'operator!='.
	 *
	 *
	 *  \section Stencil-Iterators
	 *
	 *  
	 *     In order to avoid a very large 3D array of _Iterator instances (eg. a weno stencil stored in a 3d array requires 7^3=343 stencil grid points)
	 *     we represent the stencil iterators in a 1D array. Thus, the mapping of stencil grid points to the 1D array depends
	 *     on the stencil represented by a StencilTubeIterator. The formats given here use the most convenient format wrt code reuse.
	 *     A 1D array is also the most convenient wrt cache reuse!
	 *
	 *     The first 7 entries (entries 0-6) are the same for all Stencil Formats.
	 *  
	 *
	 *   \section Stencil-Formats
	 *
	 *   The following diagrams show the mapping from stencil grid points to entries in the 1D array of _Iterator instances.
	 *   Stencil Format mappings are given by z slices in increasing z, x is horizontal (increasing to the right), y is vertical (increasing downwards).
	 *
	 *
	 *   <UL>
	 *   <LI>SF_FIRSTORDER:   ( 7 stencil grid points )
	 * 
	 \verbatim

	    - - -           - 3 -         - - -
	    - 5 -           1 0 2         - 6 -
	    - - -           - 4 -         - - -
         \endverbatim
	 *
	 *
	 *   <LI>SF_FIRSTORDER_NEIGHBORS:   ( 23 stencil grid points )
	 *  
	 *   Includes the grid points necesary for computing the upwind approximation of first order derivatives
	 *   in the center stencil grid point and its 6 immediate neighboring grid points.
	 *
	 \verbatim
			              12
	           -  8  -         13  3 14         - 20 -
	     7     9  5 10      15  1  0  2 16     21  6 22      24
	           - 11  -         17  4 18         - 23 -
				      19
         \endverbatim
	 *
	 *
	 *   <LI>SF_FIRSTORDER_CURVATURE: ( 19 stencil grid points )
	 \verbatim

	    -  7  -       11 3 12         - 15  -
	    8  5  9        1 0  2        16  6 17 
	    - 10  -       13 4 14         - 18  - 
	 \endverbatim
	 *  
	 *
	 *  <LI>SF_WENO: (used by both ENO and WENO schemes) (19 stencil grid points)
	 *
	 \verbatim

	   0-6:   (0,0,0), (-1,0,0), (1,0,0), (0,-1,0), (0,1,0), (0,0,-1), (0,0,1),  // corresponds to SF_FIRSTORDER
	   7-10:  (-3,0,0), (-2,0,0), (2,0,0), (3,0,0)                               // additional grid points in x direction
	   11-14: (0,-3,0), (0,-2,0), (0,2,0), (0,3,0)                               // additional grid points in y direction
	   15-18: (0,0,-3), (0,0,-2), (0,0,2), (0,0,3)                               // additional grid points in z direction
	 \endverbatim  
	 *
	 *  <LI>SF_WENO_CURVATURE: (used by both ENO and WENO schemes) (31 stencil grid points)
	 *
	 \verbatim

	   0-18: Same as SF_WENO
	   19-30: Same as 7-18 in SF_FIRSTORDER_CURVATURE

	                                    11 
	                                    12
	                 - 19  -         23  3 24         - 27  -
	    15    16    20  5 21    7  8  1  0  2  9 10   28 6 29    17    18
	                 - 22  -         25  4 26         - 30  - 
	                                    13
	                                    14
	 
	 \endverbatim
	 *  <LI>SF_VOXEL (8 stencil grid points). Can be used for marching cubes.
	 \verbatim

	                0  1      4  5
	                2  3      6  7
	 \endverbatim
	 *  <LI>SF_VOXEL_GRAD (32 stencil grid points). Can be used for marching cubes and supports computation of 
	 *  gradients at the corner points.
	 \verbatim
	                      12  13              20 21               
                8   9     14   0   1   15     22   4  5  23     28 29          
                10 11     16   2   3   17     24   6  7  25     30 31
                              18  19              26 27 
	 \endverbatim
	 *  <LI>SF_BOX (27 stencil grid points). 
	 \verbatim

	   19  7 20       11 3 12        23 15 24
	    8  5  9        1 0  2        16  6 17 
	   21 10 22       13 4 14        25 18 26 

	 \endverbatim
	 *
	 *
	 *  <LI>SF_NONE (1 stencil grid point). No stencil.
	 *  </UL>
	 *
	 *
	 *  \section MC Marching Cubes
	 *
	 *  The marching cubes algorithm can be implemented using a StencilTubeIterator with stencil-format 'SF_VOXEL'
	 *  or 'SF_VOXEL_GRAD' and TubeType 'BETA_TUBE'.
	 *  'SF_VOXEL_GRAD' includes the possibilty for computing the gradient at the corner grid points of the voxel
	 *  using the method <i>gradientAtVoxel()</i>.
	 *  It sequentially returns all grid points in the zero-crossing and provides
	 *  constant time access to the seven neighbors that together with the current grid point constitute the eight
	 *  corners of a voxel centered between grid points.
	 *  Please note that with a 'ZERO_CROSSING_WIDTH' equal to dx, the narrow band should be at least 3 voxels on
	 *  each side of the zero-crossing when using 'SF_VOXEL', and at least 4 voxels on each side of the zero-crossing
	 *  when using 'SF_VOXEL_GRAD'.
	 *  
	 *
	 *  Note that the stencilFormat must be used when instantiating and incrementing an iterator when 
	 *  iterating over a tube with a stencil. Iterating over the entire narrow band with the SF_FIRSTORDER stencilFormat, for example,
	 *  is done as follows:
	 *  EntireStencilTubeIterator iter = grid.beginStencilTubeIterator<ENTIRE_TUBE, SF_FIRSTORDER, false>();
	 *  EntireStencilTubeIterator iend = grid.endStencilTubeIterator<ENTIRE_TUBE, SF_FIRSTORDER, false>();
	 *  while (iter != iend)
	 *  {
	 *	// ...do the processing here...
	 *
	 *	iter.operator++<SF_FIRSTORDER, false>();
	 *  } 
	 *
	 *  Note that if just using the method call iter++ implicitly corresponds to the method call
	 *  iter.operator++<SF_NONE, false>(), i.e. a method call with a stencil consisting only of the center grid point.
	 */
	template<TubeType tt> 
	class StencilTubeIterator
	{
	public:

	    StencilTubeIterator() { }
	    virtual ~StencilTubeIterator();

	    inline void commit() { }

	    /*! \brief prefix operator++ method */
	    template <StencilFormat iterStencil, bool copyElements> inline
		StencilTubeIterator& operator++();

	    //inline StencilTubeIterator& operator++() { return operator++<SF_NONE, false>(); }
	    /*! \brief postfix operator++ method. Note that it has default template arguments StencilFormat=SF_NONE and CopyElements=false
	     *  contrary to the prefix operator. */
	    inline StencilTubeIterator& operator++(int) { return operator++<SF_NONE, false>(); }

	    inline bool operator==(const StencilTubeIterator& iter) const;
	    inline bool operator!=(const StencilTubeIterator& iter) const;
	    /*! \brief The value of the center element of the stencil. Identical to the getValue() method */
	    inline Real operator*() const;
	    /*! \brief Get the index (i,j,k) of the center element of the stencil. */
	    inline void getIndex(Index* i, Index *j, Index *k) const;
	    /*! \brief Get index i of the center element of the stencil. */
	    inline Index getI() const { return i; }
	    /*! \brief Get index j of the center element of the stencil. */
	    inline Index getJ() const { return j; }
	    /*! \brief Get index k of the center element of the stencil. */
	    inline Index getK() const { return k; }

	    /*! \brief Get index i of element m of the stencil. */
	    inline Index getI(UInt m) { return si[m].getI(); }
	    /*! \brief Get index j of element m of the stencil. */
	    inline Index getJ(UInt m) { return si[m].getJ(); }
	    /*! \brief Get index k of element m of the stencil. */
	    inline Index getK(UInt m) { return si[m].getK(); }

	    /*! \brief Get linear <i>array-index</i> of the current (stencil-center) grid point.
	     * This index can be used to index into additional value arrays decoupled
	     * from the data structure itself. */
	    inline UInt getArrayIndex() const;
	    inline UInt getLinearIndex() const { return getArrayIndex(); }
	    /*! \brief Get linear <i>array-index</i> of stencil grid point number 'i' 
	     * (See the definition of the Stencil Formats).
	     * This index can be used to index into additional value arrays decoupled
	     * from the data structure itself. */
	    inline UInt getArrayIndex(int i) const;


	    /*! \brief first order one sided difference to the left in the x direction */
	    inline Real dm1x() const;
	    /*! \brief first order one sided difference to the right in the x direction */
	    inline Real dp1x() const;
	    /*! \brief first order one sided difference to the left in the y direction */
	    inline Real dm1y() const;
	    /*! \brief first order one sided difference to the right in the y direction */
	    inline Real dp1y() const;
	    /*! \brief first order one sided difference to the left in the z direction */
	    inline Real dm1z() const;
	    /*! \brief first order one sided difference to the right in the z direction */
	    inline Real dp1z() const;

	    /*! \brief second order central difference in the x direction */
	    inline Real dc2x() const;
	    /*! \brief second order central difference in the y direction */
	    inline Real dc2y() const;
	    /*! \brief second order central difference in the z direction */
	    inline Real dc2z() const;

	    /*! \brief second order accurate second partial derivative in x direction */
	    inline Real d2xx() const;
	    /*! \brief second order accurate second partial derivative in y direction */
	    inline Real d2yy() const;
	    /*! \brief second order accurate second partial derivative in z direction */
	    inline Real d2zz() const;

	    /*! \brief second order accurate second partial derivative df^2/(dxdy) */
	    template <StencilFormat iterStencil> inline
		Real d2xy() const;
	    /*! \brief second order accurate second partial derivative df^2/(dxdz)*/
	    template <StencilFormat iterStencil> inline
		Real d2xz() const;
	    /*! \brief second order accurate second partial derivative df^2/(dydz)*/
	    template <StencilFormat iterStencil> inline
		Real d2yz() const;
            
	    /*! \brief fifth order WENO one sided finite difference to the left in the x direction */
	    inline Real dm5x() const;
	    /*! \brief fifth order WENO one sided finite difference to the right in the x direction */
	    inline Real dp5x() const;
	    /*! \brief fifth order WENO one sided finite difference to the left in the y direction */
	    inline Real dm5y() const;
	    /*! \brief fifth order WENO one sided finite difference to the right in the y direction */
	    inline Real dp5y() const;
	    /*! \brief fifth order WENO one sided finite difference to the left in the z direction */
	    inline Real dm5z() const;
	    /*! \brief fifth order WENO one sided finite difference to the right in the z direction */
	    inline Real dp5z() const;


	    template <StencilFormat iterStencil> inline
		Real meanCurvature() const;

	    inline Real laplacian() const;
	    inline void normal(Real n[3]) const;
	    inline void gradient(Real g[3]) const;
	    inline Real gradientLength() const;

	    /*! \brief Gets the gradient at grid point number <i>i</i> in the stencil. */
	    template <StencilFormat iterStencil, int i> inline
		void gradientAtVoxel(Real g[3]) const;


	    /*! \brief Sets the value 'v' directly in buffer 'i' of the DT-Grid at the position 
	     * determined by the array index of iterator number 'j' in the stencil.
	     * Default value of 'i' is the buffer being iterated over. */
	    inline void setValue(Real v, const int i=-1, const int j=0);
	    //inline void setValue(Real v, int i=0, int j=0);


	    /*! \brief Get value at corresponding position in buffer 'i' (can determined from the array-index also returned by the getArrayIndex() method.
	     * Default value of 'i' is the buffer being iterated over. */
	    inline Real getValue(const int i=-1, const int j=0);
	    //inline Real getValue(int i=0, int j=0);

	    /*! \brief Gets the value of TubeIterator instance number <i>i</i>. */
	    inline Real operator()(UInt i) const;


	    inline void getIndex(UInt i, Index *x, Index *y, Index *z) const;

	    IterationOrder getIterationOrder() { return IteratorOrder_XYZ; }

	    inline UInt getIv1D(UInt m) const { return si[m].getIv1D(); }

	    template<StencilFormat iterStencil> void retrieveStencilValues();

	    /*! returns true if iterator 'm' is in an existing grid point inside the narrow band, false otherwise */
	    inline bool isValid(int m) { return si[m].isValid(); }

	    inline void getLocator(Locator& loc) const { si[0].getLocator(loc); }

	protected:

	    template <StencilFormat iterStencil, bool copyElements>
	    	void init(DTGrid<Traits> *parent, UInt bufferId1, UInt bufferId2, UInt begin, bool initSecondBuffer);

	    template <StencilFormat iterStencil, bool copyElements>
	    	void init(DTGrid<Traits> *parent, UInt bufferId1, UInt bufferId2, UInt begin, bool initSecondBuffer, Locator *loc);

	    template<StencilFormat iterStencil>
		void updateStencil();

	    template<StencilFormat iterStencil>
		void initStencilUsingRandomAccess(UInt bufferId);

	    inline void gradient(Real g[3], UInt x1, UInt x2, UInt y1, UInt y2, UInt z1, UInt z2) const;

	    // helper method for the WENO one-sided difference methods
	    inline Real computeWENO(Real v1, Real v2, Real v3, Real v4, Real v5) const;

	protected:
	    DTGrid *parent;
	    Real *v1, *v2;
	    UInt numValues1, numValues2; // number of values in the values array of the parent
	    Index i, j, k;			 // index (i,j,k) in grid
	    TubeIterator si[StencilTraits::maxNumberOfStencilIterators];         // stencil of iterators	    
	    Real beta, gamma, dx;                // width of beta and gamma tubes and dx
	    Real dxc1, dxc2, dxc3, dxc4;         // constants used during finite differencing
	    bool fastUpdate;                     // used for optimization purposes in the updateStencil()
	                                         // method if the stencil has only moved one grid point in
	                                         // the direction of the compression (Z)

	    friend class DTGrid;
	};



	// declaration of custom iterators

	/*! \brief ZeroCrossingIterator sequentially iterates a stencil over the Zero-Crossing-Tube */
	typedef StencilTubeIterator<ZERO_CROSSING_TUBE> ZeroCrossingIterator;
	/*! \brief BetaTubeIterator sequentially iterates a stencil over the Beta-Tube */
	typedef StencilTubeIterator<BETA_TUBE> BetaTubeIterator;
	/*! \brief GammaTubeIterator sequentially iterates a stencil over the Gamma-Tube */
	typedef StencilTubeIterator<GAMMA_TUBE> GammaTubeIterator;
	/*! \brief EntireTubeIterator sequentially iterates a stencil over the entire tube (i.e. the entire narrow band) */
	typedef StencilTubeIterator<ENTIRE_TUBE> EntireTubeIterator;
	/*! \brief InsideTubeIterator sequentially iterates a stencil over the part of the narrow band that is inside */
	typedef StencilTubeIterator<INSIDE_TUBE> InsideTubeIterator;
	/*! \brief OutsideTubeIterator sequentially iterates a stencil over the part of the narrow band that is outside */
	typedef StencilTubeIterator<OUTSIDE_TUBE> OutsideTubeIterator;



	// construction and destruction

	/*! \brief Construct an empty DT-Grid with identity transformation. 
	 *  This constructor does not make assumtions on the data type stored. */
	DTGrid(InitParams initParams=InitParams());

	template<typename InitIter> DTGrid(InitIter beginIter, InitIter endIter, InitParams initParams);

	/*! \brief Construct an empty DT-Grid with bounding box [0,dim[0]-1], [0,dim[1]-1], [0,dim[2]-1] and
	 * rotation (XYZ Euler angles) (rx, ry, rz) and translation (tx, ty, tz). 
	 * This constructor does not make assumptions on the data type stored. */
	DTGrid(Index dim[3], Real rx, Real ry, Real rz, Real tx, Real ty, Real tz, InitParams initParams);

	/*! Shallow copy */
	DTGrid(const DTGrid& dtgrid2);

	virtual ~DTGrid();

	
	void init(InitParams initParams=InitParams());

	template<StencilFormat stencilFormat> Index getStencilXWidth();

	void getInitParams(InitParams *initParams) const;

	void saveSparseVolume(const std::string& fileName, std::string outputFormat="");
	void loadSparseVolume(const std::string& fileName);



	/*! \brief Push grid point (x,y,z) with value <i>val</i> into this DT-Grid.
	 *  Note that grid points must be pushed in lexicographic (x,y,z) order. 
	 *  Also note that a sequence of push() operations must be enclosed by calls
	 *  to beginSafePush() and endSafePush(). */
	void push(Index x, Index y, Index z, Data val);
	/*! \brief Begin a sequence of push operations. */
	void beginSafePush();
	/*! \brief End as sequence of push operations. */
	void endSafePush();
	template <bool safePush, bool safePushValue> void push3D(Index x, Index y, Index z, Data val);
	template <bool safePush> void push3D(Index x, Index y, Index z);
	template <bool safePush> void push2D(Index x, Index y, UInt val);
	template <bool safePush> void push1D(Index x, UInt val);
	void pushAndMend3D(Index x, Index y, Index z, Data val);



	// iterator related

	/*! \brief Return a StencilTubeIterator that points to the first element in the zero-crossing. */
	template <StencilFormat iterStencil, bool copyElements> ZeroCrossingIterator beginZeroCrossing();
	/*! \brief Return a StencilTubeIterator that points to the last element in the zero-crossing. */
	template <StencilFormat iterStencil, bool copyElements> ZeroCrossingIterator endZeroCrossing();
	/*! \brief Return a StencilTubeIterator that points to the first element in the beta-tube. */
	template <StencilFormat iterStencil, bool copyElements> BetaTubeIterator beginBetaTube();
	/*! \brief Return a StencilTubeIterator that points to the last element in the beta-tube. */
	template <StencilFormat iterStencil, bool copyElements> BetaTubeIterator endBetaTube();
	/*! \brief Return a StencilTubeIterator that points to the first element in the gamma-tube. */
	template <StencilFormat iterStencil, bool copyElements> GammaTubeIterator beginGammaTube();
	/*! \brief Return a StencilTubeIterator that points to the last element in the gamma-tube. */
	template <StencilFormat iterStencil, bool copyElements> GammaTubeIterator endGammaTube();
	/*! \brief Return a StencilTubeIterator that points to the first element in the narrow band. */
	template <StencilFormat iterStencil, bool copyElements> EntireTubeIterator beginEntireTube();
	/*! \brief Return a StencilTubeIterator that points to the last element in the narrow band. */
	template <StencilFormat iterStencil, bool copyElements> EntireTubeIterator endEntireTube();

	/*! \brief Return a StencilTubeIterator that points to the first element in the narrow band. */
	template <StencilFormat iterStencil, bool copyElements> InsideTubeIterator beginInsideTube();
	/*! \brief Return a StencilTubeIterator that points to the last element in the narrow band. */
	template <StencilFormat iterStencil, bool copyElements> InsideTubeIterator endInsideTube();

	/*! \brief Return a StencilTubeIterator that points to the first element in the narrow band. */
	template <StencilFormat iterStencil, bool copyElements> OutsideTubeIterator beginOutsideTube();
	/*! \brief Return a StencilTubeIterator that points to the last element in the narrow band. */
	template <StencilFormat iterStencil, bool copyElements> OutsideTubeIterator endOutsideTube();


	/*! \brief Return a TubeIterator that points to the element in the narrow band (no stencil) pointed to by the locator. */
	TubeIterator beginTubeIterator(const Locator& loc, bool inCorrectXYColumn, bool valid, UInt id=0) const;
	/*! \brief Return a TubeIterator that points to the first element in the narrow band (no stencil). */
	TubeIterator beginTubeIterator(UInt id=0) const;
	/*! \brief Return a TubeIterator that points to the last element in the narrow band (no stencil). */
	TubeIterator endTubeIterator(UInt id=0) const;

	/*! \brief Return a TubeIterator that points to the first element in the narrow band (no stencil). */
	TubeTopologyIterator beginTubeTopologyIterator() const;
	/*! \brief Return a TubeIterator that points to the last element in the narrow band (no stencil). */
	TubeTopologyIterator endTubeTopologyIterator() const;

	/*! \brief Return a VolumeIterator that points to the first element in narrow band unioned with the inside. */
	VolumeIterator beginVolumeIterator(UInt id=0) const;
	/*! \brief Return a VolumeIterator that points to the last element in the narrow band unioned with the inside. */
	VolumeIterator endVolumeIterator(UInt id=0) const;


	/*! \brief Generic iterator-method that returns a StencilTubeIterator that points to the first element in the given tube.
	 * Semantics:
	 * <UL>
	 *  <LI> If both buffers i1 and i2 are specified, those buffers are used (buffer i2 is allocated if not equal to i1).
	 *  <LI> If no buffers are specified, buffers 0 and 1 are used, and buffer 1 allocated. (default parameters).
	 *  <LI> If only buffer i1 is specified, buffer i2 is set to buffer 1 and is assumed allocated. 
	 * </UL>
	 */
	template<TubeType tt, StencilFormat iterStencil, bool copyElements>
	    StencilTubeIterator<tt> beginStencilTubeIterator(int i1=0, int i2=-1);
	/*! \brief Generic iterator-method that returns a StencilTubeIterator that points to the last element in the given tube. */
	template<TubeType tt, StencilFormat iterStencil, bool copyElements>
	    StencilTubeIterator<tt> endStencilTubeIterator(int i1=0);


	/*! \brief Generic iterator-method that returns a StencilTubeIterator that points to the first element in the given tube.
	 * Semantics:
	 * <UL>
	 *  <LI> If both buffers i1 and i2 are specified, those buffers are used (buffer i2 is allocated if not equal to i1).
	 *  <LI> If no buffers are specified, buffers 0 and 1 are used, and buffer 1 allocated. (default parameters).
	 *  <LI> If only buffer i1 is specified, buffer i2 is set to buffer 1 and is assumed allocated. 
	 * </UL>
	 *
	 * 'numIntervals' specifies the number of intervals for which we want individual iterators.
	 * The iterators are output in the 'iterators' array. The end-iterator of iterator i is iterator i+1.
	 * The lenght of the iterators array must be 'numIntervals+1'. Note that this array is allocated by the caller of the method.
	 */
	template<TubeType tt, StencilFormat iterStencil, bool copyElements>
	    void getStencilTubeIterators(UInt numIntervals, StencilTubeIterator<tt> *iterators, int i1=0, int i2=-1);


	// random access operators	

	/*! \brief Random access. Returns the value of grid point (x,y,z). */
	inline Data operator()(Index x, Index y, Index z) const;
	/*! \brief Random access. Returns true if grid point (x,y,z) exists in the grid.
	 * If grid point (x,y,z) exists, val is set to the value at grid point (x,y,z). */
	inline bool operator()(Index x, Index y, Index z, Data *val) const;
	/*! \brief Random access operator that determines if a given grid point is inside or outside.
	 *  inside() is specifically tailored for level sets.
	 * \return Returns true if grid point (x, y, z) is less than zero. */
	inline bool inside(Index i, Index j, Index k) const;
	/*! \brief Random access. Returns the value at position 'pos', taking into account the transformation of the grid. */
	inline Data operator()(const Matrix::Vector3<Real>& pos);

	//inline Data& operator()(UInt i) { return va3D[i]; }
	inline Data operator()(UInt i) { return va3D[i]; }
	inline void setValue(UInt i, Data val) { va3D[i] = val; }

       /*!
	* \brief closestPoint computes the closest point on the surface from a given grid coordinate.
	*
	*  closestPoint() assumes that DT-Grid stores a signed distance field.
	*
	*  Input
	*  \arg p : indices of grid-point at which closest point should be evaluated.
	*
	*  Output:
	*  \arg n : normal (points away from closest point), valid if *inNarrowBand == true
	*  \arg cd : closest distance to interface, clamped to the width of the narrow band if *inNarrowBand == false
	*  \arg inside : true if p is inside or on interface, false otherwise
	*  \arg inNarrowBand : true if p lies in the narrow band
	*
	*  In order to find the closest point, the caller of this method should do the following:
	*  
	*  cp[0] = p[0] - n[0] * (cd / dx); 
	*  cp[1] = p[1] - n[1] * (cd / dx); 
	*  cp[2] = p[2] - n[2] * (cd / dx); 
	*
	*  The narrow band should be close to a signed distance function, |n| = 1 .
	*/
	inline void closestPoint(const Index p[3], Matrix::Vector3<Real>& n, Real *cd, bool *inside, bool *inNarrowBand) const;
	
    inline void closestPoint(const Index x, const Index y, const Index z, Matrix::Vector3<Real>& n, Real *cd, bool *inside, bool *inNarrowBand) const;
	
    bool lookupDistanceAndNormal(const Real pos[3], Real& dist, Matrix::Vector3<Real>& normal) const;
    // neighbor search methods 
	// (alternate versions of these methods can easily be constructed depending on which output is required)

	/*! \brief get locator of neighbor (x-1,y,z)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x-1,y,z)
	 * \return true if neighbor Locator is valid, false otherwise */
	inline bool getXM(const Locator& loc, Locator& nbLoc) const;
	/*! \brief get locator of neighbor (x+1,y,z)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x+1,y,z)
	 * \return true if neighbor Locator is valid, false otherwise */
	inline bool getXP(const Locator& loc, Locator& nbLoc) const;
	/*! \brief get locator of neighbor (x, y-1,z)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x,y-1,z)
	 * \return true if neighbor Locator is valid, false otherwise */
	inline bool getYM(const Locator& loc, Locator& nbLoc) const;
	/*! \brief get locator of neighbor (x, y+1,z)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x,y+1,z)
	 * \return true if neighbor Locator is valid, false otherwise */
	inline bool getYP(const Locator& loc, Locator& nbLoc) const;
	/*! \brief get locator of neighbor (x, y, z-1)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x,y,z-1)
	 * \return true if neighbor Locator is valid, false otherwise */
	inline bool getZM(const Locator& loc, Locator& nbLoc) const;
	/*! \brief get locator of neighbor (x, y, z+1)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x,y,z+1)
	 * \return true if neighbor Locator is valid, false otherwise */
	inline bool getZP(const Locator& loc, Locator& nbLoc) const;

	/*! \brief get locator of neighbor (x-1,y,z), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x-1,y,z)
	 * Assumes that the neighbor exists. */
	inline void getEXM(const Locator& loc, Locator& nbLoc) const;
	/*! \brief get locator of neighbor (x+1,y,z), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x+1,y,z)
	 * Assumes that the neighbor exists. */
	inline void getEXP(const Locator& loc, Locator& nbLoc) const;
	/*! \brief get locator of neighbor (x, y-1,z), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x,y-1,z)
	 * Assumes that the neighbor exists. */
	inline void getEYM(const Locator& loc, Locator& nbLoc) const;
	/*! \brief get locator of neighbor (x, y+1,z), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x,y+1,z)
	 * Assumes that the neighbor exists. */
	inline void getEYP(const Locator& loc, Locator& nbLoc) const;
	/*! \brief get locator of neighbor (x, y, z-1), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x,y,z-1)
	 * Assumes that the neighbor exists. */
	inline void getEZM(const Locator& loc, Locator& nbLoc) const;
	/*! \brief get locator of neighbor (x, y, z+1), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg nbLoc : Output Locator of (x,y,z+1)
	 * Assumes that the neighbor exists. */
	inline void getEZP(const Locator& loc, Locator& nbLoc) const;


	/*! \brief get locator of neighbor (x-1,y,z)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline UInt getXMArrayIndex(const Locator& loc, UInt *arrayIndex) const;
	/*! \brief get locator of neighbor (x+1,y,z)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline UInt getXPArrayIndex(const Locator& loc, UInt *arrayIndex) const;
	/*! \brief get locator of neighbor (x, y-1,z)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline UInt getYMArrayIndex(const Locator& loc, UInt *arrayIndex) const;
	/*! \brief get locator of neighbor (x, y+1,z)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline UInt getYPArrayIndex(const Locator& loc, UInt *arrayIndex) const;
	/*! \brief get locator of neighbor (x, y, z-1)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline UInt getZMArrayIndex(const Locator& loc, UInt *arrayIndex) const;
	/*! \brief get locator of neighbor (x, y, z+1)
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline UInt getZPArrayIndex(const Locator& loc, UInt *arrayIndex) const;



	/*! \brief get locator of neighbor (x-1,y,z), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline void getEXMArrayIndex(const Locator& loc, UInt *arrayIndex) const;
	/*! \brief get locator of neighbor (x+1,y,z), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline void getEXPArrayIndex(const Locator& loc, UInt *arrayIndex) const;
	/*! \brief get locator of neighbor (x, y-1,z), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline void getEYMArrayIndex(const Locator& loc, UInt *arrayIndex) const;
	/*! \brief get locator of neighbor (x, y+1,z), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline void getEYPArrayIndex(const Locator& loc, UInt *arrayIndex) const;
	/*! \brief get locator of neighbor (x, y, z-1), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline void getEZMArrayIndex(const Locator& loc, UInt *arrayIndex) const;
	/*! \brief get locator of neighbor (x, y, z+1), which we assume exists
	 * \arg loc   : Input Locator of (x,y,z)
	 * \arg arrayIndex : array index of grid point (x,y,z), valid if method returns 1
	 * \return 1 if neighbor exists, 0 otherwise */
	inline void getEZPArrayIndex(const Locator& loc, UInt *arrayIndex) const;



	/*! getLocator returns the Locator of position (x,y,z) if it exists or the Locator of the lexicographically smallest element
	*  that is lexicographically larger than (x,y,z) */
	inline bool getExistingLocator(Index x, Index y, Index z, Locator *loc, bool *inCorrectXYColumn) const;

	/*! \brief getLocator computes the Locator of position (x,y,z)
	 * \return true if Locator is valid, false otherwise */
	inline bool getLocator(Index x, Index y, Index z, Locator *loc) const;


    inline bool getLocatorAndValue(Index x, Index y, Index z, Locator& loc, Data& value) const;

	/*! \brief getLocator computes the Locator of the grid point at position 'i' in the 3D value array.
	 * i must be less than numVa3D */
	inline void getLocator(UInt i, Locator *loc) const;

	/*! \brief Returns the value at the grid point specified by Locator, <i>loc</i>, in internal buffer, <i>i</i> (default is 0). */
	inline Data getValue(const Locator& loc, UInt id = 0) const { return values[id][loc.iv3D]; }
	//inline Data getValue(const Locator& loc, UInt id = 0) const { return values[id][loc.iv3D]; }

	/*! 
	 * \brief Get values at all 8 corners of voxel with lower left corner (x,y,z). 
	 *  The values in the <i>values</i> array are the values of the grid points (x,y,z),(x,y,z+1),(x,y+1,z),(x,y+1,z+1),(x+1,y,z),(x+1,y,z+1),(x+1,y+1,z),(x+1,y+1,z+1)  */
	inline void getVoxelValues(Index x, Index y, Index z, Data values[8]) const;



	// CSG


	/*! \brief This CSG method computes the CSG operation indicated
	 *  by CSGType between this DT-Grid and dtgrid2 and assigns the result to this DT-Grid.
	 *  It takes into account the transformation from this DTGrid to dtgrid2.
	 *
	 * TemplateInterpolator requirements:
	 * <OL>
	 *  <LI> Data interp(Data values[8], Vector3<Data> unitOffset)
	 * </OL>
	 *
	 * Note that the CSG operation does not change the bounding box of this DTGrid. */
	template<class TemplateInterpolator> void CSG(const DTGrid& dtgrid2, CSGType csgType, TemplateInterpolator *interpolator);


	/*! \brief This CSG method computes the CSG operation indicated
	 *  by CSGType between this DT-Grid and dtgrid2 and assigns the result to this DT-Grid.
	 *  It assumes that this DT-Grid and dtgrid2 are Grid Aligned (GA), i.e. defined in the same grid coordinate system. */
	void CSG_GA(DTGrid& dtgrid2, CSGType csgType);


	/*! \brief This CSG method computes the CSG operation indicated
	 *  by CSGType between this DT-Grid and dtgrid2 and assigns the result to this DT-Grid.
	 *  It assumes that this DT-Grid and dtgrid2 are Grid Aligned (GA), i.e. defined in the same grid coordinate system. 
	 *  func1 is a function applied to this dtgrid before doing the CSG, and func2 is a function applied to dtgrid 2.
	 *  So essentially this method is a fast way to do the following CSG operation: CSG( (*func1)(dtgrid1), (*func2)(dtgrid2) ).
	 *  func1 and func2 are function objects and must support the following method: Real operator()(Real phiValue).
	 */
	template<typename FuncT1, typename FuncT2>
	void CSG_GA(DTGrid& dtgrid2, CSGType csgType, FuncT1 *func1, FuncT2 *func2);


	/*! \brief This CSG method computes the CSG operation indicated
	 *  by CSGType between this DT-Grid and the grid point values pointed to by the TemplateIterators and assigns the result to this DT-Grid.
	 *  It assumes that the two grids are Grid Aligned (GA), i.e. defined in the same grid coordinate system.
	 *
	 * TemplateIterator requirements:
	 * <OL>
	 *  <LI> bool operator++() 
	 *  <LI> bool operator++(int)
	 *  <LI> bool operator!=()
	 *  <LI> Data getValue() 
	 *  <LI> getIndex(Index *x, Index *y, Index *z)
	 * </OL>
	 *
	 * Note that the CSG operation does not change the bounding box of this DTGrid. */
	template<class TemplateIterator> 
	void CSG_GA(TemplateIterator iter2, TemplateIterator iend2, CSGType csgType);


	/*! \brief This CSG method computes the CSG operation indicated
	 *  by CSGType between this DT-Grid and the grid point values pointed to by the TemplateIterators and assigns the result to this DT-Grid.
	 *  It assumes that the two grids are Grid Aligned (GA), i.e. defined in the same grid coordinate system.
	 *
	 * TemplateIterator requirements:
	 * <OL>
	 *  <LI> bool operator++() 
	 *  <LI> bool operator++(int)
	 *  <LI> bool operator!=()
	 *  <LI> Data getValue() 
	 *  <LI> getIndex(Index *x, Index *y, Index *z)
	 * </OL>
	 *
	 * Note that the CSG operation does not change the bounding box of this DTGrid. */
	template<class TemplateIterator, typename FuncT1, typename FuncT2> 
	void CSG_GA(TemplateIterator iter2, TemplateIterator iend2, CSGType csgType, FuncT1 *func1, FuncT2 *func2);


	// transformations

	Matrix::Matrix4x4<Real> getWorldToGrid() const { return worldToGrid; }
	Matrix::Matrix4x4<Real> getGridToWorld() const { return gridToWorld; }
	/*! \brief Changes the scale (and dx) and re-scales all distance values in the grid */ 
	void setScale(Real scale);
	void setTranslation(Real tx, Real ty, Real tz);
	/*! Angles in radians */
	void setRotationXYZ(Real rx, Real ry, Real rz);
	Real getScale();
	template<typename Real2> void getTranslation(Real2 *tx, Real2 *ty, Real2 *tz);
	/*! Angles in radians */
	template<typename Real2> void getRotationXYZ(Real2 *rx, Real2 *ry, Real2 *rz);

	void translateGrid(Index t[3]);


	// misc

	DTGrid<Traits> *copy(Index bbox[3][2]);
	void paste(DTGrid<Traits> *subset, Index bbox[3][2]);


	/*! \brief <i>dim</i> is set to the effective bounding box of this DT-Grid. 
	 * The <i>effective grid size</i> is the actual extents of the grid points contained in the DT-Grid. */
	void effectiveGridSize(Index dim[3][2]) const;
	/*! \brief <i>bbox</i> is set to the bounding box of this DT-Grid. The bounding box may be different
	 *  from the actual extents, e.g. when storing open level sets. However, the bounding box will always
	 *  be at least as big as the effective grid size */
	void boundingBox(Index bbox[3][2]) const;
	/*! \brief <i>bbox</i> is set to the bounding box of this DT-Grid. The bounding box may be different
	 *  from the actual extents, e.g. when storing open level sets. However, the bounding box will always
	 *  be at least as big as the effective grid size */
	template<typename Index2>
	void boundingBox(Index2 bbox[3][2]) const;
	/*! \brief Gets the dimensions of the bounding box along each dimension, uses boundingBox() as a sub-routine */
	void boundingBoxDim(Index bboxDim[3]);
	void setOpenLevelSetBBox(Index openLevelSetBBox[3][2]);
	/*! \brief Return the number of values (or grid points) stored in the DT-Grid. */ 
	UInt getNumValues(UInt id=0) const { return numValues[id]; }
	/*! \brief Assumes that the DTGrid stores a level set. 
	 *  For open level sets, the assumption of the rebuildNarrowBandMethod is that the narrow band was fully contained 
	 *  within the openLevelSetBoundingBox prior to the rebuild. */
	virtual void rebuildNarrowBand();
	/*! \brief Assumes that the DTGrid stores a level set */
	virtual void rebuildNarrowBand(UInt **permutationArray);
	/*! \brief Dilate this DT-Grid by <i>width</i> grid points */
	void dilate(Index width);
	/*! \brief Dilate this DT-Grid, <i>dtgrid2</i>, by <i>width</i> grid points and assign it to this DT-Grid. */
	void dilate(Index width, DTGrid *dtgrid2, bool doDelete=true);
	/*! \brief Allocate buffer of <i>size</i> values and assign it to the internal buffer with id <i>id</i> */ 
	void allocateBuffer(UInt id, UInt size);
	/*! \brief Free buffer i */
	void freeBuffer(UInt i) { freeBuffer(i, true); }
	/*! \brief Swap buffers i1 and i2 */
	void swapBuffers(UInt i1=0, UInt i2=1);
	/*! \brief Clear everything in this DT-Grid. */
	virtual void clear() { clear(true); }

	void setGamma(Real gamma) { this->gamma = gamma; }
	Real getGamma() { return gamma; }
	Real getBeta() { return beta; }
	Data getInsideConstant() { return insideConstant; }
	Data getOutsideConstant() { return outsideConstant; }
	/* Identical to getScale() */
	Real getDx() { return dx; }
	int getNarrowBandWidth() { return (int)floor( gamma/dx + (Real)0.5 ) + (Traits::includeSafeBand?1:0); }

	inline bool doesElementExist(Index x, Index y, Index z) const;

	inline UInt getMemUsed() const;
	UInt getMaxMemUsed() const;
	void resetMaxMemUsed();
	void updateMaxMemUsed(UInt additionalMemUsed=0);
	UInt getValueMemUsed() const;
	UInt getMaxValueMemUsed() const;
	UInt getIndexMemUsed() const;
	UInt getMaxIndexMemUsed() const;
	inline UInt getMemUsed(UInt id) const;
	UInt getValueMemUsed(UInt id) const;
	UInt getIndexMemUsed(UInt id) const;

	UInt getVa1D(Index x);
	UInt getVa2D(Index x, Index y);
	UInt getVa2DIndex(Index x, Index y);
	UInt getVa1DIndex(Index x);
	UInt getVa2D(UInt i);
	UInt getVa1D(UInt i);
	Data getValue(UInt i) const { return va3D[i]; }



	// These methods needed for the Volumetric DT-Grid
	// Ideally the methods should be protected, and the Volumetric DT-Grid declared a friend of this DT-Grid.
	UInt *getVa1D() { return va1D; }
	UInt *getVa2D() { return va2D; }
	Data *getVa3D() { return va3D; }

	UInt *getAA1D() { return aa1D; }
	UInt *getAA2D() { return aa2D; }
	UInt *getAA3D() { return aa3D; }

	Index *getXIndex() { return xIndex; }
	Index *getYIndex() { return yIndex; }
	Index *getZIndex() { return zIndex; }

	UInt getNumVa1D() const { return numVa1D; }
	UInt getNumVa2D() { return numVa2D; }
	UInt getNumVa3D() { return numVa3D; }
	UInt getNumXIndex() { return numXIndex; }
	UInt getNumYIndex() { return numYIndex; }
	UInt getNumZIndex() { return numZIndex; }
	UInt getNumAA1D() { return numXIndex/2; }
	UInt getNumAA2D() { return numYIndex/2; }
	UInt getNumAA3D() { return numZIndex/2; }

	void initBoundingBox();

	void clear(bool doDelete);



	/** computeGradient() assumes that the DTGrid stores a scalar field.
	 * Input: Real v[6]: is the values of the six adjacent grid points needed for computing the gradient
	 *                   using centralized second order finite differences.
	 * Output: Real n[3]: is the computed gradient (normal). */
	inline void computeGradient(Matrix::Vector3<Real>& n, Real v[6]) const;

    inline void computeGradient(const Locator& loc, Matrix::Vector3<Real>& n) const;
	/*! 
	 * Input:
	 *        n : The index (can be an x, y or z index)
	 *   nIndex : The index array (can be xIndex, yIndex or zIndex of DT-Grid, see below)
	 *  (fi,li) : The range of indices in 'nIndex' where findIndex() should search for 'n'.
	 * Output:
	 *       in : nIndex[in] is either 'n' or the smallest value in 'nIndex' larger than 'n', if 'n' is within the range of a p-column.
	 *            If 'n' is NOT within the range of a p-column, 'n' is clamped to the closest end point. If the closest end point is li (last index),
	 *            'in' is set to li+1 to maintain the invariant that it is the smallest grid point that is larger.
	 * Returns true if 'n' is inside a connected component and false otherwise */
	inline bool findIndex(Index n, Index *nIndex, int fi, int li, int *in) const;
	/** get the value at grid point (x,y,z). (fi,li) delimits the search in the 'xIndex' array, see below. */
	inline Data get(Index x, Index y, Index z, int fi, int li) const;
	/** get the value at grid point (y,z). (fi,li) delimits the search in the 'yIndex' array, see below. */
	inline Data get(Index y, Index z, int fi, int li) const;
	/** get the value at grid point (z). (fi,li) delimits the search in the 'zIndex' array, see below. */
	inline Data get(Index z, int fi, int li) const;
	/** get neighbor (x-1,y,z). ix is the index of 'x' in the xIndex array */
	inline Data getXM(Index x, Index y, Index z, int ix) const;
	/** get neighbor (x+1,y,z). ix is the index of 'x' in the xIndex array */
	inline Data getXP(Index x, Index y, Index z, int ix) const;
	/** get neighbor (y-1,z). iy is the index of 'y' in the yIndex array */
	inline Data getYM(Index y, Index z, int iy) const;
	/** get neighbor (y+1,z). iy is the index of 'y' in the yIndex array */
	inline Data getYP(Index y, Index z, int iy) const;
	/** get neighbor (z-1). iz is the index of 'z' in the zIndex array */
	inline Data getZM(Index z, int iz) const;
	/** get neighbor (z+1). iz is the index of 'z' in the zIndex array */
	inline Data getZP(Index z, int iz) const;



	// The following iterators:
	// _Iterator1D, _Iterator2D, _StencilIterator1D, _StencilIterator2D
	// are used by the dilation algorithm.

	class _Iterator1D
	{
	public:
	    _Iterator1D() { }
	    _Iterator1D(DTGrid *parent, bool begin);
	    
	    inline _Iterator1D& operator++();
	    inline _Iterator1D& operator++(int) { return operator++(); }
	    inline bool operator==(const _Iterator1D& iter) const;
	    inline bool operator!=(const _Iterator1D& iter) const;

	    inline UInt getValue() const;
	    inline void setValue(UInt value);
	    inline Index getX() const;

	protected:
	    DTGrid *parent;
	    UInt *va1D;
	    UInt value;
	    Index x;
	    UInt iv1D;
	    UInt ic1D;
	};


	class _Iterator2D
	{
	public:
	    _Iterator2D() { }
	    _Iterator2D(DTGrid *parent, bool begin);
	    
	    inline _Iterator2D& operator++();
	    inline _Iterator2D& operator++(int) { return operator++(); }
	    inline bool operator==(const _Iterator2D& iter) const;
	    inline bool operator!=(const _Iterator2D& iter) const;

	    inline UInt getArrayIndex() const { return iv2D; }
	    inline UInt getValue() const;
	    inline void setValue(UInt value);
	    inline Index getX() const;
	    inline Index getY() const;

	protected:
	    DTGrid *parent;
	    UInt *va2D;
	    UInt value;
	    Index x, y;
	    UInt iv1D, iv2D;
	    UInt ic1D, ic2D;
	};

	
	class _StencilIterator1D
	{
	public:
	    _StencilIterator1D() { valid = false; }
	    _StencilIterator1D(DTGrid *parent, Index width);
	    _StencilIterator1D(Index width, Index *xIndex, UInt iv1D, UInt *va1D, UInt sii, UInt eii);
	    _StencilIterator1D(Index width, Index *xIndex, UInt *aa1D, UInt *va1D, UInt sii, UInt eii);

	    inline bool isEndValid();
	    inline bool isStartValid();

	    inline void getEndValues(UInt *s, UInt *e);
	    inline void getStartValues(UInt *s, UInt *e);

	    inline void increment(Index xi);

	    inline Index getStartIndex() const;
	    inline Index getEndIndex() const;
	    inline void setStartAndEndIndex(Index s, Index e);

	    inline Index getCenterX() const { return (s+e)/2; }
	    inline bool isAtEnd() const { return ic1Ds > eii; }
	    inline Index getIncrement() const;

	    inline void reset();
	    inline void invalidate();
	    inline bool isValid();

	protected:
	    Index s, e;  // the start(s) and end(e) coordinates of the stencil
	    Index width;
	    UInt *va1D;
	    UInt iv1D;
	    UInt *aa1D;

	    Index *xIndex;
	    UInt eii, sii;

	    // start iterator
	    UInt iv1Ds;
	    UInt ic1Ds;
	    Index xs;

	    // end iterator
	    UInt iv1De;
	    UInt ic1De;
	    Index xe;

	    bool valid;
	};


	class _StencilIterator2D
	{
	public:
	    _StencilIterator2D() { }
	    _StencilIterator2D(DTGrid *parent, Index width);
	    ~_StencilIterator2D();

	    inline bool isEndValid(UInt k);
	    inline bool isStartValid(UInt k);

	    inline void getEndValues(UInt *s, UInt *e, UInt k);
	    inline void getStartValues(UInt *s, UInt *e, UInt k);

	    inline void increment(Index xi, Index yi);
	    inline void increment();

	    inline Index getCenterX() { return si1D.getCenterX(); }
	    inline Index getCenterY() { return (s+e)/2; }

	protected:
	    DTGrid *parent;
	    Index width;
	    Index s, e;
	    //UInt lsi;
	    Index lsi;                 // last stencil index in si2D array. The lsi index is a cyclic index.
	    _StencilIterator1D si1D;
	    //_StencilIterator1D *si2D;
	    _StencilIterator1D si2D[2*Traits::maxDilationWidth+1];

	    _Iterator2D iter2D;
	};


	_Iterator1D begin1D();
	_Iterator1D end1D();

	_Iterator2D begin2D();
	_Iterator2D end2D();

	_StencilIterator1D beginStencil1D(Index width);
	_StencilIterator2D beginStencil2D(Index width);



    protected:


	// SafeStorage is used by the safePush methods

	struct SafeStorage
	{
	    std::vector<UInt> va1D;        
	    std::vector<Index> xIndex;              
	    std::vector<UInt> aa1D;        

	    std::vector<UInt> va2D;         
	    std::vector<Index> yIndex;              
	    std::vector<UInt> aa2D;         

	    std::vector<Index> zIndex;              
	    std::vector<UInt> aa3D;
	    std::vector<Data> va3D;
	};

	

	template<class TemplateInterpolator>
	class TransformationIterator
	{
	public:

	    TransformationIterator() : grid2(NULL) { }
	    ~TransformationIterator() { }

	    TransformationIterator(
		const DTGrid *grid2, 
		Index bbox[3][2], const Matrix::Matrix4x4<Real>& trans, 
		TemplateInterpolator *interpolator, Data outsideConstant);

	    TransformationIterator(Index bbox[3][2]);

	    Data getValue() const;
	    bool operator!=(const TransformationIterator& oti) const;
	    bool operator==(const TransformationIterator& oti) const;
	    void getIndex(Index *i, Index *j, Index *k) const;
	    Index getI() const { return i; }
	    Index getJ() const { return j; }
	    Index getK() const { return k; }

	    TransformationIterator operator++();
	    TransformationIterator operator++(int) { return this->operator++(); }

	protected:
	    Data value;
	    TemplateInterpolator *interpolator;
	    Index bbox[3][2];   // all grid points inclusive
	    Matrix::Matrix4x4<Real> gridToGrid2;
	    const DTGrid *grid2;
	    Index i, j, k;
	    Data outsideConstant;
	};


	template<IteratorType it>
	    _Iterator<it> begin(UInt id=0) const;
	template<IteratorType it>
	    _Iterator<it> end(UInt id=0) const;


	template <class TemplateInterpolator> TransformationIterator<TemplateInterpolator> beginTransformationIterator(Index bbox[3][2], const Matrix::Matrix4x4<Real>& trans, TemplateInterpolator *interpolator, Data outsideConstant) const;
	template <class TemplateInterpolator> TransformationIterator<TemplateInterpolator> endTransformationIterator(Index bbox[3][2]) const;

	void dilate1D(Index width, DTGrid *dtgrid2);
	void dilate2D(Index width, DTGrid *dtgrid2);
	void dilate3D(Index width, DTGrid *dtgrid2);
	void initBuffers();
	void freeBuffer(UInt i, bool doDelete);
	void setEndMarkers();
	void assign(DTGrid *dtgrid2);
	void initTransformation();






	// data structure storage

	UInt *va1D;         // indices that point to the start of each (x)-column in yIndex
	UInt numVa1D;    
	Index *xIndex;              // start and end x indices of the connected components contained in the projection of the narrow band onto the x axis
	UInt numXIndex; 
	UInt lastXIndex;
	UInt *aa1D;         // points into va1D. gives the start (and in this case also end) index of each connected component in the x direction

	UInt *va2D;         // indices that point to the start of each (x,y)-column in zIndex
	UInt numVa2D;
	Index *yIndex;              // start and end y indices of the connected components in the y direction contained in the projection of the narrow band onto the (x,y) plane
	UInt numYIndex;
	UInt lastYIndex;
	UInt *aa2D;         // points into va2D. gives the start (and in this case also end) index of each connected component in the y direction

	Index *zIndex;              // start and end z index of each connected component in the z direction
	UInt numZIndex;
	UInt lastZIndex;
	UInt *aa3D;         // points into va3D. gives the start (and in this case also end) index of each connected component in the z direction

	UInt numVa3D;
	Data *va3D;				   // values of grid points in narrow band (will always correspond to values[0])

	Data *values[maxNumBuffers];             // the value-buffers 
	UInt numValues[maxNumBuffers];   // number of values in the value-buffers

	// prefs

	Real beta, gamma, dx;
	Data insideConstant, outsideConstant;


	// transformation, order: S(cale) R(otation) T(ranslation)

	Matrix::Matrix4x4<Real> gridToWorld, worldToGrid;
	Index bbox[3][2];                      // (min,max) for x, y and z in grid coordinates
	Index openLevelSetBBox[3][2];          // limiting bounding box used when Traits::openLevelSetsSupport = true
	bool openLevelSetBBoxValid;
	Real tx, ty, tz;
	Real rx, ry, rz;   // Angles in radians


	// misc

	SafeStorage *safeStorage;

	UInt maxMemUsed;
	UInt maxValueMemUsed;
	UInt maxIndexMemUsed;


	template<TubeType tt> friend class StencilTubeIterator;
    };


}


#include "DTGrid_Impl.h"


#endif
