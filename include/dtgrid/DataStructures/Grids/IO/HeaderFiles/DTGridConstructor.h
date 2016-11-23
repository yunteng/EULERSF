/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#ifndef _grids_dtgridconstructor_h_
#define _grids_dtgridconstructor_h_

#include <limits> 
#include "SvolSaverInterface.h"
#include "SvolLoaderInterface.h"
#include <string>


namespace Grids
{
    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    class DTGridConstructor
    {
    public:

	class SaveOutput
	{
	public:
	    UIntType storageRequirements;
	    UIntType numValues;
	};

	class SaveInput
	{
	public:

	    SaveInput() : skipValuesLargerThanGamma(true), mendSvol(true), saveSpecificSlices(false), useSpecifiedNumComponents(false) { }

	public:
	    RealType dx, beta, gamma, insideConstant, outsideConstant;
	    IndexType integralTranslation[3];
	    RealType translation[3];
	    RealType rotation[3];
	    std::string fileName;
	    bool skipValuesLargerThanGamma;
	    bool mendSvol;
	    bool saveSpecificSlices;
	    UIntType sliceBounds[2];

	    // FOR INITIALIZATION WITH SPECIFIED NUMBER OF COMPONENTS, FASTER WHEN KNOWN //
	    bool useSpecifiedNumComponents;
	    UIntType numXIndex, numYIndex, numZIndex, numVa1D, numVa2D, numVa3D;
	    IndexType bbox[3][2];
	    ///////////////////////////////////////////////////////////////////////////////
	};

	typedef typename SvolLoaderInterface<IndexType, RealType, DataType, UIntType>::LoadOutput LoadOutput;

	typedef typename SvolLoaderInterface<IndexType, RealType, DataType, UIntType>::InitParameters LoadInput;


    public:

	/* Saves only values with absolute value less than gamma.  */
	template<typename IteratorType>
	void save(const SaveInput& saveInput, SaveOutput *saveOutput, SvolSaverInterface<IndexType, RealType, DataType, UIntType> *svolSaver, IteratorType& iter);

	void load(const LoadInput& loadInput, LoadOutput *loadOutput, SvolLoaderInterface<IndexType, RealType, DataType, UIntType> *svolLoader);

    protected:

	/*! If pushToSaver is true, the value is passed onto the svol-saver */ 
	template<bool pushToSaver>
	void push3D(IndexType x, IndexType y, IndexType z, DataType val, SvolSaverInterface<IndexType, RealType, DataType, UIntType> *svolSaver);

	template<bool pushToSaver>
	void pushAndMend3D(IndexType x, IndexType y, IndexType z, DataType val, SvolSaverInterface<IndexType, RealType, DataType, UIntType> *svolSaver);

	template<bool pushToSaver>
	void push2D(IndexType x, IndexType y, UIntType val, SvolSaverInterface<IndexType, RealType, DataType, UIntType> *svolSaver);

	template<bool pushToSaver>
	void push1D(IndexType x, UIntType val, SvolSaverInterface<IndexType, RealType, DataType, UIntType> *svolSaver);

	// variables used during prediction of the number of values;
	IndexType lastXIndex, lastYIndex, lastZIndex;
	DataType lastVa3D;
	UIntType numVa3D, numVa2D, numVa1D, numXIndex, numYIndex, numZIndex;
    };




    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////



    //////////////////////////////////////////////////////////////////////////
    // LOADING
    //////////////////////////////////////////////////////////////////////////


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    void DTGridConstructor<IndexType, RealType, DataType, UIntType>::load(const LoadInput& loadInput, LoadOutput *loadOutput, SvolLoaderInterface<IndexType, RealType, DataType, UIntType> *svolLoader)
    {
	svolLoader->init(loadInput);
	svolLoader->load(loadOutput);
    }





    //////////////////////////////////////////////////////////////////////////
    // SAVING
    //////////////////////////////////////////////////////////////////////////


    template<typename IndexType, typename RealType, typename DataType, typename UIntType>
    template<typename IteratorType>
    void DTGridConstructor<IndexType, RealType, DataType, UIntType>::save(const SaveInput& saveInput, SaveOutput *saveOutput, SvolSaverInterface<IndexType, RealType, DataType, UIntType> *svolSaver, IteratorType& iter)
    {
	IndexType i, j, k;
	DataType v;
	RealType gamma;
	IndexType bbox[3][2];
	IndexType integralTranslation[3];
	IteratorType iter2;

	iter2 = iter;
	gamma = saveInput.gamma;
	integralTranslation[0] = saveInput.integralTranslation[0];
	integralTranslation[1] = saveInput.integralTranslation[1];
	integralTranslation[2] = saveInput.integralTranslation[2];

	if (saveInput.useSpecifiedNumComponents)
	{
	    numXIndex = saveInput.numXIndex;
	    numYIndex = saveInput.numYIndex;
	    numZIndex = saveInput.numZIndex;
	    numVa1D = saveInput.numVa1D;
	    numVa2D = saveInput.numVa2D;
	    numVa3D = saveInput.numVa3D;
	    bbox[0][0] = saveInput.bbox[0][0];
	    bbox[0][1] = saveInput.bbox[0][1];
	    bbox[1][0] = saveInput.bbox[1][0];
	    bbox[1][1] = saveInput.bbox[1][1];
	    bbox[2][0] = saveInput.bbox[2][0];
	    bbox[2][1] = saveInput.bbox[2][1];
	}
	else
	{

	    bbox[0][0] = std::numeric_limits<IndexType>::max();
	    bbox[1][0] = std::numeric_limits<IndexType>::max();
	    bbox[2][0] = std::numeric_limits<IndexType>::max();
	    bbox[0][1] = -std::numeric_limits<IndexType>::max();
	    bbox[1][1] = -std::numeric_limits<IndexType>::max();
	    bbox[2][1] = -std::numeric_limits<IndexType>::max();

	    numXIndex = 0;
	    numYIndex = 0;
	    numZIndex = 0;
	    numVa1D = 0;
	    numVa2D = 0;
	    numVa3D = 0;
	    // The initialization of these indices not significant
	    lastXIndex = 0;
	    lastYIndex = 0;
	    lastZIndex = 0;
	    lastVa3D = 0;

	    // FIRST PREDICT THE NUMBER OF INDICES AND VALUES THAT WILL BE ADDED TO THE DT-GRID

	    if (iter.hasNext())
	    {
		if (saveInput.saveSpecificSlices)
		{
		    UIntType currentSlice = 0;
		    IndexType lastSeenX = iter.getI();

		    while (iter.hasNext() && currentSlice<=saveInput.sliceBounds[1])
		    {
			i = (IndexType)iter.getI();
			j = (IndexType)iter.getJ();
			k = (IndexType)iter.getK();

			v = static_cast<DataType>(iter.getValue());

			if (!saveInput.skipValuesLargerThanGamma || fabs(v)<gamma)
			{
			    if (currentSlice >= saveInput.sliceBounds[0])
			    {
				if (saveInput.mendSvol)
				    pushAndMend3D<false>(i, j, k, v, svolSaver);
				else
				    push3D<false>(i, j, k, v, svolSaver);			

				if (i<bbox[0][0]) { bbox[0][0] = i; }
				if (i>bbox[0][1]) { bbox[0][1] = i; }
				if (j<bbox[1][0]) { bbox[1][0] = j; }
				if (j>bbox[1][1]) { bbox[1][1] = j; }
				if (k<bbox[2][0]) { bbox[2][0] = k; }
				if (k>bbox[2][1]) { bbox[2][1] = k; }
			    }
			}

			if (i != lastSeenX)
			{
			    currentSlice++;
			    lastSeenX = i;
			}

			iter.next();
		    }
		}
		else
		{
		    while (iter.hasNext())
		    {
			i = (IndexType)iter.getI();
			j = (IndexType)iter.getJ();
			k = (IndexType)iter.getK();

			v = static_cast<DataType>(iter.getValue());

			if (!saveInput.skipValuesLargerThanGamma || fabs(v)<gamma)
			{
			    if (saveInput.mendSvol)
				pushAndMend3D<false>(i, j, k, v, svolSaver);
			    else
				push3D<false>(i, j, k, v, svolSaver);			

			    if (i<bbox[0][0]) { bbox[0][0] = i; }
			    if (i>bbox[0][1]) { bbox[0][1] = i; }
			    if (j<bbox[1][0]) { bbox[1][0] = j; }
			    if (j>bbox[1][1]) { bbox[1][1] = j; }
			    if (k<bbox[2][0]) { bbox[2][0] = k; }
			    if (k>bbox[2][1]) { bbox[2][1] = k; }

			}
			iter.next();
		    }
		}
	    }

	}


	bbox[0][0] += integralTranslation[0];
	bbox[0][1] += integralTranslation[0];
	bbox[1][0] += integralTranslation[1];
	bbox[1][1] += integralTranslation[1];
	bbox[2][0] += integralTranslation[2];
	bbox[2][1] += integralTranslation[2];


	////////////////////////////////////////////////////////////////

	// the minimal storage requirements of a dtgrid with support for random access using linear search within each column (can be improved by using only the required number of bits to store indices and aa-indices).
	saveOutput->storageRequirements = numVa3D * sizeof(DataType) + numVa2D * 2* sizeof(UIntType) + numVa1D * 2 * sizeof(UIntType) + numXIndex* sizeof(IndexType) + numYIndex * sizeof(IndexType) + numZIndex * sizeof(IndexType); 
	saveOutput->numValues = numVa3D;

	////////////////////////////////////////////////////////////////
	typename SvolSaverInterface<IndexType, RealType, DataType, UIntType>::InitParameters svolSaverInitParams;

	svolSaverInitParams.numXIndex = numXIndex;
	svolSaverInitParams.numYIndex = numYIndex;
	svolSaverInitParams.numZIndex = numZIndex;
	svolSaverInitParams.numVa1D = numVa1D;
	svolSaverInitParams.numVa2D = numVa2D;
	svolSaverInitParams.numVa3D = numVa3D;
	svolSaverInitParams.dx = saveInput.dx;
	svolSaverInitParams.beta = saveInput.beta;
	svolSaverInitParams.gamma = saveInput.gamma;
	svolSaverInitParams.insideConstant = saveInput.insideConstant;
	svolSaverInitParams.outsideConstant = saveInput.outsideConstant;
	svolSaverInitParams.translation[0] = saveInput.translation[0];
	svolSaverInitParams.translation[1] = saveInput.translation[1];
	svolSaverInitParams.translation[2] = saveInput.translation[2];
	svolSaverInitParams.rotation[0] = saveInput.rotation[0];
	svolSaverInitParams.rotation[1] = saveInput.rotation[1];
	svolSaverInitParams.rotation[2] = saveInput.rotation[2];
	svolSaverInitParams.bbox[0][0] = bbox[0][0];
	svolSaverInitParams.bbox[0][1] = bbox[0][1];
	svolSaverInitParams.bbox[1][0] = bbox[1][0];
	svolSaverInitParams.bbox[1][1] = bbox[1][1];
	svolSaverInitParams.bbox[2][0] = bbox[2][0];
	svolSaverInitParams.bbox[2][1] = bbox[2][1];
	svolSaverInitParams.fileName = saveInput.fileName;

	//svolSaverInitParams.printASCII(std::cout);

	svolSaver->init(svolSaverInitParams);

	////////////////////////////////////////////////////////////////

	if (iter2.hasNext())
	{
	    if (saveInput.saveSpecificSlices)
	    {
		UIntType currentSlice = 0;
		IndexType lastSeenX = iter2.getI() + integralTranslation[0];

		while (iter2.hasNext() && currentSlice<=saveInput.sliceBounds[1])
		{
		    i = iter2.getI() + integralTranslation[0];
		    j = iter2.getJ() + integralTranslation[1];
		    k = iter2.getK() + integralTranslation[2];

		    v = static_cast<DataType>(iter2.getValue());

		    if (!saveInput.skipValuesLargerThanGamma || fabs(v)<gamma)
		    {
			if (currentSlice >= saveInput.sliceBounds[0])
			{
			    if (saveInput.mendSvol)
				pushAndMend3D<true>(i, j, k, v, svolSaver);
			    else
				push3D<true>(i, j, k, v, svolSaver);			
			}
		    }

		    if (i != lastSeenX)
		    {
			currentSlice++;
			lastSeenX = i;
		    }

		    iter2.next();
		}
	    }
	    else
	    {

		while (iter2.hasNext())
		{
		    i = iter2.getI() + integralTranslation[0];
		    j = iter2.getJ() + integralTranslation[1];
		    k = iter2.getK() + integralTranslation[2];

		    v = static_cast<DataType>(iter2.getValue());

		    if (!saveInput.skipValuesLargerThanGamma || fabs(v)<gamma)
		    {
			if (saveInput.mendSvol)
			    pushAndMend3D<true>(i, j, k, v, svolSaver);
			else
			    push3D<true>(i, j, k, v, svolSaver);			

		    }

		    iter2.next();

		}
	    }

	}

	iter2.commit();

	svolSaver->flushAll();

    }




    template<typename IndexType, typename RealType, typename DataType, typename UIntType> template<bool pushToSaver>
    void DTGridConstructor<IndexType, RealType, DataType, UIntType>::push1D(IndexType x, UIntType val, SvolSaverInterface<IndexType, RealType, DataType, UIntType> *svolSaver)
    {
	if ( pushToSaver ? (svolSaver->xIndexSize() > 0 && svolSaver->xIndexBack() == x-1) : (numXIndex > 0 && lastXIndex == x-1) )
	{
	    // this element is part of an existing connected component
	    if (pushToSaver)
	    {
		svolSaver->xIndexBack() = x;
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->aa1DBack() = svolSaver->va1DSize();
		}
	    }
	    else
	    {
		lastXIndex = x;
	    }
	}
	else
	{
	    // begin new connected component
	    if (pushToSaver)
	    {
		svolSaver->pushAA1D( svolSaver->va1DSize() );
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->pushAA1D( svolSaver->va1DSize() );
		}
		svolSaver->pushXIndex( x );
		svolSaver->pushXIndex( x );
	    }
	    else
	    {
		numXIndex += 2;
		lastXIndex = x;
	    }
	}

	if (pushToSaver)
	{
	    svolSaver->pushVa1D(val);
	}
	else
	{
	    numVa1D++;
	}
    }



    template<typename IndexType, typename RealType, typename DataType, typename UIntType> template<bool pushToSaver>
    void DTGridConstructor<IndexType, RealType, DataType, UIntType>::push2D(IndexType x, IndexType y, UIntType val, SvolSaverInterface<IndexType, RealType, DataType, UIntType> *svolSaver)
    {
	if ( pushToSaver ? (svolSaver->va2DSize()==0 || svolSaver->xIndexBack() != x) : (numVa2D==0 || lastXIndex != x) )
	{
	    // there are no elements in this column
	    push1D<pushToSaver>(x, (pushToSaver ? svolSaver->yIndexSize() : numYIndex), svolSaver);

	    // begin new connected component
	    if (pushToSaver)
	    {
		svolSaver->pushAA2D(svolSaver->va2DSize());
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->pushAA2D( svolSaver->va2DSize() );
		}
		svolSaver->pushYIndex(y);
		svolSaver->pushYIndex(y);
	    }
	    else
	    {
		numYIndex += 2;
		lastYIndex = y;
	    }
	}
	else if ( pushToSaver ? (svolSaver->yIndexBack() == y-1) : (lastYIndex == y-1) )
	{
	    // this element is part of an existing connected component
	    if (pushToSaver)
	    {
		svolSaver->yIndexBack() = y;
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->aa2DBack() = svolSaver->va2DSize();
		}
	    }
	    else
	    {
		lastYIndex = y;
	    }
	}
	else
	{
	    // begin new connected component
	    if (pushToSaver)
	    {
		svolSaver->pushAA2D(svolSaver->va2DSize());
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->pushAA2D( svolSaver->va2DSize() );
		}
		svolSaver->pushYIndex(y);
		svolSaver->pushYIndex(y);
	    }
	    else
	    {
		numYIndex += 2;
		lastYIndex = y;
	    }
	}

	if (pushToSaver)
	{
	    svolSaver->pushVa2D(val);
	}
	else
	{
	    numVa2D++;
	}
    }



    template<typename IndexType, typename RealType, typename DataType, typename UIntType> template<bool pushToSaver>
    void DTGridConstructor<IndexType, RealType, DataType, UIntType>::push3D(IndexType x, IndexType y, IndexType z, DataType val, SvolSaverInterface<IndexType, RealType, DataType, UIntType> *svolSaver)
    {
	if ( pushToSaver ? ( svolSaver->va3DSize() == 0 || svolSaver->xIndexBack() != x || svolSaver->yIndexBack() != y ) : ( numVa3D==0 || lastXIndex != x || lastYIndex != y) )
	{
	    // there are no elements in this column
	    push2D<pushToSaver>(x, y, ( pushToSaver ? svolSaver->zIndexSize() : numZIndex ), svolSaver);

	    // begin new connected component
	    if (pushToSaver)
	    {
		svolSaver->pushAA3D(svolSaver->va3DSize());
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->pushAA3D( svolSaver->va3DSize() );
		}
		svolSaver->pushZIndex(z);
		svolSaver->pushZIndex(z);
	    }
	    else
	    {
		numZIndex += 2;
		lastZIndex = z;
	    }
	}
	else if ( pushToSaver ? (svolSaver->zIndexBack() == z-1) : (lastZIndex == z-1) )
	{
	    // this element is part of an existing connected component
	    if (pushToSaver)
	    {
		svolSaver->zIndexBack() = z;
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->aa3DBack() = svolSaver->va3DSize();
		}
	    }
	    else
	    {
		lastZIndex = z;
	    }
	}
	else
	{
	    // begin new connected component
	    if (pushToSaver)
	    {
		svolSaver->pushAA3D(svolSaver->va3DSize());
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->pushAA3D( svolSaver->va3DSize() );
		}
		svolSaver->pushZIndex(z);
		svolSaver->pushZIndex(z);
	    }
	    else
	    {
		numZIndex += 2;
		lastZIndex = z;
	    }
	}

	if (pushToSaver)
	{
	    svolSaver->pushVa3D(val);
	}
	else
	{
	    lastVa3D = val;
	    numVa3D++;
	}
    }



    template<typename IndexType, typename RealType, typename DataType, typename UIntType> template<bool pushToSaver>
    void DTGridConstructor<IndexType, RealType, DataType, UIntType>::pushAndMend3D(IndexType x, IndexType y, IndexType z, DataType val, SvolSaverInterface<IndexType, RealType, DataType, UIntType> *svolSaver)
    {
	if ( pushToSaver ? ( svolSaver->va3DSize() == 0 || svolSaver->xIndexBack() != x || svolSaver->yIndexBack() != y ) : ( numVa3D==0 || lastXIndex != x || lastYIndex != y) )
	{
	    if ( pushToSaver && (svolSaver->va3DSize() != 0) )
	    {
		// make sure the end of last connected component was positive
		// we do not handle if the last grid point added was erroneous
		svolSaver->va3DBack() = fabs(svolSaver->va3DBack());
	    }
	    // also make sure that the start of the first connected component in p-column is positive
	    val = fabs(val);

	    // there are no elements in this column
	    push2D<pushToSaver>(x, y, ( pushToSaver ? svolSaver->zIndexSize() : numZIndex ), svolSaver);

	    if (pushToSaver)
	    {
		// begin new connected component
		svolSaver->pushAA3D(svolSaver->va3DSize());
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->pushAA3D( svolSaver->va3DSize() );
		}
		svolSaver->pushZIndex(z);
		svolSaver->pushZIndex(z);
	    }
	    else
	    {
		numZIndex += 2;
		lastZIndex = z;
	    }
	}
	else if ( pushToSaver ? (svolSaver->zIndexBack() == z-1) : (lastZIndex == z-1) )
	{
	    if (pushToSaver)
	    {
		// this element is part of an existing connected component
		svolSaver->zIndexBack() = z;
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->aa3DBack() = svolSaver->va3DSize();
		}
	    }
	    else
	    {
		lastZIndex = z;
	    }
	}
	else
	{
	    // begin new connected component
	    if (pushToSaver)
	    {
		svolSaver->pushAA3D(svolSaver->va3DSize());
		if (svolSaver->pushTwoAA())
		{
		    svolSaver->pushAA3D( svolSaver->va3DSize() );
		}
		svolSaver->pushZIndex(z);
		svolSaver->pushZIndex(z);
	    }
	    else
	    {
		numZIndex += 2;
		lastZIndex = z;
	    }

	    // make sure that the "adjacent" ends of connected components have consistent signs.
	    // favour outside
	    if ( val<0 && ( pushToSaver ? svolSaver->va3DBack()>0 : lastVa3D>0 ) )
	    {
		val = fabs(val);
	    }
	    else if ( pushToSaver && val>0 && svolSaver->va3DBack()<0 )
	    {
		svolSaver->va3DBack() = fabs(svolSaver->va3DBack());
	    }
	}

	if (pushToSaver)
	{
	    svolSaver->pushVa3D(val);
	}
	else
	{
	    numVa3D++;
	    lastVa3D = val;
	}

    }


}





#endif
