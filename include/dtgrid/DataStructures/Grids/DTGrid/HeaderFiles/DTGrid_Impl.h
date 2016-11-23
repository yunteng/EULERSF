/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/


#undef FIND_INDEX_EXISTS
#undef FIND_INDEX
#undef INCREMENT_FAST_Z
#undef INCREMENT_UNTIL


namespace Grids
{



#define INCREMENT_FAST_Z(outsideValue, zn, iterId)				    \
	if ( si[iterId].valid )							    \
	{									    \
	    if (si[iterId].z < zn)						    \
	    {									    \
		if (si[iterId].z == si[iterId].zIndex[si[iterId].ic3D])		    \
		{								    \
		    si[iterId].ic3D++;						    \
		    si[iterId].z = si[iterId].zIndex[si[iterId].ic3D];		    \
		    if (parent->va2D[si[iterId].iv2D] == si[iterId].ic3D)	    \
		    {								    \
			si[iterId].iv2D++;					    \
			if (si[iterId].y == parent->yIndex[si[iterId].ic2D])	    \
			{							    \
			    si[iterId].ic2D++;					    \
			    if (parent->va1D[si[iterId].iv1D] == si[iterId].ic2D)   \
			    {							    \
				si[iterId].iv1D++;				    \
				if (si[iterId].x == parent->xIndex[si[iterId].ic1D])\
				{						    \
				    si[iterId].ic1D++;				    \
				    si[iterId].x = parent->xIndex[si[iterId].ic1D]; \
				    si[iterId].ic1D++;				    \
				}						    \
				else						    \
				{						    \
				    si[iterId].x++;				    \
				}						    \
			    }							    \
			    si[iterId].y = parent->yIndex[si[iterId].ic2D];	    \
			    si[iterId].ic2D++;					    \
			}							    \
			else							    \
			{							    \
			    si[iterId].y++;					    \
			}							    \
			si[iterId].valid = false;				    \
		    }								    \
		    si[iterId].ic3D++;						    \
		    si[iterId].iv3D++;						    \
		    si[iterId].value = outsideValue;				    \
		}								    \
		else								    \
		{								    \
		    si[iterId].z++;						    \
		    si[iterId].iv3D++;						    \
		    si[iterId].value = si[iterId].va3D[si[iterId].iv3D];	    \
		}								    \
	    }									    \
	    else if (si[iterId].z == zn)					    \
	    {									    \
		si[iterId].value = si[iterId].va3D[si[iterId].iv3D];		    \
	    }									    \
	}








#define INCREMENT_UNTIL(outsideValue, xn, yn, zn, iterId)					\
	if (si[iterId].x < xn)									\
	{											\
	    do											\
	    {											\
		si[iterId]++;									\
	    }											\
	    while (si[iterId].x < xn);								\
	}											\
	if (xn==si[iterId].x)									\
	{											\
	    if (si[iterId].y < yn)								\
	    {											\
	        do										\
		{										\
    		    si[iterId]++;								\
		}										\
		while (si[iterId].x==xn && si[iterId].y<yn);					\
	    }											\
	    if (xn==si[iterId].x && yn==si[iterId].y)						\
	    {											\
		while (xn==si[iterId].x && yn==si[iterId].y && si[iterId].z < zn)		\
		{										\
		    si[iterId].iv3D++;								\
		    if (si[iterId].z == si[iterId].zIndex[si[iterId].ic3D])			\
		    {										\
			si[iterId].ic3D++;							\
			if (parent->va2D[si[iterId].iv2D] == si[iterId].ic3D)			\
			{									\
			    si[iterId].iv2D++;							\
			    if (si[iterId].y == parent->yIndex[si[iterId].ic2D])		\
			    {									\
				si[iterId].ic2D++;						\
				if (parent->va1D[si[iterId].iv1D] == si[iterId].ic2D)		\
				{								\
				    si[iterId].iv1D++;						\
				    if (si[iterId].x == parent->xIndex[si[iterId].ic1D])	\
				    {								\
					si[iterId].ic1D++;					\
					si[iterId].x = parent->xIndex[si[iterId].ic1D];		\
					si[iterId].ic1D++;					\
				    }								\
				    else							\
				    {								\
					si[iterId].x++;						\
				    }								\
				}								\
				si[iterId].y = parent->yIndex[si[iterId].ic2D];			\
				si[iterId].ic2D++;						\
			    }									\
			    else								\
			    {									\
				si[iterId].y++;							\
			    }									\
			}									\
			si[iterId].z = si[iterId].zIndex[si[iterId].ic3D];			\
			si[iterId].ic3D++;							\
		    }										\
		    else									\
		    {										\
			si[iterId].z++;								\
		    }										\
		}										\
		if (xn==si[iterId].x && yn==si[iterId].y)					\
		{										\
		    si[iterId].valid = true;							\
		    if (zn==si[iterId].z)							\
		    {										\
			si[iterId].value = si[iterId].va3D[si[iterId].iv3D];			\
		    }										\
		    else									\
		    {										\
			/* si[iterId].value = ( si[iterId].va3D[si[iterId].iv3D] < 0 ? si[iterId].insideConstant : si[iterId].outsideConstant ); */ \
			si[iterId].value = outsideValue;					\
		    }										\
		}										\
		else										\
		{										\
		    si[iterId].valid = false;							\
		    /* si[iterId].value = ( si[iterId].va3D[si[iterId].iv3D] < 0 ? si[iterId].insideConstant : si[iterId].outsideConstant ); */ \
		    si[iterId].value = outsideValue;						\
		}										\
	    }											\
	    else										\
	    {											\
		si[iterId].valid = false;							\
		/* si[iterId].value = ( si[iterId].va3D[si[iterId].iv3D] < 0 ? si[iterId].insideConstant : si[iterId].outsideConstant ); */ \
		si[iterId].value = outsideValue;						\
	    }											\
	}											\
	else											\
	{											\
	    si[iterId].valid = false;								\
	    si[iterId].value = outsideValue;							\
	}



#define FIND_INDEX(n, nIndex, fi, li, in, returnValue)          \
    if (Traits::randomAccessType == USE_LINEAR_SEARCH)         \
    {                                                           \
	if (n < nIndex[fi] || n > nIndex[li])                   \
	{							\
	    return returnValue;					\
	}							\
	in = fi+1;						\
	while (n > nIndex[in])					\
	{							\
	    in += 2;						\
	}							\
	in -= 1;						\
	if (n < nIndex[in])					\
	{							\
	    return returnValue;					\
	}					                \
    }								\
    else							\
    {								\
	if (n < nIndex[fi] || n > nIndex[li])                   \
	{							\
	    return returnValue;					\
	}							\
        {							\
	    int in1, in2, it;		                        \
	    in1 = 0;						\
	    in2 = (li-fi)+1;					\
	    while (in2 != in1+2)				\
	    {							\
		it = ((in1+in2)>>2)<<1;				\
		if (n < nIndex[it+fi])				\
		{						\
		    in2 = it;					\
		}						\
		else						\
		{						\
		    in1 = it;					\
		}						\
	    }							\
	    in = in1 = in1 + fi;				\
	    if ( n > nIndex[in1+1] )				\
	    {							\
		return returnValue;				\
	    }							\
	}							\
    }

// previously defined as FIND_INDEX_EXISTS_LINEAR, changed by Yun
#define FIND_INDEX_EXISTS(n, nIndex, fi, li, in)		\
    if (Traits::randomAccessType == USE_LINEAR_SEARCH)         \
    {                                                           \
	in = fi+1;						\
	while (n > nIndex[in])					\
	{							\
	    in += 2;						\
	}							\
	in -= 1;						\
    }								\
    else							\
    {								\
        {							\
	    int in1, in2, it;		                        \
	    in1 = 0;						\
	    in2 = (li-fi)+1;					\
	    while (in2 != in1+2)				\
	    {							\
		it = ((in1+in2)>>2)<<1;				\
		if (n < nIndex[it+fi])				\
		{						\
		    in2 = it;					\
		}						\
		else						\
		{						\
		    in1 = it;					\
		}						\
	    }							\
	    in = in1 = in1 + fi;				\
	}							\
    }



///////////////////////////////////////////////////////////////////////////////////
// CLOSEST POINT METHOD
///////////////////////////////////////////////////////////////////////////////////

/**
*  Input:
*           p : indices of grid-point at which closest point should be evaluated.
*
*  Output:
*           n : normal (points towards closest point), valid if *inNarrowBand == true
*           cd : closest distance to interface, clamped to the width of the narrow band if *inNarrowBand == false
*           inside : true if p is inside or on interface, false otherwise
*           inNarrowBand : true if p lies in the narrow band
*
*  The narrow band should be a signed distance function, |n| = 1.
*/
template<class Traits>
void DTGrid<Traits>::closestPoint(const Index p[3], Matrix::Vector3<Real>& n, Real *cd, bool *inside, bool *inNarrowBand) const
{
    Index x = p[0];
    Index y = p[1];
    Index z = p[2];
    closestPoint(x, y, z, n, cd, inside, inNarrowBand);
}



template<class Traits>
void DTGrid<Traits>::closestPoint(const Index x, const Index y, const Index z, Matrix::Vector3<Real>& n, Real *cd, bool *inside, bool *inNarrowBand) const
{
    int ix=0, iy=0, iz=0, fi, li, tmp;
    Real v[6];
    // Index x = p[0];
    // Index y = p[1];
    // Index z = p[2];


    // determine if 'p' exists in the narrow band
    fi = 0;
    li = lastXIndex;
    if (findIndex(x, xIndex, fi, li, &ix))
    {
		tmp = aa1D[(ix)>>1] + x - xIndex[ix];
		fi = va1D[tmp];
		li = va1D[tmp+1]-1;
		if (findIndex(y, yIndex, fi, li, &iy))
		{
		    tmp = aa2D[(iy)>>1] + y - yIndex[iy];
		    iz = fi = va2D[tmp];
		    li = va2D[tmp+1]-1;
		    if (findIndex(z, zIndex, fi, li, &iz))
		    {
				// at this point we know the grid point (x,y,z) exists

				tmp = aa3D[(iz)>>1] + z - zIndex[iz];

				*cd = va3D[tmp];
				*inside = (*cd < 0);
				*inNarrowBand = true;

				// at this point we know 'p' exists in the narrow band

				v[0] = getXM(x, y, z, ix);  // find neighbor to the left in the x direction
				v[1] = getXP(x, y, z, ix);  // find neighbor to the right in the x direction
				v[2] = getYM(y, z, iy);     // find neighbor to the left in the y direction
				v[3] = getYP(y, z, iy);     // find neighbor to the right in the y direction
				v[4] = getZM(z, iz);        // find neighbor to the left in the z direction
				v[5] = getZP(z, iz);        // find neighbor to the right in the z direction

				computeGradient(n, v);
				return;
		    }
		    else
		    {
				*inNarrowBand = false;
				*inside = ( va3D[aa3D[(iz)>>1]] < 0 ? true : false );
		    }
		}
		else
		{
		    if (!Traits::openLevelSetsSupport)
		    {
				*inside = *inNarrowBand = false;
		    }
		    else
		    {
				//OpenLS:
				UInt aa2DIndex = aa2D[iy>>1];
				UInt va2DValue = va2D[aa2DIndex];
				UInt aa3DIndex = aa3D[va2DValue>>1];
				*inside =  (va3D[aa3DIndex] < 0 ? true : false);
				*inNarrowBand = false;
		    }
		}
    }
    else
    {
		if (!Traits::openLevelSetsSupport)
		{
		    *inside = *inNarrowBand = false;
		}
		else
		{
		    //OpenLS:
		    UInt aa1DIndex = aa1D[ix>>1];
		    UInt va1DValue = va1D[aa1DIndex];
		    UInt aa2DIndex = aa2D[va1DValue>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    *inside =  (va3D[aa3DIndex] < 0 ? true : false);
		    *inNarrowBand = false;
		}
    }
}
/*
template<class Traits>
bool DTGrid<Traits>::lookupDistanceAndNormal(const Real pos[3], Real& dist, Matrix::Vector3<Real>& normal) const
{
    Real cds[8] = {0.0};
    Matrix::Vector3<Real> normals[8];
    bool gridInside[8];
    bool gridInNarrowBand[8];

    Matrix::Vector3<Real> vec3 = worldToGrid * Matrix::Vector3<Real>(pos);
    // convert vec3 to a unit offset
    Matrix::Vector3<Index> minVec = Matrix::Vector3<Real>::floor(vec3);
    vec3 = vec3 - Matrix::Vector3<Real>::floor(vec3);

    Index x = minVec[0];
    Index y = minVec[1];
    Index z = minVec[2];

    Locator gridLoc[8];

    int i = 0;
    for(i = 0; i < 8; i++){
    	gridInNarrowBand[i] = getLocatorAndValue(x + i / 4, y + (i % 4) / 2, z + i % 2, gridLoc[i], cds[i]);
    	if(gridInNarrowBand[i]){
    		break;
    	}
    }
    // none of the 8 grid points are in the narrow band
    if(i == 8){
    	dist = cds[0];
    	return false;
    }
    // not all the grid points exist in the narrow band, just use partial result
    if(i != 0){
    	computeGradient(gridLoc[i], normal);
    	dist = cds[i];
    	return true;
    }

	gridInNarrowBand[1] = getZP(gridLoc[0], gridLoc[1]);
	gridInNarrowBand[2] = getYP(gridLoc[0], gridLoc[2]);
	gridInNarrowBand[4] = getXP(gridLoc[0], gridLoc[4]);

	bool useInterpolation = true;
	if(gridInNarrowBand[1]){
		gridInNarrowBand[3] = getYP(gridLoc[1], gridLoc[3]);
		gridInNarrowBand[5] = getXP(gridLoc[1], gridLoc[5]);
		if(gridInNarrowBand[5]){
			gridInNarrowBand[7] = getYP(gridLoc[5], gridLoc[7]);
		}else{
			useInterpolation = false;
		}
		if(gridInNarrowBand[2]){
			gridInNarrowBand[6] = getXP(gridLoc[2], gridLoc[6]);
		}else{
			useInterpolation = false;
		}
	}else{
		useInterpolation =  false;
	}
	if(!useInterpolation){
		computeGradient(gridLoc[i], normal);
    	dist = cds[i];
    	return true;
	}
    	

	for(int i = 0; i < 8; i++){
		if(gridInNarrowBand[i])
			cds[i] = va3D[gridLoc[i].iv3D];   		
	}
	for(int i = 0; i < 8; i++){
		if(gridInNarrowBand[i])
			computeGradient(gridLoc[i], normals[i]);  		
	}


    

    Math::TrilinearInterpolator<Real, Real> distInterpolator = Math::TrilinearInterpolator<Real, Real>();
    dist = distInterpolator.interp(cds, vec3, gridInNarrowBand);

    Math::TrilinearInterpolator<Matrix::Vector3<Real>, Real> normalInterpolator = Math::TrilinearInterpolator<Matrix::Vector3<Real>, Real>();
    normal = normalInterpolator.interp(normals, vec3, gridInNarrowBand);
    normal.normalize();
    if(normal.length() == 0){
        for(int i = 0; i < 8; i++){
            if(gridInNarrowBand[i]){
                normal = normals[i];
                break;
            }
        }
    }

    return true;
}*/

template<class Traits>
bool DTGrid<Traits>::lookupDistanceAndNormal(const Real pos[3], Real& dist, Matrix::Vector3<Real>& normal) const
{
    Real cds[8] = {0.0};
    Matrix::Vector3<Real> normals[8];
    bool gridInside[8];
    bool gridInNarrowBand[8];

    Matrix::Vector3<Real> vec3 = worldToGrid * Matrix::Vector3<Real>(pos);
    // convert vec3 to a unit offset
    Matrix::Vector3<Index> minVec = Matrix::Vector3<Real>::floor(vec3);
    vec3 = vec3 - Matrix::Vector3<Real>::floor(vec3);

    Index x = minVec[0];
    Index y = minVec[1];
    Index z = minVec[2];

    closestPoint(x, y, z, normals[0], &cds[0], &gridInside[0], &gridInNarrowBand[0]);
    closestPoint(x, y, z + 1, normals[1], &cds[1], &gridInside[1], &gridInNarrowBand[1]);
    closestPoint(x, y + 1, z, normals[2], &cds[2], &gridInside[2], &gridInNarrowBand[2]);
    closestPoint(x, y + 1, z + 1, normals[3], &cds[3], &gridInside[3], &gridInNarrowBand[3]);
    closestPoint(x + 1, y, z, normals[4], &cds[4], &gridInside[4], &gridInNarrowBand[4]);
    closestPoint(x + 1, y, z + 1, normals[5], &cds[5], &gridInside[5], &gridInNarrowBand[5]);
    closestPoint(x + 1, y + 1, z, normals[6], &cds[6], &gridInside[6], &gridInNarrowBand[6]);
    closestPoint(x + 1, y + 1, z + 1, normals[7], &cds[7], &gridInside[7], &gridInNarrowBand[7]);

    int existInNarrowBand = 0;
    int outsideCnt = 0;
    for(int i = 0; i < 8; i++){
        existInNarrowBand += gridInNarrowBand[i] ? 1 : 0;
        outsideCnt += gridInside[i] ? 0 : 1;
    }
    if(existInNarrowBand == 0){
    	if(outsideCnt == 0){
    		dist = insideConstant;
    		return false;
    	}else if(outsideCnt == 8){
    		dist = outsideConstant;
        	return false;
    	}else{
    		Core::throwDefaultException(std::string("lookup point on the boundary but not captured by the narrow band!!!"), __FILE__, __LINE__);
    	}
    }
    

    Math::TrilinearInterpolator<Real, Real> distInterpolator = Math::TrilinearInterpolator<Real, Real>();
    dist = distInterpolator.interp(cds, vec3, gridInNarrowBand);

    Math::TrilinearInterpolator<Matrix::Vector3<Real>, Real> normalInterpolator = Math::TrilinearInterpolator<Matrix::Vector3<Real>, Real>();
    normal = normalInterpolator.interp(normals, vec3, gridInNarrowBand);
    normal.normalize();
    if(normal.length() == 0){
        for(int i = 0; i < 8; i++){
            if(gridInNarrowBand[i]){
                normal = normals[i];
                break;
            }
        }
    }

    return true;
}

template<class Traits>
bool DTGrid<Traits>::lookupDistanceAndNormalFast(const Matrix::Vector3<Real>& pos, Real& dist, Matrix::Vector3<Real>& normal) const
{
    Real cd = 0.0;
    // Matrix::Vector3<Real> normal;
    bool gridInside;
    bool gridInNarrowBand;

    Matrix::Vector3<Real> vec3 = worldToGrid * pos;
    // convert vec3 to a unit offset
    Matrix::Vector3<Index> minVec = Matrix::Vector3<Real>::round(vec3);
    vec3 = vec3 - Matrix::Vector3<Real>::floor(vec3);

    Index x = minVec[0];
    Index y = minVec[1];
    Index z = minVec[2];

    closestPoint(x, y, z, normal, &cd, &gridInside, &gridInNarrowBand);
    // closestPoint(x, y, z + 1, normals[1], &cds[1], &gridInside[1], &gridInNarrowBand[1]);
    // closestPoint(x, y + 1, z, normals[2], &cds[2], &gridInside[2], &gridInNarrowBand[2]);
    // closestPoint(x, y + 1, z + 1, normals[3], &cds[3], &gridInside[3], &gridInNarrowBand[3]);
    // closestPoint(x + 1, y, z, normals[4], &cds[4], &gridInside[4], &gridInNarrowBand[4]);
    // closestPoint(x + 1, y, z + 1, normals[5], &cds[5], &gridInside[5], &gridInNarrowBand[5]);
    // closestPoint(x + 1, y + 1, z, normals[6], &cds[6], &gridInside[6], &gridInNarrowBand[6]);
    // closestPoint(x + 1, y + 1, z + 1, normals[7], &cds[7], &gridInside[7], &gridInNarrowBand[7]);

    // int existInNarrowBand = 0;
    // int outsideCnt = 0;
    // for(int i = 0; i < 8; i++){
    //     existInNarrowBand += gridInNarrowBand[i] ? 1 : 0;
    //     outsideCnt += gridInside[i] ? 0 : 1;
    // }
    if(!gridInNarrowBand){
    	if(gridInside){
    		// Core::throwDefaultException(std::string("lookup point inside the mesh but not captured by the narrow band!!!"), __FILE__, __LINE__);
    		dist = insideConstant;
    		return false;
    	}else{
    		dist = outsideConstant;
        	return false;
    	}
    }
    
    dist = cd;
    // Math::TrilinearInterpolator<Real, Real> distInterpolator = Math::TrilinearInterpolator<Real, Real>();
    // dist = distInterpolator.interp(cds, vec3, gridInNarrowBand);

    // Math::TrilinearInterpolator<Matrix::Vector3<Real>, Real> normalInterpolator = Math::TrilinearInterpolator<Matrix::Vector3<Real>, Real>();
    // normal = normalInterpolator.interp(normals, vec3, gridInNarrowBand);
    // normal.normalize();
    // if(normal.length() == 0){
    //     for(int i = 0; i < 8; i++){
    //         if(gridInNarrowBand[i]){
    //             normal = normals[i];
    //             break;
    //         }
    //     }
    // }

    return true;
}

///////////////////////////////////////////////////////////////////////////////////
// NEIGHBOR SEARCH METHODS
///////////////////////////////////////////////////////////////////////////////////
template<class Traits>
bool DTGrid<Traits>::getLocatorAndValue(Index x, Index y, Index z, Locator& loc, Data& value) const
{
	// The Locator is defined as follows
    // struct Locator
    // {
    //     UInt iv1D, iv2D, iv3D;
    //     UInt ic1D, ic2D, ic3D;
    //     Index x, y, z;
    // };

	int fi, li, tmp, ix, iy, iz;

    // determine if 'p' exists in the narrow band
    fi = 0;
    li = lastXIndex;

    if(findIndex(x, xIndex, fi, li, &ix)){
    	loc.ic1D = ix;
    	loc.x = x;
    	loc.iv1D = tmp = aa1D[(ix)>>1] + x - xIndex[ix];
    	fi = va1D[tmp];
    	li = va1D[tmp+1]-1;
    	if(findIndex(y, yIndex, fi, li, &iy)){
    		loc.ic2D = iy;
	    	loc.y = y;
	    	loc.iv2D = tmp = aa2D[(iy)>>1] + y - yIndex[iy];
	    	fi = va2D[tmp];
	    	li = va2D[tmp+1]-1;
	    	if(findIndex(z, zIndex, fi, li, &iz)){
				// at this point we know the grid point (x,y,z) exists
				loc.z = z;
				loc.ic3D = iz;
				loc.iv3D = aa3D[(iz)>>1] + z - zIndex[iz];
				value = va3D[loc.iv3D];
				return true;
			}else{
				value = va3D[aa3D[(iz)>>1]] < 0 ? insideConstant : outsideConstant;
				return false;
			}
    	}else{
    		value = outsideConstant;
    		return false;
    	}
    }else{
    	value = outsideConstant;
    	return false;
    }

}

template<class Traits>
void DTGrid<Traits>::getVoxelValues(Index x, Index y, Index z, Data values[8]) const
{
    values[0] = (*this)(x,y,z);
    values[1] = (*this)(x,y,z+1);
    values[2] = (*this)(x,y+1,z);
    values[3] = (*this)(x,y+1,z+1);
    values[4] = (*this)(x+1,y,z);
    values[5] = (*this)(x+1,y,z+1);
    values[6] = (*this)(x+1,y+1,z);
    values[7] = (*this)(x+1,y+1,z+1);
}



// get neighbor (x-1,y,z)
template<class Traits>
bool DTGrid<Traits>::getXM(const Locator& loc, Locator& nbLoc) const
{
    int fi, li, iy, iz, tmp;

    if ( loc.x-1 < xIndex[loc.ic1D] )
    {
	return false;
    }

    // now we know that p-column x-1 exists
    // and we need to check if the grid point (x-1,y,z) exists...
    fi = va1D[loc.iv1D-1];
    li = va1D[loc.iv1D]-1;

    FIND_INDEX(loc.y, yIndex, fi, li, iy, false);

    tmp = aa2D[(iy)>>1] + loc.y - yIndex[iy];
    fi = va2D[tmp];
    li = va2D[tmp+1]-1;

    FIND_INDEX(loc.z, zIndex, fi, li, iz, false);

    nbLoc.iv1D = loc.iv1D-1;
    nbLoc.iv2D = tmp;
    nbLoc.iv3D = aa3D[(iz)>>1] + loc.z - zIndex[iz];
    nbLoc.ic1D = loc.ic1D;
    nbLoc.ic2D = iy;
    nbLoc.ic3D = iz;
    nbLoc.x = loc.x-1;
    nbLoc.y = loc.y;
    nbLoc.z = loc.z;

    return true;
}

// get neighbor (x+1,y,z)
template<class Traits>
bool DTGrid<Traits>::getXP(const Locator& loc, Locator& nbLoc) const
{
    int fi, li, iy, iz, tmp;

    if ( loc.x+1 > xIndex[loc.ic1D+1] )
    {
	return false;
    }

    // now we know that p-column x+1 exists
    // and we need to check if the grid point (x+1,y,z) exists...
    fi = va1D[loc.iv1D+1];
    li = va1D[loc.iv1D+2]-1;

    FIND_INDEX(loc.y, yIndex, fi, li, iy, false);

    tmp = aa2D[(iy)>>1] + loc.y - yIndex[iy];
    fi = va2D[tmp];
    li = va2D[tmp+1]-1;

    FIND_INDEX(loc.z, zIndex, fi, li, iz, false);

    nbLoc.iv1D = loc.iv1D+1;
    nbLoc.iv2D = tmp;
    nbLoc.iv3D = aa3D[(iz)>>1] + loc.z - zIndex[iz];
    nbLoc.ic1D = loc.ic1D;
    nbLoc.ic2D = iy;
    nbLoc.ic3D = iz;
    nbLoc.x = loc.x+1;
    nbLoc.y = loc.y;
    nbLoc.z = loc.z;

    return true;
}


// get neighbor (x, y-1,z)
template<class Traits>
bool DTGrid<Traits>::getYM(const Locator& loc, Locator& nbLoc) const
{
    int fi, li, iz;

    if ( loc.y-1 < yIndex[loc.ic2D] )
    {
	return false;
    }
    // now we know that p-column y-1 exists
    // check if the grid point (x,y-1,z) exists
    fi = va2D[loc.iv2D-1];
    li = va2D[loc.iv2D]-1;
    FIND_INDEX(loc.z, zIndex, fi, li, iz, false);

    nbLoc.iv1D = loc.iv1D;
    nbLoc.iv2D = loc.iv2D-1;
    nbLoc.iv3D = aa3D[(iz)>>1] + loc.z - zIndex[iz];
    nbLoc.ic1D = loc.ic1D;
    nbLoc.ic2D = loc.ic2D;
    nbLoc.ic3D = iz;
    nbLoc.x = loc.x;
    nbLoc.y = loc.y-1;
    nbLoc.z = loc.z;

    return true;
}

// get neighbor (x, y+1,z)
template<class Traits>
bool DTGrid<Traits>::getYP(const Locator& loc, Locator& nbLoc) const
{
    int fi, li, iz;

    if ( loc.y+1 > yIndex[loc.ic2D+1] )
    {
	return false;
    }

    // now we know that p-column y+1 exists
    // check if the grid point (x,y+1,z) exists
    fi = va2D[loc.iv2D+1];
    li = va2D[loc.iv2D+2]-1;
    FIND_INDEX(loc.z, zIndex, fi, li, iz, false);

    nbLoc.iv1D = loc.iv1D;
    nbLoc.iv2D = loc.iv2D+1;
    nbLoc.iv3D = aa3D[(iz)>>1] + loc.z - zIndex[iz];
    nbLoc.ic1D = loc.ic1D;
    nbLoc.ic2D = loc.ic2D;
    nbLoc.ic3D = iz;
    nbLoc.x = loc.x;
    nbLoc.y = loc.y+1;
    nbLoc.z = loc.z;

    return true;
}

// get neighbor (x, y, z-1)
template<class Traits>
bool DTGrid<Traits>::getZM(const Locator& loc, Locator& nbLoc) const
{
    if ( loc.z-1 >= zIndex[loc.ic3D] )
    {
	nbLoc.iv1D = loc.iv1D;
	nbLoc.iv2D = loc.iv2D;
	nbLoc.iv3D = loc.iv3D-1;
	nbLoc.ic1D = loc.ic1D;
	nbLoc.ic2D = loc.ic2D;
	nbLoc.ic3D = loc.ic3D;
	nbLoc.x = loc.x;
	nbLoc.y = loc.y;
	nbLoc.z = loc.z-1;

	return true;
    }
    return false;
}

// get neighbor (x, y, z+1)
template<class Traits>
bool DTGrid<Traits>::getZP(const Locator& loc, Locator& nbLoc) const
{
    if ( loc.z+1 <= zIndex[loc.ic3D+1] )
    {
	nbLoc.iv1D = loc.iv1D;
	nbLoc.iv2D = loc.iv2D;
	nbLoc.iv3D = loc.iv3D+1;
	nbLoc.ic1D = loc.ic1D;
	nbLoc.ic2D = loc.ic2D;
	nbLoc.ic3D = loc.ic3D;
	nbLoc.x = loc.x;
	nbLoc.y = loc.y;
	nbLoc.z = loc.z+1;

	return true;
    }
    return false;
}




template<class Traits>
void DTGrid<Traits>::getEXM(const Locator& loc, Locator& nbLoc) const
{
    int fi, li, iy, iz, tmp;

    fi = va1D[loc.iv1D-1];
    li = va1D[loc.iv1D]-1;

    FIND_INDEX_EXISTS(loc.y, yIndex, fi, li, iy);

    tmp = aa2D[(iy)>>1] + loc.y - yIndex[iy];
    fi = va2D[tmp];
    li = va2D[tmp+1]-1;

    FIND_INDEX_EXISTS(loc.z, zIndex, fi, li, iz);

    nbLoc.iv1D = loc.iv1D-1;
    nbLoc.iv2D = tmp;
    nbLoc.iv3D = aa3D[(iz)>>1] + loc.z - zIndex[iz];
    nbLoc.ic1D = loc.ic1D;
    nbLoc.ic2D = iy;
    nbLoc.ic3D = iz;
    nbLoc.x = loc.x-1;
    nbLoc.y = loc.y;
    nbLoc.z = loc.z;
}

template<class Traits>
void DTGrid<Traits>::getEXP(const Locator& loc, Locator& nbLoc) const
{
    int fi, li, iy, iz, tmp;

    // now we know that p-column x-1 exists
    // and we need to check if the grid point (x-1,y,z) exists...
    fi = va1D[loc.iv1D+1];
    li = va1D[loc.iv1D+2]-1;

    FIND_INDEX_EXISTS(loc.y, yIndex, fi, li, iy);

    tmp = aa2D[(iy)>>1] + loc.y - yIndex[iy];
    fi = va2D[tmp];
    li = va2D[tmp+1]-1;

    FIND_INDEX_EXISTS(loc.z, zIndex, fi, li, iz);

    nbLoc.iv1D = loc.iv1D+1;
    nbLoc.iv2D = tmp;
    nbLoc.iv3D = aa3D[(iz)>>1] + loc.z - zIndex[iz];
    nbLoc.ic1D = loc.ic1D;
    nbLoc.ic2D = iy;
    nbLoc.ic3D = iz;
    nbLoc.x = loc.x+1;
    nbLoc.y = loc.y;
    nbLoc.z = loc.z;
}

template<class Traits>
void DTGrid<Traits>::getEYM(const Locator& loc, Locator& nbLoc) const
{
    int fi, li, iz;

    fi = va2D[loc.iv2D-1];
    li = va2D[loc.iv2D]-1;

    FIND_INDEX_EXISTS(loc.z, zIndex, fi, li, iz);

    nbLoc.iv1D = loc.iv1D;
    nbLoc.iv2D = loc.iv2D-1;
    nbLoc.iv3D = aa3D[(iz)>>1] + loc.z - zIndex[iz];
    nbLoc.ic1D = loc.ic1D;
    nbLoc.ic2D = loc.ic2D;
    nbLoc.ic3D = iz;
    nbLoc.x = loc.x;
    nbLoc.y = loc.y-1;
    nbLoc.z = loc.z;
}

template<class Traits>
void DTGrid<Traits>::getEYP(const Locator& loc, Locator& nbLoc) const
{
    int fi, li, iz;

    fi = va2D[loc.iv2D+1];
    li = va2D[loc.iv2D+2]-1;

    FIND_INDEX_EXISTS(loc.z, zIndex, fi, li, iz);

    nbLoc.iv1D = loc.iv1D;
    nbLoc.iv2D = loc.iv2D+1;
    nbLoc.iv3D = aa3D[(iz)>>1] + loc.z - zIndex[iz];
    nbLoc.ic1D = loc.ic1D;
    nbLoc.ic2D = loc.ic2D;
    nbLoc.ic3D = iz;
    nbLoc.x = loc.x;
    nbLoc.y = loc.y+1;
    nbLoc.z = loc.z;
}

template<class Traits>
void DTGrid<Traits>::getEZM(const Locator& loc, Locator& nbLoc) const
{
    nbLoc.iv1D = loc.iv1D;
    nbLoc.iv2D = loc.iv2D;
    nbLoc.iv3D = loc.iv3D-1;
    nbLoc.ic1D = loc.ic1D;
    nbLoc.ic2D = loc.ic2D;
    nbLoc.ic3D = loc.ic3D;
    nbLoc.x = loc.x;
    nbLoc.y = loc.y;
    nbLoc.z = loc.z-1;
}

template<class Traits>
void DTGrid<Traits>::getEZP(const Locator& loc, Locator& nbLoc) const
{
    nbLoc.iv1D = loc.iv1D;
    nbLoc.iv2D = loc.iv2D;
    nbLoc.iv3D = loc.iv3D+1;
    nbLoc.ic1D = loc.ic1D;
    nbLoc.ic2D = loc.ic2D;
    nbLoc.ic3D = loc.ic3D;
    nbLoc.x = loc.x;
    nbLoc.y = loc.y;
    nbLoc.z = loc.z+1;
}



template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getXMArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    int fi, li, iy, iz, tmp;

    if ( loc.x-1 < xIndex[loc.ic1D] )
    {
	return 0;
    }

    // now we know that p-column x-1 exists
    // and we need to check if the grid point (x-1,y,z) exists...
    fi = va1D[loc.iv1D-1];
    li = va1D[loc.iv1D]-1;

    FIND_INDEX(loc.y, yIndex, fi, li, iy, 0);

    tmp = aa2D[(iy)>>1] + loc.y - yIndex[iy];
    fi = va2D[tmp];
    li = va2D[tmp+1]-1;

    FIND_INDEX(loc.z, zIndex, fi, li, iz, 0);

    *arrayIndex = aa3D[(iz)>>1] + loc.z - zIndex[iz];

    return 1U;
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getXPArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    int fi, li, iy, iz, tmp;

    if ( loc.x+1 > xIndex[loc.ic1D+1] )
    {
	return 0;
    }

    // now we know that p-column x-1 exists
    // and we need to check if the grid point (x-1,y,z) exists...
    fi = va1D[loc.iv1D+1];
    li = va1D[loc.iv1D+2]-1;

    FIND_INDEX(loc.y, yIndex, fi, li, iy, 0);

    tmp = aa2D[(iy)>>1] + loc.y - yIndex[iy];
    fi = va2D[tmp];
    li = va2D[tmp+1]-1;

    FIND_INDEX(loc.z, zIndex, fi, li, iz, 0);

    *arrayIndex = aa3D[(iz)>>1] + loc.z - zIndex[iz];

    return 1U;
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getYMArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    int fi, li, iz;

    if ( loc.y-1 < yIndex[loc.ic2D] )
    {
	return 0;
    }
    // now we know that p-column y-1 exists
    // check if the grid point (x,y-1,z) exists
    fi = va2D[loc.iv2D-1];
    li = va2D[loc.iv2D]-1;

    FIND_INDEX(loc.z, zIndex, fi, li, iz, 0);

    *arrayIndex = aa3D[(iz)>>1] + loc.z - zIndex[iz];

    return 1U;
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getYPArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    int fi, li, iz;

    if ( loc.y+1 > yIndex[loc.ic2D+1] )
    {
	return 0;
    }

    // now we know that p-column y+1 exists
    // check if the grid point (x,y+1,z) exists
    fi = va2D[loc.iv2D+1];
    li = va2D[loc.iv2D+2]-1;
    FIND_INDEX(loc.z, zIndex, fi, li, iz, 0);

    *arrayIndex = aa3D[(iz)>>1] + loc.z - zIndex[iz];

    return 1U;
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getZMArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    if ( loc.z-1 >= zIndex[loc.ic3D] )
    {
	*arrayIndex = loc.iv3D-1;
	return 1U;
    }
    return 0;
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getZPArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    if ( loc.z+1 <= zIndex[loc.ic3D+1] )
    {
	*arrayIndex = loc.iv3D+1;
	return 1U;
    }
    return 0;
}




template<class Traits>
void DTGrid<Traits>::getEXMArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    int fi, li, iy, iz, tmp;

    fi = va1D[loc.iv1D-1];
    li = va1D[loc.iv1D]-1;

    FIND_INDEX_EXISTS(loc.y, yIndex, fi, li, iy);

    tmp = aa2D[(iy)>>1] + loc.y - yIndex[iy];
    fi = va2D[tmp];
    li = va2D[tmp+1]-1;

    FIND_INDEX_EXISTS(loc.z, zIndex, fi, li, iz);

    *arrayIndex = aa3D[(iz)>>1] + loc.z - zIndex[iz];
}

template<class Traits>
void DTGrid<Traits>::getEXPArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    int fi, li, iy, iz, tmp;

    fi = va1D[loc.iv1D+1];
    li = va1D[loc.iv1D+2]-1;

    FIND_INDEX_EXISTS(loc.y, yIndex, fi, li, iy);

    tmp = aa2D[(iy)>>1] + loc.y - yIndex[iy];
    fi = va2D[tmp];
    li = va2D[tmp+1]-1;

    FIND_INDEX_EXISTS(loc.z, zIndex, fi, li, iz);

    *arrayIndex = aa3D[(iz)>>1] + loc.z - zIndex[iz];
}

template<class Traits>
void DTGrid<Traits>::getEYMArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    int fi, li, iz;

    fi = va2D[loc.iv2D-1];
    li = va2D[loc.iv2D]-1;

    FIND_INDEX_EXISTS(loc.z, zIndex, fi, li, iz);

    *arrayIndex = aa3D[(iz)>>1] + loc.z - zIndex[iz];
}

template<class Traits>
void DTGrid<Traits>::getEYPArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    int fi, li, iz;

    fi = va2D[loc.iv2D+1];
    li = va2D[loc.iv2D+2]-1;

    FIND_INDEX_EXISTS(loc.z, zIndex, fi, li, iz);

    *arrayIndex = aa3D[(iz)>>1] + loc.z - zIndex[iz];
}

template<class Traits>
void DTGrid<Traits>::getEZMArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    *arrayIndex = loc.iv3D-1;
}

template<class Traits>
void DTGrid<Traits>::getEZPArrayIndex(const Locator& loc, UInt *arrayIndex) const
{
    *arrayIndex = loc.iv3D+1;
}





/*! \brief getLocator computes the Locator of the grid point at position 'i' in the 3D value array.
* i must be less than numVa3D */
template<class Traits>
void DTGrid<Traits>::getLocator(UInt i, Locator *loc) const
{
    // 3D
    loc->iv3D = i;
    loc->ic3D = ((Math::binarySearch(aa3D, 0, numZIndex>>1, i)));
    // special case when only storing aa indices for starts of connected components
    if ( i<aa3D[loc->ic3D] )
    {
	loc->ic3D -= 1;
    }
    loc->ic3D = (loc->ic3D<<1);
    loc->z = zIndex[loc->ic3D] + (Index)( i - aa3D[(loc->ic3D)>>1] );
    loc->iv2D = Math::binarySearch(va2D, 0, numVa2D, loc->ic3D);
    // the point returned by the binary search is larger! We need it to be lower!
    if (va2D[loc->iv2D]!=loc->ic3D) { loc->iv2D--; }

    // 2D
    loc->ic2D = ((Math::binarySearch(aa2D, 0, numYIndex>>1, loc->iv2D)));
    // special case when only storing aa indices for starts of connected components
    if ( loc->iv2D<aa2D[loc->ic2D] )
    {
	loc->ic2D -= 1;
    }
    loc->ic2D = (loc->ic2D<<1);
    loc->y = yIndex[loc->ic2D] + (Index)( loc->iv2D - aa2D[(loc->ic2D)>>1] );
    loc->iv1D = Math::binarySearch(va1D, 0, numVa1D, loc->ic2D);
    // the point returned by the binary search is larger! We need it to be lower!
    if (va1D[loc->iv1D]!=loc->ic2D) { loc->iv1D--; }

    // 1D
    loc->ic1D = ((Math::binarySearch(aa1D, 0, numXIndex>>1, loc->iv1D)>>1)<<1);
    // special case when only storing aa indices for starts of connected components
    if ( loc->iv1D<aa1D[loc->ic1D] )
    {
	loc->ic1D -= 1;
    }
    loc->ic1D = (loc->ic1D<<1);
    loc->x = xIndex[loc->ic1D] + (Index)( loc->iv1D - aa1D[(loc->ic1D)>>1] );
}

template<class Traits>
bool DTGrid<Traits>::getLocator(const Matrix::Vector3<Real>& pos, Locator *loc) const
{
	Matrix::Vector3<Real> vec3 = worldToGrid * pos;
	Matrix::Vector3<Index> minVec = Matrix::Vector3<Real>::floor(vec3);
	return getLocator(minVec[0], minVec[1], minVec[2], loc);
}


// getLocator returns the Locator of position (x,y,z)
template<class Traits>
bool DTGrid<Traits>::getLocator(Index x, Index y, Index z, Locator *loc) const
{
    // The Locator is defined as follows
    // struct Locator
    // {
    //     UInt iv1D, iv2D, iv3D;
    //     UInt ic1D, ic2D, ic3D;
    //     Index x, y, z;
    // };

    int fi, li, tmp, ix, iy, iz;

    // determine if 'p' exists in the narrow band
    fi = 0;
    li = lastXIndex;

    FIND_INDEX(x, xIndex, fi, li, ix, false);

    loc->ic1D = ix;
    loc->x = x;
    loc->iv1D = tmp = aa1D[(ix)>>1] + x - xIndex[ix];
    fi = va1D[tmp];
    li = va1D[tmp+1]-1;

    FIND_INDEX(y, yIndex, fi, li, iy, false);

    loc->ic2D = iy;
    loc->y = y;
    loc->iv2D = tmp = aa2D[(iy)>>1] + y - yIndex[iy];
    fi = va2D[tmp];
    li = va2D[tmp+1]-1;

    FIND_INDEX(z, zIndex, fi, li, iz, false);

    // at this point we know the grid point (x,y,z) exists
    loc->z = z;
    loc->ic3D = iz;
    loc->iv3D = aa3D[(iz)>>1] + z - zIndex[iz];

    return true;
}



/*! getLocator returns the Locator of position (x,y,z) if it exists or the Locator of the lexicographically smallest element
*  that is lexicographically larger than (x,y,z) */
template<class Traits>
bool DTGrid<Traits>::getExistingLocator(Index x, Index y, Index z, Locator *loc, bool *inCorrectXYColumn) const
{
    // The Locator is defined as follows
    // struct Locator
    // {
    //     UInt iv1D, iv2D, iv3D;
    //     UInt ic1D, ic2D, ic3D;
    //     Index x, y, z;
    // };


    int fi, li, tmp, ix, iy, iz;

    // determine if 'p' exists in the narrow band
    fi = 0;
    li = lastXIndex;
    if (findIndex(x, xIndex, fi, li, &ix))
    {
	loc->ic1D = ix;
	loc->x = x;
	loc->iv1D = tmp = aa1D[(ix)>>1] + x - xIndex[ix];
	fi = va1D[tmp];
	li = va1D[tmp+1]-1;
	if (findIndex(y, yIndex, fi, li, &iy))
	{
	    loc->ic2D = iy;
	    loc->y = y;
	    loc->iv2D = tmp = aa2D[(iy)>>1] + y - yIndex[iy];
	    fi = va2D[tmp];
	    li = va2D[tmp+1]-1;

	    *inCorrectXYColumn = true;

	    if (findIndex(z, zIndex, fi, li, &iz))
	    {
		// at this point we know the grid point (x,y,z) exists
		loc->z = z;
		loc->ic3D = iz;
		loc->iv3D = aa3D[(iz)>>1] + z - zIndex[iz];

		return true;
	    }
	    else
	    {
		// here we do not need to check if we moved to a new x-p-column as this will be taken care of in
		// the following incrementUntil operation, but we do it anyway to maintain the invariant that the
		// grid point returned is the smallest grid point that is larger if the grid point we look for does not
		// exist
		loc->ic3D = iz;  
		loc->z = zIndex[loc->ic3D];
		loc->iv3D = aa3D[(loc->ic3D)>>1];

		if ( loc->ic3D == va2D[loc->iv2D+1] )
		{
		    *inCorrectXYColumn = false;
		    loc->iv2D++;

		    if ( loc->y == yIndex[loc->ic2D+1] )
		    {
			loc->ic2D += 2;
			loc->y = yIndex[loc->ic2D];

			if ( loc->ic2D == va1D[loc->iv1D+1] )
			{
			    loc->iv1D++;

			    if ( loc->x == xIndex[loc->ic1D+1] )
			    {
				loc->ic1D += 2;
				loc->x = xIndex[loc->ic1D];
			    }
			    else
			    {
				loc->x++;
			    }
			}
		    }
		    else
		    {
			loc->y++;
		    }
		}
	    }
	}
	else
	{
	    // here we do not need to check if we moved to a new x-p-column as this will be taken care of in
	    // the following incrementUntil operation, but we do it anyway to maintain the invariant that the
	    // grid point returned is the smallest grid point that is larger if the grid point we look for does not
	    // exist
	    loc->ic2D = iy;  
	    loc->y = yIndex[loc->ic2D];
	    loc->iv2D = aa2D[(loc->ic2D)>>1];

	    if ( loc->ic2D == va1D[loc->iv1D+1] )
	    {
		loc->iv1D++;

		if ( loc->x == xIndex[loc->ic1D+1] )
		{
		    loc->ic1D += 2;
		    loc->x = xIndex[loc->ic1D];
		}
		else
		{
		    loc->x++;
		}
	    }

	    loc->ic3D = va2D[loc->iv2D];
	    loc->z = zIndex[loc->ic3D];
	    loc->iv3D = aa3D[(loc->ic3D)>>1];
	    *inCorrectXYColumn = false;
	}
    }
    else
    {
	loc->ic1D = ix;
	loc->x = xIndex[loc->ic1D];
	loc->iv1D = aa1D[(loc->ic1D)>>1];

	loc->ic2D = va1D[loc->iv1D];  
	loc->y = yIndex[loc->ic2D];
	loc->iv2D = aa2D[(loc->ic2D)>>1];	    

	loc->ic3D = va2D[loc->iv2D];
	loc->z = zIndex[loc->ic3D];
	loc->iv3D = aa3D[(loc->ic3D)>>1];

	*inCorrectXYColumn = false;
    }


    return false;
}


template<class Traits>
bool DTGrid<Traits>::doesElementExist(Index x, Index y, Index z) const
{
    int fi, li, tmp, ix, iy, iz;

    // determine if '(x,y,z)' exists in the narrow band
    fi = 0;
    li = lastXIndex;
    if (findIndex(x, xIndex, fi, li, &ix))
    {
	tmp = aa1D[(ix)>>1] + x - xIndex[ix];
	fi = va1D[tmp];
	li = va1D[tmp+1]-1;
	if (findIndex(y, yIndex, fi, li, &iy))
	{
	    tmp = aa2D[(iy)>>1] + y - yIndex[iy];
	    fi = va2D[tmp];
	    li = va2D[tmp+1]-1;
	    if (findIndex(z, zIndex, fi, li, &iz))
	    {
		return true;
	    }
	}
    }

    return false;
}



///////////////////////////////////////////////////////////////////////////////////
// CLOSEST POINT HELPERS
///////////////////////////////////////////////////////////////////////////////////


template<class Traits>
void DTGrid<Traits>::computeGradient(Matrix::Vector3<Real>& n, Real v[6]) const
{
    n[0] = (v[1] - v[0]) / ( 2*dx );
    n[1] = (v[3] - v[2]) / ( 2*dx );
    n[2] = (v[5] - v[4]) / ( 2*dx );
    n.normalize();
}

template<class Traits>
void DTGrid<Traits>::computeGradient(const Locator& loc, Matrix::Vector3<Real>& n) const
{
	Real v[6];
	v[0] = getXM(loc.x, loc.y, loc.z, loc.ic1D);  // find neighbor to the left in the x direction
	v[1] = getXP(loc.x, loc.y, loc.z, loc.ic1D);  // find neighbor to the right in the x direction
	v[2] = getYM(loc.y, loc.z, loc.ic2D);     // find neighbor to the left in the y direction
	v[3] = getYP(loc.y, loc.z, loc.ic2D);     // find neighbor to the right in the y direction
	v[4] = getZM(loc.z, loc.ic3D);        // find neighbor to the left in the z direction
	v[5] = getZP(loc.z, loc.ic3D);        // find neighbor to the right in the z direction
	computeGradient(n, v);
}

// ix points to start of connected component
template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::getXM(Index x, Index y, Index z, int ix) const
{
    Index tmpI;
    int tmp;

    x--;
    // does the neighbor exist?
    if (x >= (tmpI=xIndex[ix]))
    {
	// construct index 'tmp' into 'va1D' array
	tmp = aa1D[(ix)>>1] + x - tmpI;
	return get(y, z, va1D[tmp], va1D[tmp+1]-1);
    }

    if (!Traits::openLevelSetsSupport)
    {
	return outsideConstant;
    }
    else
    {
	//OpenLS:
	UInt aa1DIndex = aa1D[ix>>1];
	UInt va1DValue = va1D[aa1DIndex];
	UInt aa2DIndex = aa2D[va1DValue>>1];
	UInt va2DValue = va2D[aa2DIndex];
	UInt aa3DIndex = aa3D[va2DValue>>1];
	return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
    }
}



// ix points to start of connected component
template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::getXP(Index x, Index y, Index z, int ix) const
{
    int tmp;

    x++;
    if (x <= xIndex[ix+1])
    {
	tmp = aa1D[(ix)>>1] + x - xIndex[ix];
	return get(y, z, va1D[tmp], va1D[tmp+1]-1);
    }

    if (!Traits::openLevelSetsSupport)
    {
	return outsideConstant;
    }
    else
    {
	//OpenLS:
	UInt aa1DIndex = aa1D[ix>>1] + xIndex[ix+1] - xIndex[ix];
	UInt va1DValue = va1D[aa1DIndex];
	UInt aa2DIndex = aa2D[va1DValue>>1];
	UInt va2DValue = va2D[aa2DIndex];
	UInt aa3DIndex = aa3D[va2DValue>>1];
	return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
    }
    return outsideConstant;
}


// iy points to start of connected component
template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::getYM(Index y, Index z, int iy) const
{
    Index tmpI;
    int tmp;

    y--;
    if (y >= (tmpI=yIndex[iy]))
    {
	tmp = aa2D[(iy)>>1] + y - tmpI;

	return get(z, va2D[tmp], va2D[tmp+1]-1);
    }

    if (!Traits::openLevelSetsSupport)
    {
	return outsideConstant;
    }
    else
    {
	//OpenLS:
	UInt aa2DIndex = aa2D[iy>>1];
	UInt va2DValue = va2D[aa2DIndex];
	UInt aa3DIndex = aa3D[va2DValue>>1];
	return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
    }
}



// iy points to start of connected component
template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::getYP(Index y, Index z, int iy) const
{
    int tmp;

    y++;
    if (y <= yIndex[iy+1])
    {
	tmp = aa2D[(iy)>>1] + y - yIndex[iy];

	return get(z, va2D[tmp], va2D[tmp+1]-1);
    }

    if (!Traits::openLevelSetsSupport)
    {
	return outsideConstant;
    }
    else
    {
	//OpenLS:
	UInt aa2DIndex = aa2D[iy>>1] + yIndex[iy+1] - yIndex[iy];
	UInt va2DValue = va2D[aa2DIndex];
	UInt aa3DIndex = aa3D[va2DValue>>1];
	return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
    }
}



// iz points to start of connected component
template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::getZM(Index z, int iz) const
{
    Index tmpI;

    z--;
    if (z >= (tmpI=zIndex[iz]))
    {
	return va3D[aa3D[(iz)>>1] + z - tmpI];
    }

    return va3D[aa3D[(iz)>>1]] < 0 ? insideConstant : outsideConstant;
}



// iy points to start of connected component
template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::getZP(Index z, int iz) const
{
    z++;
    if (z <= zIndex[iz+1])
    {
	return va3D[aa3D[(iz)>>1] + z - zIndex[iz]];
    }

    return va3D[aa3D[(iz+1)>>1]+zIndex[iz+1]-zIndex[iz]] < 0 ? insideConstant : outsideConstant;
}


template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::get(Index z, int fi, int li) const
{
    if (Traits::randomAccessType == USE_BINARY_SEARCH)
    {
	int iz, iz2, it;

	// inside min, max range of z indices in this column?
	if (z < zIndex[fi])
	{
	    return va3D[ aa3D[(fi)>>1] ] < 0 ? insideConstant : outsideConstant;
	}
	else if (z > zIndex[li])
	{
	    return va3D[ aa3D[(li)>>1] + (zIndex[li] - zIndex[li-1]) ] < 0 ? insideConstant : outsideConstant; 
	}
	else
	{
	    iz = 0;
	    iz2 = (li-fi)+1;

	    // binary search
	    while (iz2 != iz+2)
	    {
		it = ((iz+iz2)>>2)<<1;  // make even

		if (z < zIndex[it+fi])
		{
		    iz2 = it;
		}
		else
		{
		    iz = it;
		}
	    }

	    // now iz points to the correct closest even index
	    iz += fi;
	    if ( z <= zIndex[iz+1] )
	    {
		return va3D[ aa3D[(iz)>>1] + z - zIndex[iz] ];
	    }
	    else
	    {
		// note that the sign will always be taken from the same p-column
		return va3D[ aa3D[(iz+2)>>1] ] < 0 ? insideConstant : outsideConstant;
	    }
	}
    }
    else
    {
	int iz;
	Index tmpI;

	if ( z >= zIndex[fi] && z <= zIndex[li] )
	{
	    iz = fi+1;
	    while (z > zIndex[iz])
	    {
		iz += 2;
	    }
	    iz -= 1;
	    if (z >= (tmpI=zIndex[iz]))
	    {
		return va3D[ aa3D[(iz)>>1] + z - tmpI ];
	    }
	    else
	    {
		return va3D[ aa3D[(iz)>>1] ] < 0 ? insideConstant : outsideConstant;
	    }
	}
	else
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (z < zIndex[fi])
		{
		    UInt aa3DIndex = aa3D[fi>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa3DIndex = aa3D[li>>1] + (zIndex[li] - zIndex[li-1]);
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	}
    }
}




template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::get(Index y, Index z, int fi, int li) const
{
    if (Traits::randomAccessType == USE_BINARY_SEARCH)
    {
	int iy=0, tmp, it;
	int iy2;

	if ( y >= yIndex[fi] && y <= yIndex[li] )
	{
	    iy = 0;
	    iy2 = (li-fi)+1;

	    // binary search
	    while (iy2 != iy+2)
	    {
		it = ((iy+iy2)>>2)<<1;  // make even

		if (y < yIndex[it+fi])
		{
		    iy2 = it;
		}
		else
		{
		    iy = it;
		}
	    }

	    // now iy points to the correct closest even index
	    iy += fi;
	    if ( y <= yIndex[iy+1] )
	    {
		tmp = aa2D[(iy)>>1] + y - yIndex[iy];
		fi = va2D[tmp];
		li = va2D[tmp+1]-1;
		return get(z, fi, li);
	    }
	}
	else
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (y < yIndex[fi])
		{
		    UInt aa2DIndex = aa2D[fi>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa2DIndex = aa2D[li>>1] + (yIndex[li] - yIndex[li-1]);
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	}

	if (!Traits::openLevelSetsSupport)
	{
	    return outsideConstant;
	}
	else
	{
	    //OpenLS:
	    UInt aa2DIndex = aa2D[iy>>1] + (yIndex[iy+1] - yIndex[iy]);
	    UInt va2DValue = va2D[aa2DIndex];
	    UInt aa3DIndex = aa3D[va2DValue>>1];
	    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	}

	return outsideConstant;
    }
    else
    {
	int iy=0, tmp;
	Index tmpI;

	if ( y >= yIndex[fi] && y <= yIndex[li] )
	{
	    iy = fi+1;
	    while (y > yIndex[iy])
	    {
		iy += 2;
	    }
	    iy -= 1;
	    if (y >= (tmpI=yIndex[iy]))
	    {
		tmp = aa2D[(iy)>>1] + y - tmpI;
		fi = va2D[tmp];
		li = va2D[tmp+1]-1;

		return get(z, fi, li);
	    }
	}
	else
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (y < yIndex[fi])
		{
		    UInt aa2DIndex = aa2D[fi>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa2DIndex = aa2D[li>>1] + (yIndex[li] - yIndex[li-1]);
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	}

	if (!Traits::openLevelSetsSupport)
	{
	    return outsideConstant;
	}
	else
	{
	    //OpenLS:
	    UInt aa2DIndex = aa2D[iy>>1];
	    UInt va2DValue = va2D[aa2DIndex];
	    UInt aa3DIndex = aa3D[va2DValue>>1];
	    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	}

	return outsideConstant;
    }
}





template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::get(Index x, Index y, Index z, int fi, int li) const
{
    if (Traits::randomAccessType == USE_BINARY_SEARCH)
    {

	int ix=0, tmp, it;
	int ix2;


	if ( x >= xIndex[fi] && x <= xIndex[li] )
	{
	    ix = 0;
	    ix2 = (li-fi)+1;

	    // binary search
	    while (ix2 != ix+2)
	    {
		it = ((ix+ix2)>>2)<<1;  // make even

		if (x < xIndex[it+fi])
		{
		    ix2 = it;
		}
		else
		{
		    ix = it;
		}
	    }

	    // now iy points to the correct closest even index
	    ix += fi;
	    if ( x <= xIndex[ix+1] )
	    {
		tmp = aa1D[(ix)>>1] + x - xIndex[ix];
		fi = va1D[tmp];
		li = va1D[tmp+1]-1;
		return get(y, z, fi, li);
	    }
	}
	else
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (x < xIndex[fi])
		{
		    //OpenLS:
		    UInt aa1DIndex = aa1D[fi>>1];
		    UInt va1DValue = va1D[aa1DIndex];
		    UInt aa2DIndex = aa2D[va1DValue>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    //OpenLS:
		    UInt aa1DIndex = aa1D[li>>1] + (xIndex[li] - xIndex[li-1]);
		    UInt va1DValue = va1D[aa1DIndex];
		    UInt aa2DIndex = aa2D[va1DValue>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	}


	if (!Traits::openLevelSetsSupport)
	{
	    return outsideConstant;
	}
	else
	{
	    //OpenLS: x > xIndex[ix+1]
	    UInt aa1DIndex = aa1D[ix>>1] + (xIndex[ix+1] - xIndex[ix]);
	    UInt va1DValue = va1D[aa1DIndex];
	    UInt aa2DIndex = aa2D[va1DValue>>1];
	    UInt va2DValue = va2D[aa2DIndex];
	    UInt aa3DIndex = aa3D[va2DValue>>1];
	    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	}

    }
    else
    {
	int ix=0, tmp;
	Index tmpI;

	if ( x >= xIndex[fi] && x <= xIndex[li] )
	{
	    ix = fi+1;
	    while (x > xIndex[ix])
	    {
		ix += 2;
	    }
	    ix -= 1;
	    if (x >= (tmpI=xIndex[ix]))
	    {
		tmp = aa1D[(ix)>>1] + x - tmpI;
		fi = va1D[tmp];
		li = va1D[tmp+1]-1;

		return get(y, z, fi, li);
	    }
	}
	else
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (x < xIndex[fi])
		{
		    //OpenLS:
		    UInt aa1DIndex = aa1D[fi>>1];
		    UInt va1DValue = va1D[aa1DIndex];
		    UInt aa2DIndex = aa2D[va1DValue>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    //OpenLS:
		    UInt aa1DIndex = aa1D[li>>1] + (xIndex[li] - xIndex[li-1]);
		    UInt va1DValue = va1D[aa1DIndex];
		    UInt aa2DIndex = aa2D[va1DValue>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	}

	if (!Traits::openLevelSetsSupport)
	{
	    return outsideConstant;
	}
	else
	{
	    //OpenLS: x < xIndex[ix]
	    UInt aa1DIndex = aa1D[ix>>1];
	    UInt va1DValue = va1D[aa1DIndex];
	    UInt aa2DIndex = aa2D[va1DValue>>1];
	    UInt va2DValue = va2D[aa2DIndex];
	    UInt aa3DIndex = aa3D[va2DValue>>1];
	    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	}

    }
}




template<class Traits>
bool DTGrid<Traits>::findIndex(Index n, Index *nIndex, int fi, int li, int *in) const
{
    if (Traits::randomAccessType == USE_BINARY_SEARCH)
    {
	int in2, in1, it;

	if (n < nIndex[fi])
	{ 
	    *in = 0;
	    return false;
	}

	if (n > nIndex[li])
	{
	    *in = li+1;
	    return false;
	}


	in1 = 0;
	in2 = (li-fi)+1;

	// binary search
	// at all times, nIndex[in1] <= n < nIndex[in2]
	while (in2 != in1+2)
	{
	    it = ((in1+in2)>>2)<<1;   // make even

	    if (n < nIndex[it+fi])
	    {
		in2 = it;
	    }
	    else
	    {
		in1 = it;
	    }
	}

	// now in1 points to the closest even index, i.e. the start of a connected component

	*in = in1 = in1 + fi;

	if ( n <= nIndex[in1+1] )
	{
	    return true;
	}

	*in += 2;

	return false;
    }
    else
    {
	if (n < nIndex[fi])
	{ 
	    *in = 0;
	    return false;
	}

	if (n > nIndex[li])
	{
	    *in = li+1;
	    return false;
	}


	*in = fi+1;
	// now in points to the end of a connected component
	while (n > nIndex[*in])
	{ 
	    *in += 2;
	} 
	*in -= 1;
	// now in points to the start of a connected component
	if (n < nIndex[*in])
	{ 
	    return false;
	}
	return true;
    }
}





///////////////////////////////////////////////////////////////////////////////////

















///////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF _ITERATOR
///////////////////////////////////////////////////////////////////////////////////

template<class Traits> template<IteratorType it>
DTGrid<Traits>::_Iterator<it>::_Iterator(const DTGrid *parent, Data *va3D, UInt numVa3D, const Locator& loc, bool valid, bool inCorrectXYColumn)
: parent(parent), va3D(va3D), numVa3D(numVa3D), zIndex(parent==NULL?NULL:parent->zIndex)
{
    if (numVa3D > 0)
    {
	insideConstant = parent->insideConstant;
	outsideConstant = parent->outsideConstant;
	ic1D = loc.ic1D + 1;
	ic2D = loc.ic2D + 1;
	ic3D = loc.ic3D + 1;
	iv1D = loc.iv1D + 1;
	iv2D = loc.iv2D + 1;
	iv3D = loc.iv3D;
	x = loc.x;
	y = loc.y;
	z = loc.z;
	if (valid)
	{
	    value = va3D[iv3D];
	}
	else
	{
	    value = va3D[iv3D] < 0 ? insideConstant : outsideConstant; 
	}
	this->inCorrectPColumn = inCorrectXYColumn;
	this->valid = valid;
	zn = zIndex[ic3D];
    }
    else
    {
	iv3D = 0;
	this->valid = false;
    }
}

template<class Traits> template<IteratorType it>
DTGrid<Traits>::_Iterator<it>::_Iterator(const DTGrid *parent, Data *va3D, UInt numVa3D, bool begin)
: parent(parent), va3D(va3D), numVa3D(numVa3D), zIndex(parent==NULL?NULL:parent->zIndex)
{
    if (begin)
    {
	if (numVa3D > 0)
	{
	    insideConstant = parent->insideConstant;
	    outsideConstant = parent->outsideConstant;
	    ic1D = ic2D = ic3D = 1;
	    iv1D = iv2D = 1;
	    iv3D = 0;
	    x = parent->xIndex[0];
	    y = parent->yIndex[0];
	    z = zIndex[0];
	    value = va3D[iv3D];
	    valid = true;
	    inCorrectPColumn = true;
	    zn = zIndex[ic3D];  // used with IT_VOLUME
	}
	else
	{
	    iv3D = 0;
	    valid = false;
	}
    }
    else
    {
	iv3D = numVa3D;
	valid = false;
    }
}

template<class Traits> template<IteratorType it>
void DTGrid<Traits>::_Iterator<it>::getLocator(Locator& loc) const
{
    loc.ic1D = ic1D - 1;
    loc.ic2D = ic2D - 1;
    loc.ic3D = ic3D - 1;
    loc.iv1D = iv1D - 1;
    loc.iv2D = iv2D - 1;
    loc.iv3D = iv3D;
    loc.x = x;
    loc.y = y;
    loc.z = z;
}


template<class Traits> template<IteratorType it>
typename DTGrid<Traits>::Data DTGrid<Traits>::_Iterator<it>::operator*() const
{
    return value;
}


template<class Traits> template<IteratorType it>
bool DTGrid<Traits>::_Iterator<it>::operator==(const _Iterator& iter) const
{
    return iter.iv3D == iv3D;
}

template<class Traits> template<IteratorType it>
bool DTGrid<Traits>::_Iterator<it>::operator!=(const _Iterator& iter) const
{
    return iter.iv3D != iv3D;
}

template<class Traits> template<IteratorType it>
#ifdef WIN32
typename DTGrid<Traits>::_Iterator<it> &DTGrid<Traits>::_Iterator<it>::operator++()
#else
typename DTGrid<Traits>::template _Iterator<it> &DTGrid<Traits>::_Iterator<it>::operator++()
#endif
{
    // we assume that iv3D < numVa3D when this method is entered

    if (it == IT_TUBE)
    {

	iv3D++;

	// set value
	value = va3D[iv3D];

	// reached end of connected component in z direction?
	if (z == zIndex[ic3D])
	{
	    ic3D++;

	    // moved to a new (x,y) column?
	    if (parent->va2D[iv2D] == ic3D)
	    {
		iv2D++;

		// reached end of connected component in y direction?
		if (y == parent->yIndex[ic2D])                         // // VTUNE caution: LCP (Length Changing Prefix), but alternative is to use long instead of short at the cost of more memory usage
		{
		    ic2D++;

		    // moved to a new (x) column?
		    if (parent->va1D[iv1D] == ic2D)
		    {
			iv1D++;

			if (x == parent->xIndex[ic1D])
			{
			    ic1D++;
			    x = parent->xIndex[ic1D];
			    ic1D++;
			}
			else
			{
			    x++;
			}
		    }

		    y = parent->yIndex[ic2D];
		    ic2D++;

		}
		else
		{
		    y++;
		}

	    }

	    z = zIndex[ic3D];
	    ic3D++;
	}
	else
	{
	    z++;
	}

    }
    else if (it == IT_VOLUME)
    {


	// reached end of connected component in z direction?
	if (z == zn)
	{
	    ic3D++;

	    if (value < 0)
	    {
		if (ic3D % 2 != 0)
		{
		    z++;
		    iv3D++;
		    value = va3D[iv3D];
		    zn = zIndex[ic3D];
		}
		else
		{
		    z++;
		    value = insideConstant;
		    zn = zIndex[ic3D]-1;
		}
	    }
	    else
	    {
		// moved to a new (x,y) column?
		if (parent->va2D[iv2D] == ic3D)
		{
		    iv2D++;

		    // reached end of connected component in y direction?
		    if (y == parent->yIndex[ic2D])
		    {
			ic2D++;

			// moved to a new (x) column?
			if (parent->va1D[iv1D] == ic2D)
			{
			    iv1D++;

			    if (x == parent->xIndex[ic1D])
			    {
				ic1D++;
				x = parent->xIndex[ic1D];
				ic1D++;
			    }
			    else
			    {
				x++;
			    }
			}

			y = parent->yIndex[ic2D];
			ic2D++;

		    }
		    else
		    {
			y++;
		    }
		}

		z = zIndex[ic3D];
		ic3D++;
		zn = zIndex[ic3D];

		iv3D++;
		value = va3D[iv3D];

	    }
	}
	else
	{
	    z++;

	    // we should only increment iv3D and set value if we are inside a connected component
	    if ( ic3D % 2 != 0 )
	    {
		// if ic3D is odd we are inside a connected component
		iv3D++;
		value = va3D[iv3D];
	    }

	}

    }


    return *this;
}



template<class Traits> template<IteratorType it>
void DTGrid<Traits>::_Iterator<it>::incrementFast()
{
    iv3D++;
    value = va3D[iv3D];
    z++;
}

template<class Traits> template<IteratorType it>
void DTGrid<Traits>::_Iterator<it>::incrementFastZ(const Data& outsideValue, Index zn)
{
    // incrementFastZ must only increment the iterator one point forward (therefore, if "z == zIndex[ic3D]" as below and we are exiting a connected componen, 
    // the value should be set to +/- gamma since the iterator will under no circumstances end up in the correct grid point
    // note that it must be possible to continue iteration with incrementUntil following a sequence of incrementFastZ operations

    // first check if this iterator is valid in this (x,y)'th p-column of grid points
    // we have to make this check since a point in the stencil may pass beyond its last point before
    // the center stencil point does! (note that incrementFastZ is only called when the center stencil grid point moves one grid point forward, but
    // we are on the border of the narrow band.

    if ( inCorrectPColumn )   // valid, ie in correct p-column?
    {
	// yes, in correct p-column

	if (z < zn)   
	{
	    // yes, we need to increment

	    if (z == zIndex[ic3D])
	    {

		// we jump to next connected component, ie
		// this stencil iterator jumps more than one grid point, ie.
		// the corresponding grid point lies outside of the narrow band

		valid = false;

		ic3D++;  // next connected component
		z = zIndex[ic3D];

		// we have to make sure that the iterator stays in the same (x,y)'th p-column
		if (parent->va2D[iv2D] == ic3D)
		{
		    // we move to next p-column
		    // adjust 1D and 2D pointers

		    iv2D++;

		    // reached end of connected component in y direction?
		    if (y == parent->yIndex[ic2D])
		    {
			ic2D++;

			// moved to a new (x) column?
			if (parent->va1D[iv1D] == ic2D)
			{
			    iv1D++;

			    if (x == parent->xIndex[ic1D])
			    {
				ic1D++;
				x = parent->xIndex[ic1D];
				ic1D++;
			    }
			    else
			    {
				x++;  
			    }
			}

			y = parent->yIndex[ic2D];
			ic2D++;

		    }
		    else
		    {
			y++;
		    }

		    ic3D++;
		    iv3D++;
		    inCorrectPColumn = false;
		}
		else
		{
		    ic3D++;
		    iv3D++;
		}

		value = outsideValue;
	    }
	    else
	    {
		// we stay within the same connected component
		z++;
		iv3D++;
		// no need to set valid to true here, it will already be true
		value = va3D[iv3D];
	    }
	}
	else if (z==zn)
	{
	    valid = true;
	    value = va3D[iv3D];
	}
    }
}

template<class Traits> template<IteratorType it>
void DTGrid<Traits>::_Iterator<it>::incrementUntil(const Data& outsideValue, Index xn, Index yn, Index zn)
{
    if (x < xn)
    {
	do
	{
	    iv1D++;

	    if (x == parent->xIndex[ic1D])
	    {
		ic1D++;
		x = parent->xIndex[ic1D];
		ic1D++;
	    }
	    else
	    {
		x++;
	    }		
	}
	while (x < xn);

	ic2D = parent->va1D[iv1D-1];  // we know iv1D >= 1
	y = parent->yIndex[ic2D];
	iv2D = parent->aa2D[(ic2D)>>1];	    
	ic2D++;

	ic3D = parent->va2D[iv2D];
	iv2D++;
	z = zIndex[ic3D];
	iv3D = parent->aa3D[(ic3D)>>1];
	ic3D++;
    }
    if (xn==x)
    {

	if (y < yn)  
	{
	    do
	    {
		iv2D++;

		if (y == parent->yIndex[ic2D])
		{
		    ic2D++;

		    // moved to a new (x) column?
		    if (parent->va1D[iv1D] == ic2D)
		    {
			iv1D++;

			if (x == parent->xIndex[ic1D])
			{
			    ic1D++;
			    x = parent->xIndex[ic1D];
			    ic1D++;
			}
			else
			{
			    x++;
			}
		    }

		    y = parent->yIndex[ic2D];
		    ic2D++;
		}
		else
		{
		    y++;
		}		
	    }
	    while (x==xn && y<yn);

	    ic3D = parent->va2D[iv2D-1];
	    z = zIndex[ic3D];
	    iv3D = parent->aa3D[(ic3D)>>1];
	    ic3D++;
	}

	if (xn==x && yn==y)
	{

	    while (xn==x && yn==y && z < zn)  
	    {
		iv3D++;
		if (z == zIndex[ic3D])
		{
		    ic3D++;

		    // moved to a new (x,y) column?
		    if (parent->va2D[iv2D] == ic3D)
		    {
			iv2D++;

			// reached end of connected component in y direction?
			if (y == parent->yIndex[ic2D])
			{
			    ic2D++;

			    // moved to a new (x) column?
			    if (parent->va1D[iv1D] == ic2D)
			    {
				iv1D++;

				if (x == parent->xIndex[ic1D])
				{
				    ic1D++;
				    x = parent->xIndex[ic1D];
				    ic1D++;
				}
				else
				{
				    x++;
				}
			    }

			    y = parent->yIndex[ic2D];
			    ic2D++;

			}
			else
			{
			    y++;
			}

		    }

		    z = zIndex[ic3D];
		    ic3D++;
		}
		else
		{
		    z++;  
		}
	    }

	    if (xn==x && yn==y)
	    {
		inCorrectPColumn = true;

		if (zn==z)
		{
		    value = va3D[iv3D];
		    valid = true;
		}
		else
		{
		    value = outsideValue;
		    valid = false;
		}
	    }
	    else
	    {
		value = outsideValue;
		inCorrectPColumn = false;
		valid = false;
	    }
	}
	else
	{
	    value = outsideValue;
	    inCorrectPColumn = false;
	    valid = false;
	}
    }
    else
    {
	value = outsideValue;
	inCorrectPColumn = false;
	valid = false;
    }
}



template<class Traits> template<IteratorType it>
bool DTGrid<Traits>::_Iterator<it>::hasNext() const
{
    return iv3D != parent->numVa3D;
}




template<class Traits> template<IteratorType it>
void DTGrid<Traits>::_Iterator<it>::getIndex(Index *x, Index *y, Index *z) const
{
    *x = this->x;
    *y = this->y;
    *z = this->z;
}



template<class Traits> template<IteratorType it> 
#ifdef WIN32
typename DTGrid<Traits>::_Iterator<it> DTGrid<Traits>::begin(UInt id) const
#else
typename DTGrid<Traits>::template _Iterator<it> DTGrid<Traits>::begin(UInt id) const
#endif
{
    return _Iterator<it>(this, &values[id], values[id].getNumValues(), true);
}



template<class Traits> template<IteratorType it>
#ifdef WIN32
typename DTGrid<Traits>::_Iterator<it> DTGrid<Traits>::end(UInt id) const
#else
typename DTGrid<Traits>::template _Iterator<it> DTGrid<Traits>::end(UInt id) const
#endif
{
    return _Iterator<it>(this, &values[id], values[id].getNumValues(), false);
}


template<class Traits>
typename DTGrid<Traits>::TubeIterator DTGrid<Traits>::beginTubeIterator(const Locator& loc, bool inCorrectXYColumn, bool valid, UInt id) const
{
    return TubeIterator(this, values[id], numValues[id], loc, valid, inCorrectXYColumn);
}


template<class Traits>
typename DTGrid<Traits>::TubeIterator DTGrid<Traits>::beginTubeIterator(UInt id) const
{
    return TubeIterator(this, values[id], numValues[id], true);
}

template<class Traits>
typename DTGrid<Traits>::TubeIterator DTGrid<Traits>::endTubeIterator(UInt id) const
{
    return TubeIterator(this, values[id], numValues[id], false);
}



template<class Traits>
typename DTGrid<Traits>::TubeTopologyIterator DTGrid<Traits>::beginTubeTopologyIterator() const
{
    return TubeTopologyIterator(this, true);
}

template<class Traits>
typename DTGrid<Traits>::TubeTopologyIterator DTGrid<Traits>::endTubeTopologyIterator() const
{
    return TubeTopologyIterator(this, false);
}


template<class Traits>
typename DTGrid<Traits>::VolumeIterator DTGrid<Traits>::beginVolumeIterator(UInt id) const
{
    return VolumeIterator(this, values[id], numValues[id], true);
}

template<class Traits>
typename DTGrid<Traits>::VolumeIterator DTGrid<Traits>::endVolumeIterator(UInt id) const
{
    return VolumeIterator(this, values[id], numValues[id], false);
}


///////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF TRANSFORMATIONITERATOR
///////////////////////////////////////////////////////////////////////////////////


template<class Traits>template <class TemplateInterpolator> 
#ifdef WIN32
typename DTGrid<Traits>::TransformationIterator<TemplateInterpolator> DTGrid<Traits>::beginTransformationIterator(Index bbox[3][2], const Matrix::Matrix4x4<Real>& trans, TemplateInterpolator *interpolator, Data outsideConstant) const
#else
typename DTGrid<Traits>::template TransformationIterator<TemplateInterpolator> DTGrid<Traits>::beginTransformationIterator(Index bbox[3][2], const Matrix::Matrix4x4<Real>& trans, TemplateInterpolator *interpolator, Data outsideConstant) const
#endif
{
    return TransformationIterator<TemplateInterpolator>(this, bbox, trans, interpolator, outsideConstant);	
}

template<class Traits> template <class TemplateInterpolator> 
#ifdef WIN32
typename DTGrid<Traits>::TransformationIterator<TemplateInterpolator> DTGrid<Traits>::endTransformationIterator(Index bbox[3][2]) const
#else
typename DTGrid<Traits>::template TransformationIterator<TemplateInterpolator> DTGrid<Traits>::endTransformationIterator(Index bbox[3][2]) const
#endif
{
    return TransformationIterator<TemplateInterpolator>(bbox);
}




template<class Traits>template <class TemplateInterpolator> 
DTGrid<Traits>::TransformationIterator<TemplateInterpolator>::TransformationIterator(const DTGrid *grid2, Index bbox[3][2], const Matrix::Matrix4x4<Real>& trans, TemplateInterpolator *interpolator, Data outsideConstant)
: grid2(grid2)
{
    // begin iterator
    gridToGrid2 = trans;
    this->interpolator = interpolator;
    this->outsideConstant = outsideConstant;
    i = bbox[0][0];
    j = bbox[1][0];
    k = bbox[2][0]-1;
    this->bbox[0][0] = bbox[0][0];
    this->bbox[0][1] = bbox[0][1];
    this->bbox[1][0] = bbox[1][0];
    this->bbox[1][1] = bbox[1][1];
    this->bbox[2][0] = bbox[2][0];
    this->bbox[2][1] = bbox[2][1];
    operator++();
}


template<class Traits>template <class TemplateInterpolator> 
DTGrid<Traits>::TransformationIterator<TemplateInterpolator>::TransformationIterator(Index bbox[3][2])
{
    // end iterator
    i = bbox[0][1]+1;
    j = bbox[1][0];
    k = bbox[2][0];	
}



template<class Traits>template <class TemplateInterpolator> 
typename DTGrid<Traits>::Data DTGrid<Traits>::TransformationIterator<TemplateInterpolator>::getValue() const
{
    return value;
}

template<class Traits>template <class TemplateInterpolator> 
bool DTGrid<Traits>::TransformationIterator<TemplateInterpolator>::operator!=(const TransformationIterator& oti) const
{
    return oti.i != i || oti.j != j || oti.k != k;
}

template<class Traits>template <class TemplateInterpolator> 
bool DTGrid<Traits>::TransformationIterator<TemplateInterpolator>::operator==(const TransformationIterator& oti) const
{
    return oti.i == i && oti.j == j && oti.k == k;
}

template<class Traits>template <class TemplateInterpolator> 
void DTGrid<Traits>::TransformationIterator<TemplateInterpolator>::getIndex(Index *i, Index *j, Index *k) const
{
    *i = this->i;
    *j = this->j;
    *k = this->k;
}


template<class Traits>template <class TemplateInterpolator> 
#ifdef WIN32
typename DTGrid<Traits>::TransformationIterator<TemplateInterpolator> DTGrid<Traits>::TransformationIterator<TemplateInterpolator>::operator++()
#else
typename DTGrid<Traits>::template TransformationIterator<TemplateInterpolator> DTGrid<Traits>::TransformationIterator<TemplateInterpolator>::operator++()
#endif
{
    // assume that the iterator is valid

    do
    {

	if (k == bbox[2][1])
	{
	    k = bbox[2][0];
	    if (j == bbox[1][1])
	    {
		j = bbox[1][0];
		i++;
		if (i > bbox[0][1])
		{
		    // we moved past the last element, just return!
		    return *this;
		}
	    }
	    else
	    {
		j++;
	    }
	}
	else
	{
	    k++;
	}

	// at this point we know that the iterator is valid, ie. we have not
	// moved past the last element

	Data values[8];  // eight values of corners of voxel

	// transform the grid coordinate (i,j,k) to the grid coordinate system
	// of the other grid
	Matrix::Vector3<Real> vec3 = gridToGrid2 * Matrix::Vector3<Real>((Real)i, (Real)j, (Real)k);
	// convert vec3 to a unit offset

	Matrix::Vector3<Index> minVec = Matrix::Vector3<Real>::floor(vec3);

	vec3 = vec3 - Matrix::Vector3<Real>::floor(vec3);

	// Ordering of corner values:
	//
	//     X:         X + 1:
	// (and Y to the right, Z up) 
	//
	//  1 -- 3       5 -- 7
	//  |    |       |    |
	//  0 -- 2       4 -- 6
	//


	// Get values at positions (x,y,z),(x,y,z+1),(x,y+1,z),(x,y+1,z+1),(x+1,y,z),(x+1,y,z+1),(x+1,y+1,z),(x+1,y+1,z+1)
	// instantaneously
	grid2->getVoxelValues(minVec[0], minVec[1], minVec[2], values);

	// interpolate to the final value
	value = interpolator->interp(values, vec3);

    } 
    while ( fabs(value) >= outsideConstant );

    return *this;
}


// END TRANSFORMATION ITERATOR IMPLEMENTATION



template<class Traits>template<StencilFormat stencilFormat>
typename DTGrid<Traits>::Index DTGrid<Traits>::getStencilXWidth()
{
    switch (stencilFormat)
    {
    case SF_NONE:
	return 1;
    case SF_FIRSTORDER:
    case SF_FIRSTORDER_CURVATURE:
    case SF_BOX:
	return 3;
    case SF_WENO:
    case SF_WENO_CURVATURE:
	return 7;
    case SF_VOXEL_GRAD:
	return 4;
    case SF_VOXEL:
	return 2;
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////////
// DTGRID INITIALIZATION AND SERIALIZATION METHODS
///////////////////////////////////////////////////////////////////////////////////

template<class Traits>
void DTGrid<Traits>::initTransformation()
{
    tx = 0;
    ty = 0;
    tz = 0;
    rx = 0;
    ry = 0;
    rz = 0;

    gridToWorld = Matrix::Matrix4x4<Real>::identity();
    worldToGrid = Matrix::Matrix4x4<Real>::identity();
}



template<class Traits>
void DTGrid<Traits>::initBoundingBox()
{
    effectiveGridSize(bbox);
}


template<class Traits>
void DTGrid<Traits>::getInitParams(InitParams *initParams) const
{
    initParams->beta = beta;
    initParams->gamma = gamma;
    initParams->dx = dx;
    initParams->insideConstant = insideConstant;
    initParams->outsideConstant = outsideConstant;
}



template<class Traits>
void DTGrid<Traits>::saveSparseVolume(const std::string& fileName, std::string outputFormat)
{
    // Default format
    if (outputFormat == "")
    {
	outputFormat = SvolSaverDTTop1<Index, Real, Data, UInt>::idString;
    }

    
    // bounding box MUST be valid in saved svol
    initBoundingBox();

    //cout << "Saving DTGrid...";
    SvolSaverInterface<Index, Real, Data, UInt> *svolSaver;
    svolSaver = Grids::SvolSaverFactory::getInstance<Index, Real, Data, UInt>(outputFormat);
    if (svolSaver == NULL)
    {
	Core::throwDefaultException(std::string("Unknown output format ")+outputFormat, __FILE__, __LINE__);
    }
    DTGridConstructor<Index, Real, Data, UInt> dtgridConstructor;

    typename DTGridConstructor<Index, Real, Data, UInt>::SaveOutput saveOutput;
    typename DTGridConstructor<Index, Real, Data, UInt>::SaveInput saveInput;

    saveInput.dx = dx;
    saveInput.beta = beta;
    saveInput.gamma = gamma;
    saveInput.insideConstant = insideConstant;
    saveInput.outsideConstant = outsideConstant;
    saveInput.integralTranslation[0] = 0;
    saveInput.integralTranslation[1] = 0;
    saveInput.integralTranslation[2] = 0;
    saveInput.translation[0] = tx;
    saveInput.translation[1] = ty;
    saveInput.translation[2] = tz;
    saveInput.rotation[0] = rx;
    saveInput.rotation[1] = ry;
    saveInput.rotation[2] = rz;
    saveInput.fileName = fileName;
    saveInput.skipValuesLargerThanGamma = false;
    saveInput.mendSvol = !Traits::openLevelSetsSupport;

    TubeIterator iter = beginTubeIterator();
    dtgridConstructor.save(saveInput, &saveOutput, svolSaver, iter);
    delete svolSaver;
    //cout << "done!" << endl;
}


template<class Traits>
void DTGrid<Traits>::loadSparseVolume(const std::string& fileName)
{
    std::string header;	

    clear();

    SvolLoaderInterface<Index, Real, Data, UInt> *svolLoader;
    svolLoader = Grids::SvolLoaderFactory::getInstance<Index, Real, Data, UInt>(fileName, &header);
    if (svolLoader == NULL)
    {
	Core::throwDefaultException(std::string("Unknown header format: ")+header, __FILE__, __LINE__);
    }
    DTGridConstructor<Index, Real, Data, UInt> dtgridConstructor;

    typename DTGridConstructor<Index, Real, Data, UInt>::LoadOutput loadOutput;
    typename DTGridConstructor<Index, Real, Data, UInt>::LoadInput loadInput;

    loadInput.numIndexEndMarkers = numIndexEndMarkers;
    loadInput.numAAEndMarkers = numAAEndMarkers;
    loadInput.numValueEndMarkers = numValueEndMarkers;
    loadInput.fileName = fileName;

    dtgridConstructor.load(loadInput, &loadOutput, svolLoader);
    delete svolLoader;

    va1D = loadOutput.va1D;
    va2D = loadOutput.va2D;
    va3D = loadOutput.va3D;

    aa1D = loadOutput.aa1D;
    aa2D = loadOutput.aa2D;
    aa3D = loadOutput.aa3D;

    numVa1D = loadOutput.numVa1D;
    numVa2D = loadOutput.numVa2D;
    numVa3D = loadOutput.numVa3D;

    numXIndex = loadOutput.numXIndex;
    numYIndex = loadOutput.numYIndex;
    numZIndex = loadOutput.numZIndex;

    xIndex = loadOutput.xIndex;
    yIndex = loadOutput.yIndex;
    zIndex = loadOutput.zIndex;

    lastXIndex = loadOutput.lastXIndex;
    lastYIndex = loadOutput.lastYIndex;
    lastZIndex = loadOutput.lastZIndex;

    dx = loadOutput.dx;
    beta = loadOutput.beta;
    gamma = loadOutput.gamma;

    insideConstant = loadOutput.insideConstant;
    outsideConstant = loadOutput.outsideConstant;

    bbox[0][0] = loadOutput.bbox[0][0];
    bbox[0][1] = loadOutput.bbox[0][1];
    bbox[1][0] = loadOutput.bbox[1][0];
    bbox[1][1] = loadOutput.bbox[1][1];
    bbox[2][0] = loadOutput.bbox[2][0];
    bbox[2][1] = loadOutput.bbox[2][1];

    tx = loadOutput.translation[0];
    ty = loadOutput.translation[1];
    tz = loadOutput.translation[2];

    rx = loadOutput.rotation[0];
    ry = loadOutput.rotation[1];
    rz = loadOutput.rotation[2];

    setEndMarkers();
    initBuffers();

    gridToWorld = Matrix::Matrix4x4<Real>::translation(tx,ty,tz) * Matrix::Matrix4x4<Real>::rotationXYZ(rx,ry,rz) * Matrix::Matrix4x4<Real>::scale(dx);
    worldToGrid = Matrix::Matrix4x4<Real>::scale(1/dx) * Matrix::Matrix4x4<Real>::rotationXYZ(rx,ry,rz).transpose() * Matrix::Matrix4x4<Real>::translation(-tx,-ty,-tz);
}



template<class Traits>
void DTGrid<Traits>::setEndMarkers()
{
    va1D[numVa1D] = numYIndex;
    va2D[numVa2D] = numZIndex;
    if (va3D!=NULL)
    {
	va3D[numVa3D] = (Data)0;
    }
    // this requirement for Index end markers is used in incrementUntil
    xIndex[numXIndex] = std::numeric_limits<Index>::max();
    yIndex[numYIndex] = std::numeric_limits<Index>::max();
    zIndex[numZIndex] = std::numeric_limits<Index>::max();
    // this requirement for aa end markers is unsed in incrementUntil since
    // we look up in the aa array eg at position numXIndex (iv1D etc may point 2 beyond the last element) 
    aa1D[(numXIndex)>>1] = numVa1D;
    aa2D[(numYIndex)>>1] = numVa2D;
    aa3D[(numZIndex)>>1] = numVa3D;
}




template<class Traits>template<typename InitIter>
DTGrid<Traits>::DTGrid(InitIter beginIter, InitIter endIter, InitParams initParams)
{
    Index i, j, k;
    Data v;
    InitIter beginIter2, endIter2;
    IterationOrder iterationOrder = beginIter.getIterationOrder();

    beta = initParams.beta;
    gamma = initParams.gamma;
    outsideConstant = initParams.outsideConstant;
    insideConstant = initParams.insideConstant;
    dx = initParams.dx;

    beginIter2 = beginIter;
    endIter2 = endIter;

    // FIRST PREDICT THE NUMBER OF INDICES AND VALUES THAT WILL BE ADDED TO THE SPARSE GRID

    UInt numIndices3D = 0, numIndices2D = 0, numIndices1D = 0;
    UInt numValues3D = 0, numValues2D = 0, numValues1D = 0;
    UInt numAccs3D, numAccs2D, numAccs1D;
    Index gp3[3], gp2[2], gp1;  // previous grid point
    gp3[0] = gp3[1] = gp3[2] = -2;
    gp2[0] = gp2[1] = -2;
    gp1 = -2;

    while (beginIter != endIter)
    {
	Index i, j, k;

	if (iterationOrder == IteratorOrder_ZYX)
	{
	    i = beginIter.getK();
	    j = beginIter.getJ();
	    k = beginIter.getI();
	}
	else
	{
	    i = beginIter.getI();
	    j = beginIter.getJ();
	    k = beginIter.getK();
	}

	v = beginIter.getValue();
	if (fabs(v)<gamma)
	{
	    numValues3D++;

	    if (i != gp1)
	    {
		numValues1D++;
		if (i != gp1+1)
		{
		    numIndices1D++;
		}
		gp1 = i;
	    }

	    if (i != gp2[0] || j != gp2[1])
	    {
		numValues2D++;
		if (i != gp2[0] || j != gp2[1]+1)
		{
		    numIndices2D++;
		}
		gp2[0] = i;
		gp2[1] = j;
	    }

	    if (i != gp3[0] || j != gp3[1] || k != gp3[2]+1)
	    {
		numIndices3D++;
	    }

	    gp3[0] = i;
	    gp3[1] = j;
	    gp3[2] = k;

	}

	++beginIter;
    }


    numAccs1D = numIndices1D * 2;
    numIndices1D = numIndices1D * 2;
    numValues1D += 1;

    numAccs2D = numIndices2D * 2;
    numIndices2D = numIndices2D * 2;
    numValues2D += 1;

    numAccs3D = numIndices3D * 2;
    numIndices3D = numIndices3D * 2;

    va1D = new UInt[numValues1D];
    numVa1D = 0;
    xIndex = new Index[numIndices1D];
    numXIndex = 0;
    lastXIndex = -1;
    aa1D = new UInt[(numIndices1D>>1)+numAAEndMarkers];

    va2D = new UInt[numValues2D];
    numVa2D = 0;
    yIndex = new Index[numIndices2D];
    numYIndex = 0;
    lastYIndex = -1;
    aa2D = new UInt[numIndices2D>>1];

    va3D = new Data[numValues3D];
    numVa3D = 0;
    zIndex = new Index[numIndices3D];
    numZIndex = 0;
    lastZIndex = -1;
    aa3D = new UInt[numIndices3D>>1];

    ////////////////////////////////////////////////////////////////

    while (beginIter2 != endIter2)
    {
	if (iterationOrder == IteratorOrder_ZYX)
	{
	    i = beginIter2.getK();
	    j = beginIter2.getJ();
	    k = beginIter2.getI();
	}
	else
	{
	    i = beginIter2.getI();
	    j = beginIter2.getJ();
	    k = beginIter2.getK();
	}

	v = beginIter2.getValue();
	if (fabs(v)<gamma)
	{
	    push3D<false,false>(i, j, k, v);			
	}

	++beginIter2;
    }

    // set end markers
    va1D[numVa1D] = numYIndex;
    va2D[numVa2D] = numZIndex;

    initBuffers();
    initTransformation();
    initBoundingBox();

    openLevelSetBBoxValid = false;
}	



// Input: 
//        Index dim[3]: implies a bounding box of [0,dim[0]-1], [0,dim[1]-1], [0,dim[2]-1]
template<class Traits>
DTGrid<Traits>::DTGrid(Index dim[3], Real rx, Real ry, Real rz, Real tx, Real ty, Real tz, InitParams initParams)
{
    beta = initParams.beta;
    gamma = initParams.gamma;
    outsideConstant = initParams.outsideConstant;
    insideConstant = initParams.insideConstant;
    dx = initParams.dx;

    this->bbox[0][0] = 0;
    this->bbox[0][1] = dim[0]-1;
    this->bbox[1][0] = 0;
    this->bbox[1][1] = dim[1]-1;
    this->bbox[2][0] = 0;
    this->bbox[2][1] = dim[2]-1;

    this->rx = rx;
    this->ry = ry;
    this->rz = rz;
    this->tx = tx;
    this->ty = ty;
    this->tz = tz;

    gridToWorld = Matrix::Matrix4x4<Real>::translation(tx,ty,tz) * Matrix::Matrix4x4<Real>::rotationXYZ(rx,ry,rz) * Matrix::Matrix4x4<Real>::scale(dx);
    worldToGrid = Matrix::Matrix4x4<Real>::scale(1/dx) * Matrix::Matrix4x4<Real>::rotationXYZ(rx,ry,rz).transpose() * Matrix::Matrix4x4<Real>::translation(-tx,-ty,-tz);

    clear(false);

    openLevelSetBBoxValid = false;
}


template<class Traits>
DTGrid<Traits>::DTGrid(const DTGrid& dtgrid2)
{
    beta = dtgrid2.beta;
    gamma = dtgrid2.gamma;
    outsideConstant = dtgrid2.outsideConstant;
    insideConstant = dtgrid2.insideConstant;
    dx = dtgrid2.dx;
    bbox[0][0] = dtgrid2.bbox[0][0];
    bbox[0][1] = dtgrid2.bbox[0][1];
    bbox[1][0] = dtgrid2.bbox[1][0];
    bbox[1][1] = dtgrid2.bbox[1][1];
    bbox[2][0] = dtgrid2.bbox[2][0];
    bbox[2][1] = dtgrid2.bbox[2][1];

    worldToGrid = dtgrid2.worldToGrid;
    gridToWorld = dtgrid2.gridToWorld;
    rx = dtgrid2.rx;
    ry = dtgrid2.ry;
    rz = dtgrid2.rz;
    tx = dtgrid2.tx;
    ty = dtgrid2.ty;
    tz = dtgrid2.tz;

    va1D = dtgrid2.va1D;        
    numVa1D = dtgrid2.numVa1D;    
    xIndex = dtgrid2.xIndex;            
    numXIndex = dtgrid2.numXIndex; 
    lastXIndex = dtgrid2.lastXIndex;
    aa1D = dtgrid2.aa1D;         

    va2D = dtgrid2.va2D;        
    numVa2D = dtgrid2.numVa2D;
    yIndex = dtgrid2.yIndex;             
    numYIndex = dtgrid2.numYIndex;
    lastYIndex = dtgrid2.lastYIndex;
    aa2D = dtgrid2.aa2D;         

    zIndex = dtgrid2.zIndex;              
    numZIndex = dtgrid2.numZIndex;
    lastZIndex = dtgrid2.lastZIndex;
    aa3D = dtgrid2.aa3D;         

    numVa3D = dtgrid2.numVa3D;
    va3D = dtgrid2.va3D;                

    openLevelSetBBoxValid = dtgrid2.openLevelSetBBoxValid;

    openLevelSetBBox[0][0] = dtgrid2.openLevelSetBBox[0][0];
    openLevelSetBBox[1][0] = dtgrid2.openLevelSetBBox[1][0];
    openLevelSetBBox[2][0] = dtgrid2.openLevelSetBBox[2][0];
    openLevelSetBBox[0][1] = dtgrid2.openLevelSetBBox[0][1];
    openLevelSetBBox[1][1] = dtgrid2.openLevelSetBBox[1][1];
    openLevelSetBBox[2][1] = dtgrid2.openLevelSetBBox[2][1];

    initBuffers();
}



template<class Traits>
void DTGrid<Traits>::init(InitParams initParams)
{
    beta = initParams.beta;
    gamma = initParams.gamma;
    outsideConstant = initParams.outsideConstant;
    insideConstant = initParams.insideConstant;
    dx = initParams.dx;

    bbox[0][0] = 0;
    bbox[0][1] = 0;
    bbox[1][0] = 0;
    bbox[1][1] = 0;
    bbox[2][0] = 0;
    bbox[2][1] = 0;

    va1D = NULL;        
    numVa1D = 0;    
    xIndex = NULL;            
    numXIndex = 0; 
    lastXIndex = 0;
    aa1D = NULL;         

    va2D = NULL;        
    numVa2D = 0;
    yIndex = NULL;             
    numYIndex = 0;
    lastYIndex = 0;
    aa2D = NULL;         

    zIndex = NULL;              
    numZIndex = 0;
    lastZIndex = 0;
    aa3D = NULL;         

    numVa3D = 0;
    va3D = NULL;                

    openLevelSetBBoxValid = false;
    	
    initBuffers();
    initTransformation();
}


template<class Traits>
DTGrid<Traits>::DTGrid(InitParams initParams)
{
    init(initParams);
}







template<class Traits>
DTGrid<Traits>::~DTGrid()
{
    delete[] va1D;
    delete[] va2D;
    delete[] xIndex;
    delete[] yIndex;
    delete[] zIndex;
    delete[] aa1D;
    delete[] aa2D;
    delete[] aa3D;

    UInt i;
    for (i=0; i<maxNumBuffers; i++)
    {
	freeBuffer(i);
    }
}



template<class Traits>template <bool safePush>
void DTGrid<Traits>::push1D(Index x, UInt val)
{
    if (safePush)
    {
	if ( safeStorage->xIndex.size() > 0 && safeStorage->xIndex.back() == x-1 )
	{
	    // this element is part of an existing connected component
	    safeStorage->xIndex.back() = x;
	}
	else
	{
	    // begin new connected component
	    safeStorage->aa1D.push_back(static_cast<UInt>(safeStorage->va1D.size()));
	    safeStorage->xIndex.push_back(x);
	    safeStorage->xIndex.push_back(x);
	}

	safeStorage->va1D.push_back(val);
    }
    else
    {
	va1D[numVa1D] = val;

	if ( numXIndex > 0 && xIndex[lastXIndex] == x-1 )
	{
	    // this element is part of an existing connected component
	    xIndex[lastXIndex] = x;
	    //aa1D[lastXIndex] += 1;
	}
	else
	{
	    // begin new connected component
	    xIndex[numXIndex] = x;
	    aa1D[numXIndex>>1] = numVa1D;
	    numXIndex++;	
	    xIndex[numXIndex] = x;
	    numXIndex++;	
	    lastXIndex = numXIndex-1;
	}
	numVa1D++;
    }
}






template<class Traits>template <bool safePush>
void DTGrid<Traits>::push2D(Index x, Index y, UInt val)
{
    if (safePush)
    {
	if (safeStorage->va2D.size()==0 || safeStorage->xIndex.back() != x )
	{
	    // there are no elements in this column
	    push1D<true>(x, static_cast<UInt>(safeStorage->yIndex.size()));
	    // begin new connected component
	    safeStorage->aa2D.push_back(static_cast<UInt>(safeStorage->va2D.size()));
	    safeStorage->yIndex.push_back(y);
	    safeStorage->yIndex.push_back(y);
	}
	else if ( safeStorage->yIndex.back() == y-1 )
	{
	    // this element is part of an existing connected component
	    safeStorage->yIndex.back() = y;
	}
	else
	{
	    // begin new connected component
	    safeStorage->aa2D.push_back(static_cast<UInt>(safeStorage->va2D.size()));
	    safeStorage->yIndex.push_back(y);
	    safeStorage->yIndex.push_back(y);
	}

	safeStorage->va2D.push_back(val);
    }
    else
    {
	va2D[numVa2D] = val;

	if (numVa2D==0 || xIndex[lastXIndex] != x )
	{
	    // there are no elements in this column
	    push1D<false>(x, numYIndex);
	    // begin new connected component
	    yIndex[numYIndex] = y;
	    aa2D[numYIndex>>1] = numVa2D;
	    numYIndex++;	
	    yIndex[numYIndex] = y;
	    numYIndex++;
	    lastYIndex = numYIndex-1;
	}
	else if ( yIndex[lastYIndex] == y-1 )
	{
	    // this element is part of an existing connected component
	    yIndex[lastYIndex] = y;
	}
	else
	{
	    // begin new connected component
	    yIndex[numYIndex] = y;
	    aa2D[numYIndex>>1] = numVa2D;
	    numYIndex++;	
	    yIndex[numYIndex] = y;
	    numYIndex++;
	    lastYIndex = numYIndex-1;
	}
	numVa2D++;
    }
}



template<class Traits>
void DTGrid<Traits>::push(Index x, Index y, Index z, Data val)
{
    push3D<true, true>(x,y,z,val);
}


template<class Traits> template<bool safePush, bool safePushValue>
void DTGrid<Traits>::push3D(Index x, Index y, Index z, Data val)
{
    if (safePush)
    {
	if (!safePushValue)
	{
	    va3D[numVa3D] = val;
	}
	else
	{
	    safeStorage->va3D.push_back(val);
	}


	if (numVa3D==0 || safeStorage->xIndex.back() != x || safeStorage->yIndex.back() != y )
	{
	    // there are no elements in this column
	    push2D<true>(x, y, static_cast<UInt>(safeStorage->zIndex.size()));
	    // begin new connected component
	    safeStorage->aa3D.push_back(numVa3D);
	    safeStorage->zIndex.push_back(z);
	    safeStorage->zIndex.push_back(z);
	}
	else if ( safeStorage->zIndex.back() == z-1 )
	{
	    // this element is part of an existing connected component
	    safeStorage->zIndex.back() = z;
	}
	else
	{
	    // begin new connected component
	    safeStorage->aa3D.push_back(numVa3D);
	    safeStorage->zIndex.push_back(z);
	    safeStorage->zIndex.push_back(z);
	}

	numVa3D++;
    }
    else
    {
	va3D[numVa3D] = val;

	if (numVa3D==0 || xIndex[lastXIndex] != x || yIndex[lastYIndex] != y )
	{
	    // there are no elements in this column
	    push2D<false>(x, y, numZIndex);
	    // begin new connected component
	    zIndex[numZIndex] = z;
	    aa3D[numZIndex>>1] = numVa3D;
	    numZIndex++;	
	    zIndex[numZIndex] = z;
	    numZIndex++;
	    lastZIndex = numZIndex-1;
	}
	else if ( zIndex[lastZIndex] == z-1 )
	{
	    // this element is part of an existing connected component
	    zIndex[lastZIndex] = z;
	}
	else
	{
	    // begin new connected component
	    zIndex[numZIndex] = z;
	    aa3D[numZIndex>>1] = numVa3D;
	    numZIndex++;	
	    zIndex[numZIndex] = z;
	    numZIndex++;	
	    lastZIndex = numZIndex-1;
	}
	numVa3D++;
    }
}




template<class Traits> template <bool safePush>
void DTGrid<Traits>::push3D(Index x, Index y, Index z)
{
    if (safePush)
    {
	if (numVa3D==0 || safeStorage->xIndex.back() != x || safeStorage->yIndex.back() != y )
	{
	    // there are no elements in this column
	    push2D<true>(x, y, safeStorage->zIndex.size());
	    // begin new connected component
	    safeStorage->aa3D.push_back(numVa3D);
	    safeStorage->zIndex.push_back(z);
	    safeStorage->zIndex.push_back(z);
	}
	else if ( safeStorage->zIndex.back() == z-1 )
	{
	    // this element is part of an existing connected component
	    safeStorage->zIndex.back() = z;
	}
	else
	{
	    // begin new connected component
	    safeStorage->aa3D.push_back(numVa3D);
	    safeStorage->zIndex.push_back(z);
	    safeStorage->zIndex.push_back(z);
	}

	numVa3D++;
    }
    else
    {
	if (numVa3D==0 || xIndex[lastXIndex] != x || yIndex[lastYIndex] != y )
	{
	    // there are no elements in this column
	    push2D<false>(x, y, numZIndex);
	    // begin new connected component
	    zIndex[numZIndex] = z;
	    aa3D[numZIndex>>1] = numVa3D;
	    numZIndex++;	
	    zIndex[numZIndex] = z;
	    numZIndex++;
	    lastZIndex = numZIndex-1;
	}
	else if ( zIndex[lastZIndex] == z-1 )
	{
	    // this element is part of an existing connected component
	    zIndex[lastZIndex] = z;
	}
	else
	{
	    // begin new connected component
	    zIndex[numZIndex] = z;
	    aa3D[numZIndex>>1] = numVa3D;
	    numZIndex++;	
	    zIndex[numZIndex] = z;
	    numZIndex++;	
	    lastZIndex = numZIndex-1;
	}
	numVa3D++;
    }
}



template<class Traits>
void DTGrid<Traits>::pushAndMend3D(Index x, Index y, Index z, Data val)
{
    va3D[numVa3D] = val;

    if (numVa3D==0 || xIndex[lastXIndex] != x || yIndex[lastYIndex] != y )
    {
	if (numVa3D != 0)
	{
	    // make sure the end of last connected component was positive
	    // we do not handle if the last grid point added was erroneous
	    va3D[numVa3D-1] = fabs(va3D[numVa3D-1]);
	}
	// also make sure that the start of the first connected component in p-column is positive
	va3D[numVa3D] = fabs(va3D[numVa3D]);

	// there are no elements in this column
	push2D<false>(x, y, numZIndex);
	// begin new connected component
	zIndex[numZIndex] = z;
	aa3D[numZIndex>>1] = numVa3D;
	numZIndex++;	
	zIndex[numZIndex] = z;
	numZIndex++;
	lastZIndex = numZIndex-1;
    }
    else if ( zIndex[lastZIndex] == z-1 )
    {
	// this element is part of an existing connected component
	zIndex[lastZIndex] = z;
    }
    else
    {
	// begin new connected component

	// make sure that the "adjacent" ends of connected components have consistent signs.
	// favour outside
	if ( va3D[numVa3D]<0 && va3D[numVa3D-1]>0 )
	{
	    va3D[numVa3D] = fabs(va3D[numVa3D]);
	}
	else if ( va3D[numVa3D]>0 && va3D[numVa3D-1]<0 )
	{
	    va3D[numVa3D-1] = fabs(va3D[numVa3D-1]);
	}

	zIndex[numZIndex] = z;
	aa3D[numZIndex>>1] = numVa3D;
	numZIndex++;	
	zIndex[numZIndex] = z;
	numZIndex++;	
	lastZIndex = numZIndex-1;
    }
    numVa3D++;
}


template<class Traits>
void DTGrid<Traits>::beginSafePush()
{
    numVa3D = 0;
    safeStorage = new SafeStorage();
}

template<class Traits>
void DTGrid<Traits>::endSafePush()
{
    UInt i;
    UInt size;

    typename vector<UInt>::iterator uiIter;
    typename vector<Index>::iterator iIter;
    typename vector<Data>::iterator rIter;

    numVa1D = static_cast<UInt>(safeStorage->va1D.size());
    va1D = new UInt[numVa1D + numValueEndMarkers];
    for (i=0, uiIter=safeStorage->va1D.begin(); i<numVa1D; i++, uiIter++)
    {
	va1D[i] = *uiIter;
    }
    safeStorage->va1D.clear();
    numXIndex = static_cast<UInt>(safeStorage->xIndex.size());
    xIndex = new Index[numXIndex + numIndexEndMarkers];
    for (i=0, iIter=safeStorage->xIndex.begin(); i<numXIndex; i++, iIter++)
    {
	xIndex[i] = *iIter;
    }
    safeStorage->xIndex.clear();
    aa1D = new UInt[(numXIndex>>1) + numAAEndMarkers];
    for (i=0, size=(numXIndex>>1), uiIter=safeStorage->aa1D.begin(); i<size; i++, uiIter++)
    {
	aa1D[i] = *uiIter;
    }
    safeStorage->aa1D.clear();

    numVa2D = static_cast<UInt>(safeStorage->va2D.size());
    va2D = new UInt[numVa2D + numValueEndMarkers];
    for (i=0, uiIter=safeStorage->va2D.begin(); i<numVa2D; i++, uiIter++)
    {
	va2D[i] = *uiIter;
    }
    safeStorage->va2D.clear();
    numYIndex = static_cast<UInt>(safeStorage->yIndex.size());
    yIndex = new Index[numYIndex + numIndexEndMarkers];
    for (i=0, iIter=safeStorage->yIndex.begin(); i<numYIndex; i++, iIter++)
    {
	yIndex[i] = *iIter;
    }
    safeStorage->yIndex.clear();
    aa2D = new UInt[(numYIndex>>1) + numAAEndMarkers];
    for (i=0, size=(numYIndex>>1), uiIter=safeStorage->aa2D.begin(); i<size; i++, uiIter++)
    {
	aa2D[i] = *uiIter;
    }
    safeStorage->aa2D.clear();

    if (safeStorage->va3D.size() > 0)
    {
	va3D = new Data[numVa3D + numValueEndMarkers];
	for (i=0, rIter=safeStorage->va3D.begin(); i<numVa3D; i++, rIter++)
	{
	    va3D[i] = *rIter;
	}
	safeStorage->va3D.clear();
    }
    numZIndex = static_cast<UInt>(safeStorage->zIndex.size());
    zIndex = new Index[numZIndex + numIndexEndMarkers];
    for (i=0, iIter=safeStorage->zIndex.begin(); i<numZIndex; i++, iIter++)
    {
	zIndex[i] = *iIter;
    }
    safeStorage->zIndex.clear();
    aa3D = new UInt[(numZIndex>>1) + numAAEndMarkers];
    for (i=0, size=(numZIndex>>1), uiIter=safeStorage->aa3D.begin(); i<size; i++, uiIter++)
    {
	aa3D[i] = *uiIter;
    }
    safeStorage->aa3D.clear();

    delete safeStorage;
    initBuffers();
    setEndMarkers();

    lastXIndex = numXIndex-1;
    lastYIndex = numYIndex-1;
    lastZIndex = numZIndex-1;

}





///////////////////////////////////////////////////////////////////////////////////
// REBUILD AND DILATION METHODS
///////////////////////////////////////////////////////////////////////////////////


template<class Traits>
void DTGrid<Traits>::dilate1D(Index width, DTGrid *dtgrid2)
{	
    CCUnion_Public<Index> ccu = CCUnion_Public<Index>(width, dtgrid2->xIndex);
    vector<Index> xIndexSafe;              
    vector<UInt> aa1DSafe;
    typename vector<UInt>::iterator aaIter, aaEnd;
    typename vector<Index>::iterator iIter;

    xIndexSafe.reserve(dtgrid2->getNumXIndex());
    aa1DSafe.reserve(dtgrid2->getNumAA1D());

    Index e, s;
    UInt i;

    ccu.addCCI(0, dtgrid2->numXIndex);

    ccu.startCCUnion();

    while (!ccu.unionCCDone())
    {
	ccu.getUnionCC(&s, &e);

	// The assumption is that the narrow band was fully contained within the bounding box prior to the rebuild 
	if (Traits::openLevelSetsSupport)
	{
	    // In this case we have to also cut the dilated level set against the openLevelSetBBox
	    s = std::max(s, openLevelSetBBox[0][0]);
	    s = std::min(s, openLevelSetBBox[0][1]);
	    e = std::max(e, openLevelSetBBox[0][0]);
	    e = std::min(e, openLevelSetBBox[0][1]);

	    xIndexSafe.push_back(s);
	    aa1DSafe.push_back(numVa1D);
	    numVa1D += e-s;
	    xIndexSafe.push_back(e);
	    numVa1D++;
	}
	else
	{
	    xIndexSafe.push_back(s);
	    aa1DSafe.push_back(numVa1D);
	    numVa1D += e-s;
	    xIndexSafe.push_back(e);
	    numVa1D++;
	}
    }

    numXIndex = static_cast<UInt>(xIndexSafe.size());
    
    xIndex = new Index[numXIndex + numIndexEndMarkers];
    aa1D = new UInt[(numXIndex>>1) + numAAEndMarkers];
    aaIter = aa1DSafe.begin();
    aaEnd = aa1DSafe.end();
    iIter = xIndexSafe.begin();
    i = 0;
    while (aaIter != aaEnd)
    {
	xIndex[i] = *iIter;
	aa1D[i>>1] = *aaIter;
	iIter++;
	i++;
	xIndex[i] = *iIter;
	iIter++;
	i++;
	aaIter++;
    }

    aa1DSafe.clear();
    xIndexSafe.clear();

    va1D = new UInt[numVa1D + numValueEndMarkers];

    // set end marker
    xIndex[numXIndex] = std::numeric_limits<Index>::max();
    aa1D[numXIndex>>1] = numVa1D;

    lastXIndex = numXIndex-1;
}



template<class Traits>
void DTGrid<Traits>::dilate2D(Index width, DTGrid *dtgrid2)
{
    CCUnion_Public<Index> ccu = CCUnion_Public<Index>(width, dtgrid2->yIndex);
    vector<Index> yIndexSafe;
    vector<UInt> aa2DSafe;
    typename vector<UInt>::iterator aaIter, aaEnd;
    typename vector<Index>::iterator iIter;

    Index e, s;
    UInt ei, si;
    UInt i;
    _Iterator1D iter, iend;
    _StencilIterator1D siter;
    Index x, x_tmp; 


    dilate1D(width, dtgrid2);

    yIndexSafe.reserve(dtgrid2->getNumYIndex());
    aa2DSafe.reserve(dtgrid2->getNumAA2D());

    // iterators over dilated 1D component
    iter = begin1D();
    iend = end1D();
    // stencil iterator over original 1D component
    siter = dtgrid2->beginStencil1D(width);

    x = iter.getX();

    while (iter != iend)
    {
	// set proper value in 2D sparse grid constituent:
	// pointer into the yIndex array 
	iter.setValue(numYIndex);


	// now add new intervals from front of stencil 
	if (siter.isEndValid())
	{
	    siter.getEndValues(&si, &ei);
	    ccu.addCCI(si, ei);
	}


	// now compute union of connected components
	ccu.startCCUnion();

	while (!ccu.unionCCDone())
	{
	    ccu.getUnionCC(&s, &e);

	    if (Traits::openLevelSetsSupport)
	    {
		// In this case we have to also cut the dilated level set against the openLevelSetBBox
		if (x == siter.getCenterX())
		{
		    s = std::max(s, openLevelSetBBox[1][0]);
		    s = std::min(s, openLevelSetBBox[1][1]);
		    e = std::max(e, openLevelSetBBox[1][0]);
		    e = std::min(e, openLevelSetBBox[1][1]);

		    yIndexSafe.push_back(s);
		    aa2DSafe.push_back(numVa2D);
		    numVa2D += e-s;
		    yIndexSafe.push_back(e);
		    numVa2D++;
		    numYIndex+=2;
		}
	    }
	    else
	    {
		yIndexSafe.push_back(s);
		aa2DSafe.push_back(numVa2D);
		numVa2D += e-s;
		yIndexSafe.push_back(e);
		numVa2D++;
		numYIndex+=2;
	    }
	}

	// now remove old intervals from back of stencil
	if (siter.isStartValid())
	{
	    siter.getStartValues(&si, &ei);

	    ccu.removeCCI();
	}

	if (Traits::openLevelSetsSupport)
	{
	    if (siter.getCenterX()==x)
	    {
		iter++;
		x_tmp = iter.getX();
		siter.increment(x_tmp-x);
		x = x_tmp;
	    }
	    else
	    {
		siter.increment(1);
	    }
	}
	else
	{
	    iter++;
	    x_tmp = iter.getX();
	    siter.increment(x_tmp-x);
	    x = x_tmp;
	}
    }

    yIndex = new Index[numYIndex + numIndexEndMarkers];
    aa2D = new UInt[(numYIndex>>1) + numAAEndMarkers];
    aaIter = aa2DSafe.begin();
    aaEnd = aa2DSafe.end();
    iIter = yIndexSafe.begin();
    i = 0;
    while (aaIter != aaEnd)
    {
	aa2D[i>>1] = *aaIter;
	aaIter++;
	yIndex[i] = *iIter;
	iIter++;
	i++;
	yIndex[i] = *iIter;
	iIter++;
	i++;
    }


    aa2DSafe.clear();
    yIndexSafe.clear();

    va2D = new UInt[numVa2D + numValueEndMarkers];

    // set end marker
    yIndex[numYIndex] = std::numeric_limits<Index>::max();
    aa2D[numYIndex>>1] = numVa2D;
    va1D[numVa1D] = numYIndex;

    lastYIndex = numYIndex-1;
}



template<class Traits>
void DTGrid<Traits>::dilate3D(Index width, DTGrid *dtgrid2)
{
    CCUnion_Public<Index> ccu = CCUnion_Public<Index>(width, dtgrid2->zIndex);
    vector<Index> zIndexSafe;              
    vector<UInt> aa3DSafe;
    typename vector<UInt>::iterator aaIter, aaEnd;
    typename vector<Index>::iterator iIter;
    Index e, s;
    UInt ei, si;
    UInt i;
    _Iterator2D iter, iend;
    Index x, x_tmp, y, y_tmp;
    UInt stencilWidth = width*2+1;

    dilate2D(width, dtgrid2);

    zIndexSafe.reserve(dtgrid2->getNumZIndex());
    aa3DSafe.reserve(dtgrid2->getNumAA3D());

    // iterators over dilated 2D component
    iter = begin2D();
    iend = end2D();
    // stencil iterator over original 2D component
    _StencilIterator2D siter = dtgrid2->beginStencil2D(width);

    x = iter.getX();
    y = iter.getY();


    while (iter != iend)
    {
	// set proper value in 2D sparse grid constituent:
	// pointer into the yIndex array 
	iter.setValue(numZIndex);


	// now add new intervals from front of stencil
	for (i=0; i<stencilWidth; i++)
	{
	    if (siter.isEndValid(i))
	    {
		siter.getEndValues(&si, &ei, i);
		ccu.addCCI(si, ei);
	    }
	}

	// now compute union of connected components
	ccu.startCCUnion();

	while (!ccu.unionCCDone())
	{
	    ccu.getUnionCC(&s, &e);

	    if (Traits::openLevelSetsSupport)
	    {
		// In this case we have to also cut the dilated level set against the openLevelSetBBox
		if (siter.getCenterX()==x && siter.getCenterY()==y)
		{
		    s = std::max(s, openLevelSetBBox[2][0]);
		    s = std::min(s, openLevelSetBBox[2][1]);
		    e = std::max(e, openLevelSetBBox[2][0]);	
		    e = std::min(e, openLevelSetBBox[2][1]);

		    zIndexSafe.push_back(s);
		    aa3DSafe.push_back(numVa3D);
		    numVa3D += e-s;
		    zIndexSafe.push_back(e);
		    numVa3D++;
		    numZIndex+=2;
		}
	    }
	    else
	    {
		zIndexSafe.push_back(s);
		aa3DSafe.push_back(numVa3D);
		numVa3D += e-s;
		zIndexSafe.push_back(e);
		numVa3D++;
		numZIndex+=2;
	    }
	}

	// now remove old intervals from back of stencil
	for (i=0; i<stencilWidth; i++)
	{
	    if (siter.isStartValid(i))
	    {
		siter.getStartValues(&si, &ei, i);
		ccu.removeCCI();
	    }
	}

	if (Traits::openLevelSetsSupport)
	{
	    if (siter.getCenterX()==x && siter.getCenterY()==y)
	    {
		iter++;
		x = iter.getX();
		y = iter.getY();
	    }
	    siter.increment();
	}
	else
	{
	    iter++;
	    x_tmp = iter.getX();
	    y_tmp = iter.getY();
	    siter.increment(x_tmp-x, y_tmp-y);
	    x = x_tmp;
	    y = y_tmp;
	}

    }


    zIndex = new Index[numZIndex + numIndexEndMarkers];


    aa3D = new UInt[(numZIndex>>1)+numAAEndMarkers];
    aaIter = aa3DSafe.begin();
    aaEnd = aa3DSafe.end();
    iIter = zIndexSafe.begin();
    i = 0;
    while (aaIter != aaEnd)
    {
	aa3D[i>>1] = *aaIter;
	aaIter++;
	zIndex[i] = *iIter;
	iIter++;
	i++;
	zIndex[i] = *iIter;
	iIter++;
	i++;
    }
    aa3DSafe.clear();
    zIndexSafe.clear();

    va3D = new Data[numVa3D + numValueEndMarkers];


    // set end markers
    zIndex[numZIndex] = std::numeric_limits<Index>::max();
    va2D[numVa2D] = numZIndex;
    aa3D[numZIndex>>1] = numVa3D;
    va3D[numVa3D] = (Data)0;

    lastZIndex = numZIndex-1;

    initBuffers();
}



template<class Traits>
void DTGrid<Traits>::dilate(Index width, DTGrid *dtgrid2, bool doDelete)
{
    clear(doDelete);	    
    dilate3D(width, dtgrid2);
}


template<class Traits>
void DTGrid<Traits>::dilate(Index width)
{
    DTGrid dtgrid2(*this);
    clear(false);
    dilate3D(width, &dtgrid2);

    // COPY ELEMENTS IN GAMMA TUBE TO NEW BUFFER
    TubeIterator iter, iend, iter2, iend2;
    Data val, vg;
    Index x, y, z, x2, y2, z2;

    // When open level sets are supported, the value assignment becomes a bit more complicated
    if (Traits::openLevelSetsSupport)
    {
	Data prevVal;
	Index prevX, prevY;

	iter2 = dtgrid2.beginTubeIterator();
	iend2 = dtgrid2.endTubeIterator();
	iter = beginTubeIterator();
	iend = endTubeIterator();

	if (iter2!=iend2)
	{
	    val = prevVal = iter2.getValue();
	}

	while (iter2 != iend2)
	{
	    iter2.getIndex(&x,&y,&z);
	    iter.getIndex(&x2,&y2,&z2);

	    if (x==x2 && y==y2 && z==z2)
	    {
		iter.setValue(val);
		iter++;
		iter2++;
		prevVal = val;
		prevX = x2;
		prevY = y2;
		val = iter2.getValue();
	    }
	    else
	    {
		// check if iter is in the same z-column as iter2.
		// if this is the case we should use val as padding value,
		// otherwise we should use prevVal as padding value until iter
		// changes column.
		if (x==x2 && y==y2)
		{
		    // iter and iter2 are in the same column.
		    // in this case we just use val as padding value.
		    vg = ( val<0 ? -gamma : gamma);
		    while (z!=z2)  // We can do this optimization here since we know that dtgrid2 is contained in this dtgrid
		    {
			// pad with +/- gamma
			iter.setValue(vg);
			iter++;
			z2 = iter.getK();
		    }

		    // iter and iter2 are now equal
		    // set value
		    iter.setValue(val);
		    iter++;
		    iter2++;
		    prevVal = val;
		    prevX = x2;
		    prevY = y2;
		    val = iter2.getValue();
		}
		else if (x==x2)
		{
		    // iter and iter2 are not in the same column.
		    // we should therefore use prevVal as padding value until iter
		    // changes column.

		    // we only use prevVal if iter did not change column itself
		    if (iter.getI()==prevX && iter.getJ()==prevY)
		    {
			vg = ( prevVal<0 ? -gamma : gamma);
			while (iter.getI()==x2 && iter.getJ()==y2)
			{
			    // pad with +/- gamma
			    iter.setValue(vg);
			    iter++;
			}
		    }
		    vg = ( val<0 ? -gamma : gamma);
		    iter.getIndex(&x2, &y2, &z2); 			
		    while (x!=x2 || y!=y2 || z!=z2)
		    {
			// pad with +/- gamma
			iter.setValue(vg);
			iter++;
			iter.getIndex(&x2, &y2, &z2); 
		    }

		    // iter and iter2 are now equal
		    // set value
		    iter.setValue(val);
		    iter++;
		    iter2++;
		    prevVal = val;
		    prevX = x2;
		    prevY = y2;
		    val = iter2.getValue();
		}
		else
		{
		    // iter and iter2 are not in the same column.
		    // we should therefore use prevVal as padding value until iter
		    // changes column.
		    vg = ( prevVal<0 ? -gamma : gamma);
		    while (iter.getI()==x2 )
		    {
			// pad with +/- gamma
			iter.setValue(vg);
			iter++;
		    }
		    vg = ( val<0 ? -gamma : gamma);
		    iter.getIndex(&x2, &y2, &z2); 			
		    while (x!=x2 || y!=y2 || z!=z2)
		    {
			// pad with +/- gamma
			iter.setValue(vg);
			iter++;
			iter.getIndex(&x2, &y2, &z2); 
		    }

		    // iter and iter2 are now equal
		    // set value
		    iter.setValue(val);
		    iter++;
		    iter2++;
		    prevVal = val;
		    prevX = x2;
		    prevY = y2;
		    val = iter2.getValue();
		}
	    }
	}	    
	while (iter != iend)
	{
	    iter.setValue(vg);
	    iter++;
	}
    }
    else
    {
	iter2 = dtgrid2.beginTubeIterator();
	iend2 = dtgrid2.endTubeIterator();
	iter = beginTubeIterator();
	iend = endTubeIterator();
	while (iter2 != iend2)
	{
	    iter2.getIndex(&x,&y,&z);
	    val = iter2.getValue();

	    vg = ( val<0 ? insideConstant : outsideConstant);

	    iter.getIndex(&x2,&y2,&z2);
	    while (x!=x2 || y!=y2 || z!=z2)
	    {
		// pad with +/- outsideConstant
		iter.setValue(vg);
		iter++;
		iter.getIndex(&x2,&y2,&z2);
	    }
	    // set value
	    iter.setValue(val);
	    iter++;  
	    iter2++;
	}	    
	while (iter != iend)
	{
	    iter.setValue(vg);
	    iter++;
	}
    }
}


template<class Traits>
void DTGrid<Traits>::rebuildNarrowBand()
{
    TubeIterator iter, iend, iter2, iend2;
    Real val, vg;
    Index x, y, z, x2, y2, z2;


    // REBUILD NARROW BAND USING THE DILATION ALGORITHM

    InitParams initParams = InitParams(dx, beta, gamma, insideConstant, outsideConstant);
    DTGrid<Traits> dtgrid2(initParams);
    dtgrid2.allocateBuffer(0, numValues[0]);


    // 1: COPY ALL ELEMENTS LESS THAN GAMMA


    dtgrid2.beginSafePush();
    dtgrid2.safeStorage->va1D.reserve(getNumVa1D());
    dtgrid2.safeStorage->va2D.reserve(getNumVa2D());
    dtgrid2.safeStorage->xIndex.reserve(getNumXIndex());
    dtgrid2.safeStorage->yIndex.reserve(getNumYIndex());
    dtgrid2.safeStorage->zIndex.reserve(getNumZIndex());
    dtgrid2.safeStorage->aa1D.reserve(getNumAA1D());
    dtgrid2.safeStorage->aa2D.reserve(getNumAA2D());
    dtgrid2.safeStorage->aa3D.reserve(getNumAA3D());

    iter = beginTubeIterator();
    iend = endTubeIterator();

    while (iter != iend)
    {
	val = iter.getValue();
	if (fabs(val) < gamma)
	{
	    iter.getIndex(&x,&y,&z);
	    dtgrid2.template push3D<true,false>(x, y, z, val);
	}
	iter++;
    }

    dtgrid2.endSafePush();


    if (dtgrid2.getNumValues() != 0)
    {
	// 2: DILATE 

	dilate(Traits::rebuildDilationWidth, &dtgrid2);


	// 3: COPY ELEMENTS IN GAMMA TUBE TO NEW BUFFER

	// When open level sets are supported, the value assignment becomes a bit more complicated
	if (Traits::openLevelSetsSupport)
	{
	    Data prevVal;
	    Index prevX, prevY;

	    iter2 = dtgrid2.beginTubeIterator();
	    iend2 = dtgrid2.endTubeIterator();
	    iter = beginTubeIterator();
	    iend = endTubeIterator();

	    if (iter2!=iend2)
	    {
		val = prevVal = iter2.getValue();
	    }

	    while (iter2 != iend2)
	    {
		iter2.getIndex(&x,&y,&z);
		iter.getIndex(&x2,&y2,&z2);

		if (x==x2 && y==y2 && z==z2)
		{
		    iter.setValue(val);
		    iter++;
		    iter2++;
		    prevVal = val;
		    prevX = x2;
		    prevY = y2;
		    val = iter2.getValue();
		}
		else
		{
		    // check if iter is in the same z-column as iter2.
		    // if this is the case we should use val as padding value,
		    // otherwise we should use prevVal as padding value until iter
		    // changes column.
		    if (x==x2 && y==y2)
		    {
			// iter and iter2 are in the same column.
			// in this case we just use val as padding value.
			vg = ( val<0 ? -gamma : gamma);
			while (z!=z2)  // We can do this optimization here since we know that dtgrid2 is contained in this dtgrid
			{
			    // pad with +/- gamma
			    iter.setValue(vg);
			    iter++;
			    z2 = iter.getK();
			}

			// iter and iter2 are now equal
			// set value
			iter.setValue(val);
			iter++;
			iter2++;
			prevVal = val;
			prevX = x2;
			prevY = y2;
			val = iter2.getValue();
		    }
		    else if (x==x2)
		    {
			// iter and iter2 are not in the same column.
			// we should therefore use prevVal as padding value until iter
			// changes column.

			// we only use prevVal if iter did not change column itself
			if (iter.getI()==prevX && iter.getJ()==prevY)
			{
			    vg = ( prevVal<0 ? -gamma : gamma);
			    while (iter.getI()==x2 && iter.getJ()==y2)
			    {
				// pad with +/- gamma
				iter.setValue(vg);
				iter++;
			    }
			}
			vg = ( val<0 ? -gamma : gamma);
			iter.getIndex(&x2, &y2, &z2); 			
			while (x!=x2 || y!=y2 || z!=z2)
			{
			    // pad with +/- gamma
			    iter.setValue(vg);
			    iter++;
			    iter.getIndex(&x2, &y2, &z2); 
			}

			// iter and iter2 are now equal
			// set value
			iter.setValue(val);
			iter++;
			iter2++;
			prevVal = val;
			prevX = x2;
			prevY = y2;
			val = iter2.getValue();
		    }
		    else
		    {
			// iter and iter2 are not in the same column.
			// we should therefore use prevVal as padding value until iter
			// changes column.
			vg = ( prevVal<0 ? -gamma : gamma);
			while (iter.getI()==x2 )
			{
			    // pad with +/- gamma
			    iter.setValue(vg);
			    iter++;
			}
			vg = ( val<0 ? -gamma : gamma);
			iter.getIndex(&x2, &y2, &z2); 			
			while (x!=x2 || y!=y2 || z!=z2)
			{
			    // pad with +/- gamma
			    iter.setValue(vg);
			    iter++;
			    iter.getIndex(&x2, &y2, &z2); 
			}

			// iter and iter2 are now equal
			// set value
			iter.setValue(val);
			iter++;
			iter2++;
			prevVal = val;
			prevX = x2;
			prevY = y2;
			val = iter2.getValue();
		    }
		}
	    }	    
	    while (iter != iend)
	    {
		iter.setValue(vg);
		iter++;
	    }
	}
	else
	{
	    iter2 = dtgrid2.beginTubeIterator();
	    iend2 = dtgrid2.endTubeIterator();
	    iter = beginTubeIterator();
	    iend = endTubeIterator();
	    while (iter2 != iend2)
	    {
		iter2.getIndex(&x,&y,&z);
		val = iter2.getValue();

		vg = ( val<0 ? -gamma : gamma);

		x2 = iter.getI();
		while (x!=x2)
		{
		    // pad with +/- gamma
		    iter.setValue(vg);
		    iter++;
		    x2 = iter.getI();
		}
		y2 = iter.getJ();
		while (y!=y2)
		{
		    // pad with +/- gamma
		    iter.setValue(vg);
		    iter++;
		    y2 = iter.getJ();
		}
		z2 = iter.getK();
		while (z!=z2)
		{
		    // pad with +/- gamma
		    iter.setValue(vg);
		    iter++;
		    z2 = iter.getK();
		}


		// set value
		iter.setValue(val);
		iter++;  
		iter2++;
	    }	    
	    while (iter != iend)
	    {
		iter.setValue(vg);
		iter++;
	    }
	}
    }
    else
    {
	clear();
    }
}



template<class Traits>
void DTGrid<Traits>::rebuildNarrowBand(UInt **permutationArray)
{
    TubeIterator iter, iend, iter2, iend2;
    Real val, vg;
    Index x, y, z, x2, y2, z2;
    std::vector<UInt> tmpPermutationArray;
    UInt permutationInside, permutationOutside, permutationCurrent;

    permutationInside = getNumValues();
    permutationOutside = permutationInside + 1;


    // REBUILD NARROW BAND USING THE DILATION ALGORITHM

    InitParams initParams = InitParams(dx, beta, gamma, insideConstant, outsideConstant);
    DTGrid<Traits> dtgrid2(initParams);
    dtgrid2.allocateBuffer(0, numValues[0]);


    // 1: COPY ALL ELEMENTS LESS THAN GAMMA


    dtgrid2.beginSafePush();
    dtgrid2.safeStorage->va1D.reserve(getNumVa1D());
    dtgrid2.safeStorage->va2D.reserve(getNumVa2D());
    dtgrid2.safeStorage->xIndex.reserve(getNumXIndex());
    dtgrid2.safeStorage->yIndex.reserve(getNumYIndex());
    dtgrid2.safeStorage->zIndex.reserve(getNumZIndex());
    dtgrid2.safeStorage->aa1D.reserve(getNumAA1D());
    dtgrid2.safeStorage->aa2D.reserve(getNumAA2D());
    dtgrid2.safeStorage->aa3D.reserve(getNumAA3D());

    iter = beginTubeIterator();
    iend = endTubeIterator();

    while (iter != iend)
    {
	val = iter.getValue();
	if (fabs(val) < gamma)
	{
	    iter.getIndex(&x,&y,&z);
	    dtgrid2.template push3D<true,false>(x, y, z, val);
	    tmpPermutationArray.push_back(iter.getArrayIndex());
	}
	iter++;
    }

    dtgrid2.endSafePush();


    if (dtgrid2.getNumValues() != 0)
    {
	// 2: DILATE 

	dilate(Traits::rebuildDilationWidth, &dtgrid2);


	// 3: COPY ELEMENTS IN GAMMA TUBE TO NEW BUFFER

	*permutationArray = new UInt[getNumValues()];

	// When open level sets are supported, the value assignment becomes a bit more complicated
	if (Traits::openLevelSetsSupport)
	{
	    Data prevVal;
	    Index prevX, prevY;

	    iter2 = dtgrid2.beginTubeIterator();
	    iend2 = dtgrid2.endTubeIterator();
	    iter = beginTubeIterator();
	    iend = endTubeIterator();

	    if (iter2!=iend2)
	    {
		val = prevVal = iter2.getValue();
	    }

	    while (iter2 != iend2)
	    {
		iter2.getIndex(&x,&y,&z);
		iter.getIndex(&x2,&y2,&z2);

		if (x==x2 && y==y2 && z==z2)
		{
		    iter.setValue(val);
		    (*permutationArray)[iter.getArrayIndex()] = tmpPermutationArray[iter2.getArrayIndex()];
		    iter++;
		    iter2++;
		    prevVal = val;
		    prevX = x2;
		    prevY = y2;
		    val = iter2.getValue();
		}
		else
		{
		    // check if iter is in the same z-column as iter2.
		    // if this is the case we should use val as padding value,
		    // otherwise we should use prevVal as padding value until iter
		    // changes column.
		    if (x==x2 && y==y2)
		    {
			// iter and iter2 are in the same column.
			// in this case we just use val as padding value.
			if (val < 0)
			{
			    vg = -gamma;
			    permutationCurrent = permutationInside;
			}
			else
			{
			    vg = gamma;
			    permutationCurrent = permutationOutside;
			}
			while (z!=z2)  // We can do this optimization here since we know that dtgrid2 is contained in this dtgrid
			{
			    // pad with +/- gamma
			    iter.setValue(vg);
			    (*permutationArray)[iter.getArrayIndex()] = permutationCurrent;
			    iter++;
			    z2 = iter.getK();
			}

			// iter and iter2 are now equal
			// set value
			iter.setValue(val);
			(*permutationArray)[iter.getArrayIndex()] = tmpPermutationArray[iter2.getArrayIndex()];
			iter++;
			iter2++;
			prevVal = val;
			prevX = x2;
			prevY = y2;
			val = iter2.getValue();
		    }
		    else if (x==x2)
		    {
			// iter and iter2 are not in the same column.
			// we should therefore use prevVal as padding value until iter
			// changes column.

			// we only use prevVal if iter did not change column itself
			if (iter.getI()==prevX && iter.getJ()==prevY)
			{
			    if (val < 0)
			    {
				vg = -gamma;
				permutationCurrent = permutationInside;
			    }
			    else
			    {
				vg = gamma;
				permutationCurrent = permutationOutside;
			    }
			    while (iter.getI()==x2 && iter.getJ()==y2)
			    {
				// pad with +/- gamma
				iter.setValue(vg);
				(*permutationArray)[iter.getArrayIndex()] = permutationCurrent;
				iter++;
			    }
			}
			if (val < 0)
			{
			    vg = -gamma;
			    permutationCurrent = permutationInside;
			}
			else
			{
			    vg = gamma;
			    permutationCurrent = permutationOutside;
			}
			iter.getIndex(&x2, &y2, &z2); 			
			while (x!=x2 || y!=y2 || z!=z2)
			{
			    // pad with +/- gamma
			    iter.setValue(vg);
			    (*permutationArray)[iter.getArrayIndex()] = permutationCurrent;
			    iter++;
			    iter.getIndex(&x2, &y2, &z2); 
			}

			// iter and iter2 are now equal
			// set value
			iter.setValue(val);
			(*permutationArray)[iter.getArrayIndex()] = tmpPermutationArray[iter2.getArrayIndex()];
			iter++;
			iter2++;
			prevVal = val;
			prevX = x2;
			prevY = y2;
			val = iter2.getValue();
		    }
		    else
		    {
			// iter and iter2 are not in the same column.
			// we should therefore use prevVal as padding value until iter
			// changes column.
			if (val < 0)
			{
			    vg = -gamma;
			    permutationCurrent = permutationInside;
			}
			else
			{
			    vg = gamma;
			    permutationCurrent = permutationOutside;
			}
			while (iter.getI()==x2 )
			{
			    // pad with +/- gamma
			    iter.setValue(vg);
			    (*permutationArray)[iter.getArrayIndex()] = permutationCurrent;
			    iter++;
			}
			if (val < 0)
			{
			    vg = -gamma;
			    permutationCurrent = permutationInside;
			}
			else
			{
			    vg = gamma;
			    permutationCurrent = permutationOutside;
			}
			iter.getIndex(&x2, &y2, &z2); 			
			while (x!=x2 || y!=y2 || z!=z2)
			{
			    // pad with +/- gamma
			    iter.setValue(vg);
			    (*permutationArray)[iter.getArrayIndex()] = permutationCurrent;
			    iter++;
			    iter.getIndex(&x2, &y2, &z2); 
			}

			// iter and iter2 are now equal
			// set value
			iter.setValue(val);
			(*permutationArray)[iter.getArrayIndex()] = tmpPermutationArray[iter2.getArrayIndex()];
			iter++;
			iter2++;
			prevVal = val;
			prevX = x2;
			prevY = y2;
			val = iter2.getValue();
		    }
		}
	    }	    
	    while (iter != iend)
	    {
		iter.setValue(vg);
		(*permutationArray)[iter.getArrayIndex()] = permutationCurrent;
		iter++;
	    }
	}
	else
	{
	    iter2 = dtgrid2.beginTubeIterator();
	    iend2 = dtgrid2.endTubeIterator();
	    iter = beginTubeIterator();
	    iend = endTubeIterator();
	    while (iter2 != iend2)
	    {
		iter2.getIndex(&x,&y,&z);
		val = iter2.getValue();

		if (val < 0)
		{
		    vg = -gamma;
		    permutationCurrent = permutationInside;
		}
		else
		{
		    vg = gamma;
		    permutationCurrent = permutationOutside;
		}

		x2 = iter.getI();
		while (x!=x2)
		{
		    // pad with +/- gamma
		    iter.setValue(vg);
		    (*permutationArray)[iter.getArrayIndex()] = permutationCurrent;
		    iter++;
		    x2 = iter.getI();
		}
		y2 = iter.getJ();
		while (y!=y2)
		{
		    // pad with +/- gamma
		    iter.setValue(vg);
		    (*permutationArray)[iter.getArrayIndex()] = permutationCurrent;
		    iter++;
		    y2 = iter.getJ();
		}
		z2 = iter.getK();
		while (z!=z2)
		{
		    // pad with +/- gamma
		    iter.setValue(vg);
		    (*permutationArray)[iter.getArrayIndex()] = permutationCurrent;
		    iter++;
		    z2 = iter.getK();
		}

		// set value
		iter.setValue(val);
		(*permutationArray)[iter.getArrayIndex()] = tmpPermutationArray[iter2.getArrayIndex()];
		iter++;  
		iter2++;
	    }	    
	    while (iter != iend)
	    {
		iter.setValue(vg);
		(*permutationArray)[iter.getArrayIndex()] = permutationCurrent;
		iter++;
	    }
	}
    }
    else
    {
	clear();
    }
}





///////////////////////////////////////////////////////////////////////////////////
// TRANSFORMATION OPERATIONS
///////////////////////////////////////////////////////////////////////////////////



template<class Traits>
void DTGrid<Traits>::setScale(Real scale)
{
    TubeIterator iter, iend;
    Real rescale = scale / dx;

    dx = scale;
    gamma *= rescale;
    beta *= rescale;
    insideConstant *= rescale;
    outsideConstant *= rescale;

    iter = beginTubeIterator();
    iend = endTubeIterator();

    while (iter != iend)
    {
	iter.setValue( iter.getValue() * rescale );
	iter++;
    }

    gridToWorld = Matrix::Matrix4x4<Real>::translation(tx,ty,tz) * Matrix::Matrix4x4<Real>::rotationXYZ(rx,ry,rz) * Matrix::Matrix4x4<Real>::scale(dx);
    worldToGrid = Matrix::Matrix4x4<Real>::scale(1/dx) * Matrix::Matrix4x4<Real>::rotationXYZ(rx,ry,rz).transpose() * Matrix::Matrix4x4<Real>::translation(-tx,-ty,-tz);
}

template<class Traits>
void DTGrid<Traits>::setTranslation(Real tx, Real ty, Real tz)
{
    this->tx = tx;
    this->ty = ty;
    this->tz = tz;

    gridToWorld = Matrix::Matrix4x4<Real>::translation(tx,ty,tz) * Matrix::Matrix4x4<Real>::rotationXYZ(rx,ry,rz) * Matrix::Matrix4x4<Real>::scale(dx);
    worldToGrid = Matrix::Matrix4x4<Real>::scale(1/dx) * Matrix::Matrix4x4<Real>::rotationXYZ(rx,ry,rz).transpose() * Matrix::Matrix4x4<Real>::translation(-tx,-ty,-tz);
}


template<class Traits>
void DTGrid<Traits>::translateGrid(Index t[3])
{
    UInt i;

    for (i=0; i<numXIndex; i++)
    {
	xIndex[i] += t[0];
    }

    for (i=0; i<numYIndex; i++)
    {
	yIndex[i] += t[1];
    }

    for (i=0; i<numZIndex; i++)
    {
	zIndex[i] += t[2];
    }

    bbox[0][0] += t[0];
    bbox[0][1] += t[0];
    bbox[1][0] += t[1];
    bbox[1][1] += t[1];
    bbox[2][0] += t[2];
    bbox[2][1] += t[2];

    if (Traits::openLevelSetsSupport)
    {
	if (openLevelSetBBoxValid)
	{
	    openLevelSetBBox[0][0] += t[0];
	    openLevelSetBBox[1][0] += t[1];
	    openLevelSetBBox[2][0] += t[2];
	    openLevelSetBBox[0][1] += t[0];
	    openLevelSetBBox[1][1] += t[1];
	    openLevelSetBBox[2][1] += t[2];
	}
    }
}


template<class Traits>
void DTGrid<Traits>::setRotationXYZ(Real rx, Real ry, Real rz)
{
    this->rx = rx;
    this->ry = ry;
    this->rz = rz;

    gridToWorld = Matrix::Matrix4x4<Real>::translation(tx,ty,tz) * Matrix::Matrix4x4<Real>::rotationXYZ(rx,ry,rz) * Matrix::Matrix4x4<Real>::scale(dx);
    worldToGrid = Matrix::Matrix4x4<Real>::scale(1/dx) * Matrix::Matrix4x4<Real>::rotationXYZ(rx,ry,rz).transpose() * Matrix::Matrix4x4<Real>::translation(-tx,-ty,-tz);
}

template<class Traits>
typename DTGrid<Traits>::Real DTGrid<Traits>::getScale()
{
    return dx;
}

template<class Traits>
template<typename Real2> 
void DTGrid<Traits>::getTranslation(Real2 *tx, Real2 *ty, Real2 *tz)
{
    *tx = this->tx;
    *ty = this->ty;
    *tz = this->tz;
}

template<class Traits>
template<typename Real2> 
void DTGrid<Traits>::getRotationXYZ(Real2 *rx, Real2 *ry, Real2 *rz)
{
    *rx = this->rx;
    *ry = this->ry;
    *rz = this->rz;
}




///////////////////////////////////////////////////////////////////////////////////
// CSG OPERATIONS
///////////////////////////////////////////////////////////////////////////////////



// This CSG operation takes a DTGrid, dtgrid2, as input and takes into account the transformation from
// this DTGrid to dtgrid2. 
// 
// TemplateInterpolator requirements:
// Data interp(Data values[8], Matrix::Vector3<Real> unitOffset)
template<class Traits> template <class TemplateInterpolator>
void DTGrid<Traits>::CSG(const DTGrid& dtgrid2, CSGType csgType, TemplateInterpolator *interpolator)
{
    TransformationIterator< TemplateInterpolator > tIter, tIterEnd;
    Index bbox2[3][2];
    Index bbox2local[3][2];
    Matrix::Matrix4x4<Real> grid2ToThis = worldToGrid * dtgrid2.getGridToWorld();
    Matrix::Matrix4x4<Real> thisToGrid2 = dtgrid2.getWorldToGrid() * gridToWorld;
    Matrix::Vector3<Real> vec3, minVec, maxVec;
    bool initialized = false;
    UInt i,j,k,m;

    if (Traits::openLevelSetsSupport)
    {
	if (openLevelSetBBoxValid)
	{
	    Core::throwDefaultException("CSG operations including transformations and open level sets not well-defined!", __FILE__, __LINE__);
	}
    }


    // 1: Compute the axis-aligned bounding box of dtgrid2 in the grid coordinate system of this DTGrid

    //dtgrid2.boundingBox(bbox2);
    // this is safer if bounding boxes are not valid
    dtgrid2.effectiveGridSize(bbox2);
    initBoundingBox();


    for (i=0; i<2; i++)
    {
	for (j=0; j<2; j++)
	{
	    for (k=0; k<2; k++)
	    {
		if (!initialized)
		{
		    maxVec = minVec = grid2ToThis * Matrix::Vector3<Real>(bbox2[0][0], bbox2[1][0], bbox2[2][0]);
		    initialized = true;
		}
		else
		{
		    vec3 = grid2ToThis * Matrix::Vector3<Real>(bbox2[0][i], bbox2[1][j], bbox2[2][k]);

		    for (m=0; m<3; m++)
		    {
			if (vec3[m] < minVec[m])
			{
			    minVec[m] = vec3[m];
			}
			else if (vec3[m] > maxVec[m])
			{
			    maxVec[m] = vec3[m];
			}
		    }

		}
	    }
	}
    }

    // clamp to grid coordinates
    minVec = Matrix::Vector3<Real>::floor(minVec);
    maxVec = Matrix::Vector3<Real>::ceil(maxVec);

    // construct bounding box in grid coordinates
    bbox2local[0][0] = (Index) minVec[0];
    bbox2local[0][1] = (Index) maxVec[0];
    bbox2local[1][0] = (Index) minVec[1];
    bbox2local[1][1] = (Index) maxVec[1];
    bbox2local[2][0] = (Index) minVec[2];
    bbox2local[2][1] = (Index) maxVec[2];


    // check if the bounding boxes do not overlap!
    for (i=0; i<3; i++)
    {
	if (bbox2local[i][0] > bbox[i][1] || bbox2local[i][1] < bbox[i][0])
	{
	    if ( csgType == CSG_DIFFERENCE )
	    {
		// the bounding boxes do not overlap, so the difference is equal to this dtgrid
		return;
	    }
	    else if ( csgType == CSG_INTERSECTION )
	    {
		// the bounding boxes do not overlap, so the intersection is equal to the empty grid
		clear();
		return;
	    }
	    // if csgType == CSG_UNION, the result will be both grids combined even though the bounding boxes
	    // do not overlap, so we do not return here!!
	}
    }

    // 2: Initialize transformation iterators

    tIter = dtgrid2.template beginTransformationIterator< TemplateInterpolator >(bbox2local, thisToGrid2, interpolator, outsideConstant);
    tIterEnd = dtgrid2.template endTransformationIterator< TemplateInterpolator >(bbox2local);


    // 3: Do CSG operation

    CSG_GA(tIter, tIterEnd, csgType);
}


template<class Traits>
void DTGrid<Traits>::CSG_GA(DTGrid& dtgrid2, CSGType csgType)
{
    StencilTubeIterator<ENTIRE_TUBE> iter, iend;

    iter = dtgrid2.beginStencilTubeIterator<ENTIRE_TUBE, SF_NONE, false>();
    iend = dtgrid2.endStencilTubeIterator<ENTIRE_TUBE, SF_NONE, false>();

    if (Traits::openLevelSetsSupport)
    {
	if (openLevelSetBBoxValid)
	{
	    if ( 
		(dtgrid2.openLevelSetBBoxValid &&
		(openLevelSetBBox[0][0]!=dtgrid2.openLevelSetBBox[0][0] ||
		openLevelSetBBox[1][0]!=dtgrid2.openLevelSetBBox[1][0] ||
		openLevelSetBBox[2][0]!=dtgrid2.openLevelSetBBox[2][0] ||
		openLevelSetBBox[0][1]!=dtgrid2.openLevelSetBBox[0][1] ||
		openLevelSetBBox[1][1]!=dtgrid2.openLevelSetBBox[1][1] ||
		openLevelSetBBox[2][1]!=dtgrid2.openLevelSetBBox[2][1])) ||
		(!dtgrid2.openLevelSetBBoxValid &&
		(openLevelSetBBox[0][0]>dtgrid2.openLevelSetBBox[0][0] ||
		openLevelSetBBox[1][0]>dtgrid2.openLevelSetBBox[1][0] ||
		openLevelSetBBox[2][0]>dtgrid2.openLevelSetBBox[2][0] ||
		openLevelSetBBox[0][1]<dtgrid2.openLevelSetBBox[0][1] ||
		openLevelSetBBox[1][1]<dtgrid2.openLevelSetBBox[1][1] ||
		openLevelSetBBox[2][1]<dtgrid2.openLevelSetBBox[2][1]))
		)
	    {
		Core::throwDefaultException("CSG operations including open level sets require openLevelSetBBox to agree, otherwise not well-defined!", __FILE__, __LINE__);
	    }
	}
    }


    CSG_GA(iter, iend, csgType);
}



// This CSG operation takes a TemplateIterator pair (iter2 and iend2) that iterates over a grid
// in the same grid coordinate system as this DTGrid.
// 
// TemplateIterator requirements:  
// bool operator++(), bool operator++(int), bool operator!=(), Data getValue(), getIndex(Index *x, Index *y, Index *z)
// 
template<class Traits> template <class TemplateIterator> 
void DTGrid<Traits>::CSG_GA(TemplateIterator iter2, TemplateIterator iend2, CSGType csgType)
{
    InitParams initParams = InitParams(dx, beta, gamma, insideConstant, outsideConstant);
    DTGrid<Traits> dtgridRes = DTGrid<Traits>(initParams);
    dtgridRes.setTranslation(tx,ty,tz);
    dtgridRes.setRotationXYZ(rx,ry,rz);
    StencilTubeIterator<ENTIRE_TUBE> iter1, iend1;


    dtgridRes.beginSafePush();

    iter1 = beginStencilTubeIterator<ENTIRE_TUBE, SF_NONE, false>();
    iend1 = endStencilTubeIterator<ENTIRE_TUBE, SF_NONE, false>();


    if (Traits::openLevelSetsSupport)
    {
	/*
	There are the following cases:
	(We assume that iter1 is lexicographically smaller than iter2)

	- Iter2 is in the same column and is lexicographically larger: 
	-> Use the sign of the value at iter2.
	- Iter2 is not in the same column and is lexicographically larger:
	-> If iter2 is not in the same zy slice, use iter2's prevval for the
	rest of this zy slice (since ls2's interface is not intersected by a curve).
	And start to use iter2's currentval when iter1 changes slice because then we know
	that ls2 is not intersected by a curve from iter1's current position to the position of iter2.
	-> If iter2 is in the same zy slice, use iter2's prevval for the rest of the column.
	And start to use iter2's currentval when iter1 changes column (same argument as above).
	- Iter2 has no more elements: 
	-> Use the sign of the last element referenced by iter2  
	*/

	Data prevVal1, prevVal2;    
	Data v1, v2;

	if (iter1 != iend1)
	{
	    v1 = iter1.getValue();
	}
	else
	{
	    // default is outside (positive) if no elements
	    v1 = static_cast<Data>(1);
	}


	if (iter2 != iend2)
	{
	    v2 = iter2.getValue();
	}
	else
	{
	    // default is outside (positive) if no elements
	    v2 = static_cast<Data>(1);
	}


	while (iter1 != iend1 || iter2 != iend2)
	{
	    Index i1[3], i2[3];

	    if (iter1 == iend1)
	    {
		v2 = iter2.getValue();

		// continue until iter2 == iend2
		while (iter2 != iend2)
		{
		    iter2.getIndex(&i2[0], &i2[1], &i2[2]);    
		    v2 = iter2.getValue();

		    if (csgType == CSG_UNION)
		    {
			if ( v1 > 0 )
			{
			    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
			}
		    }
		    else if (csgType == CSG_INTERSECTION)
		    {
			if ( v1 < 0 )
			{
			    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
			}
		    }
		    else if (csgType == CSG_DIFFERENCE)
		    {
			if ( v1 < 0 )
			{
			    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], -v2 );
			}
		    }

		    iter2++;
		}
	    }
	    else if (iter2 == iend2)
	    {
		v1 = iter1.getValue();

		// continue until iter1 == iend1
		while (iter1 != iend1)
		{
		    iter1.getIndex(&i1[0], &i1[1], &i1[2]);    
		    v1 = iter1.getValue();

		    if (csgType == CSG_UNION)
		    {
			if ( v2 > 0 )
			{
			    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
			}
		    }
		    else if (csgType == CSG_INTERSECTION)
		    {
			if ( v2 < 0 )
			{
			    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
			}
		    }
		    else if (csgType == CSG_DIFFERENCE)
		    {
			if ( v2 > 0 )
			{
			    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
			}
		    }

		    iter1++;
		    v1 = iter1.getValue();
		}
	    }
	    else
	    {	
		iter1.getIndex(&i1[0], &i1[1], &i1[2]);    
		iter2.getIndex(&i2[0], &i2[1], &i2[2]);
		prevVal1 = v1;
		prevVal2 = v2;
		v1 = iter1.getValue();
		v2 = iter2.getValue();


		// determine cases:

		if (i1[0] == i2[0])
		{
		    if (i1[1] == i2[1])
		    {
			if (i1[2] == i2[2])
			{
			    // equal
			    if (csgType == CSG_UNION)
			    {
				dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], min( v1, v2 ));
			    }
			    else if (csgType == CSG_INTERSECTION)
			    {
				dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], max( v1, v2 ));
			    }
			    else if (csgType == CSG_DIFFERENCE)
			    {
				dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], max( v1, -v2 ));
			    }

			    iter1++;
			    iter2++;
			}
			else
			{
			    // not equal, but in same column
			    if (i1[2] < i2[2])
			    {
				// iter1 smallest
				while (iter1!=iend1 && (i1[0] == i2[0] && i1[1] == i2[1] && i1[2] < i2[2]))
				{
				    v1 = iter1.getValue();
				    if (csgType == CSG_UNION)
				    {
					if ( v2 > 0 )
					{
					    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
					}
				    }
				    else if (csgType == CSG_INTERSECTION)
				    {
					if ( v2 < 0 )
					{
					    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
					}
				    }
				    else if (csgType == CSG_DIFFERENCE)
				    {
					if ( v2 > 0 )
					{
					    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
					}
				    }

				    iter1++;
				    if (iter1 != iend1)
				    {
					iter1.getIndex(&i1[0], &i1[1], &i1[2]);    
				    }
				}
			    }
			    else
			    {
				// iter2 smallest
				while (iter2!=iend2 && (i1[0] == i2[0] && i1[1] == i2[1] && i2[2] < i1[2]))
				{
				    v2 = iter2.getValue();
				    if (csgType == CSG_UNION)
				    {
					if ( v1 > 0 )
					{
					    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
					}
				    }
				    else if (csgType == CSG_INTERSECTION)
				    {
					if ( v1 < 0 )
					{
					    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
					}
				    }
				    else if (csgType == CSG_DIFFERENCE)
				    {
					if ( v1 < 0 )
					{
					    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], -v2 );
					}
				    }

				    iter2++;
				    if (iter2!=iend2)
				    {
					iter2.getIndex(&i2[0], &i2[1], &i2[2]);    
				    }
				}
			    }
			}
		    }
		    else
		    {
			// not equal, and not in same column
			// but in same slice 
			if (i1[1] < i2[1])
			{
			    // Use prevVal2 for the rest of the column and then switch to v2
			    Data v2Active = prevVal2;

			    // iter1 smallest
			    while (iter1!=iend1 && (i1[0] == i2[0] && (i1[1] < i2[1] || (i1[1] == i2[1] && i1[2] < i2[2]))))
			    {
				v1 = iter1.getValue();
				if (csgType == CSG_UNION)
				{
				    if ( v2Active > 0 )
				    {
					dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
				    }
				}
				else if (csgType == CSG_INTERSECTION)
				{
				    if ( v2Active < 0 )
				    {
					dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
				    }
				}
				else if (csgType == CSG_DIFFERENCE)
				{
				    if ( v2Active > 0 )
				    {
					dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
				    }
				}

				iter1++;
				if (iter1 != iend1)
				{
				    // did we switch column?
				    if (iter1.getJ() != i1[1])
				    {
					v2Active = v2;
				    }
				    iter1.getIndex(&i1[0], &i1[1], &i1[2]);
				}
			    }
			}
			else
			{
			    // Use prevVal1 for the rest of the column and then switch to v1
			    Data v1Active = prevVal1;

			    // iter2 smallest
			    while (iter2!=iend2 && (i1[0] == i2[0] && (i2[1] < i1[1] || (i1[1] == i2[1] && i2[2] < i1[2]))))
			    {
				v2 = iter2.getValue();
				if (csgType == CSG_UNION)
				{
				    if ( v1Active > 0 )
				    {
					dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
				    }
				}
				else if (csgType == CSG_INTERSECTION)
				{
				    if ( v1Active < 0 )
				    {
					dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
				    }
				}
				else if (csgType == CSG_DIFFERENCE)
				{
				    if ( v1Active < 0 )
				    {
					dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], -v2 );
				    }
				}

				iter2++;
				if (iter2!=iend2)
				{
				    // did we switch column?
				    if (iter2.getJ() != i2[1])
				    {
					v1Active = v1;
				    }
				    iter2.getIndex(&i2[0], &i2[1], &i2[2]);
				}
			    }
			}
		    }
		}
		else
		{
		    // not equal, and not in same slice
		    if (i1[0] < i2[0])
		    {
			// iter1 smallest

			// Use prevVal2 for the rest of the slice and then switch to v2
			Data v2Active = prevVal2;

			// iter1 smallest
			while (iter1!=iend1 && (i1[0] < i2[0] || (i1[0] == i2[0] && (i1[1] < i2[1] || (i1[1] == i2[1] && i1[2] < i2[2])))))
			{
			    v1 = iter1.getValue();
			    if (csgType == CSG_UNION)
			    {
				if ( v2Active > 0 )
				{
				    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
				}
			    }
			    else if (csgType == CSG_INTERSECTION)
			    {
				if ( v2Active < 0 )
				{
				    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
				}
			    }
			    else if (csgType == CSG_DIFFERENCE)
			    {
				if ( v2Active > 0 )
				{
				    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
				}
			    }

			    iter1++;
			    if (iter1 != iend1)
			    {
				// did we switch slice?
				if (iter1.getI() != i1[0])
				{
				    v2Active = v2;
				}
				iter1.getIndex(&i1[0], &i1[1], &i1[2]);
			    }
			}
		    }
		    else
		    {
			// iter2 smallest

			// Use prevVal1 for the rest of the slice and then switch to v1
			Data v1Active = prevVal1;

			while (iter2 != iend2 && (i2[0] < i1[0] || (i1[0] == i2[0] && (i2[1] < i1[1] || (i1[1] == i2[1] && i2[2] < i1[2])))))
			{
			    v2 = iter2.getValue();
			    if (csgType == CSG_UNION)
			    {
				if ( v1Active > 0 )
				{
				    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
				}
			    }
			    else if (csgType == CSG_INTERSECTION)
			    {
				if ( v1Active < 0 )
				{
				    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
				}
			    }
			    else if (csgType == CSG_DIFFERENCE)
			    {
				if ( v1Active < 0 )
				{
				    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], -v2 );
				}
			    }

			    iter2++;
			    if (iter2 != iend2)
			    {
				// did we switch slice?
				if (iter2.getI() != i2[0])
				{
				    v1Active = v1;
				}
				iter2.getIndex(&i2[0], &i2[1], &i2[2]);
			    }
			}
		    }
		}
	    }
	}
    }
    else
    {
	while (iter1 != iend1 || iter2 != iend2)
	{
	    Index i1[3], i2[3];
	    UInt i;
	    Data v1, v2;

	    if (iter1 == iend1)
	    {
		iter2.getIndex(&i2[0], &i2[1], &i2[2]);
		i1[0] = i2[0]+1;  // force i1 to be lexicographically larger
		v2 = iter2.getValue();
		v1 = outsideConstant;		
	    }
	    else if (iter2 == iend2)
	    {
		iter1.getIndex(&i1[0], &i1[1], &i1[2]);
		i2[0] = i1[0]+1;  // force i2 to be lexicographically larger
		v1 = iter1.getValue();
		v2 = outsideConstant;
	    }
	    else
	    {	
		iter1.getIndex(&i1[0], &i1[1], &i1[2]);    
		iter2.getIndex(&i2[0], &i2[1], &i2[2]);    
		v1 = iter1.getValue();
		v2 = iter2.getValue();
	    }



	    // next determine the lexicographically smallest element of the two iterators

	    for (i=0; i<3; i++)
	    {
		if ( i1[i] < i2[i] )
		{
		    // i1 smallest
		    if (csgType == CSG_UNION)
		    {
			if ( v2 > 0 )
			{
			    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
			}
		    }
		    else if (csgType == CSG_INTERSECTION)
		    {
			if ( v2 < 0 )
			{
			    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
			}
		    }
		    else if (csgType == CSG_DIFFERENCE)
		    {
			if ( v2 > 0 )
			{
			    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
			}
		    }
		    iter1++;
		    break;
		}
		else if ( i1[i] > i2[i] )
		{
		    // i2 smallest
		    if (csgType == CSG_UNION)
		    {
			if ( v1 > 0 )
			{
			    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
			}
		    }
		    else if (csgType == CSG_INTERSECTION)
		    {
			if ( v1 < 0 )
			{
			    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
			}
		    }
		    else if (csgType == CSG_DIFFERENCE)
		    {
			if ( v1 < 0 )
			{
			    dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], -v2 );
			}
		    }
		    iter2++;
		    break;
		}
	    }



	    if ( i==3 )
	    {
		// the two grid points were equal
		if (csgType == CSG_UNION)
		{
		    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], min( v1, v2 ));
		}
		else if (csgType == CSG_INTERSECTION)
		{
		    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], max( v1, v2 ));
		}
		else if (csgType == CSG_DIFFERENCE)
		{
		    dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], max( v1, -v2 ));
		}

		iter1++;
		iter2++;
	    }
	}
    }

    dtgridRes.endSafePush();
    clear();
    assign(&dtgridRes);
    dtgridRes.clear(false);
}





template<class Traits>
template<typename FuncT1, typename FuncT2>
void DTGrid<Traits>::CSG_GA(DTGrid& dtgrid2, CSGType csgType, FuncT1 *func1, FuncT2 *func2)
{
    StencilTubeIterator<ENTIRE_TUBE> iter, iend;

    if (Traits::openLevelSetsSupport)
    {
	if (openLevelSetBBoxValid)
	{
	    if ( 
		(dtgrid2.openLevelSetBBoxValid &&
		(openLevelSetBBox[0][0]!=dtgrid2.openLevelSetBBox[0][0] ||
		openLevelSetBBox[1][0]!=dtgrid2.openLevelSetBBox[1][0] ||
		openLevelSetBBox[2][0]!=dtgrid2.openLevelSetBBox[2][0] ||
		openLevelSetBBox[0][1]!=dtgrid2.openLevelSetBBox[0][1] ||
		openLevelSetBBox[1][1]!=dtgrid2.openLevelSetBBox[1][1] ||
		openLevelSetBBox[2][1]!=dtgrid2.openLevelSetBBox[2][1])) ||
		(!dtgrid2.openLevelSetBBoxValid &&
		(openLevelSetBBox[0][0]>dtgrid2.openLevelSetBBox[0][0] ||
		openLevelSetBBox[1][0]>dtgrid2.openLevelSetBBox[1][0] ||
		openLevelSetBBox[2][0]>dtgrid2.openLevelSetBBox[2][0] ||
		openLevelSetBBox[0][1]<dtgrid2.openLevelSetBBox[0][1] ||
		openLevelSetBBox[1][1]<dtgrid2.openLevelSetBBox[1][1] ||
		openLevelSetBBox[2][1]<dtgrid2.openLevelSetBBox[2][1]))
		)
	    {
		Core::throwDefaultException("CSG operations including open level sets require openLevelSetBBox to agree, otherwise not well-defined!", __FILE__, __LINE__);
	    }
	}
    }

    iter = dtgrid2.beginStencilTubeIterator<ENTIRE_TUBE, SF_NONE, false>();
    iend = dtgrid2.endStencilTubeIterator<ENTIRE_TUBE, SF_NONE, false>();

    CSG_GA(iter, iend, csgType, func1, func2);
}


template<class Traits> template <class TemplateIterator, typename FuncT1, typename FuncT2> 
void DTGrid<Traits>::CSG_GA(TemplateIterator iter2, TemplateIterator iend2, CSGType csgType, FuncT1 *func1, FuncT2 *func2)
{
    InitParams initParams = InitParams(dx, beta, gamma, insideConstant, outsideConstant);
    DTGrid<Traits> dtgridRes = DTGrid<Traits>(initParams);
    dtgridRes.setTransformation(tx,ty,tz);
    dtgridRes.setRotationXYZ(rx,ry,rz);
    StencilTubeIterator<ENTIRE_TUBE> iter1, iend1;


    dtgridRes.beginSafePush();

    iter1 = beginStencilTubeIterator<ENTIRE_TUBE, SF_NONE, false>();
    iend1 = endStencilTubeIterator<ENTIRE_TUBE, SF_NONE, false>();


    while (iter1 != iend1 || iter2 != iend2)
    {
	Index i1[3], i2[3];
	UInt i;
	Data v1, v2;

	if (iter1 == iend1)
	{
	    iter2.getIndex(&i2[0], &i2[1], &i2[2]);
	    i1[0] = i2[0]+1;  // force i2 to be lexicographically larger
	    v2 = (*func2)(iter2.getValue());
	    v1 = (*func1)(outsideConstant);		
	}
	else if (iter2 == iend2)
	{
	    iter1.getIndex(&i1[0], &i1[1], &i1[2]);
	    i2[0] = i1[0]+1;  // force i1 to be lexicographically larger
	    v1 = (*func1)(iter1.getValue());
	    v2 = (*func2)(outsideConstant);
	}
	else
	{	
	    iter1.getIndex(&i1[0], &i1[1], &i1[2]);    
	    iter2.getIndex(&i2[0], &i2[1], &i2[2]);    
	    v1 = (*func1)(iter1.getValue());
	    v2 = (*func2)(iter2.getValue());
	}



	// next determine the lexicographically smallest element of the two iterators

	for (i=0; i<3; i++)
	{
	    if ( i1[i] < i2[i] )
	    {
		// i1 smallest
		if (csgType == CSG_UNION)
		{
		    if ( v2 > 0 )
		    {
			dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
		    }
		}
		else if (csgType == CSG_INTERSECTION)
		{
		    if ( v2 < 0 )
		    {
			dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
		    }
		}
		else if (csgType == CSG_DIFFERENCE)
		{
		    if ( v2 > 0 )
		    {
			dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], v1 );
		    }
		}
		iter1++;
		break;
	    }
	    else if ( i1[i] > i2[i] )
	    {
		// i2 smallest
		if (csgType == CSG_UNION)
		{
		    if ( v1 > 0 )
		    {
			dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
		    }
		}
		else if (csgType == CSG_INTERSECTION)
		{
		    if ( v1 < 0 )
		    {
			dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], v2 );
		    }
		}
		else if (csgType == CSG_DIFFERENCE)
		{
		    if ( v1 < 0 )
		    {
			dtgridRes.template push3D<true, true>(i2[0], i2[1], i2[2], -v2 );
		    }
		}
		iter2++;
		break;
	    }
	}



	if ( i==3 )
	{
	    // the two grid points were equal
	    if (csgType == CSG_UNION)
	    {
		dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], min( v1, v2 ));
	    }
	    else if (csgType == CSG_INTERSECTION)
	    {
		dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], max( v1, v2 ));
	    }
	    else if (csgType == CSG_DIFFERENCE)
	    {
		dtgridRes.template push3D<true, true>(i1[0], i1[1], i1[2], max( v1, -v2 ));
	    }

	    iter1++;
	    iter2++;
	}



    }

    dtgridRes.endSafePush();
    clear();
    assign(&dtgridRes);
    dtgridRes.clear(false);
}



///////////////////////////////////////////////////////////////////////////////////
// MISC METHODS
///////////////////////////////////////////////////////////////////////////////////

// random access

template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::operator()(const Matrix::Vector3<Real>& pos)
{
    Data values[8];
    Data value;
    Matrix::Vector3<Real> vec3 = worldToGrid * pos;
    Math::TrilinearInterpolator<Data, Real> interpolator = Math::TrilinearInterpolator<Data, Real>();

    // convert vec3 to a unit offset
    Matrix::Vector3<Index> minVec = Matrix::Vector3<Real>::floor(vec3);
    vec3 = vec3 - Matrix::Vector3<Real>::floor(vec3);

    // Ordering of corner values:
    //
    //     X:         X + 1:
    // (and Y to the right, Z up) 
    //
    //  1 -- 3       5 -- 7
    //  |    |       |    |
    //  0 -- 2       4 -- 6
    //

    // Get values at positions (x,y,z),(x,y,z+1),(x,y+1,z),(x,y+1,z+1),(x+1,y,z),(x+1,y,z+1),(x+1,y+1,z),(x+1,y+1,z+1) instantaneously
    getVoxelValues(minVec[0], minVec[1], minVec[2], values);

    // interpolate to the final value
    value = interpolator.interp(values, vec3);
    
    return value;
}

template<class Traits>
typename DTGrid<Traits>::Data DTGrid<Traits>::operator()(Index x, Index y, Index z) const
{
    if (Traits::randomAccessType == USE_LINEAR_SEARCH)
    {
	UInt ic, tmp, fi, li;
	Index tmpI;


	// SEARCH IN X
	fi = 0;
	li = lastXIndex;
	if ( x < xIndex[fi] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		UInt aa1DIndex = aa1D[fi>>1];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	}
	else if ( x > xIndex[li] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		UInt aa1DIndex = aa1D[li>>1] + xIndex[li] - xIndex[li-1];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	}
	ic = fi+1;
	while (x > xIndex[ic])
	{
	    ic += 2;
	}
	ic -= 1;
	if (x < (tmpI=xIndex[ic]))
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		//OpenLS:
		UInt aa1DIndex = aa1D[ic>>1];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	}
	tmp = aa1D[(ic)>>1] + x - tmpI;
	fi = va1D[tmp];
	li = va1D[tmp+1]-1;


	// SEARCH IN Y
	if ( y < yIndex[fi] || y > yIndex[li] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (y < yIndex[fi])
		{
		    UInt aa2DIndex = aa2D[fi>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa2DIndex = aa2D[li>>1] + (yIndex[li] - yIndex[li-1]);
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	}
	ic = fi+1;
	while (y > yIndex[ic])
	{
	    ic += 2;
	}
	ic -= 1;
	if (y < (tmpI=yIndex[ic]))
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		//OpenLS:
		UInt aa2DIndex = aa2D[ic>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	}
	tmp = aa2D[(ic)>>1] + y - tmpI;
	fi = va2D[tmp];
	li = va2D[tmp+1]-1;


	// SEARCH IN Z
	if ( z < zIndex[fi] || z > zIndex[li] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (z < zIndex[fi])
		{
		    UInt aa3DIndex = aa3D[fi>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa3DIndex = aa3D[li>>1] + (zIndex[li] - zIndex[li-1]);
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	}
	ic = fi+1;
	while (z > zIndex[ic])
	{
	    ic += 2;
	}
	ic -= 1;
	if (z >= (tmpI=zIndex[ic]))
	{
	    return va3D[ aa3D[(ic)>>1] + z - tmpI ];
	}
	else
	{
	    return va3D[ aa3D[(ic)>>1] ] < 0 ? insideConstant : outsideConstant;
	}
    }
    else //if (Traits::randomAccessType == USE_BINARY_SEARCH)
    {
	UInt ix, ix2, tmp, it, li, fi;


	// SEARCH IN X
	fi = 0;
	li = lastXIndex;
	if ( x < xIndex[fi] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		UInt aa1DIndex = aa1D[fi>>1];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	}
	else if ( x > xIndex[li] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		UInt aa1DIndex = aa1D[li>>1] + xIndex[li] - xIndex[li-1];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	}

	ix = 0;
	ix2 = (li-fi)+1;

	// binary search
	while (ix2 != ix+2)
	{
	    it = ((ix+ix2)>>2)<<1;  // make even

	    if (x < xIndex[it+fi])
	    {
		ix2 = it;
	    }
	    else
	    {
		ix = it;
	    }
	}

	// now iy points to the correct closest even index
	ix += fi;
	if ( x > xIndex[ix+1] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		//OpenLS:
		UInt aa1DIndex = aa1D[ix>>1] + xIndex[ix+1] - xIndex[ix];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	}

	tmp = aa1D[(ix)>>1] + x - xIndex[ix];
	fi = va1D[tmp];
	li = va1D[tmp+1]-1;



	// SEARCH IN Y
	if ( y < yIndex[fi] || y > yIndex[li] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (y < yIndex[fi])
		{
		    UInt aa2DIndex = aa2D[fi>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa2DIndex = aa2D[li>>1] + (yIndex[li] - yIndex[li-1]);
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	}

	ix = 0;
	ix2 = (li-fi)+1;

	// binary search
	while (ix2 != ix+2)
	{
	    it = ((ix+ix2)>>2)<<1;  // make even

	    if (y < yIndex[it+fi])
	    {
		ix2 = it;
	    }
	    else
	    {
		ix = it;
	    }
	}

	// now ix points to the correct closest even index
	ix += fi;
	if ( y > yIndex[ix+1] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		//OpenLS:
		UInt aa2DIndex = aa2D[ix>>1] + yIndex[ix+1] - yIndex[ix];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	}

	tmp = aa2D[(ix)>>1] + y - yIndex[ix];
	fi = va2D[tmp];
	li = va2D[tmp+1]-1;



	// SEARCH IN Z
	// inside min, max range of z indices in this column?
	if (z < zIndex[fi] || z > zIndex[li])
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		return outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (z < zIndex[fi])
		{
		    UInt aa3DIndex = aa3D[fi>>1];
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa3DIndex = aa3D[li>>1] + (zIndex[li] - zIndex[li-1]);
		    return (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	}
	else
	{
	    ix = 0;
	    ix2 = (li-fi)+1;

	    // binary search
	    while (ix2 != ix+2)
	    {
		it = ((ix+ix2)>>2)<<1;  // make even

		if (z < zIndex[it+fi])
		{
		    ix2 = it;
		}
		else
		{
		    ix = it;
		}
	    }

	    // now ix points to the correct closest even index
	    ix += fi;
	    if ( z <= zIndex[ix+1] )
	    {
		return va3D[ aa3D[(ix)>>1] + z - zIndex[ix] ];
	    }
	    else
	    {
		// note that the sign will always be taken from the same p-column
		return va3D[ aa3D[(ix+2)>>1] ] < 0 ? insideConstant : outsideConstant;
	    }
	}

    }
}



template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getVa2DIndex(Index x, Index y)
{
    int ix, iy, fi, li, tmp;

    fi = 0;
    li = lastXIndex;
    findIndex(x, xIndex, fi, li, &ix);
    tmp = aa1D[(ix)>>1] + x - xIndex[ix];
    fi = va1D[tmp];
    li = va1D[tmp+1]-1;
    findIndex(y, yIndex, fi, li, &iy);
    return aa2D[(iy)>>1] + y - yIndex[iy];
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getVa1DIndex(Index x)
{
    int ix, fi, li;

    fi = 0;
    li = lastXIndex;
    findIndex(x, xIndex, fi, li, &ix);
    return aa1D[(ix)>>1] + x - xIndex[ix];
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getVa2D(UInt i)
{
    return va2D[i];
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getVa1D(UInt i)
{
    return va1D[i];
}


template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getVa1D(Index x)
{
    int ix, fi, li;

    fi = 0;
    li = lastXIndex;
    findIndex(x, xIndex, fi, li, &ix);
    return va1D[ aa1D[(ix)>>1] + x - xIndex[ix] ];
}


/*! We assume that (x,y) exists in the grid */
template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getVa2D(Index x, Index y)
{
    int ix, iy, fi, li, tmp;

    fi = 0;
    li = lastXIndex;
    findIndex(x, xIndex, fi, li, &ix);
    tmp = aa1D[(ix)>>1] + x - xIndex[ix];
    fi = va1D[tmp];
    li = va1D[tmp+1]-1;
    findIndex(y, yIndex, fi, li, &iy);
    return va2D[ aa2D[(iy)>>1] + y - yIndex[iy] ];
}


// random access
template<class Traits>
bool DTGrid<Traits>::operator()(Index x, Index y, Index z, Data *val) const
{
    if (Traits::randomAccessType == USE_LINEAR_SEARCH)
    {
	UInt ic, tmp, fi, li;
	Index tmpI;

	// SEARCH IN X
	fi = 0;
	li = lastXIndex;

	if ( x < xIndex[fi] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		UInt aa1DIndex = aa1D[fi>>1];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		*val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	    return false;
	}
	else if ( x > xIndex[li] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		UInt aa1DIndex = aa1D[li>>1] + xIndex[li] - xIndex[li-1];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		*val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	    return false;
	}

	ic = fi+1;
	while (x > xIndex[ic])
	{
	    ic += 2;
	}
	ic -= 1;
	if (x < (tmpI=xIndex[ic]))
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
		return false;
	    }
	    else
	    {
		//OpenLS:
		UInt aa1DIndex = aa1D[ic>>1];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		*val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		return false;
	    }
	}
	tmp = aa1D[(ic)>>1] + x - tmpI;
	fi = va1D[tmp];
	li = va1D[tmp+1]-1;


	// SEARCH IN Y
	if ( y < yIndex[fi] || y > yIndex[li] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
		return false;
	    }
	    else
	    {
		// OpenLS:
		if (y < yIndex[fi])
		{
		    UInt aa2DIndex = aa2D[fi>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    *val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa2DIndex = aa2D[li>>1] + (yIndex[li] - yIndex[li-1]);
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    *val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		return false;
	    }
	}
	ic = fi+1;
	while (y > yIndex[ic])
	{
	    ic += 2;
	}
	ic -= 1;
	if (y < (tmpI=yIndex[ic]))
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
		return false;
	    }
	    else
	    {
		//OpenLS:
		UInt aa2DIndex = aa2D[ic>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		*val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		return false;
	    }
	}
	tmp = aa2D[(ic)>>1] + y - tmpI;
	fi = va2D[tmp];
	li = va2D[tmp+1]-1;


	// SEARCH IN Z
	if ( z < zIndex[fi] || z > zIndex[li] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (z < zIndex[fi])
		{
		    UInt aa3DIndex = aa3D[fi>>1];
		    *val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa3DIndex = aa3D[li>>1] + (zIndex[li] - zIndex[li-1]);
		    *val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	    return false;
	}
	ic = fi+1;
	while (z > zIndex[ic])
	{
	    ic += 2;
	}
	ic -= 1;
	if (z >= (tmpI=zIndex[ic]))
	{
	    *val = va3D[ aa3D[(ic)>>1] + z - tmpI ];
	    return true;
	}
	else
	{
	    *val = va3D[ aa3D[(ic)>>1] ] < 0 ? insideConstant : outsideConstant;
	    return false;
	}

    }
    else //if(Traits::randomAccessType == USE_BINARY_SEARCH)
    {
	UInt ix, ix2, tmp, it, li, fi;


	// SEARCH IN X
	fi = 0;
	li = lastXIndex;

	if ( x < xIndex[fi] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		UInt aa1DIndex = aa1D[fi>>1];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		*val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	    return false;
	}
	else if ( x > xIndex[li] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		UInt aa1DIndex = aa1D[li>>1] + xIndex[li] - xIndex[li-1];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		*val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	    return false;
	}

	ix = 0;
	ix2 = (li-fi)+1;

	// binary search
	while (ix2 != ix+2)
	{
	    it = ((ix+ix2)>>2)<<1;  // make even

	    if (x < xIndex[it+fi])
	    {
		ix2 = it;
	    }
	    else
	    {
		ix = it;
	    }
	}

	// now ix points to the correct closest even index
	ix += fi;
	if ( x > xIndex[ix+1] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
	    }
	    else
	    {
		//OpenLS:
		UInt aa1DIndex = aa1D[ix>>1] + xIndex[ix+1] - xIndex[ix];
		UInt va1DValue = va1D[aa1DIndex];
		UInt aa2DIndex = aa2D[va1DValue>>1];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		*val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	    return false;
	}

	tmp = aa1D[(ix)>>1] + x - xIndex[ix];
	fi = va1D[tmp];
	li = va1D[tmp+1]-1;



	// SEARCH IN Y
	if ( y < yIndex[fi] || y > yIndex[li] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (y < yIndex[fi])
		{
		    UInt aa2DIndex = aa2D[fi>>1];
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    *val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa2DIndex = aa2D[li>>1] + (yIndex[li] - yIndex[li-1]);
		    UInt va2DValue = va2D[aa2DIndex];
		    UInt aa3DIndex = aa3D[va2DValue>>1];
		    *val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	    return false;
	}

	ix = 0;
	ix2 = (li-fi)+1;

	// binary search
	while (ix2 != ix+2)
	{
	    it = ((ix+ix2)>>2)<<1;  // make even

	    if (y < yIndex[it+fi])
	    {
		ix2 = it;
	    }
	    else
	    {
		ix = it;
	    }
	}

	// now ix points to the correct closest even index
	ix += fi;
	if ( y > yIndex[ix+1] )
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
	    }
	    else
	    {
		//OpenLS:
		UInt aa2DIndex = aa2D[ix>>1] + yIndex[ix+1] - yIndex[ix];
		UInt va2DValue = va2D[aa2DIndex];
		UInt aa3DIndex = aa3D[va2DValue>>1];
		*val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
	    }
	    return false;
	}

	tmp = aa2D[(ix)>>1] + y - yIndex[ix];
	fi = va2D[tmp];
	li = va2D[tmp+1]-1;



	// SEARCH IN Z
	// inside min, max range of z indices in this column?
	if (z < zIndex[fi] || z > zIndex[li])
	{
	    if (!Traits::openLevelSetsSupport)
	    {
		*val = outsideConstant;
	    }
	    else
	    {
		// OpenLS:
		if (z < zIndex[fi])
		{
		    UInt aa3DIndex = aa3D[fi>>1];
		    *val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
		else
		{
		    UInt aa3DIndex = aa3D[li>>1] + (zIndex[li] - zIndex[li-1]);
		    *val = (va3D[aa3DIndex] < 0 ? insideConstant : outsideConstant);
		}
	    }
	    return false;
	}
	else
	{
	    ix = 0;
	    ix2 = (li-fi)+1;

	    // binary search
	    while (ix2 != ix+2)
	    {
		it = ((ix+ix2)>>2)<<1;  // make even

		if (z < zIndex[it+fi])
		{
		    ix2 = it;
		}
		else
		{
		    ix = it;
		}
	    }

	    // now ix points to the correct closest even index
	    ix += fi;
	    if ( z <= zIndex[ix+1] )
	    {
		*val = va3D[ aa3D[(ix)>>1] + z - zIndex[ix] ];
		return true;
	    }
	    else
	    {
		// note that the sign will always be taken from the same p-column
		*val = va3D[ aa3D[(ix+2)>>1] ] < 0 ? insideConstant : outsideConstant;
		return false;
	    }
	}
    }
}

template<class Traits>
void DTGrid<Traits>::boundingBoxIndices(const Matrix::Vector3<Real>& mins, const Matrix::Vector3<Real>& maxs, Index iMins[3], Index iMaxs[3])
{
	Matrix::Vector3<Real> fMins = worldToGrid * mins;
	Matrix::Vector3<Real> fMaxs = worldToGrid * maxs;

	Matrix::Vector3<Index> minVec = Matrix::Vector3<Real>::floor(fMins);
	Matrix::Vector3<Index> maxVec = Matrix::Vector3<Real>::ceil(fMaxs);

	for(int x = 0; x < 3; x++){
		iMins[x] = minVec[x] - 1;
		iMaxs[x] = maxVec[x] + 1;
	}
}

template<class Traits>
void DTGrid<Traits>::setOpenLevelSetBBox(Index openLevelSetBBox[3][2])
{
    openLevelSetBBoxValid = true;
    this->openLevelSetBBox[0][0] = openLevelSetBBox[0][0];
    this->openLevelSetBBox[0][1] = openLevelSetBBox[0][1];
    this->openLevelSetBBox[1][0] = openLevelSetBBox[1][0];
    this->openLevelSetBBox[1][1] = openLevelSetBBox[1][1];
    this->openLevelSetBBox[2][0] = openLevelSetBBox[2][0];
    this->openLevelSetBBox[2][1] = openLevelSetBBox[2][1];
}


template<class Traits> template<typename Index2>
void DTGrid<Traits>::boundingBox(Index2 bbox[3][2]) const
{
    bbox[0][0] = static_cast<Index2>(this->bbox[0][0]);
    bbox[0][1] = static_cast<Index2>(this->bbox[0][1]);
    bbox[1][0] = static_cast<Index2>(this->bbox[1][0]);
    bbox[1][1] = static_cast<Index2>(this->bbox[1][1]);
    bbox[2][0] = static_cast<Index2>(this->bbox[2][0]);
    bbox[2][1] = static_cast<Index2>(this->bbox[2][1]);
}

template<class Traits>
void DTGrid<Traits>::boundingBox(Index bbox[3][2]) const
{
    bbox[0][0] = this->bbox[0][0];
    bbox[0][1] = this->bbox[0][1];
    bbox[1][0] = this->bbox[1][0];
    bbox[1][1] = this->bbox[1][1];
    bbox[2][0] = this->bbox[2][0];
    bbox[2][1] = this->bbox[2][1];
}


template<class Traits>
void DTGrid<Traits>::boundingBoxDim(Index bboxDim[3])
{
    Index bbox[3][2];
    boundingBox(bbox);
    bboxDim[0] = bbox[0][1] - bbox[0][0] + 1;
    bboxDim[1] = bbox[1][1] - bbox[1][0] + 1;
    bboxDim[2] = bbox[2][1] - bbox[2][0] + 1;
}


template<class Traits>
void DTGrid<Traits>::effectiveGridSize(Index dim[3][2]) const
{
    // X
    if (numXIndex > 0)
    {
	dim[0][0] = xIndex[0];
	dim[0][1] = xIndex[lastXIndex];
    }
    else
    {
	dim[0][0] = 0;
	dim[0][1] = 0;
    }


    // Y
    if (numYIndex > 0)
    {
	// find effective grid size in Y
	UInt i;
	Index tmp;
	Index minY = yIndex[0];
	Index maxY = yIndex[0];
	for (i=1; i<numYIndex; i++)
	{
	    tmp = yIndex[i];
	    if (tmp > maxY)
	    {
		maxY = tmp;
	    }
	    else if (tmp < minY)
	    {   
		minY = tmp;
	    }   
	}
	dim[1][0] = minY;
	dim[1][1] = maxY;
    }
    else
    {
	dim[1][0] = 0;
	dim[1][1] = 0;
    }

    // Z
    if (numZIndex > 0)
    {
	// find effective grid size in Z
	UInt i;
	Index tmp;
	Index minZ = zIndex[0];
	Index maxZ = zIndex[0];
	for (i=1; i<numZIndex; i++)
	{
	    tmp = zIndex[i];
	    if (tmp > maxZ)
	    {
		maxZ = tmp;
	    }
	    else if (tmp < minZ)
	    {   
		minZ = tmp;
	    }   
	}
	dim[2][0] = minZ;
	dim[2][1] = maxZ;
    }
    else
    {
	dim[2][0] = 0;
	dim[2][1] = 0;
    }
}



template<class Traits>
DTGrid<Traits> *DTGrid<Traits>::copy(Index bbox[3][2])
{
    // TODO: Currently this is just a naive implementation that does a full traversal.
    //       Can be optimized by combining searches in the 1D, 2D and 3D components whenever
    //       the traversal reaches a boundary of the bbox

    Index i, j, k;
    TubeIterator iter, iend;


    InitParams initParams = InitParams(dx, beta, gamma, insideConstant, outsideConstant);
    DTGrid<Traits> *copy = new DTGrid<Traits>(initParams);
    copy->beginSafePush();

    iter = beginTubeIterator();
    iend = endTubeIterator();

    while (iter != iend)
    {
	iter.getIndex(&i, &j, &k);
	if (i>=bbox[0][0] && i<=bbox[0][1] && j>=bbox[1][0] && j<=bbox[1][1] && k>=bbox[2][0] && k<=bbox[2][1])
	{
	    copy->push3D<true, true>(i, j, k, iter.getValue());
	}
	iter++;
    }

    copy->endSafePush();

    return copy;
}


template<class Traits>
void DTGrid<Traits>::paste(DTGrid<Traits> *subset, Index bbox[3][2])
{
    // TODO: Can be done in less memory if one first runs a trial traversal to compute the size of all the arrays

    InitParams initParams = InitParams(dx, beta, gamma, insideConstant, outsideConstant);
    DTGrid<Traits> *res = new DTGrid<Traits>(initParams);
    res->setTranslation(tx, ty, tz);
    res->setRotationXYZ(rx, ry, rz);

    TubeIterator iter1, iend1, iter2, iend2;

    iter1 = beginTubeIterator();
    iend1 = endTubeIterator();
    iter2 = subset->beginTubeIterator();
    iend2 = subset->endTubeIterator();

    res->beginSafePush();

    while (iter1 != iend1 || iter2 != iend2)
    {
	Index i1[3], i2[3];
	UInt i;

	if (iter1 == iend1)
	{
	    iter2.getIndex(&i2[0], &i2[1], &i2[2]);
	    i1[0] = i2[0]+1;  // force i2 to be lexicographically larger
	}
	else if (iter2 == iend2)
	{
	    iter1.getIndex(&i1[0], &i1[1], &i1[2]);
	    i2[0] = i1[0]+1;  // force i1 to be lexicographically larger
	}
	else
	{	
	    iter1.getIndex(&i1[0], &i1[1], &i1[2]);    
	    iter2.getIndex(&i2[0], &i2[1], &i2[2]);    
	}



	// next determine the lexicographically smallest element of the two iterators

	for (i=0; i<3; i++)
	{
	    if ( i1[i] < i2[i] )
	    {
		// i1 smallest

		// now check if i1 is inside the bbox
		bool insideBBox = (i1[0]>=bbox[0][0] && i1[0]<=bbox[0][1] && i1[1]>=bbox[1][0] && i1[1]<=bbox[1][1] && i1[2]>=bbox[2][0] && i1[2]<=bbox[2][1]);

		// only do something if we are outside the bbox
		if (!insideBBox)
		{
		    res->push3D<true, true>(i1[0], i1[1], i1[2], iter1.getValue());
		}

		iter1++;
		break;
	    }
	    else if ( i1[i] > i2[i] )
	    {
		// i2 smallest

		// now check if i2 is inside the bbox
		bool insideBBox = (i2[0]>=bbox[0][0] && i2[0]<=bbox[0][1] && i2[1]>=bbox[1][0] && i2[1]<=bbox[1][1] && i2[2]>=bbox[2][0] && i2[2]<=bbox[2][1]);

		// only do something if we are inside the bbox
		if (insideBBox)
		{
		    res->push3D<true, true>(i2[0], i2[1], i2[2], iter2.getValue());
		}

		iter2++;
		break;
	    }
	}



	if ( i==3 )
	{
	    // the two grid points were equal

	    // now check if the point is inside the bbox
	    bool insideBBox = (i1[0]>=bbox[0][0] && i1[0]<=bbox[0][1] && i1[1]>=bbox[1][0] && i1[1]<=bbox[1][1] && i1[2]>=bbox[2][0] && i1[2]<=bbox[2][1]);

	    if (insideBBox)
	    {
		res->push3D<true, true>(i2[0], i2[1], i2[2], iter2.getValue());
	    }
	    else
	    {
		res->push3D<true, true>(i1[0], i1[1], i1[2], iter1.getValue());
	    }

	    iter1++;
	    iter2++;
	}

    }
    
    res->endSafePush();

    clear();
    assign(res);
    res->clear(false);
}



// Generic iterator methods

// Semantics:
// - If both buffers i1 and i2 are specified, those buffers are used (buffer i2 is allocated if not equal to i1).
// - If no buffers are specified, buffers 0 and 1 are used, and buffer 1 allocated. (default parameters).
// - If only buffer i1 is specified, buffer i2 is set to buffer 1 and is assumed allocated.
template<class Traits>
template<TubeType tt, StencilFormat iterStencil, bool copyElements>
#ifdef WIN32
typename DTGrid<Traits>::StencilTubeIterator<tt> DTGrid<Traits>::beginStencilTubeIterator(int i1, int i2)
#else
typename DTGrid<Traits>::template StencilTubeIterator<tt> DTGrid<Traits>::beginStencilTubeIterator(int i1, int i2)
#endif
{
    bool initSecondBuffer;
    if (i2 == -1)
    {
	i2 = 1;
	initSecondBuffer = false;
    }
    else if (i2 == i1)
    {
	initSecondBuffer = false;
    }
    else
    {
	initSecondBuffer = true;
    }
    StencilTubeIterator<tt> iter;
    iter.template init<iterStencil, copyElements>(this, i1, i2, 1, initSecondBuffer);
    return iter;
}


template<class Traits>
template<TubeType tt, StencilFormat iterStencil, bool copyElements>
#ifdef WIN32
typename DTGrid<Traits>::StencilTubeIterator<tt> DTGrid<Traits>::endStencilTubeIterator(int i1)
#else
typename DTGrid<Traits>::template StencilTubeIterator<tt> DTGrid<Traits>::endStencilTubeIterator(int i1)
#endif
{
    StencilTubeIterator<tt> iter;
    iter.template init<iterStencil, copyElements>(this, i1, 0, 0, false);
    return iter;
}

template<class Traits>
template<TubeType tt, StencilFormat iterStencil, bool copyElements>
void DTGrid<Traits>::getStencilTubeIterators(UInt numIntervals, StencilTubeIterator<tt> *iterators, int i1, int i2)
{
    UInt delta = getNumValues(i1) / numIntervals;
    UInt currentPos = delta;
    UInt i;
    Locator loc;

    bool initSecondBuffer;
    if (i2 == -1)
    {
	i2 = 1;
	initSecondBuffer = false;
    }
    else if (i2 == i1)
    {
	initSecondBuffer = false;
    }
    else
    {
	initSecondBuffer = true;
    }

    iterators[0].template init<iterStencil, copyElements>(this, i1, i2, 1, initSecondBuffer);

    for (i=1; i<numIntervals; i++)
    {
	getLocator(currentPos, &loc);
	currentPos += delta;
	iterators[i].template init<iterStencil, copyElements>(this, i1, i2, 1, initSecondBuffer, &loc);
    }

    iterators[numIntervals] = endStencilTubeIterator<tt, iterStencil, copyElements>(i1);
}


template<class Traits> template<StencilFormat iterStencil, bool copyElements>
typename DTGrid<Traits>::ZeroCrossingIterator DTGrid<Traits>::beginZeroCrossing()
{
    return beginStencilTubeIterator<ZERO_CROSSING_TUBE, iterStencil, copyElements>(0);
}


template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::ZeroCrossingIterator DTGrid<Traits>::endZeroCrossing()
{
    return endStencilTubeIterator<ZERO_CROSSING_TUBE, iterStencil, copyElements>(0);
}

template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::BetaTubeIterator DTGrid<Traits>::beginBetaTube()
{
    return beginStencilTubeIterator<BETA_TUBE, iterStencil, copyElements>(0);
}


template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::BetaTubeIterator DTGrid<Traits>::endBetaTube()
{
    return endStencilTubeIterator<BETA_TUBE, iterStencil, copyElements>(0);
}


template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::GammaTubeIterator DTGrid<Traits>::beginGammaTube()
{
    return beginStencilTubeIterator<GAMMA_TUBE, iterStencil, copyElements>(0);
}


template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::GammaTubeIterator DTGrid<Traits>::endGammaTube()
{
    return endStencilTubeIterator<GAMMA_TUBE, iterStencil, copyElements>(0);
}

template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::EntireTubeIterator DTGrid<Traits>::beginEntireTube()
{
    return beginStencilTubeIterator<ENTIRE_TUBE, iterStencil, copyElements>(0);
}

template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::EntireTubeIterator DTGrid<Traits>::endEntireTube()
{
    return endStencilTubeIterator<ENTIRE_TUBE, iterStencil, copyElements>(0);
}


template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::InsideTubeIterator DTGrid<Traits>::beginInsideTube()
{
    return beginStencilTubeIterator<INSIDE_TUBE, iterStencil, copyElements>(0);
}

template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::InsideTubeIterator DTGrid<Traits>::endInsideTube()
{
    return endStencilTubeIterator<INSIDE_TUBE, iterStencil, copyElements>(0);
}


template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::OutsideTubeIterator DTGrid<Traits>::beginOutsideTube()
{
    return beginStencilTubeIterator<OUTSIDE_TUBE, iterStencil, copyElements>(0);
}

template<class Traits> template<StencilFormat iterStencil, bool copyElements> 
typename DTGrid<Traits>::OutsideTubeIterator DTGrid<Traits>::endOutsideTube()
{
    return endStencilTubeIterator<OUTSIDE_TUBE, iterStencil, copyElements>(0);
}



template<class Traits>
void DTGrid<Traits>::initBuffers()
{
    UInt i;
    for (i=1; i<maxNumBuffers; i++)
    {
	values[i] = NULL;
	numValues[i] = 0;
    }

    values[0] = va3D;
    numValues[0] = numVa3D;

    safeStorage = NULL;
}


template<class Traits>
void DTGrid<Traits>::freeBuffer(UInt i, bool doDelete)
{
    if (doDelete)
    {
	delete values[i];
    }

    values[i] = NULL;
    numValues[i] = 0;

    va3D = values[0];
    numVa3D = numValues[0];
}

template<class Traits>
void DTGrid<Traits>::swapBuffers(UInt i1, UInt i2)
{
    UInt tmpInt;
    Data *tmpValues;

    tmpInt = numValues[i1];
    numValues[i1] = numValues[i2];
    numValues[i2] = tmpInt;

    tmpValues = values[i1];
    values[i1] = values[i2];
    values[i2] = tmpValues;

    va3D = values[0];
    numVa3D = numValues[0];

    setEndMarkers();
}


template<class Traits>
void DTGrid<Traits>::allocateBuffer(UInt id, UInt size)
{
    delete values[id];
    values[id] = new Data[size+numValueEndMarkers];
    numValues[id] = size;

    va3D = values[0];
    numVa3D = numValues[0];
}



template<class Traits>
void DTGrid<Traits>::assign(DTGrid *dtgrid2)
{
    beta = dtgrid2->beta;
    gamma = dtgrid2->gamma;
    outsideConstant = dtgrid2->outsideConstant;
    insideConstant = dtgrid2->insideConstant;
    dx = dtgrid2->dx;
    bbox[0][0] = dtgrid2->bbox[0][0];
    bbox[0][1] = dtgrid2->bbox[0][1];
    bbox[1][0] = dtgrid2->bbox[1][0];
    bbox[1][1] = dtgrid2->bbox[1][1];
    bbox[2][0] = dtgrid2->bbox[2][0];
    bbox[2][1] = dtgrid2->bbox[2][1];

    va1D = dtgrid2->va1D;        
    numVa1D = dtgrid2->numVa1D;    
    xIndex = dtgrid2->xIndex;            
    numXIndex = dtgrid2->numXIndex; 
    lastXIndex = dtgrid2->lastXIndex;
    aa1D = dtgrid2->aa1D;         

    va2D = dtgrid2->va2D;        
    numVa2D = dtgrid2->numVa2D;
    yIndex = dtgrid2->yIndex;             
    numYIndex = dtgrid2->numYIndex;
    lastYIndex = dtgrid2->lastYIndex;
    aa2D = dtgrid2->aa2D;         

    zIndex = dtgrid2->zIndex;              
    numZIndex = dtgrid2->numZIndex;
    lastZIndex = dtgrid2->lastZIndex;
    aa3D = dtgrid2->aa3D;         

    numVa3D = dtgrid2->numVa3D;
    va3D = dtgrid2->va3D;                

    tx = dtgrid2->tx;
    ty = dtgrid2->ty;
    tz = dtgrid2->tz;

    rx = dtgrid2->rx;
    ry = dtgrid2->ry;
    rz = dtgrid2->rz;

    openLevelSetBBoxValid = dtgrid2->openLevelSetBBoxValid;

    openLevelSetBBox[0][0] = dtgrid2->openLevelSetBBox[0][0];
    openLevelSetBBox[1][0] = dtgrid2->openLevelSetBBox[1][0];
    openLevelSetBBox[2][0] = dtgrid2->openLevelSetBBox[2][0];
    openLevelSetBBox[0][1] = dtgrid2->openLevelSetBBox[0][1];
    openLevelSetBBox[1][1] = dtgrid2->openLevelSetBBox[1][1];
    openLevelSetBBox[2][1] = dtgrid2->openLevelSetBBox[2][1];

    initBuffers();
}


template<class Traits>
void DTGrid<Traits>::clear(bool doDelete)
{
    UInt i;
    for (i=0; i<maxNumBuffers; i++)
    {
	freeBuffer(i, doDelete);
    }

    if (doDelete)
    {
	delete[] va1D;
	delete[] xIndex;
	delete[] aa1D;

	delete[] va2D;
	delete[] yIndex;
	delete[] aa2D;

	delete[] zIndex;
	delete[] aa3D;
    }

    va1D = NULL;        
    numVa1D = 0; 
    xIndex = NULL;            
    numXIndex = 0; 
    lastXIndex = 0;
    aa1D = NULL;         

    va2D = NULL;        
    numVa2D = 0;
    yIndex = NULL;             
    numYIndex = 0;
    lastYIndex = 0;
    aa2D = NULL;         

    zIndex = NULL;              
    numZIndex = 0;
    lastZIndex = 0;
    aa3D = NULL;         

    numVa3D = 0;
    va3D = NULL;                
}



template<class Traits>
bool DTGrid<Traits>::inside(Index i, Index j, Index k) const
{
    return (*this)(i,j,k) <= 0;
}



///////////////////////////////////////////////////////////////////////////////////
// ITERATOR METHODS
///////////////////////////////////////////////////////////////////////////////////


template<class Traits>  template<TubeType tt>  template<StencilFormat iterStencil, bool copyElements>
void DTGrid<Traits>::StencilTubeIterator<tt>::init(DTGrid<Traits> *parent, UInt bufferId1, UInt bufferId2, UInt begin, bool initSecondBuffer, Locator *loc)
{
    if (begin)
    {
	// begin
	this->parent = parent;
	this->beta = parent->beta;
	this->gamma = parent->gamma;
	this->dx = parent->dx;
	dxc1 = 1.0 / dx;
	dxc2 = 1.0 / (2.0*dx);
	dxc3 = 1.0 / (dx*dx);
	dxc4 = 1.0 / (4.0*dx*dx);
	v1 = parent->values[bufferId1];
	numValues1 = parent->numValues[bufferId1];

	if (numValues1 != 0)
	{
	    if (initSecondBuffer)
	    {
		parent->allocateBuffer(bufferId2, numValues1);		    
	    }

	    v2 = parent->values[bufferId2];
	    numValues2 = parent->numValues[bufferId2];

	    TubeIterator iter = parent->beginTubeIterator(*loc, true, true, bufferId1);
	    TubeIterator iend = parent->endTubeIterator(bufferId1);

	    if (tt == ZERO_CROSSING_TUBE)
	    {
		while (fabs(*iter) > (Traits::zeroCrossingWidth*dx) && iter!=iend)
		{
		    if (copyElements)
		    {
			v2[iter.getArrayIndex()] = iter.getValue();
		    }
		    ++iter;
		}
	    }

	    if (tt == BETA_TUBE)
	    {
		while (fabs(*iter) >=beta && iter!=iend)
		{
		    if (copyElements)
		    {
			v2[iter.getArrayIndex()] = iter.getValue();
		    }
		    ++iter;
		}
	    }

	    if (tt == GAMMA_TUBE)
	    {
		while (fabs(*iter) >=gamma && iter!=iend)
		{
		    if (copyElements)
		    {
			v2[iter.getArrayIndex()] = iter.getValue();
		    }
		    ++iter;
		}
	    }

	    if (tt == INSIDE_TUBE)
	    {
		while ( *iter > 0 && iter!=iend)
		{
		    if (copyElements)
		    {
			v2[iter.getArrayIndex()] = iter.getValue();
		    }
		    ++iter;
		}
	    }

	    // For tt == OUTSIDE_TUBE no action is needed since the first element in the dtgrid is always outside!


	    if (iter != iend)
	    {
		i = iter.getI();
		j = iter.getJ();
		k = iter.getK();

		si[0] = iter;

		initStencilUsingRandomAccess<iterStencil>(bufferId1);
	    }
	    else
	    {
		si[0] = parent->endTubeIterator();
	    }
	}
	else
	{
	    si[0] = parent->endTubeIterator();
	}
    }
    else
    {
	// end
	if (parent == NULL)
	{
	}
	else
	{
	    si[0] = parent->endTubeIterator();
	}
    }
}


template<class Traits>  template<TubeType tt>  template<StencilFormat iterStencil, bool copyElements>
void DTGrid<Traits>::StencilTubeIterator<tt>::init(DTGrid<Traits> *parent, UInt bufferId1, UInt bufferId2, UInt begin, bool initSecondBuffer)
{
    if (begin)
    {
	// begin
	this->parent = parent;
	this->beta = parent->beta;
	this->gamma = parent->gamma;
	this->dx = parent->dx;
	dxc1 = static_cast<Real>(1.0) / dx;
	dxc2 = static_cast<Real>(1.0) / (static_cast<Real>(2.0)*dx);
	dxc3 = static_cast<Real>(1.0) / (dx*dx);
	dxc4 = static_cast<Real>(1.0) / (static_cast<Real>(4.0)*dx*dx);
	v1 = parent->values[bufferId1];
	numValues1 = parent->numValues[bufferId1];

	if (numValues1 != 0)
	{
	    if (initSecondBuffer)
	    {
		parent->allocateBuffer(bufferId2, numValues1);		    
	    }

	    v2 = parent->values[bufferId2];
	    numValues2 = parent->numValues[bufferId2];

	    TubeIterator iter = parent->beginTubeIterator(bufferId1);
	    TubeIterator iend = parent->endTubeIterator(bufferId1);

	    if (tt == ZERO_CROSSING_TUBE)
	    {
		while (fabs(*iter) > Traits::zeroCrossingWidth*dx && iter!=iend)
		{
		    if (copyElements)
		    {
			v2[iter.getArrayIndex()] = iter.getValue();
		    }
		    ++iter;
		}
	    }

	    if (tt == BETA_TUBE)
	    {
		while (fabs(*iter) >=beta && iter!=iend)
		{
		    if (copyElements)
		    {
			v2[iter.getArrayIndex()] = iter.getValue();
		    }
		    ++iter;
		}
	    }

	    if (tt == GAMMA_TUBE)
	    {
		while (fabs(*iter) >=gamma && iter!=iend)
		{
		    if (copyElements)
		    {
			v2[iter.getArrayIndex()] = iter.getValue();
		    }
		    ++iter;
		}
	    }

	    if (tt == INSIDE_TUBE)
	    {
		while ( *iter > 0 && iter!=iend)
		{
		    if (copyElements)
		    {
			v2[iter.getArrayIndex()] = iter.getValue();
		    }
		    ++iter;
		}
	    }

	    // For tt == OUTSIDE_TUBE no action is needed since the first element in the dtgrid is always outside!


	    if (iter != iend)
	    {
		int m;

		i = iter.getI();
		j = iter.getJ();
		k = iter.getK();

		si[0] = iter;

		for (m=1; m<StencilTraits::getStencilLength<iterStencil>(); m++)
		{
		    si[m] = parent->beginTubeIterator(bufferId1);
		}

		fastUpdate = false;

		updateStencil<iterStencil>();
	    }
	    else
	    {
		si[0] = parent->endTubeIterator();
	    }
	}
	else
	{
	    si[0] = parent->endTubeIterator();
	}
    }
    else
    {
	// end
	if (parent == NULL)
	{
	}
	else
	{
	    si[0] = parent->endTubeIterator();
	}
    }
}



template<class Traits>  template<TubeType tt> 
DTGrid<Traits>::StencilTubeIterator<tt>::~StencilTubeIterator()
{
}



template<class Traits>  template<TubeType tt> 
template<StencilFormat iterStencil> 
void DTGrid<Traits>::StencilTubeIterator<tt>::retrieveStencilValues()
{
    // now set values for all stencil grid points
    // use +/- gamma for values not inside the narrow band
    int m;
    int size = (int)(StencilTraits::getStencilLength<iterStencil>());
    for (m=0; m<size; m++)
    {
	if (si[m].isValid())
	{
	    si[m].retrieveValue();
	}
    }
}



template<class Traits>  template<TubeType tt>  template<StencilFormat iterStencil, bool copyElements>
#ifdef WIN32
typename DTGrid<Traits>::StencilTubeIterator<tt>& DTGrid<Traits>::StencilTubeIterator<tt>::operator++()
#else
typename DTGrid<Traits>::template StencilTubeIterator<tt>& DTGrid<Traits>::StencilTubeIterator<tt>::operator++()
#endif
{	

    // increment is allowed as we do not point to the last element
    si[0]++;

    if (tt == ZERO_CROSSING_TUBE)
    {
	while (fabs(si[0].getValue()) > Traits::zeroCrossingWidth*dx)
	{
	    if (copyElements)
	    {
		v2[si[0].getArrayIndex()] = si[0].getValue();
	    }
	    si[0]++;
	}
    }

    if (tt == BETA_TUBE)
    {
	while (fabs(si[0].getValue()) >= parent->beta /*&& si[0].hasNext()*/)
	{
	    if (copyElements)
	    {
		v2[si[0].getArrayIndex()] = si[0].getValue();
	    }
	    si[0]++;
	}
    }

    if (tt == GAMMA_TUBE)
    {
	while (fabs(si[0].getValue()) >= parent->gamma /*&& si[0].hasNext()*/)
	{
	    if (copyElements)
	    {
		v2[si[0].getArrayIndex()] = si[0].getValue();
	    }
	    si[0]++;
	}
    }


    if (tt == INSIDE_TUBE)
    {
	while ( si[0].getValue() > 0 )   
	{
	    if (copyElements)
	    {
		v2[si[0].getArrayIndex()] = si[0].getValue();
	    }
	    si[0]++;
	}
    }


    if (tt == OUTSIDE_TUBE)
    {
	while ( si[0].getValue() < 0 )
	{
	    if (copyElements)
	    {
		v2[si[0].getArrayIndex()] = si[0].getValue();
	    }
	    si[0]++;
	}
    }


    if (!si[0].hasNext())
    {
	return *this;
    }

    if ( i == si[0].getI() && j == si[0].getJ() && si[0].getK() == k+1)  // // VTUNE caution: LCP (Length Changing Prefix), but alternative is to use long instead of short at the cost of more memory usage
    {
	// A fast update is possible if the stencil only has moved one grid point in the z direction
	fastUpdate = true;
	k++;                        // // VTUNE caution: LCP (Length Changing Prefix), but alternative is to use long instead of short at the cost of more memory usage
    }
    else
    {
	fastUpdate = false;
	i = si[0].getI();
	j = si[0].getJ();
	k = si[0].getK();    // // VTUNE caution: LCP (Length Changing Prefix), but alternative is to use long instead of short at the cost of more memory usage
    }

    // in the case of SF_NONE, this call is optimized away by the compiler
    updateStencil<iterStencil>();

    return *this;
}




template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::getValue(const int i, const int j)
{
    if (i==-1)
    {
	return si[j].getValue();
    }
    else
    {
	return parent->values[i][si[j].getArrayIndex()];
    }
}

template<class Traits>  template<TubeType tt> 
void DTGrid<Traits>::StencilTubeIterator<tt>::setValue(Real v, const int i, const int j)
{
    if (i==-1)
    {
	v1[si[j].getArrayIndex()] = v;
	si[j].setCachedValue(v);
    }
    else
    {
	parent->values[i][si[j].getArrayIndex()] = v;
    }
}


template<class Traits>  template<TubeType tt> 
bool DTGrid<Traits>::StencilTubeIterator<tt>::operator==(const StencilTubeIterator& iter) const
{
    return si[0] == iter.si[0];
}

template<class Traits>  template<TubeType tt> 
bool DTGrid<Traits>::StencilTubeIterator<tt>::operator!=(const StencilTubeIterator& iter) const
{
    return si[0] != iter.si[0];
}

template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::operator*() const
{
    return si[0].getValue();
}

template<class Traits>  template<TubeType tt> 
void DTGrid<Traits>::StencilTubeIterator<tt>::getIndex(Index* i, Index *j, Index *k) const
{
    *i = this->i;
    *j = this->j;
    *k = this->k;
}


template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::UInt DTGrid<Traits>::StencilTubeIterator<tt>::getArrayIndex() const
{
    return si[0].getArrayIndex();
}

template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::UInt DTGrid<Traits>::StencilTubeIterator<tt>::getArrayIndex(int i) const
{
    return si[i].getArrayIndex();
}

// first order one sided difference to the left in the z direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dm1z() const
{
    return (si[0].getValue() - si[5].getValue()) * dxc1;
}


// first order one sided difference to the right in the z direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dp1z() const
{
    return (si[6].getValue() - si[0].getValue()) * dxc1;
}


// first order one sided difference to the left in the x direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dm1x() const
{
    return (si[0].getValue() - si[1].getValue()) * dxc1;
}

// first order one sided difference to the right in the x direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dp1x() const
{
    return (si[2].getValue() - si[0].getValue()) * dxc1;
}


// first order one sided difference to the left in the y direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dm1y() const
{
    return (si[0].getValue() - si[3].getValue()) * dxc1;
}


// first order one sided difference to the right in the y direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dp1y() const
{
    return (si[4].getValue() - si[0].getValue()) * dxc1;
}


template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dc2x() const
{
    return (si[2].getValue() - si[1].getValue()) * dxc2;
}

template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dc2y() const
{
    return (si[4].getValue() - si[3].getValue()) * dxc2;
}


template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dc2z() const
{
    return (si[6].getValue() - si[5].getValue()) * dxc2;
}


// second order accurate second partial derivative in x direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::d2xx() const
{
    return ( si[2].getValue() - 2*si[0].getValue() + si[1].getValue() ) * dxc3;	
}

// second order accurate second partial derivative in y direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::d2yy() const
{
    return ( si[4].getValue() - 2*si[0].getValue() + si[3].getValue() ) * dxc3;	
}

// second order accurate second partial derivative in z direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::d2zz() const
{
    return ( si[6].getValue() - 2*si[0].getValue() + si[5].getValue() ) * dxc3;	
}

// second order accurate second partial derivative
template<class Traits>  template<TubeType tt>  template<StencilFormat iterStencil>
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::d2xy() const
{
    if (iterStencil == SF_FIRSTORDER_CURVATURE || iterStencil == SF_BOX)
    {
	return ( si[11].getValue() + si[14].getValue() - si[13].getValue() - si[12].getValue()) * dxc4; // / (4 * Math::pow2(dx));
    }
    else if (iterStencil == SF_WENO_CURVATURE)
    {
	return ( si[23].getValue() + si[26].getValue() - si[24].getValue() - si[25].getValue() ) * dxc4; // / (4 * Math::pow2(dx));	    
    }
    else
    {
	return 0;
    }
}

// second order accurate second partial derivative
template<class Traits>  template<TubeType tt>  template<StencilFormat iterStencil>
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::d2xz() const
{
    if (iterStencil == SF_FIRSTORDER_CURVATURE || iterStencil == SF_BOX)
    {
	return ( si[8].getValue() + si[17].getValue() - si[16].getValue() - si[9].getValue() ) * dxc4; // / (4 * Math::pow2(dx));
    }
    else if (iterStencil == SF_WENO_CURVATURE)
    {
	return ( si[20].getValue() + si[29].getValue() - si[21].getValue() - si[28].getValue() ) * dxc4; // / (4 * Math::pow2(dx));	    
    }
    else
    {
	return 0;
    }
}

// second order accurate second partial derivative
template<class Traits>  template<TubeType tt>  template<StencilFormat iterStencil>
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::d2yz() const
{
    if (iterStencil == SF_FIRSTORDER_CURVATURE || iterStencil == SF_BOX)
    {
	return ( si[7].getValue() + si[18].getValue() - si[15].getValue() - si[10].getValue()) * dxc4; // / (4 * Math::pow2(dx));
    }
    else if (iterStencil == SF_WENO_CURVATURE)
    {
	return ( si[19].getValue() + si[30].getValue() - si[22].getValue() - si[27].getValue() ) * dxc4; // / (4 * Math::pow2(dx));	    
    }
    else
    {
	return 0;
    }
}



template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::computeWENO(Real v1, Real v2, Real v3, Real v4, Real v5) const
{
    const Real a_dxth = this->dxc1;

    Real v[5];
    v[0] = v1 * a_dxth;
    v[1] = v2 * a_dxth;
    v[2] = v3 * a_dxth;
    v[3] = v4 * a_dxth;
    v[4] = v5 * a_dxth;

    // calculate divs once.
    static const Real a_fourth = Real(1.0)/Real(4.0);
    static const Real a_sixth = Real(1.0)/Real(6.0);
    static const Real a_twelvth = Real(1.0)/Real(12.0);

    // calculate phi1 to phi3
    Real phi[3];
    phi[0] = (11*v[2] + 2*v[0] - 7*v[1]) * a_sixth; /// 6.0;
    phi[1] = (5*v[2] + 2*v[3] - v[1]) * a_sixth; /// 6.0;
    phi[2] = (5*v[3] + 2*v[2] - v[4]) * a_sixth; /// 6.0;

    // Calculate smoothness estimates
    Real s[3];
    s[0] = (static_cast<Real>(13.0) * (Math::pow2(v[0]-static_cast<Real>(2.0)*v[1]+v[2])) ) * a_twelvth +  Math::pow2(v[0]+3*v[2]-4*v[1]) * a_fourth;
    s[1] = (static_cast<Real>(13.0) * (Math::pow2(v[1]-static_cast<Real>(2.0)*v[2]+v[3])) ) * a_twelvth +  Math::pow2(v[1]-v[3]) * a_fourth;
    s[2] = (static_cast<Real>(13.0) * (Math::pow2(v[2]-static_cast<Real>(2.0)*v[3]+v[4])) ) * a_twelvth +  Math::pow2(v[4]+3*v[2]-4*v[3]) * a_fourth;

    // Calculate epsilon and alphas
    Real alpha[3];

    Real epsilon = static_cast<Real>(1e-6);

    alpha[0] = static_cast<Real>(0.1) / Math::pow2(s[0]+epsilon);
    alpha[1] = static_cast<Real>(0.6) / Math::pow2(s[1]+epsilon);
    alpha[2] = static_cast<Real>(0.3) / Math::pow2(s[2]+epsilon);

    // Finally calculate and return phi
    Real norm = (alpha[0] + alpha[1] + alpha[2]);

    return (alpha[0]*phi[0] + alpha[1]*phi[1] + alpha[2]*phi[2])/norm;
}



// fifth order WENO one sided finite difference to the left in the x direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dm5x() const
{
    Real v1, v2, v3, v4, v5;

    v1 = si[8].getValue() - si[7].getValue();
    v2 = si[1].getValue() - si[8].getValue();
    v3 = si[0].getValue() - si[1].getValue();
    v4 = si[2].getValue() - si[0].getValue();
    v5 = si[9].getValue() - si[2].getValue();

    return computeWENO(v1, v2, v3, v4, v5);
}

// fifth order WENO one sided finite difference to the right in the x direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dp5x() const
{
    Real v1, v2, v3, v4, v5;

    v1 = si[10].getValue() - si[9].getValue();
    v2 = si[9].getValue() - si[2].getValue();
    v3 = si[2].getValue() - si[0].getValue();
    v4 = si[0].getValue() - si[1].getValue();
    v5 = si[1].getValue() - si[8].getValue();

    return computeWENO(v1, v2, v3, v4, v5);
}

// fifth order WENO one sided finite difference to the left in the y direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dm5y() const
{
    Real v1, v2, v3, v4, v5;

    v1 = si[12].getValue() - si[11].getValue();
    v2 = si[3].getValue() - si[12].getValue();
    v3 = si[0].getValue() - si[3].getValue();
    v4 = si[4].getValue() - si[0].getValue();
    v5 = si[13].getValue() - si[4].getValue();

    return computeWENO(v1, v2, v3, v4, v5);
}

// fifth order WENO one sided finite difference to the right in the y direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dp5y() const
{
    Real v1, v2, v3, v4, v5;

    v1 = si[14].getValue() - si[13].getValue();
    v2 = si[13].getValue() - si[4].getValue();
    v3 = si[4].getValue() - si[0].getValue();
    v4 = si[0].getValue() - si[3].getValue();
    v5 = si[3].getValue() - si[12].getValue();

    return computeWENO(v1, v2, v3, v4, v5);
}

// fifth order WENO one sided finite difference to the left in the z direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dm5z() const
{
    Real v1, v2, v3, v4, v5;

    v1 = si[16].getValue() - si[15].getValue();
    v2 = si[5].getValue() - si[16].getValue();
    v3 = si[0].getValue() - si[5].getValue();
    v4 = si[6].getValue() - si[0].getValue();
    v5 = si[17].getValue() - si[6].getValue();

    return computeWENO(v1, v2, v3, v4, v5);
}

// fifth order WENO one sided finite difference to the right in the z direction
template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::dp5z() const
{
    Real v1, v2, v3, v4, v5;

    v1 = si[18].getValue() - si[17].getValue();
    v2 = si[17].getValue() - si[6].getValue();
    v3 = si[6].getValue() - si[0].getValue();
    v4 = si[0].getValue() - si[5].getValue();
    v5 = si[5].getValue() - si[16].getValue();

    return computeWENO(v1, v2, v3, v4, v5);
}



template<class Traits>  template<TubeType tt> 
void DTGrid<Traits>::StencilTubeIterator<tt>::gradient(Real g[3]) const
{
    g[0] = dc2x();
    g[1] = dc2y();
    g[2] = dc2z();
}

template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::gradientLength() const
{
    return sqrt( Math::pow2(dc2x()) + Math::pow2(dc2y()) + Math::pow2(dc2z()) );
}

template<class Traits>  template<TubeType tt> 
void DTGrid<Traits>::StencilTubeIterator<tt>::normal(Real n[3]) const
{
    Real len;

    n[0] = dc2x();
    n[1] = dc2y();
    n[2] = dc2z();

    len = sqrt( Math::pow2(n[0])+Math::pow2(n[1])+Math::pow2(n[2]) );

    if (len != 0)
    {
	n[0] = n[0] / len;
	n[1] = n[1] / len;
	n[2] = n[2] / len;
    }
}

template<class Traits>  template<TubeType tt>
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::laplacian() const
{
    return d2xx() + d2yy() + d2zz();
}

template<class Traits>  template<TubeType tt>  template<StencilFormat iterStencil> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::meanCurvature() const
{
    // Curvature is computed using the traditional finite difference method
    Real Dx = dc2x();
    Real Dy = dc2y();
    Real Dz = dc2z();

    Real Dxx = d2xx();
    Real Dyy = d2yy();
    Real Dzz = d2zz();

    Real Dxy = d2xy<iterStencil>();
    Real Dxz = d2xz<iterStencil>();
    Real Dyz = d2yz<iterStencil>();

    Real denom = Math::pow3( sqrt( Math::pow2(Dx) + Math::pow2(Dy) + Math::pow2(Dz) ) );

    Real curvature = ( Math::pow2(Dx)*Dyy - 2*Dx*Dy*Dxy + Math::pow2(Dy)*Dxx + Math::pow2(Dx)*Dzz - 2*Dx*Dz*Dxz + Math::pow2(Dz)*Dxx + Math::pow2(Dy)*Dzz - 2*Dy*Dz*Dyz + Math::pow2(Dz)*Dyy) / denom;

    // In 3D, the divergence of the normal must be multiplied by 0.5 to obtain mean curvature
    return static_cast<Real>(0.5) * curvature;
}



template<class Traits>  template<TubeType tt> 
typename DTGrid<Traits>::Real DTGrid<Traits>::StencilTubeIterator<tt>::operator()(UInt i) const
{
    return si[i].getValue();
}


template<class Traits>  template<TubeType tt>  template<StencilFormat iterStencil>
void DTGrid<Traits>::StencilTubeIterator<tt>::updateStencil()
{	
    // THE INCREMENT METHODS ENSURE:
    // the invariant on the column index is that
    // 1) it is the correct index
    // 2) it is the smallest index that is larger than the correct index

    Data outsideValue = ( si[0].getValue() < 0 ? parent->insideConstant : parent->outsideConstant );


    if (iterStencil == SF_FIRSTORDER)
    {
	if ( ( tt == BETA_TUBE || tt == ZERO_CROSSING_TUBE ) 
	    || (Traits::includeSafeBand && tt == GAMMA_TUBE) // this is only valid if the level set includes a safe band!
	    // It is always safe to use this optimization when includeSafeBand is enabled. This is because of the additional safeband layer of grid points.
	    || (!Traits::includeSafeBand && tt == GAMMA_TUBE && fabs(si[0].getValue())<dx)
	    )
	{
	    if (fastUpdate)
	    {
		// In this case we know that the elements of the stencil ALWAYS exist in the narrow band
		si[1].incrementFast();	
		si[2].incrementFast();	
		si[3].incrementFast();	
		si[4].incrementFast();
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1,j,k,1);
		    INCREMENT_UNTIL(outsideValue, i+1,j,k,2);
		    INCREMENT_UNTIL(outsideValue, i,j-1,k,3);
		    INCREMENT_UNTIL(outsideValue, i,j+1,k,4);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);
		}
	    }

	    // note that iterators at the elements plus/minus Z are not maintained in this case, only the value is set!
	    si[5].getValue() = v1[si[0].getArrayIndex()-1];	
	    si[6].getValue() = v1[si[0].getArrayIndex()+1];

	    if (!Traits::includeSafeBand && tt == GAMMA_TUBE)
	    {
		si[5].incrementArrayIndexAndZ();	
		si[6].incrementArrayIndexAndZ();
	    }
	    else if (Traits::maintainArrayIndexInStencil)
	    {
		si[5].getArrayIndex() = si[0].getArrayIndex()-1;	
		si[6].getArrayIndex() = si[0].getArrayIndex()+1;
	    }
	}
	else 
	{
	    if (fastUpdate)
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_FAST_Z(outsideValue, k,1);
		    INCREMENT_FAST_Z(outsideValue, k,2);
		    INCREMENT_FAST_Z(outsideValue, k,3);
		    INCREMENT_FAST_Z(outsideValue, k,4);
		    INCREMENT_FAST_Z(outsideValue, k-1,5);
		    INCREMENT_FAST_Z(outsideValue, k+1,6);
		}
		else
		{
		    si[1].incrementFastZ(outsideValue, k);	
		    si[2].incrementFastZ(outsideValue, k);	
		    si[3].incrementFastZ(outsideValue, k);	
		    si[4].incrementFastZ(outsideValue, k);	
		    si[5].incrementFastZ(outsideValue, k-1);	
		    si[6].incrementFastZ(outsideValue, k+1);
		}
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1,j,k,1);
		    INCREMENT_UNTIL(outsideValue, i+1,j,k,2);
		    INCREMENT_UNTIL(outsideValue, i,j-1,k,3);
		    INCREMENT_UNTIL(outsideValue, i,j+1,k,4);
		    INCREMENT_UNTIL(outsideValue, i,j,k-1,5);
		    INCREMENT_UNTIL(outsideValue, i,j,k+1,6);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[5].incrementUntil(outsideValue, i, j, k-1);	
		    si[6].incrementUntil(outsideValue, i, j, k+1);
		}
	    }
	}
    } 
    else if (iterStencil == SF_FIRSTORDER_CURVATURE)
    {
	if ( ( tt == BETA_TUBE || tt == ZERO_CROSSING_TUBE )
	    || (Traits::includeSafeBand && tt == GAMMA_TUBE) // this is only valid if the level set includes a safe band!
	    )
	{
	    if (fastUpdate)
	    {
		// In this case we know that the elements of the stencil ALWAYS exist in the narrow band
		si[1].incrementFast();	
		si[2].incrementFast();	
		si[3].incrementFast();	
		si[4].incrementFast();
		si[11].incrementFast();	
		si[12].incrementFast();	
		si[13].incrementFast();	
		si[14].incrementFast();
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1,j,k,1);
		    INCREMENT_UNTIL(outsideValue, i+1,j,k,2);
		    INCREMENT_UNTIL(outsideValue, i,j-1,k,3);
		    INCREMENT_UNTIL(outsideValue, i,j+1,k,4);
		    INCREMENT_UNTIL(outsideValue, i-1,j-1,k,11);
		    INCREMENT_UNTIL(outsideValue, i+1,j-1,k,12);
		    INCREMENT_UNTIL(outsideValue, i-1,j+1,k,13);
		    INCREMENT_UNTIL(outsideValue, i+1,j+1,k,14);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[11].incrementUntil(outsideValue, i-1, j-1, k);
		    si[12].incrementUntil(outsideValue, i+1, j-1, k);
		    si[13].incrementUntil(outsideValue, i-1, j+1, k);
		    si[14].incrementUntil(outsideValue, i+1, j+1, k);
		}
	    }
	    // note that iterators at the elements plus/minus Z are not maintained in this case, only the value is set!

	    si[5].getValue() = v1[si[0].getArrayIndex()-1];	
	    si[7].getValue() = v1[si[3].getArrayIndex()-1];	
	    si[8].getValue() = v1[si[1].getArrayIndex()-1];	
	    si[9].getValue() = v1[si[2].getArrayIndex()-1];	
	    si[10].getValue() = v1[si[4].getArrayIndex()-1];	
	    si[6].getValue() = v1[si[0].getArrayIndex()+1];	
	    si[15].getValue() = v1[si[3].getArrayIndex()+1];	
	    si[16].getValue() = v1[si[1].getArrayIndex()+1];	
	    si[17].getValue() = v1[si[2].getArrayIndex()+1];	
	    si[18].getValue() = v1[si[4].getArrayIndex()+1];

	    if (Traits::maintainArrayIndexInStencil)
	    {
		si[5].getArrayIndex() = si[0].getArrayIndex()-1;	
		si[7].getArrayIndex() = si[3].getArrayIndex()-1;	
		si[8].getArrayIndex() = si[1].getArrayIndex()-1;	
		si[9].getArrayIndex() = si[2].getArrayIndex()-1;	
		si[10].getArrayIndex() = si[4].getArrayIndex()-1;	
		si[6].getArrayIndex() = si[0].getArrayIndex()+1;	
		si[15].getArrayIndex() = si[3].getArrayIndex()+1;	
		si[16].getArrayIndex() = si[1].getArrayIndex()+1;	
		si[17].getArrayIndex() = si[2].getArrayIndex()+1;	
		si[18].getArrayIndex() = si[4].getArrayIndex()+1;
	    }
	}
	else 
	{
	    if (fastUpdate)
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_FAST_Z(outsideValue, k,1);
		    INCREMENT_FAST_Z(outsideValue, k,2);
		    INCREMENT_FAST_Z(outsideValue, k,3);
		    INCREMENT_FAST_Z(outsideValue, k,4);
		    INCREMENT_FAST_Z(outsideValue, k-1,5);
		    INCREMENT_FAST_Z(outsideValue, k+1,6);
		    INCREMENT_FAST_Z(outsideValue, k-1,7);
		    INCREMENT_FAST_Z(outsideValue, k-1,8);
		    INCREMENT_FAST_Z(outsideValue, k-1,9);
		    INCREMENT_FAST_Z(outsideValue, k-1,10);
		    INCREMENT_FAST_Z(outsideValue, k,11);
		    INCREMENT_FAST_Z(outsideValue, k,12);
		    INCREMENT_FAST_Z(outsideValue, k,13);
		    INCREMENT_FAST_Z(outsideValue, k,14);
		    INCREMENT_FAST_Z(outsideValue, k+1,15);
		    INCREMENT_FAST_Z(outsideValue, k+1,16);
		    INCREMENT_FAST_Z(outsideValue, k+1,17);
		    INCREMENT_FAST_Z(outsideValue, k+1,18);
		}
		else
		{
		    si[1].incrementFastZ(outsideValue, k);	
		    si[2].incrementFastZ(outsideValue, k);	
		    si[3].incrementFastZ(outsideValue, k);	
		    si[4].incrementFastZ(outsideValue, k);	
		    si[5].incrementFastZ(outsideValue, k-1);	
		    si[6].incrementFastZ(outsideValue, k+1);
		    si[7].incrementFastZ(outsideValue, k-1);
		    si[8].incrementFastZ(outsideValue, k-1);
		    si[9].incrementFastZ(outsideValue, k-1);
		    si[10].incrementFastZ(outsideValue, k-1);
		    si[11].incrementFastZ(outsideValue, k);
		    si[12].incrementFastZ(outsideValue, k);
		    si[13].incrementFastZ(outsideValue, k);
		    si[14].incrementFastZ(outsideValue, k);
		    si[15].incrementFastZ(outsideValue, k+1);
		    si[16].incrementFastZ(outsideValue, k+1);
		    si[17].incrementFastZ(outsideValue, k+1);
		    si[18].incrementFastZ(outsideValue, k+1);
		}
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1,j,k,1);
		    INCREMENT_UNTIL(outsideValue, i+1,j,k,2);
		    INCREMENT_UNTIL(outsideValue, i,j-1,k,3);
		    INCREMENT_UNTIL(outsideValue, i,j+1,k,4);
		    INCREMENT_UNTIL(outsideValue, i,j,k-1,5);
		    INCREMENT_UNTIL(outsideValue, i,j,k+1,6);
		    INCREMENT_UNTIL(outsideValue, i,j-1,k-1,7);
		    INCREMENT_UNTIL(outsideValue, i-1,j,k-1,8);
		    INCREMENT_UNTIL(outsideValue, i+1,j,k-1,9);
		    INCREMENT_UNTIL(outsideValue, i,j+1,k-1,10);
		    INCREMENT_UNTIL(outsideValue, i-1,j-1,k,11);
		    INCREMENT_UNTIL(outsideValue, i+1,j-1,k,12);
		    INCREMENT_UNTIL(outsideValue, i-1,j+1,k,13);
		    INCREMENT_UNTIL(outsideValue, i+1,j+1,k,14);
		    INCREMENT_UNTIL(outsideValue, i,j-1,k+1,15);
		    INCREMENT_UNTIL(outsideValue, i-1,j,k+1,16);
		    INCREMENT_UNTIL(outsideValue, i+1,j,k+1,17);
		    INCREMENT_UNTIL(outsideValue, i,j+1,k+1,18);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[5].incrementUntil(outsideValue, i, j, k-1);	
		    si[6].incrementUntil(outsideValue, i, j, k+1);
		    si[7].incrementUntil(outsideValue, i, j-1, k-1);
		    si[8].incrementUntil(outsideValue, i-1, j, k-1);
		    si[9].incrementUntil(outsideValue, i+1, j, k-1);
		    si[10].incrementUntil(outsideValue, i, j+1, k-1);
		    si[11].incrementUntil(outsideValue, i-1, j-1, k);
		    si[12].incrementUntil(outsideValue, i+1, j-1, k);
		    si[13].incrementUntil(outsideValue, i-1, j+1, k);
		    si[14].incrementUntil(outsideValue, i+1, j+1, k);
		    si[15].incrementUntil(outsideValue, i, j-1, k+1);
		    si[16].incrementUntil(outsideValue, i-1, j, k+1);
		    si[17].incrementUntil(outsideValue, i+1, j, k+1);
		    si[18].incrementUntil(outsideValue, i, j+1, k+1);
		}
	    }
	}
    }
    else if (iterStencil == SF_BOX)
    {
	if ( ( tt == BETA_TUBE || tt == ZERO_CROSSING_TUBE )
	    || (Traits::includeSafeBand && tt == GAMMA_TUBE) // this is only valid if the level set includes a safe band!
	    )
	{
	    if (fastUpdate)
	    {
		// In this case we know that the elements of the stencil ALWAYS exist in the narrow band
		si[1].incrementFast();	
		si[2].incrementFast();	
		si[3].incrementFast();	
		si[4].incrementFast();
		si[11].incrementFast();	
		si[12].incrementFast();	
		si[13].incrementFast();	
		si[14].incrementFast();
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1, j, k, 1);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k, 2);	
		    INCREMENT_UNTIL(outsideValue, i, j-1, k, 3);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k, 4);	
		    INCREMENT_UNTIL(outsideValue, i-1, j-1, k, 11);
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k, 12);
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k, 13);
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k, 14);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[11].incrementUntil(outsideValue, i-1, j-1, k);
		    si[12].incrementUntil(outsideValue, i+1, j-1, k);
		    si[13].incrementUntil(outsideValue, i-1, j+1, k);
		    si[14].incrementUntil(outsideValue, i+1, j+1, k);
		}
	    }
	    // note that iterators at the elements plus/minus Z are not maintained in this case, only the value is set!
	    si[5].getValue() = v1[si[0].getArrayIndex()-1];	
	    si[7].getValue() = v1[si[3].getArrayIndex()-1];	
	    si[8].getValue() = v1[si[1].getArrayIndex()-1];	
	    si[9].getValue() = v1[si[2].getArrayIndex()-1];	
	    si[10].getValue() = v1[si[4].getArrayIndex()-1];	
	    si[6].getValue() = v1[si[0].getArrayIndex()+1];	
	    si[15].getValue() = v1[si[3].getArrayIndex()+1];	
	    si[16].getValue() = v1[si[1].getArrayIndex()+1];	
	    si[17].getValue() = v1[si[2].getArrayIndex()+1];	
	    si[18].getValue() = v1[si[4].getArrayIndex()+1];	
	    si[19].getValue() = v1[si[11].getArrayIndex()-1];	
	    si[20].getValue() = v1[si[12].getArrayIndex()-1];	
	    si[21].getValue() = v1[si[13].getArrayIndex()-1];	
	    si[22].getValue() = v1[si[14].getArrayIndex()-1];	
	    si[23].getValue() = v1[si[11].getArrayIndex()+1];	
	    si[24].getValue() = v1[si[12].getArrayIndex()+1];	
	    si[25].getValue() = v1[si[13].getArrayIndex()+1];	
	    si[26].getValue() = v1[si[14].getArrayIndex()+1];	

	    if (Traits::maintainArrayIndexInStencil)
	    {
		si[5].getArrayIndex() = si[0].getArrayIndex()-1;	
		si[7].getArrayIndex() = si[3].getArrayIndex()-1;	
		si[8].getArrayIndex() = si[1].getArrayIndex()-1;	
		si[9].getArrayIndex() = si[2].getArrayIndex()-1;	
		si[10].getArrayIndex() = si[4].getArrayIndex()-1;	
		si[6].getArrayIndex() = si[0].getArrayIndex()+1;	
		si[15].getArrayIndex() = si[3].getArrayIndex()+1;	
		si[16].getArrayIndex() = si[1].getArrayIndex()+1;	
		si[17].getArrayIndex() = si[2].getArrayIndex()+1;	
		si[18].getArrayIndex() = si[4].getArrayIndex()+1;	
		si[19].getArrayIndex() = si[11].getArrayIndex()-1;	
		si[20].getArrayIndex() = si[12].getArrayIndex()-1;	
		si[21].getArrayIndex() = si[13].getArrayIndex()-1;	
		si[22].getArrayIndex() = si[14].getArrayIndex()-1;	
		si[23].getArrayIndex() = si[11].getArrayIndex()+1;	
		si[24].getArrayIndex() = si[12].getArrayIndex()+1;	
		si[25].getArrayIndex() = si[13].getArrayIndex()+1;	
		si[26].getArrayIndex() = si[14].getArrayIndex()+1;
	    }
	}
	else 
	{
	    if (fastUpdate)
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_FAST_Z(outsideValue, k, 1);	
		    INCREMENT_FAST_Z(outsideValue, k, 2);	
		    INCREMENT_FAST_Z(outsideValue, k, 3);	
		    INCREMENT_FAST_Z(outsideValue, k, 4);	
		    INCREMENT_FAST_Z(outsideValue, k-1, 5);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 6);
		    INCREMENT_FAST_Z(outsideValue, k-1, 7);
		    INCREMENT_FAST_Z(outsideValue, k-1, 8);
		    INCREMENT_FAST_Z(outsideValue, k-1, 9);
		    INCREMENT_FAST_Z(outsideValue, k-1, 10);
		    INCREMENT_FAST_Z(outsideValue, k, 11);
		    INCREMENT_FAST_Z(outsideValue, k, 12);
		    INCREMENT_FAST_Z(outsideValue, k, 13);
		    INCREMENT_FAST_Z(outsideValue, k, 14);
		    INCREMENT_FAST_Z(outsideValue, k+1, 15);
		    INCREMENT_FAST_Z(outsideValue, k+1, 16);
		    INCREMENT_FAST_Z(outsideValue, k+1, 17);
		    INCREMENT_FAST_Z(outsideValue, k+1, 18);
		    INCREMENT_FAST_Z(outsideValue, k-1, 19);
		    INCREMENT_FAST_Z(outsideValue, k-1, 20);
		    INCREMENT_FAST_Z(outsideValue, k-1, 21);
		    INCREMENT_FAST_Z(outsideValue, k-1, 22);
		    INCREMENT_FAST_Z(outsideValue, k+1, 23);
		    INCREMENT_FAST_Z(outsideValue, k+1, 24);
		    INCREMENT_FAST_Z(outsideValue, k+1, 25);
		    INCREMENT_FAST_Z(outsideValue, k+1, 26);
		}
		else
		{
		    si[1].incrementFastZ(outsideValue, k);	
		    si[2].incrementFastZ(outsideValue, k);	
		    si[3].incrementFastZ(outsideValue, k);	
		    si[4].incrementFastZ(outsideValue, k);	
		    si[5].incrementFastZ(outsideValue, k-1);	
		    si[6].incrementFastZ(outsideValue, k+1);
		    si[7].incrementFastZ(outsideValue, k-1);
		    si[8].incrementFastZ(outsideValue, k-1);
		    si[9].incrementFastZ(outsideValue, k-1);
		    si[10].incrementFastZ(outsideValue, k-1);
		    si[11].incrementFastZ(outsideValue, k);
		    si[12].incrementFastZ(outsideValue, k);
		    si[13].incrementFastZ(outsideValue, k);
		    si[14].incrementFastZ(outsideValue, k);
		    si[15].incrementFastZ(outsideValue, k+1);
		    si[16].incrementFastZ(outsideValue, k+1);
		    si[17].incrementFastZ(outsideValue, k+1);
		    si[18].incrementFastZ(outsideValue, k+1);
		    si[19].incrementFastZ(outsideValue, k-1);
		    si[20].incrementFastZ(outsideValue, k-1);
		    si[21].incrementFastZ(outsideValue, k-1);
		    si[22].incrementFastZ(outsideValue, k-1);
		    si[23].incrementFastZ(outsideValue, k+1);
		    si[24].incrementFastZ(outsideValue, k+1);
		    si[25].incrementFastZ(outsideValue, k+1);
		    si[26].incrementFastZ(outsideValue, k+1);
		}
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1, j, k, 1);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k, 2);	
		    INCREMENT_UNTIL(outsideValue, i, j-1, k, 3);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k, 4);	
		    INCREMENT_UNTIL(outsideValue, i, j, k-1, 5);	
		    INCREMENT_UNTIL(outsideValue, i, j, k+1, 6);
		    INCREMENT_UNTIL(outsideValue, i, j-1, k-1, 7);
		    INCREMENT_UNTIL(outsideValue, i-1, j, k-1, 8);
		    INCREMENT_UNTIL(outsideValue, i+1, j, k-1, 9);
		    INCREMENT_UNTIL(outsideValue, i, j+1, k-1, 10);
		    INCREMENT_UNTIL(outsideValue, i-1, j-1, k, 11);
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k, 12);
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k, 13);
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k, 14);
		    INCREMENT_UNTIL(outsideValue, i, j-1, k+1, 15);
		    INCREMENT_UNTIL(outsideValue, i-1, j, k+1, 16);
		    INCREMENT_UNTIL(outsideValue, i+1, j, k+1, 17);
		    INCREMENT_UNTIL(outsideValue, i, j+1, k+1, 18);
		    INCREMENT_UNTIL(outsideValue, i-1, j-1, k-1, 19);
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k-1, 20);
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k-1, 21);
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k-1, 22);
		    INCREMENT_UNTIL(outsideValue, i-1, j-1, k+1, 23);
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k+1, 24);
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k+1, 25);
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k+1, 26);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[5].incrementUntil(outsideValue, i, j, k-1);	
		    si[6].incrementUntil(outsideValue, i, j, k+1);
		    si[7].incrementUntil(outsideValue, i, j-1, k-1);
		    si[8].incrementUntil(outsideValue, i-1, j, k-1);
		    si[9].incrementUntil(outsideValue, i+1, j, k-1);
		    si[10].incrementUntil(outsideValue, i, j+1, k-1);
		    si[11].incrementUntil(outsideValue, i-1, j-1, k);
		    si[12].incrementUntil(outsideValue, i+1, j-1, k);
		    si[13].incrementUntil(outsideValue, i-1, j+1, k);
		    si[14].incrementUntil(outsideValue, i+1, j+1, k);
		    si[15].incrementUntil(outsideValue, i, j-1, k+1);
		    si[16].incrementUntil(outsideValue, i-1, j, k+1);
		    si[17].incrementUntil(outsideValue, i+1, j, k+1);
		    si[18].incrementUntil(outsideValue, i, j+1, k+1);
		    si[19].incrementUntil(outsideValue, i-1, j-1, k-1);
		    si[20].incrementUntil(outsideValue, i+1, j-1, k-1);
		    si[21].incrementUntil(outsideValue, i-1, j+1, k-1);
		    si[22].incrementUntil(outsideValue, i+1, j+1, k-1);
		    si[23].incrementUntil(outsideValue, i-1, j-1, k+1);
		    si[24].incrementUntil(outsideValue, i+1, j-1, k+1);
		    si[25].incrementUntil(outsideValue, i-1, j+1, k+1);
		    si[26].incrementUntil(outsideValue, i+1, j+1, k+1);
		}
	    }
	}
    }
    else if (iterStencil == SF_VOXEL)
    {    
	if ( ( tt == BETA_TUBE || tt == ZERO_CROSSING_TUBE )
	    || (Traits::includeSafeBand && tt == GAMMA_TUBE) // this is only valid if the level set includes a safe band!
	    )
	{
	    if (fastUpdate)
	    {
		// In this case we know that the elements of the stencil ALWAYS exist in the narrow band
		si[1].incrementFast();	
		si[2].incrementFast();	
		si[3].incrementFast();	
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i+1, j, k, 1);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k, 2);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k, 3);	
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i+1, j, k);	
		    si[2].incrementUntil(outsideValue, i, j+1, k);	
		    si[3].incrementUntil(outsideValue, i+1, j+1, k);
		}
	    }

	    si[4].getValue() = v1[si[0].getArrayIndex()+1];	
	    si[5].getValue() = v1[si[1].getArrayIndex()+1];
	    si[6].getValue() = v1[si[2].getArrayIndex()+1];
	    si[7].getValue() = v1[si[3].getArrayIndex()+1];

	    if (Traits::maintainArrayIndexInStencil)
	    {
		si[4].getArrayIndex() = si[0].getArrayIndex()+1;	
		si[5].getArrayIndex() = si[1].getArrayIndex()+1;
		si[6].getArrayIndex() = si[2].getArrayIndex()+1;
		si[7].getArrayIndex() = si[3].getArrayIndex()+1;
	    }

	}
	else 
	{
	    if (fastUpdate)
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_FAST_Z(outsideValue, k, 1);	
		    INCREMENT_FAST_Z(outsideValue, k, 2);	
		    INCREMENT_FAST_Z(outsideValue, k, 3);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 4);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 5);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 6);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 7);
		}
		else
		{
		    si[1].incrementFastZ(outsideValue, k);	
		    si[2].incrementFastZ(outsideValue, k);	
		    si[3].incrementFastZ(outsideValue, k);	
		    si[4].incrementFastZ(outsideValue, k+1);	
		    si[5].incrementFastZ(outsideValue, k+1);	
		    si[6].incrementFastZ(outsideValue, k+1);	
		    si[7].incrementFastZ(outsideValue, k+1);	
		}
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i+1, j, k, 1);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k, 2);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k, 3);	
		    INCREMENT_UNTIL(outsideValue, i, j, k+1, 4);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k+1, 5);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k+1, 6);
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k+1, 7);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i+1, j, k);	
		    si[2].incrementUntil(outsideValue, i, j+1, k);	
		    si[3].incrementUntil(outsideValue, i+1, j+1, k);	
		    si[4].incrementUntil(outsideValue, i, j, k+1);	
		    si[5].incrementUntil(outsideValue, i+1, j, k+1);	
		    si[6].incrementUntil(outsideValue, i, j+1, k+1);
		    si[7].incrementUntil(outsideValue, i+1, j+1, k+1);
		}
	    }
	}
    }
    else if (iterStencil == SF_WENO)
    {
	// Note that here we cannot use incrementFast() inside GammaTube
	if ( ( tt == BETA_TUBE || tt == ZERO_CROSSING_TUBE )
	    || (Traits::includeSafeBand && tt == GAMMA_TUBE) // this is only valid if the level set includes a safe band!
	    // It is always safe to use this optimization when includeSafeBand is enabled. This is because of the additional safeband layer of grid points.
	    || (!Traits::includeSafeBand && tt == GAMMA_TUBE && fabs(si[0].getValue())<dx)
	    )
	{
	    if (fastUpdate)
	    {
		// In this case we know that the elements of the stencil ALWAYS exist in the narrow band
		si[1].incrementFast();	
		si[2].incrementFast();	
		si[3].incrementFast();	
		si[4].incrementFast();
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_FAST_Z(outsideValue, k,7);
		    INCREMENT_FAST_Z(outsideValue, k,8);
		    INCREMENT_FAST_Z(outsideValue, k,9);
		    INCREMENT_FAST_Z(outsideValue, k,10);
		    INCREMENT_FAST_Z(outsideValue, k,11);
		    INCREMENT_FAST_Z(outsideValue, k,12);
		    INCREMENT_FAST_Z(outsideValue, k,13);
		    INCREMENT_FAST_Z(outsideValue, k,14);
		    INCREMENT_FAST_Z(outsideValue, k-3,15);
		    INCREMENT_FAST_Z(outsideValue, k-2,16);
		    INCREMENT_FAST_Z(outsideValue, k+2,17);
		    INCREMENT_FAST_Z(outsideValue, k+3,18);
		}
		else
		{
		    si[7].incrementFastZ(outsideValue, k);
		    si[8].incrementFastZ(outsideValue, k);
		    si[9].incrementFastZ(outsideValue, k);
		    si[10].incrementFastZ(outsideValue, k);
		    si[11].incrementFastZ(outsideValue, k);
		    si[12].incrementFastZ(outsideValue, k);
		    si[13].incrementFastZ(outsideValue, k);
		    si[14].incrementFastZ(outsideValue, k);
		    si[15].incrementFastZ(outsideValue, k-3);
		    si[16].incrementFastZ(outsideValue, k-2);
		    si[17].incrementFastZ(outsideValue, k+2);
		    si[18].incrementFastZ(outsideValue, k+3);
		}
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1,j,k,1);
		    INCREMENT_UNTIL(outsideValue, i+1,j,k,2);
		    INCREMENT_UNTIL(outsideValue, i,j-1,k,3);
		    INCREMENT_UNTIL(outsideValue, i,j+1,k,4);
		    INCREMENT_UNTIL(outsideValue, i-3,j,k,7);
		    INCREMENT_UNTIL(outsideValue, i-2,j,k,8);
		    INCREMENT_UNTIL(outsideValue, i+2,j,k,9);
		    INCREMENT_UNTIL(outsideValue, i+3,j,k,10);
		    INCREMENT_UNTIL(outsideValue, i,j-3,k,11);
		    INCREMENT_UNTIL(outsideValue, i,j-2,k,12);
		    INCREMENT_UNTIL(outsideValue, i,j+2,k,13);
		    INCREMENT_UNTIL(outsideValue, i,j+3,k,14);
		    INCREMENT_UNTIL(outsideValue, i,j,k-3,15);
		    INCREMENT_UNTIL(outsideValue, i,j,k-2,16);
		    INCREMENT_UNTIL(outsideValue, i,j,k+2,17);
		    INCREMENT_UNTIL(outsideValue, i,j,k+3,18);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[7].incrementUntil(outsideValue, i-3, j, k);
		    si[8].incrementUntil(outsideValue, i-2, j, k);
		    si[9].incrementUntil(outsideValue, i+2, j, k);
		    si[10].incrementUntil(outsideValue, i+3, j, k);
		    si[11].incrementUntil(outsideValue, i, j-3, k);
		    si[12].incrementUntil(outsideValue, i, j-2, k);
		    si[13].incrementUntil(outsideValue, i, j+2, k);
		    si[14].incrementUntil(outsideValue, i, j+3, k);
		    si[15].incrementUntil(outsideValue, i, j, k-3);
		    si[16].incrementUntil(outsideValue, i, j, k-2);
		    si[17].incrementUntil(outsideValue, i, j, k+2);
		    si[18].incrementUntil(outsideValue, i, j, k+3);
		}
	    }
	    // note that iterators at the elements plus/minus Z are not maintained in this case, only the value is set!
	    si[5].getValue() = v1[si[0].getArrayIndex()-1];	
	    si[6].getValue() = v1[si[0].getArrayIndex()+1];	

	    if (!Traits::includeSafeBand && tt == GAMMA_TUBE)
	    {
		si[5].incrementArrayIndexAndZ();
		si[6].incrementArrayIndexAndZ();
	    }
	    else if (Traits::maintainArrayIndexInStencil)
	    {
		si[5].getArrayIndex() = si[0].getArrayIndex()-1;	
		si[6].getArrayIndex() = si[0].getArrayIndex()+1;	
	    }
	}
	else 
	{
	    if (fastUpdate)
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_FAST_Z(outsideValue, k,1);
		    INCREMENT_FAST_Z(outsideValue, k,2);
		    INCREMENT_FAST_Z(outsideValue, k,3);
		    INCREMENT_FAST_Z(outsideValue, k,4);
		    INCREMENT_FAST_Z(outsideValue, k-1,5);
		    INCREMENT_FAST_Z(outsideValue, k+1,6);
		    INCREMENT_FAST_Z(outsideValue, k,7);
		    INCREMENT_FAST_Z(outsideValue, k,8);
		    INCREMENT_FAST_Z(outsideValue, k,9);
		    INCREMENT_FAST_Z(outsideValue, k,10);
		    INCREMENT_FAST_Z(outsideValue, k,11);
		    INCREMENT_FAST_Z(outsideValue, k,12);
		    INCREMENT_FAST_Z(outsideValue, k,13);
		    INCREMENT_FAST_Z(outsideValue, k,14);
		    INCREMENT_FAST_Z(outsideValue, k-3,15);
		    INCREMENT_FAST_Z(outsideValue, k-2,16);
		    INCREMENT_FAST_Z(outsideValue, k+2,17);
		    INCREMENT_FAST_Z(outsideValue, k+3,18);
		}
		else
		{
		    si[1].incrementFastZ(outsideValue, k);	
		    si[2].incrementFastZ(outsideValue, k);	
		    si[3].incrementFastZ(outsideValue, k);	
		    si[4].incrementFastZ(outsideValue, k);	
		    si[5].incrementFastZ(outsideValue, k-1);	
		    si[6].incrementFastZ(outsideValue, k+1);
		    si[7].incrementFastZ(outsideValue, k);
		    si[8].incrementFastZ(outsideValue, k);
		    si[9].incrementFastZ(outsideValue, k);
		    si[10].incrementFastZ(outsideValue, k);
		    si[11].incrementFastZ(outsideValue, k);
		    si[12].incrementFastZ(outsideValue, k);
		    si[13].incrementFastZ(outsideValue, k);
		    si[14].incrementFastZ(outsideValue, k);
		    si[15].incrementFastZ(outsideValue, k-3);
		    si[16].incrementFastZ(outsideValue, k-2);
		    si[17].incrementFastZ(outsideValue, k+2);
		    si[18].incrementFastZ(outsideValue, k+3);
		}
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1,j,k,1);
		    INCREMENT_UNTIL(outsideValue, i+1,j,k,2);
		    INCREMENT_UNTIL(outsideValue, i,j-1,k,3);
		    INCREMENT_UNTIL(outsideValue, i,j+1,k,4);
		    INCREMENT_UNTIL(outsideValue, i,j,k-1,5);
		    INCREMENT_UNTIL(outsideValue, i,j,k+1,6);
		    INCREMENT_UNTIL(outsideValue, i-3,j,k,7);
		    INCREMENT_UNTIL(outsideValue, i-2,j,k,8);
		    INCREMENT_UNTIL(outsideValue, i+2,j,k,9);
		    INCREMENT_UNTIL(outsideValue, i+3,j,k,10);
		    INCREMENT_UNTIL(outsideValue, i,j-3,k,11);
		    INCREMENT_UNTIL(outsideValue, i,j-2,k,12);
		    INCREMENT_UNTIL(outsideValue, i,j+2,k,13);
		    INCREMENT_UNTIL(outsideValue, i,j+3,k,14);
		    INCREMENT_UNTIL(outsideValue, i,j,k-3,15);
		    INCREMENT_UNTIL(outsideValue, i,j,k-2,16);
		    INCREMENT_UNTIL(outsideValue, i,j,k+2,17);
		    INCREMENT_UNTIL(outsideValue, i,j,k+3,18);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[5].incrementUntil(outsideValue, i, j, k-1);	
		    si[6].incrementUntil(outsideValue, i, j, k+1);
		    si[7].incrementUntil(outsideValue, i-3, j, k);
		    si[8].incrementUntil(outsideValue, i-2, j, k);
		    si[9].incrementUntil(outsideValue, i+2, j, k);
		    si[10].incrementUntil(outsideValue, i+3, j, k);
		    si[11].incrementUntil(outsideValue, i, j-3, k);
		    si[12].incrementUntil(outsideValue, i, j-2, k);
		    si[13].incrementUntil(outsideValue, i, j+2, k);
		    si[14].incrementUntil(outsideValue, i, j+3, k);
		    si[15].incrementUntil(outsideValue, i, j, k-3);
		    si[16].incrementUntil(outsideValue, i, j, k-2);
		    si[17].incrementUntil(outsideValue, i, j, k+2);
		    si[18].incrementUntil(outsideValue, i, j, k+3);
		}
	    }
	}
    }
    else if (iterStencil == SF_WENO_CURVATURE)
    {
	if ( ( tt == BETA_TUBE || tt == ZERO_CROSSING_TUBE )
	    || (Traits::includeSafeBand && tt == GAMMA_TUBE) // this is only valid if the level set includes a safe band!
	    )
	{
	    if (fastUpdate)
	    {
		// In this case we know that the elements of the stencil ALWAYS exist in the narrow band
		si[1].incrementFast();	
		si[2].incrementFast();	
		si[3].incrementFast();	
		si[4].incrementFast();
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_FAST_Z(outsideValue, k,7);
		    INCREMENT_FAST_Z(outsideValue, k,8);
		    INCREMENT_FAST_Z(outsideValue, k,9);
		    INCREMENT_FAST_Z(outsideValue, k,10);
		    INCREMENT_FAST_Z(outsideValue, k,11);
		    INCREMENT_FAST_Z(outsideValue, k,12);
		    INCREMENT_FAST_Z(outsideValue, k,13);
		    INCREMENT_FAST_Z(outsideValue, k,14);
		    INCREMENT_FAST_Z(outsideValue, k-3,15);
		    INCREMENT_FAST_Z(outsideValue, k-2,16);
		    INCREMENT_FAST_Z(outsideValue, k+2,17);
		    INCREMENT_FAST_Z(outsideValue, k+3,18);
		}
		else
		{
		    si[7].incrementFastZ(outsideValue, k);
		    si[8].incrementFastZ(outsideValue, k);
		    si[9].incrementFastZ(outsideValue, k);
		    si[10].incrementFastZ(outsideValue, k);
		    si[11].incrementFastZ(outsideValue, k);
		    si[12].incrementFastZ(outsideValue, k);
		    si[13].incrementFastZ(outsideValue, k);
		    si[14].incrementFastZ(outsideValue, k);
		    si[15].incrementFastZ(outsideValue, k-3);
		    si[16].incrementFastZ(outsideValue, k-2);
		    si[17].incrementFastZ(outsideValue, k+2);
		    si[18].incrementFastZ(outsideValue, k+3);
		}
		si[23].incrementFast();
		si[24].incrementFast();
		si[25].incrementFast();
		si[26].incrementFast();
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1, j, k,1);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k,2);	
		    INCREMENT_UNTIL(outsideValue, i, j-1, k,3);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k,4);	
		    INCREMENT_UNTIL(outsideValue, i-3, j, k,7);
		    INCREMENT_UNTIL(outsideValue, i-2, j, k,8);
		    INCREMENT_UNTIL(outsideValue, i+2, j, k,9);
		    INCREMENT_UNTIL(outsideValue, i+3, j, k,10);
		    INCREMENT_UNTIL(outsideValue, i, j-3, k,11);
		    INCREMENT_UNTIL(outsideValue, i, j-2, k,12);
		    INCREMENT_UNTIL(outsideValue, i, j+2, k,13);
		    INCREMENT_UNTIL(outsideValue, i, j+3, k,14);
		    INCREMENT_UNTIL(outsideValue, i, j, k-3,15);
		    INCREMENT_UNTIL(outsideValue, i, j, k-2,16);
		    INCREMENT_UNTIL(outsideValue, i, j, k+2,17);
		    INCREMENT_UNTIL(outsideValue, i, j, k+3,18);
		    INCREMENT_UNTIL(outsideValue, i-1, j-1, k,23);
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k,24);
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k,25);
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k,26);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[7].incrementUntil(outsideValue, i-3, j, k);
		    si[8].incrementUntil(outsideValue, i-2, j, k);
		    si[9].incrementUntil(outsideValue, i+2, j, k);
		    si[10].incrementUntil(outsideValue, i+3, j, k);
		    si[11].incrementUntil(outsideValue, i, j-3, k);
		    si[12].incrementUntil(outsideValue, i, j-2, k);
		    si[13].incrementUntil(outsideValue, i, j+2, k);
		    si[14].incrementUntil(outsideValue, i, j+3, k);
		    si[15].incrementUntil(outsideValue, i, j, k-3);
		    si[16].incrementUntil(outsideValue, i, j, k-2);
		    si[17].incrementUntil(outsideValue, i, j, k+2);
		    si[18].incrementUntil(outsideValue, i, j, k+3);
		    si[23].incrementUntil(outsideValue, i-1, j-1, k);
		    si[24].incrementUntil(outsideValue, i+1, j-1, k);
		    si[25].incrementUntil(outsideValue, i-1, j+1, k);
		    si[26].incrementUntil(outsideValue, i+1, j+1, k);
		}
	    }
	    // note that iterators at the elements plus/minus Z are not maintained in this case, only the value is set!
	    si[5].getValue() = v1[si[0].getArrayIndex()-1];	
	    si[6].getValue() = v1[si[0].getArrayIndex()+1];	
	    si[19].getValue() = v1[si[3].getArrayIndex()-1];	
	    si[20].getValue() = v1[si[1].getArrayIndex()-1];	
	    si[21].getValue() = v1[si[2].getArrayIndex()-1];	
	    si[22].getValue() = v1[si[4].getArrayIndex()-1];	
	    si[27].getValue() = v1[si[3].getArrayIndex()+1];	
	    si[28].getValue() = v1[si[1].getArrayIndex()+1];	
	    si[29].getValue() = v1[si[2].getArrayIndex()+1];	
	    si[30].getValue() = v1[si[4].getArrayIndex()+1];	

	    if (Traits::maintainArrayIndexInStencil)
	    {
		si[5].getArrayIndex() = si[0].getArrayIndex()-1;	
		si[6].getArrayIndex() = si[0].getArrayIndex()+1;	
		si[19].getArrayIndex() = si[3].getArrayIndex()-1;	
		si[20].getArrayIndex() = si[1].getArrayIndex()-1;	
		si[21].getArrayIndex() = si[2].getArrayIndex()-1;	
		si[22].getArrayIndex() = si[4].getArrayIndex()-1;	
		si[27].getArrayIndex() = si[3].getArrayIndex()+1;	
		si[28].getArrayIndex() = si[1].getArrayIndex()+1;	
		si[29].getArrayIndex() = si[2].getArrayIndex()+1;	
		si[30].getArrayIndex() = si[4].getArrayIndex()+1;	
	    }
	}
	else 
	{
	    if (fastUpdate)
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_FAST_Z(outsideValue, k,1);	
		    INCREMENT_FAST_Z(outsideValue, k,2);	
		    INCREMENT_FAST_Z(outsideValue, k,3);	
		    INCREMENT_FAST_Z(outsideValue, k,4);	
		    INCREMENT_FAST_Z(outsideValue, k-1,5);	
		    INCREMENT_FAST_Z(outsideValue, k+1,6);
		    INCREMENT_FAST_Z(outsideValue, k,7);
		    INCREMENT_FAST_Z(outsideValue, k,8);
		    INCREMENT_FAST_Z(outsideValue, k,9);
		    INCREMENT_FAST_Z(outsideValue, k,10);
		    INCREMENT_FAST_Z(outsideValue, k,11);
		    INCREMENT_FAST_Z(outsideValue, k,12);
		    INCREMENT_FAST_Z(outsideValue, k,13);
		    INCREMENT_FAST_Z(outsideValue, k,14);
		    INCREMENT_FAST_Z(outsideValue, k-3,15);
		    INCREMENT_FAST_Z(outsideValue, k-2,16);
		    INCREMENT_FAST_Z(outsideValue, k+2,17);
		    INCREMENT_FAST_Z(outsideValue, k+3,18);
		    INCREMENT_FAST_Z(outsideValue, k-1,19);
		    INCREMENT_FAST_Z(outsideValue, k-1,20);
		    INCREMENT_FAST_Z(outsideValue, k-1,21);
		    INCREMENT_FAST_Z(outsideValue, k-1,22);
		    INCREMENT_FAST_Z(outsideValue, k,23);
		    INCREMENT_FAST_Z(outsideValue, k,24);
		    INCREMENT_FAST_Z(outsideValue, k,25);
		    INCREMENT_FAST_Z(outsideValue, k,26);
		    INCREMENT_FAST_Z(outsideValue, k+1,27);
		    INCREMENT_FAST_Z(outsideValue, k+1,28);
		    INCREMENT_FAST_Z(outsideValue, k+1,29);
		    INCREMENT_FAST_Z(outsideValue, k+1,30);
		}
		else
		{
		    si[1].incrementFastZ(outsideValue, k);	
		    si[2].incrementFastZ(outsideValue, k);	
		    si[3].incrementFastZ(outsideValue, k);	
		    si[4].incrementFastZ(outsideValue, k);	
		    si[5].incrementFastZ(outsideValue, k-1);	
		    si[6].incrementFastZ(outsideValue, k+1);
		    si[7].incrementFastZ(outsideValue, k);
		    si[8].incrementFastZ(outsideValue, k);
		    si[9].incrementFastZ(outsideValue, k);
		    si[10].incrementFastZ(outsideValue, k);
		    si[11].incrementFastZ(outsideValue, k);
		    si[12].incrementFastZ(outsideValue, k);
		    si[13].incrementFastZ(outsideValue, k);
		    si[14].incrementFastZ(outsideValue, k);
		    si[15].incrementFastZ(outsideValue, k-3);
		    si[16].incrementFastZ(outsideValue, k-2);
		    si[17].incrementFastZ(outsideValue, k+2);
		    si[18].incrementFastZ(outsideValue, k+3);
		    si[19].incrementFastZ(outsideValue, k-1);
		    si[20].incrementFastZ(outsideValue, k-1);
		    si[21].incrementFastZ(outsideValue, k-1);
		    si[22].incrementFastZ(outsideValue, k-1);
		    si[23].incrementFastZ(outsideValue, k);
		    si[24].incrementFastZ(outsideValue, k);
		    si[25].incrementFastZ(outsideValue, k);
		    si[26].incrementFastZ(outsideValue, k);
		    si[27].incrementFastZ(outsideValue, k+1);
		    si[28].incrementFastZ(outsideValue, k+1);
		    si[29].incrementFastZ(outsideValue, k+1);
		    si[30].incrementFastZ(outsideValue, k+1);
		}
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1, j, k,1);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k,2);	
		    INCREMENT_UNTIL(outsideValue, i, j-1, k,3);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k,4);	
		    INCREMENT_UNTIL(outsideValue, i, j, k-1,5);	
		    INCREMENT_UNTIL(outsideValue, i, j, k+1,6);
		    INCREMENT_UNTIL(outsideValue, i-3, j, k,7);
		    INCREMENT_UNTIL(outsideValue, i-2, j, k,8);
		    INCREMENT_UNTIL(outsideValue, i+2, j, k,9);
		    INCREMENT_UNTIL(outsideValue, i+3, j, k,10);
		    INCREMENT_UNTIL(outsideValue, i, j-3, k,11);
		    INCREMENT_UNTIL(outsideValue, i, j-2, k,12);
		    INCREMENT_UNTIL(outsideValue, i, j+2, k,13);
		    INCREMENT_UNTIL(outsideValue, i, j+3, k,14);
		    INCREMENT_UNTIL(outsideValue, i, j, k-3,15);
		    INCREMENT_UNTIL(outsideValue, i, j, k-2,16);
		    INCREMENT_UNTIL(outsideValue, i, j, k+2,17);
		    INCREMENT_UNTIL(outsideValue, i, j, k+3,18);
		    INCREMENT_UNTIL(outsideValue, i, j-1, k-1,19);
		    INCREMENT_UNTIL(outsideValue, i-1, j, k-1,20);
		    INCREMENT_UNTIL(outsideValue, i+1, j, k-1,21);
		    INCREMENT_UNTIL(outsideValue, i, j+1, k-1,22);
		    INCREMENT_UNTIL(outsideValue, i-1, j-1, k,23);
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k,24);
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k,25);
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k,26);
		    INCREMENT_UNTIL(outsideValue, i, j-1, k+1,27);
		    INCREMENT_UNTIL(outsideValue, i-1, j, k+1,28);
		    INCREMENT_UNTIL(outsideValue, i+1, j, k+1,29);
		    INCREMENT_UNTIL(outsideValue, i, j+1, k+1,30);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[5].incrementUntil(outsideValue, i, j, k-1);	
		    si[6].incrementUntil(outsideValue, i, j, k+1);
		    si[7].incrementUntil(outsideValue, i-3, j, k);
		    si[8].incrementUntil(outsideValue, i-2, j, k);
		    si[9].incrementUntil(outsideValue, i+2, j, k);
		    si[10].incrementUntil(outsideValue, i+3, j, k);
		    si[11].incrementUntil(outsideValue, i, j-3, k);
		    si[12].incrementUntil(outsideValue, i, j-2, k);
		    si[13].incrementUntil(outsideValue, i, j+2, k);
		    si[14].incrementUntil(outsideValue, i, j+3, k);
		    si[15].incrementUntil(outsideValue, i, j, k-3);
		    si[16].incrementUntil(outsideValue, i, j, k-2);
		    si[17].incrementUntil(outsideValue, i, j, k+2);
		    si[18].incrementUntil(outsideValue, i, j, k+3);
		    si[19].incrementUntil(outsideValue, i, j-1, k-1);
		    si[20].incrementUntil(outsideValue, i-1, j, k-1);
		    si[21].incrementUntil(outsideValue, i+1, j, k-1);
		    si[22].incrementUntil(outsideValue, i, j+1, k-1);
		    si[23].incrementUntil(outsideValue, i-1, j-1, k);
		    si[24].incrementUntil(outsideValue, i+1, j-1, k);
		    si[25].incrementUntil(outsideValue, i-1, j+1, k);
		    si[26].incrementUntil(outsideValue, i+1, j+1, k);
		    si[27].incrementUntil(outsideValue, i, j-1, k+1);
		    si[28].incrementUntil(outsideValue, i-1, j, k+1);
		    si[29].incrementUntil(outsideValue, i+1, j, k+1);
		    si[30].incrementUntil(outsideValue, i, j+1, k+1);
		}
	    }
	}
    }
    else if (iterStencil == SF_FIRSTORDER_NEIGHBORS)
    {    
	if (Traits::includeSafeBand && ( tt == BETA_TUBE || tt == ZERO_CROSSING_TUBE ) )
	{
	    if (fastUpdate)
	    {
		// In this case we know that the elements of the stencil ALWAYS exist in the narrow band
		si[1].incrementFast();	
		si[2].incrementFast();	
		si[3].incrementFast();	
		si[4].incrementFast();	
		si[5].incrementFast();	
		si[6].incrementFast();	
		si[7].incrementFast();	
		si[8].incrementFast();	
		si[9].incrementFast();	
		si[10].incrementFast();	
		si[11].incrementFast();	
		si[12].incrementFast();	
		si[13].incrementFast();	
		si[14].incrementFast();	
		si[15].incrementFast();	
		si[16].incrementFast();	
		si[17].incrementFast();	
		si[18].incrementFast();	
		si[19].incrementFast();	
		si[20].incrementFast();	
		si[21].incrementFast();	
		si[22].incrementFast();	
		si[23].incrementFast();	
		si[24].incrementFast();	
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1, j, k, 1);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k, 2);	
		    INCREMENT_UNTIL(outsideValue, i, j-1, k, 3);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k, 4);	
		    INCREMENT_UNTIL(outsideValue, i, j, k-1, 5);	
		    INCREMENT_UNTIL(outsideValue, i, j, k+1, 6);	
		    INCREMENT_UNTIL(outsideValue, i, j, k-2, 7);	
		    INCREMENT_UNTIL(outsideValue, i, j-1, k-1, 8);	
		    INCREMENT_UNTIL(outsideValue, i-1, j, k-1, 9);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k-1, 10);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k-1, 11);
		    INCREMENT_UNTIL(outsideValue, i, j-2, k, 12);	
		    INCREMENT_UNTIL(outsideValue, i-1, j-1, k, 13);	
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k, 14);	
		    INCREMENT_UNTIL(outsideValue, i-2, j, k, 15);	
		    INCREMENT_UNTIL(outsideValue, i+2, j, k, 16);	
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k, 17);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k, 18);	
		    INCREMENT_UNTIL(outsideValue, i, j+2, k, 19);
		    INCREMENT_UNTIL(outsideValue, i, j-1, k+1, 20);	
		    INCREMENT_UNTIL(outsideValue, i-1, j, k+1, 21);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k+1, 22);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k+1, 23);
		    INCREMENT_UNTIL(outsideValue, i, j, k+2, 24);	
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[5].incrementUntil(outsideValue, i, j, k-1);	
		    si[6].incrementUntil(outsideValue, i, j, k+1);	
		    si[7].incrementUntil(outsideValue, i, j, k-2);	
		    si[8].incrementUntil(outsideValue, i, j-1, k-1);	
		    si[9].incrementUntil(outsideValue, i-1, j, k-1);	
		    si[10].incrementUntil(outsideValue, i+1, j, k-1);	
		    si[11].incrementUntil(outsideValue, i, j+1, k-1);
		    si[12].incrementUntil(outsideValue, i, j-2, k);	
		    si[13].incrementUntil(outsideValue, i-1, j-1, k);	
		    si[14].incrementUntil(outsideValue, i+1, j-1, k);	
		    si[15].incrementUntil(outsideValue, i-2, j, k);	
		    si[16].incrementUntil(outsideValue, i+2, j, k);	
		    si[17].incrementUntil(outsideValue, i-1, j+1, k);	
		    si[18].incrementUntil(outsideValue, i+1, j+1, k);	
		    si[19].incrementUntil(outsideValue, i, j+2, k);
		    si[20].incrementUntil(outsideValue, i, j-1, k+1);	
		    si[21].incrementUntil(outsideValue, i-1, j, k+1);	
		    si[22].incrementUntil(outsideValue, i+1, j, k+1);	
		    si[23].incrementUntil(outsideValue, i, j+1, k+1);
		    si[24].incrementUntil(outsideValue, i, j, k+2);	
		}
	    }
	}
	else 
	{
	    if (fastUpdate)
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_FAST_Z(outsideValue, k, 1);	
		    INCREMENT_FAST_Z(outsideValue, k, 2);	
		    INCREMENT_FAST_Z(outsideValue, k, 3);	
		    INCREMENT_FAST_Z(outsideValue, k, 4);	
		    INCREMENT_FAST_Z(outsideValue, k-1, 5);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 6);	
		    INCREMENT_FAST_Z(outsideValue, k-2, 7);	
		    INCREMENT_FAST_Z(outsideValue, k-1, 8);	
		    INCREMENT_FAST_Z(outsideValue, k-1, 9);	
		    INCREMENT_FAST_Z(outsideValue, k-1, 10);	
		    INCREMENT_FAST_Z(outsideValue, k-1, 11);
		    INCREMENT_FAST_Z(outsideValue, k, 12);	
		    INCREMENT_FAST_Z(outsideValue, k, 13);	
		    INCREMENT_FAST_Z(outsideValue, k, 14);	
		    INCREMENT_FAST_Z(outsideValue, k, 15);	
		    INCREMENT_FAST_Z(outsideValue, k, 16);	
		    INCREMENT_FAST_Z(outsideValue, k, 17);	
		    INCREMENT_FAST_Z(outsideValue, k, 18);	
		    INCREMENT_FAST_Z(outsideValue, k, 19);
		    INCREMENT_FAST_Z(outsideValue, k+1, 20);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 21);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 22);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 23);	
		    INCREMENT_FAST_Z(outsideValue, k+2, 24);	
		}
		else
		{
		    si[1].incrementFastZ(outsideValue, k);	
		    si[2].incrementFastZ(outsideValue, k);	
		    si[3].incrementFastZ(outsideValue, k);	
		    si[4].incrementFastZ(outsideValue, k);	
		    si[5].incrementFastZ(outsideValue, k-1);	
		    si[6].incrementFastZ(outsideValue, k+1);	
		    si[7].incrementFastZ(outsideValue, k-2);	
		    si[8].incrementFastZ(outsideValue, k-1);	
		    si[9].incrementFastZ(outsideValue, k-1);	
		    si[10].incrementFastZ(outsideValue, k-1);	
		    si[11].incrementFastZ(outsideValue, k-1);
		    si[12].incrementFastZ(outsideValue, k);	
		    si[13].incrementFastZ(outsideValue, k);	
		    si[14].incrementFastZ(outsideValue, k);	
		    si[15].incrementFastZ(outsideValue, k);	
		    si[16].incrementFastZ(outsideValue, k);	
		    si[17].incrementFastZ(outsideValue, k);	
		    si[18].incrementFastZ(outsideValue, k);	
		    si[19].incrementFastZ(outsideValue, k);
		    si[20].incrementFastZ(outsideValue, k+1);	
		    si[21].incrementFastZ(outsideValue, k+1);	
		    si[22].incrementFastZ(outsideValue, k+1);	
		    si[23].incrementFastZ(outsideValue, k+1);	
		    si[24].incrementFastZ(outsideValue, k+2);	
		}
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i-1, j, k, 1);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k, 2);	
		    INCREMENT_UNTIL(outsideValue, i, j-1, k, 3);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k, 4);	
		    INCREMENT_UNTIL(outsideValue, i, j, k-1, 5);	
		    INCREMENT_UNTIL(outsideValue, i, j, k+1, 6);	
		    INCREMENT_UNTIL(outsideValue, i, j, k-2, 7);	
		    INCREMENT_UNTIL(outsideValue, i, j-1, k-1, 8);	
		    INCREMENT_UNTIL(outsideValue, i-1, j, k-1, 9);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k-1, 10);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k-1, 11);
		    INCREMENT_UNTIL(outsideValue, i, j-2, k, 12);	
		    INCREMENT_UNTIL(outsideValue, i-1, j-1, k, 13);	
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k, 14);	
		    INCREMENT_UNTIL(outsideValue, i-2, j, k, 15);	
		    INCREMENT_UNTIL(outsideValue, i+2, j, k, 16);	
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k, 17);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k, 18);	
		    INCREMENT_UNTIL(outsideValue, i, j+2, k, 19);
		    INCREMENT_UNTIL(outsideValue, i, j-1, k+1, 20);	
		    INCREMENT_UNTIL(outsideValue, i-1, j, k+1, 21);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k+1, 22);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k+1, 23);
		    INCREMENT_UNTIL(outsideValue, i, j, k+2, 24);	
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i-1, j, k);	
		    si[2].incrementUntil(outsideValue, i+1, j, k);	
		    si[3].incrementUntil(outsideValue, i, j-1, k);	
		    si[4].incrementUntil(outsideValue, i, j+1, k);	
		    si[5].incrementUntil(outsideValue, i, j, k-1);	
		    si[6].incrementUntil(outsideValue, i, j, k+1);	
		    si[7].incrementUntil(outsideValue, i, j, k-2);	
		    si[8].incrementUntil(outsideValue, i, j-1, k-1);	
		    si[9].incrementUntil(outsideValue, i-1, j, k-1);	
		    si[10].incrementUntil(outsideValue, i+1, j, k-1);	
		    si[11].incrementUntil(outsideValue, i, j+1, k-1);
		    si[12].incrementUntil(outsideValue, i, j-2, k);	
		    si[13].incrementUntil(outsideValue, i-1, j-1, k);	
		    si[14].incrementUntil(outsideValue, i+1, j-1, k);	
		    si[15].incrementUntil(outsideValue, i-2, j, k);	
		    si[16].incrementUntil(outsideValue, i+2, j, k);	
		    si[17].incrementUntil(outsideValue, i-1, j+1, k);	
		    si[18].incrementUntil(outsideValue, i+1, j+1, k);	
		    si[19].incrementUntil(outsideValue, i, j+2, k);
		    si[20].incrementUntil(outsideValue, i, j-1, k+1);	
		    si[21].incrementUntil(outsideValue, i-1, j, k+1);	
		    si[22].incrementUntil(outsideValue, i+1, j, k+1);	
		    si[23].incrementUntil(outsideValue, i, j+1, k+1);
		    si[24].incrementUntil(outsideValue, i, j, k+2);	
		}
	    }
	}
    }
    else if (iterStencil == SF_VOXEL_GRAD)
    {    
	if ( Traits::includeSafeBand && ( tt == BETA_TUBE || tt == ZERO_CROSSING_TUBE ) )
	{
	    if (fastUpdate)
	    {
		// In this case we know that the elements of the stencil ALWAYS exist in the narrow band
		si[1].incrementFast();	
		si[2].incrementFast();	
		si[3].incrementFast();	
		si[4].incrementFast();	
		si[5].incrementFast();	
		si[6].incrementFast();	
		si[7].incrementFast();	
		si[8].incrementFast();	
		si[9].incrementFast();	
		si[10].incrementFast();	
		si[11].incrementFast();	
		si[12].incrementFast();	
		si[13].incrementFast();	
		si[14].incrementFast();	
		si[15].incrementFast();	
		si[16].incrementFast();	
		si[17].incrementFast();	
		si[18].incrementFast();	
		si[19].incrementFast();	
		si[20].incrementFast();	
		si[21].incrementFast();	
		si[22].incrementFast();	
		si[23].incrementFast();	
		si[24].incrementFast();	
		si[25].incrementFast();	
		si[26].incrementFast();	
		si[27].incrementFast();	
		si[28].incrementFast();	
		si[29].incrementFast();	
		si[30].incrementFast();	
		si[31].incrementFast();	
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i+1, j, k, 1);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k, 2);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k, 3);	
		    INCREMENT_UNTIL(outsideValue, i, j, k+1, 4);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k+1, 5);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k+1, 6);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k+1, 7);	
		    INCREMENT_UNTIL(outsideValue, i, j, k-1, 8);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k-1, 9);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k-1, 10);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k-1, 11);
		    INCREMENT_UNTIL(outsideValue, i, j-1, k, 12);	
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k, 13);	
		    INCREMENT_UNTIL(outsideValue, i-1, j, k, 14);	
		    INCREMENT_UNTIL(outsideValue, i+2, j, k, 15);	
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k, 16);	
		    INCREMENT_UNTIL(outsideValue, i+2, j+1, k, 17);	
		    INCREMENT_UNTIL(outsideValue, i, j+2, k, 18);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+2, k, 19);
		    INCREMENT_UNTIL(outsideValue, i, j-1, k+1, 20);	
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k+1, 21);	
		    INCREMENT_UNTIL(outsideValue, i-1, j, k+1, 22);	
		    INCREMENT_UNTIL(outsideValue, i+2, j, k+1, 23);	
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k+1, 24);	
		    INCREMENT_UNTIL(outsideValue, i+2, j+1, k+1, 25);	
		    INCREMENT_UNTIL(outsideValue, i, j+2, k+1, 26);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+2, k+1, 27);
		    INCREMENT_UNTIL(outsideValue, i, j, k+2, 28);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k+2, 29);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k+2, 30);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k+2, 31);
		}
		else
		{
		    si[1].incrementUntil(outsideValue, i+1, j, k);	
		    si[2].incrementUntil(outsideValue, i, j+1, k);	
		    si[3].incrementUntil(outsideValue, i+1, j+1, k);	
		    si[4].incrementUntil(outsideValue, i, j, k+1);	
		    si[5].incrementUntil(outsideValue, i+1, j, k+1);	
		    si[6].incrementUntil(outsideValue, i, j+1, k+1);	
		    si[7].incrementUntil(outsideValue, i+1, j+1, k+1);	
		    si[8].incrementUntil(outsideValue, i, j, k-1);	
		    si[9].incrementUntil(outsideValue, i+1, j, k-1);	
		    si[10].incrementUntil(outsideValue, i, j+1, k-1);	
		    si[11].incrementUntil(outsideValue, i+1, j+1, k-1);
		    si[12].incrementUntil(outsideValue, i, j-1, k);	
		    si[13].incrementUntil(outsideValue, i+1, j-1, k);	
		    si[14].incrementUntil(outsideValue, i-1, j, k);	
		    si[15].incrementUntil(outsideValue, i+2, j, k);	
		    si[16].incrementUntil(outsideValue, i-1, j+1, k);	
		    si[17].incrementUntil(outsideValue, i+2, j+1, k);	
		    si[18].incrementUntil(outsideValue, i, j+2, k);	
		    si[19].incrementUntil(outsideValue, i+1, j+2, k);
		    si[20].incrementUntil(outsideValue, i, j-1, k+1);	
		    si[21].incrementUntil(outsideValue, i+1, j-1, k+1);	
		    si[22].incrementUntil(outsideValue, i-1, j, k+1);	
		    si[23].incrementUntil(outsideValue, i+2, j, k+1);	
		    si[24].incrementUntil(outsideValue, i-1, j+1, k+1);	
		    si[25].incrementUntil(outsideValue, i+2, j+1, k+1);	
		    si[26].incrementUntil(outsideValue, i, j+2, k+1);	
		    si[27].incrementUntil(outsideValue, i+1, j+2, k+1);
		    si[28].incrementUntil(outsideValue, i, j, k+2);	
		    si[29].incrementUntil(outsideValue, i+1, j, k+2);	
		    si[30].incrementUntil(outsideValue, i, j+1, k+2);	
		    si[31].incrementUntil(outsideValue, i+1, j+1, k+2);
		}
	    }
	}
	else 
	{
	    if (fastUpdate)
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_FAST_Z(outsideValue, k, 1);	
		    INCREMENT_FAST_Z(outsideValue, k, 2);	
		    INCREMENT_FAST_Z(outsideValue, k, 3);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 4);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 5);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 6);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 7);	
		    INCREMENT_FAST_Z(outsideValue, k-1, 8);	
		    INCREMENT_FAST_Z(outsideValue, k-1, 9);	
		    INCREMENT_FAST_Z(outsideValue, k-1, 10);	
		    INCREMENT_FAST_Z(outsideValue, k-1, 11);
		    INCREMENT_FAST_Z(outsideValue, k, 12);	
		    INCREMENT_FAST_Z(outsideValue, k, 13);	
		    INCREMENT_FAST_Z(outsideValue, k, 14);	
		    INCREMENT_FAST_Z(outsideValue, k, 15);	
		    INCREMENT_FAST_Z(outsideValue, k, 16);	
		    INCREMENT_FAST_Z(outsideValue, k, 17);	
		    INCREMENT_FAST_Z(outsideValue, k, 18);	
		    INCREMENT_FAST_Z(outsideValue, k, 19);
		    INCREMENT_FAST_Z(outsideValue, k+1, 20);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 21);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 22);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 23);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 24);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 25);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 26);	
		    INCREMENT_FAST_Z(outsideValue, k+1, 27);
		    INCREMENT_FAST_Z(outsideValue, k+2, 28);	
		    INCREMENT_FAST_Z(outsideValue, k+2, 29);	
		    INCREMENT_FAST_Z(outsideValue, k+2, 30);	
		    INCREMENT_FAST_Z(outsideValue, k+2, 31);	
		}
		else
		{
		    si[1].incrementFastZ(outsideValue, k);	
		    si[2].incrementFastZ(outsideValue, k);	
		    si[3].incrementFastZ(outsideValue, k);	
		    si[4].incrementFastZ(outsideValue, k+1);	
		    si[5].incrementFastZ(outsideValue, k+1);	
		    si[6].incrementFastZ(outsideValue, k+1);	
		    si[7].incrementFastZ(outsideValue, k+1);	
		    si[8].incrementFastZ(outsideValue, k-1);	
		    si[9].incrementFastZ(outsideValue, k-1);	
		    si[10].incrementFastZ(outsideValue, k-1);	
		    si[11].incrementFastZ(outsideValue, k-1);
		    si[12].incrementFastZ(outsideValue, k);	
		    si[13].incrementFastZ(outsideValue, k);	
		    si[14].incrementFastZ(outsideValue, k);	
		    si[15].incrementFastZ(outsideValue, k);	
		    si[16].incrementFastZ(outsideValue, k);	
		    si[17].incrementFastZ(outsideValue, k);	
		    si[18].incrementFastZ(outsideValue, k);	
		    si[19].incrementFastZ(outsideValue, k);
		    si[20].incrementFastZ(outsideValue, k+1);	
		    si[21].incrementFastZ(outsideValue, k+1);	
		    si[22].incrementFastZ(outsideValue, k+1);	
		    si[23].incrementFastZ(outsideValue, k+1);	
		    si[24].incrementFastZ(outsideValue, k+1);	
		    si[25].incrementFastZ(outsideValue, k+1);	
		    si[26].incrementFastZ(outsideValue, k+1);	
		    si[27].incrementFastZ(outsideValue, k+1);
		    si[28].incrementFastZ(outsideValue, k+2);	
		    si[29].incrementFastZ(outsideValue, k+2);	
		    si[30].incrementFastZ(outsideValue, k+2);	
		    si[31].incrementFastZ(outsideValue, k+2);	
		}
	    }
	    else
	    {
		if (Traits::inlineUpdateStencilCalls)
		{
		    INCREMENT_UNTIL(outsideValue, i+1, j, k, 1);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k, 2);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k, 3);	
		    INCREMENT_UNTIL(outsideValue, i, j, k+1, 4);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k+1, 5);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k+1, 6);
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k+1, 7);	    
		    INCREMENT_UNTIL(outsideValue, i, j, k-1, 8);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k-1, 9);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k-1, 10);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k-1, 11);
		    INCREMENT_UNTIL(outsideValue, i, j-1, k, 12);	
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k, 13);	
		    INCREMENT_UNTIL(outsideValue, i-1, j, k, 14);	
		    INCREMENT_UNTIL(outsideValue, i+2, j, k, 15);	
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k, 16);	
		    INCREMENT_UNTIL(outsideValue, i+2, j+1, k, 17);	
		    INCREMENT_UNTIL(outsideValue, i, j+2, k, 18);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+2, k, 19);
		    INCREMENT_UNTIL(outsideValue, i, j-1, k+1, 20);	
		    INCREMENT_UNTIL(outsideValue, i+1, j-1, k+1, 21);	
		    INCREMENT_UNTIL(outsideValue, i-1, j, k+1, 22);	
		    INCREMENT_UNTIL(outsideValue, i+2, j, k+1, 23);	
		    INCREMENT_UNTIL(outsideValue, i-1, j+1, k+1, 24);	
		    INCREMENT_UNTIL(outsideValue, i+2, j+1, k+1, 25);	
		    INCREMENT_UNTIL(outsideValue, i, j+2, k+1, 26);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+2, k+1, 27);
		    INCREMENT_UNTIL(outsideValue, i, j, k+2, 28);	
		    INCREMENT_UNTIL(outsideValue, i+1, j, k+2, 29);	
		    INCREMENT_UNTIL(outsideValue, i, j+1, k+2, 30);	
		    INCREMENT_UNTIL(outsideValue, i+1, j+1, k+2, 31);
		}
		else
		{

		    //    Z slices increasing to right, X increasing to right, Y increasing upwards
		    //
		    //                              18  19              26 27 
		    //                10 11     16   2   3   17     24   6  7  25     30 31
		    //                8   9     14   0   1   15     22   4  5  23     28 29          
		    //      	                    12  13              20 21               

		    si[1].incrementUntil(outsideValue, i+1, j, k);	
		    si[2].incrementUntil(outsideValue, i, j+1, k);	
		    si[3].incrementUntil(outsideValue, i+1, j+1, k);	
		    si[4].incrementUntil(outsideValue, i, j, k+1);	
		    si[5].incrementUntil(outsideValue, i+1, j, k+1);	
		    si[6].incrementUntil(outsideValue, i, j+1, k+1);
		    si[7].incrementUntil(outsideValue, i+1, j+1, k+1);	    
		    si[8].incrementUntil(outsideValue, i, j, k-1);	
		    si[9].incrementUntil(outsideValue, i+1, j, k-1);	
		    si[10].incrementUntil(outsideValue, i, j+1, k-1);	
		    si[11].incrementUntil(outsideValue, i+1, j+1, k-1);
		    si[12].incrementUntil(outsideValue, i, j-1, k);	
		    si[13].incrementUntil(outsideValue, i+1, j-1, k);	
		    si[14].incrementUntil(outsideValue, i-1, j, k);	
		    si[15].incrementUntil(outsideValue, i+2, j, k);	
		    si[16].incrementUntil(outsideValue, i-1, j+1, k);	
		    si[17].incrementUntil(outsideValue, i+2, j+1, k);	
		    si[18].incrementUntil(outsideValue, i, j+2, k);	
		    si[19].incrementUntil(outsideValue, i+1, j+2, k);
		    si[20].incrementUntil(outsideValue, i, j-1, k+1);	
		    si[21].incrementUntil(outsideValue, i+1, j-1, k+1);	
		    si[22].incrementUntil(outsideValue, i-1, j, k+1);	
		    si[23].incrementUntil(outsideValue, i+2, j, k+1);	
		    si[24].incrementUntil(outsideValue, i-1, j+1, k+1);	
		    si[25].incrementUntil(outsideValue, i+2, j+1, k+1);	
		    si[26].incrementUntil(outsideValue, i, j+2, k+1);	
		    si[27].incrementUntil(outsideValue, i+1, j+2, k+1);
		    si[28].incrementUntil(outsideValue, i, j, k+2);	
		    si[29].incrementUntil(outsideValue, i+1, j, k+2);	
		    si[30].incrementUntil(outsideValue, i, j+1, k+2);	
		    si[31].incrementUntil(outsideValue, i+1, j+1, k+2);
		}
	    }
	}
    }
}










template<class Traits>  template<TubeType tt>  template<StencilFormat iterStencil>
void DTGrid<Traits>::StencilTubeIterator<tt>::initStencilUsingRandomAccess(UInt bufferId)
{	
    Locator loc;
    bool valid, inCorrectXYColumn;

    if (iterStencil == SF_FIRSTORDER)
    {
	valid = parent->getExistingLocator(i-1, j, k, &loc, &inCorrectXYColumn); si[1] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);	
	valid = parent->getExistingLocator(i+1, j, k, &loc, &inCorrectXYColumn); si[2] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);	
	valid = parent->getExistingLocator(i, j-1, k, &loc, &inCorrectXYColumn); si[3] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);	
	valid = parent->getExistingLocator(i, j+1, k, &loc, &inCorrectXYColumn); si[4] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);	
	valid = parent->getExistingLocator(i, j, k-1, &loc, &inCorrectXYColumn); si[5] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);	
	valid = parent->getExistingLocator(i, j, k+1, &loc, &inCorrectXYColumn); si[6] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
    } 
    else if (iterStencil == SF_FIRSTORDER_CURVATURE)
    {
	valid = parent->getExistingLocator(i-1, j, k, &loc, &inCorrectXYColumn); si[1] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);	
	valid = parent->getExistingLocator(i+1, j, k, &loc, &inCorrectXYColumn);	si[2] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-1, k, &loc, &inCorrectXYColumn);	si[3] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k, &loc, &inCorrectXYColumn);	si[4] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k-1, &loc, &inCorrectXYColumn);	si[5] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k+1, &loc, &inCorrectXYColumn);  si[6] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-1, k-1, &loc, &inCorrectXYColumn);  si[7] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j, k-1, &loc, &inCorrectXYColumn);  si[8] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j, k-1, &loc, &inCorrectXYColumn);  si[9] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k-1, &loc, &inCorrectXYColumn);  si[10] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j-1, k, &loc, &inCorrectXYColumn);  si[11] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j-1, k, &loc, &inCorrectXYColumn);  si[12] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j+1, k, &loc, &inCorrectXYColumn);  si[13] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j+1, k, &loc, &inCorrectXYColumn);  si[14] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-1, k+1, &loc, &inCorrectXYColumn);  si[15] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j, k+1, &loc, &inCorrectXYColumn);  si[16] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j, k+1, &loc, &inCorrectXYColumn);  si[17] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k+1, &loc, &inCorrectXYColumn);  si[18] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
    }
    else if (iterStencil == SF_BOX)
    {
	valid = parent->getExistingLocator(i-1, j, k, &loc, &inCorrectXYColumn);	si[1] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j, k, &loc, &inCorrectXYColumn);	si[2] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-1, k, &loc, &inCorrectXYColumn);	si[3] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k, &loc, &inCorrectXYColumn);	si[4] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k-1, &loc, &inCorrectXYColumn);	si[5] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k+1, &loc, &inCorrectXYColumn);  si[6] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-1, k-1, &loc, &inCorrectXYColumn);  si[7] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j, k-1, &loc, &inCorrectXYColumn);  si[8] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j, k-1, &loc, &inCorrectXYColumn);  si[9] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k-1, &loc, &inCorrectXYColumn);  si[10] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j-1, k, &loc, &inCorrectXYColumn);  si[11] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j-1, k, &loc, &inCorrectXYColumn);  si[12] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j+1, k, &loc, &inCorrectXYColumn);  si[13] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j+1, k, &loc, &inCorrectXYColumn);  si[14] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-1, k+1, &loc, &inCorrectXYColumn);  si[15] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j, k+1, &loc, &inCorrectXYColumn);  si[16] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j, k+1, &loc, &inCorrectXYColumn);  si[17] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k+1, &loc, &inCorrectXYColumn);  si[18] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j-1, k-1, &loc, &inCorrectXYColumn);  si[19] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j-1, k-1, &loc, &inCorrectXYColumn);  si[20] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j+1, k-1, &loc, &inCorrectXYColumn);  si[21] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j+1, k-1, &loc, &inCorrectXYColumn);  si[22] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j-1, k+1, &loc, &inCorrectXYColumn);  si[23] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j-1, k+1, &loc, &inCorrectXYColumn);  si[24] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j+1, k+1, &loc, &inCorrectXYColumn);  si[25] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j+1, k+1, &loc, &inCorrectXYColumn);  si[26] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
    }
    else if (iterStencil == SF_VOXEL)
    {    
	valid = parent->getExistingLocator(i-1, j, k, &loc, &inCorrectXYColumn);	si[1] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j, k, &loc, &inCorrectXYColumn);	si[2] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-1, k, &loc, &inCorrectXYColumn);	si[3] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k, &loc, &inCorrectXYColumn);	si[4] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k-1, &loc, &inCorrectXYColumn);	si[5] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k+1, &loc, &inCorrectXYColumn);  si[6] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-3, j, k, &loc, &inCorrectXYColumn);  si[7] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-2, j, k, &loc, &inCorrectXYColumn);  si[8] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+2, j, k, &loc, &inCorrectXYColumn);  si[9] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+3, j, k, &loc, &inCorrectXYColumn);  si[10] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-3, k, &loc, &inCorrectXYColumn);  si[11] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-2, k, &loc, &inCorrectXYColumn);  si[12] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+2, k, &loc, &inCorrectXYColumn);  si[13] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+3, k, &loc, &inCorrectXYColumn);  si[14] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k-3, &loc, &inCorrectXYColumn);  si[15] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k-2, &loc, &inCorrectXYColumn);  si[16] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k+2, &loc, &inCorrectXYColumn);  si[17] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k+3, &loc, &inCorrectXYColumn);  si[18] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
    }
   else if (iterStencil == SF_VOXEL_GRAD)
    {    
	valid = parent->getExistingLocator(i+1, j, k, &loc, &inCorrectXYColumn);	si[1] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k, &loc, &inCorrectXYColumn);	si[2] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j+1, k, &loc, &inCorrectXYColumn);	si[3] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k+1, &loc, &inCorrectXYColumn);	si[4] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j, k+1, &loc, &inCorrectXYColumn);	si[5] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k+1, &loc, &inCorrectXYColumn);  si[6] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j+1, k+1, &loc, &inCorrectXYColumn);  si[7] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k-1, &loc, &inCorrectXYColumn);	si[8] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j, k-1, &loc, &inCorrectXYColumn);	si[9] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k-1, &loc, &inCorrectXYColumn);	si[10] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j+1, k-1, &loc, &inCorrectXYColumn);  si[11] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-1, k, &loc, &inCorrectXYColumn);	si[12] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j-1, k, &loc, &inCorrectXYColumn);	si[13] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j, k, &loc, &inCorrectXYColumn);	si[14] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+2, j, k, &loc, &inCorrectXYColumn);	si[15] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j+1, k, &loc, &inCorrectXYColumn);	si[16] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+2, j+1, k, &loc, &inCorrectXYColumn);	si[17] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+2, k, &loc, &inCorrectXYColumn);	si[18] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j+2, k, &loc, &inCorrectXYColumn);  si[19] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j-1, k+1, &loc, &inCorrectXYColumn);	 si[20] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j-1, k+1, &loc, &inCorrectXYColumn);	si[21] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j, k+1, &loc, &inCorrectXYColumn); 	si[22] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+2, j, k+1, &loc, &inCorrectXYColumn); 	si[23] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i-1, j+1, k+1, &loc, &inCorrectXYColumn);	si[24] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+2, j+1, k+1, &loc, &inCorrectXYColumn);	si[25] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+2, k+1, &loc, &inCorrectXYColumn);	 si[26] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j+2, k+1, &loc, &inCorrectXYColumn);  si[27] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j, k+2, &loc, &inCorrectXYColumn);	si[28] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j, k+2, &loc, &inCorrectXYColumn);	 si[29] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i, j+1, k+2, &loc, &inCorrectXYColumn);	 si[30] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
	valid = parent->getExistingLocator(i+1, j+1, k+2, &loc, &inCorrectXYColumn);	si[31] = parent->beginTubeIterator(loc, inCorrectXYColumn, valid, bufferId);
    }
}




template<class Traits>  template<TubeType tt>
void DTGrid<Traits>::StencilTubeIterator<tt>::gradient(Real g[3], UInt x1, UInt x2, UInt y1, UInt y2, UInt z1, UInt z2) const
{	
    g[0] = ( si[x2].getValue() - si[x1].getValue() ) * dxc2;	
    g[1] = ( si[y2].getValue() - si[y1].getValue() ) * dxc2;	
    g[2] = ( si[z2].getValue() - si[z1].getValue() ) * dxc2;	
}


template<class Traits>  template<TubeType tt>  template <StencilFormat iterStencil, int i>
void DTGrid<Traits>::StencilTubeIterator<tt>::gradientAtVoxel(Real g[3]) const
{
    if (iterStencil == SF_VOXEL_GRAD)
    {
	switch (i)
	{
	case 0:
	    gradient(g, 14, 1, 12, 2, 8, 4);
	    break;
	case 1:
	    gradient(g, 0, 15, 13, 3, 9, 5);
	    break;
	case 2:
	    gradient(g, 16, 3, 0, 18, 10, 6);
	    break;
	case 3:
	    gradient(g, 2, 17, 1, 19, 11, 7);
	    break;
	case 4:
	    gradient(g, 22, 5, 20, 6, 0, 28);
	    break;
	case 5:
	    gradient(g, 4, 23, 21, 7, 1, 29);
	    break;
	case 6:
	    gradient(g, 24, 7, 4, 26, 2, 30);
	    break;
	case 7:
	    gradient(g, 6, 25, 5, 27, 3, 31);
	    break;
	}
    }
}






template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getMaxMemUsed() const
{
    return maxMemUsed;
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getMaxValueMemUsed() const
{
    return maxValueMemUsed;
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getMaxIndexMemUsed() const
{
    return maxIndexMemUsed;
}


template<class Traits>
void DTGrid<Traits>::resetMaxMemUsed()
{
    maxMemUsed = getMemUsed();
    maxValueMemUsed = getValueMemUsed();
    maxIndexMemUsed = getIndexMemUsed();
}

template<class Traits>
void DTGrid<Traits>::updateMaxMemUsed(UInt additionalMemUsed)
{
    UInt currentMemUsed = getMemUsed()+additionalMemUsed;
    if (currentMemUsed > maxMemUsed)
    {
	maxMemUsed = currentMemUsed;
	maxValueMemUsed = getValueMemUsed();
	maxIndexMemUsed = getIndexMemUsed();
    }
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getMemUsed(UInt id) const
{
    return getValueMemUsed(id) + getIndexMemUsed(id);
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getValueMemUsed(UInt id) const
{
    return numValues[id] * sizeof(Data);
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getIndexMemUsed(UInt id) const
{
    return 
	(numXIndex+numIndexEndMarkers)*sizeof(Index) +
	(numYIndex+numIndexEndMarkers)*sizeof(Index) +
	(numZIndex+numIndexEndMarkers)*sizeof(Index) +
	(numVa1D+numValueEndMarkers)*sizeof(UInt) +
	(numVa2D+numValueEndMarkers)*sizeof(UInt) +
	(numXIndex/2+numAAEndMarkers)*sizeof(UInt) +
	(numYIndex/2+numAAEndMarkers)*sizeof(UInt) +
	(numZIndex/2+numAAEndMarkers)*sizeof(UInt);
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getMemUsed() const
{
    return getIndexMemUsed() + getValueMemUsed();
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getValueMemUsed() const
{
    UInt i;
    UInt memUsed = 0;
    for (i=0; i<maxNumBuffers; i++)
    {
	memUsed += numValues[i] * sizeof(Data);
    }

    return memUsed;
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::getIndexMemUsed() const
{
    return 
	(numXIndex+numIndexEndMarkers)*sizeof(Index) +
	(numYIndex+numIndexEndMarkers)*sizeof(Index) +
	(numZIndex+numIndexEndMarkers)*sizeof(Index) +
	(numVa1D+numValueEndMarkers)*sizeof(UInt) +
	(numVa2D+numValueEndMarkers)*sizeof(UInt) +
	(numXIndex/2)*sizeof(UInt) +
	(numYIndex/2)*sizeof(UInt) +
	(numZIndex/2)*sizeof(UInt);
}





///////////////////////////////////////////////////////////////////////////////////
// PROTECTED ITERATORS
///////////////////////////////////////////////////////////////////////////////////

// Initiation methods	


template<class Traits>
typename DTGrid<Traits>::_Iterator1D DTGrid<Traits>::begin1D()
{
    return _Iterator1D(this, true);
}

template<class Traits>
typename DTGrid<Traits>::_Iterator1D DTGrid<Traits>::end1D()
{
    return _Iterator1D(this, false);
}


template<class Traits>
typename DTGrid<Traits>::_Iterator2D DTGrid<Traits>::begin2D()
{
    return _Iterator2D(this, true);
}

template<class Traits>
typename DTGrid<Traits>::_Iterator2D DTGrid<Traits>::end2D()
{
    return _Iterator2D(this, false);
}


template<class Traits>
typename DTGrid<Traits>::_StencilIterator1D DTGrid<Traits>::beginStencil1D(Index width)
{
    return _StencilIterator1D(this, width);
}


template<class Traits>
typename DTGrid<Traits>::_StencilIterator2D DTGrid<Traits>::beginStencil2D(Index width)
{
    return _StencilIterator2D(this, width);
}


// _Iterator1D


template<class Traits>
DTGrid<Traits>::_Iterator1D::_Iterator1D(DTGrid *parent, bool begin)
{
    this->parent = parent;

    if (begin)
    {
	if (parent->numVa1D > 0)
	{
	    va1D = parent->va1D;
	    ic1D = 1;
	    iv1D = 0;
	    x = parent->xIndex[0];
	    value = va1D[iv1D];
	}
	else
	{
	    iv1D = 0;
	}
    }
    else
    {
	iv1D = parent->numVa1D;
    }
}

template<class Traits>
typename DTGrid<Traits>::_Iterator1D& DTGrid<Traits>::_Iterator1D::operator++()
{
    // we assume that iv1D < numVa1D when this method is entered

    iv1D++;

    // set value
    value = va1D[iv1D];

    // reached end of connected component in z direction?
    if (x == parent->xIndex[ic1D])
    {
	ic1D++;
	x = parent->xIndex[ic1D];
	ic1D++;
    }
    else
    {
	x++;
    }

    return *this;
}

template<class Traits>
bool DTGrid<Traits>::_Iterator1D::operator==(const _Iterator1D& iter) const
{
    return iv1D == iter.iv1D;
}

template<class Traits>
bool DTGrid<Traits>::_Iterator1D::operator!=(const _Iterator1D& iter) const
{
    return iv1D != iter.iv1D;
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::_Iterator1D::getValue() const
{
    return value;
}

template<class Traits>
void DTGrid<Traits>::_Iterator1D::setValue(UInt value)
{
    va1D[iv1D] = value;
}

template<class Traits>
typename DTGrid<Traits>::Index DTGrid<Traits>::_Iterator1D::getX() const
{
    return x;
}




// _Iterator2D


template<class Traits>
DTGrid<Traits>::_Iterator2D::_Iterator2D(DTGrid *parent, bool begin)
{
    this->parent = parent;

    if (begin)
    {
	if (parent->numVa2D > 0)
	{
	    va2D = parent->va2D;
	    ic1D = ic2D = 1;
	    iv1D = 1;
	    iv2D = 0;
	    x = parent->xIndex[0];
	    y = parent->yIndex[0];
	    value = va2D[iv2D];
	}
	else
	{
	    iv2D = 0;
	}
    }
    else
    {
	iv2D = parent->numVa2D;
    }
}

template<class Traits>
typename DTGrid<Traits>::_Iterator2D& DTGrid<Traits>::_Iterator2D::operator++()
{
    // we assume that iv2D < numVa2D when this method is entered

    iv2D++;

    // set value
    value = va2D[iv2D];

    // reached end of connected component in z direction?
    if (y == parent->yIndex[ic2D])
    {
	ic2D++;

	// moved to a new (x) column?
	if (parent->va1D[iv1D] == ic2D)
	{
	    iv1D++;

	    if (x == parent->xIndex[ic1D])
	    {
		ic1D++;
		x = parent->xIndex[ic1D];
		ic1D++;
	    }
	    else
	    {
		x++;
	    }
	}

	y = parent->yIndex[ic2D];
	ic2D++;
    }
    else
    {
	y++;
    }

    return *this;
}

template<class Traits>
bool DTGrid<Traits>::_Iterator2D::operator==(const _Iterator2D& iter) const
{
    return iv2D == iter.iv2D;
}

template<class Traits>
bool DTGrid<Traits>::_Iterator2D::operator!=(const _Iterator2D& iter) const
{
    return iv2D != iter.iv2D;
}

template<class Traits>
typename DTGrid<Traits>::UInt DTGrid<Traits>::_Iterator2D::getValue() const
{
    return value;
}

template<class Traits>
void DTGrid<Traits>::_Iterator2D::setValue(UInt value)
{
    va2D[iv2D] = value;
}

template<class Traits>
typename DTGrid<Traits>::Index DTGrid<Traits>::_Iterator2D::getX() const
{
    return x;
}

template<class Traits>
typename DTGrid<Traits>::Index DTGrid<Traits>::_Iterator2D::getY() const
{
    return y;
}



// Stencil Iterator 1D

template<class Traits>
DTGrid<Traits>::_StencilIterator1D::_StencilIterator1D(DTGrid *parent, Index width)
{
    this->width = width;
    this->sii = 0;
    this->eii = parent->numXIndex;
    this->xIndex = parent->xIndex;
    this->aa1D = parent->aa1D;
    this->va1D = parent->va1D;
    valid = true;

    reset();
}

template<class Traits>
DTGrid<Traits>::_StencilIterator1D::_StencilIterator1D(Index width, Index *xIndex, UInt *aa1D, UInt *va1D, UInt sii, UInt eii)
{
    this->width = width;
    this->xIndex = xIndex;
    this->aa1D = aa1D;
    this->va1D = va1D;
    this->sii = sii;
    this->eii = eii;
    valid = true;

    reset();
}

template<class Traits>
bool DTGrid<Traits>::_StencilIterator1D::isEndValid()
{
    return ic1De <= eii  &&  xe == e;
}

template<class Traits>
bool DTGrid<Traits>::_StencilIterator1D::isStartValid()
{
    return ic1Ds <= eii  &&  xs == s;
}


template<class Traits>
void DTGrid<Traits>::_StencilIterator1D::getEndValues(UInt *s, UInt *e)
{
    *s = va1D[iv1De];

    // increment end iterator
    iv1De++;
    if (xe == xIndex[ic1De])
    {
	ic1De++;
	xe = xIndex[ic1De];
	ic1De++;
    }
    else
    {
	xe++;
    }

    *e = va1D[iv1De];
}

template<class Traits>
void DTGrid<Traits>::_StencilIterator1D::getStartValues(UInt *s, UInt *e)
{
    *s = va1D[iv1Ds];

    // increment start iterator
    iv1Ds++;
    if (xs == xIndex[ic1Ds])
    {
	ic1Ds++;
	xs = xIndex[ic1Ds];
	ic1Ds++;
    }
    else
    {
	xs++;
    }

    *e = va1D[iv1Ds];
}

template<class Traits>
void DTGrid<Traits>::_StencilIterator1D::increment(Index xi)
{
    s += xi;
    e += xi;
}

template<class Traits>
typename DTGrid<Traits>::Index DTGrid<Traits>::_StencilIterator1D::getStartIndex() const
{
    return s;
}

template<class Traits>
typename DTGrid<Traits>::Index DTGrid<Traits>::_StencilIterator1D::getEndIndex() const
{
    return e;
}


template<class Traits>
void DTGrid<Traits>::_StencilIterator1D::setStartAndEndIndex(Index s, Index e)
{ 
    this->s=s; 
    this->e=e; 
}

template<class Traits>
void DTGrid<Traits>::_StencilIterator1D::reset()
{
    if (valid)
    {
	// initialize start iterator
	iv1Ds = aa1D[(sii)>>1];
	ic1Ds = sii+1;
	xs = xIndex[sii];

	// initialize end iterator
	iv1De = aa1D[(sii)>>1];
	ic1De = sii+1;
	xe = xIndex[sii];

	s = xs-2*width;
	e = xe;
    }
}

template<class Traits>
void DTGrid<Traits>::_StencilIterator1D::invalidate()
{
    valid = false;
}

template<class Traits>
bool DTGrid<Traits>::_StencilIterator1D::isValid()
{
    return valid;
}

template<class Traits>
typename DTGrid<Traits>::Index DTGrid<Traits>::_StencilIterator1D::getIncrement() const
{
    Index xsDist;
    Index xeDist;

    if (xs-s == 0)
    {
	if (xs == xIndex[ic1Ds])
	{
	    if (xIndex[ic1Ds]!=xIndex[ic1Ds+1])
		xsDist = xIndex[ic1Ds+1]-s;
	    else
		xsDist = xIndex[ic1Ds+2]-s;
	}
	else
	{
	    xsDist = 1;
	}
    }
    else
    {
	if (xs<s)
	{
	    xsDist = std::numeric_limits<Index>::max();
	}
	else
	{
	    if (e<=xs)
	    {
		xsDist = 1;
	    }
	    else
	    {
		xsDist = xs-s; 
	    }
	}
    }

    if (ic1De <= eii)
    {
	if (xe-e == 0)
	{
	    if (xe == xIndex[ic1De])
	    {
		if (xIndex[ic1De]!=xIndex[ic1De+1])
		    xeDist = xIndex[ic1De+1]-e;
		else
		    xeDist = xIndex[ic1De+2]-e;
	    }
	    else
	    {
		xeDist = 1;
	    }
	}
	else
	{
	    if (xe<e)
	    {
		xeDist = std::numeric_limits<Index>::max();
	    }
	    else
	    {
		xeDist = xe-e; 
	    }
	}
    }
    else
    {
	xeDist = std::numeric_limits<Index>::max();
    }

    Index res = std::min(xeDist, xsDist);

    return res;
}

// Stencil Iterator 2D

template<class Traits>
DTGrid<Traits>::_StencilIterator2D::_StencilIterator2D(DTGrid *parent, Index width)
{
    UInt si, ei;

    this->parent = parent;
    this->width = width;

    si1D = parent->beginStencil1D(width);
    lsi = 0;

    // assume si1D is valid
    si1D.getEndValues(&si, &ei);
    si2D[0] = _StencilIterator1D(width, parent->yIndex, parent->aa2D, parent->va2D, si, ei);
    s = si2D[0].getStartIndex();
    e = si2D[0].getEndIndex();

}


template<class Traits>
DTGrid<Traits>::_StencilIterator2D::~_StencilIterator2D()
{
}


template<class Traits>
bool DTGrid<Traits>::_StencilIterator2D::isEndValid(UInt k)
{
    return si2D[k].isValid() && si2D[k].isEndValid();
}

template<class Traits>
bool DTGrid<Traits>::_StencilIterator2D::isStartValid(UInt k)
{
    return si2D[k].isValid() && si2D[k].isStartValid();
}

template<class Traits>
void DTGrid<Traits>::_StencilIterator2D::getEndValues(UInt *s, UInt *e, UInt k)
{
    si2D[k].getEndValues(s, e);
}

template<class Traits>
void DTGrid<Traits>::_StencilIterator2D::getStartValues(UInt *s, UInt *e, UInt k)
{
    si2D[k].getStartValues(s, e);
}

template<class Traits>
void DTGrid<Traits>::_StencilIterator2D::increment(Index xi, Index yi)
{
    if (xi == 0)
    {
	// only increment in the y direction

	Index k;

	e += yi;
	s += yi;

	for (k=0; k<(2*width+1); k++)
	{
	    si2D[k].increment(yi);
	}
    }
    else
    {
	// increment in x direction	
	// and initialize in y direction

	Index minY;
	UInt si, ei;
	Index k, m;

	si1D.increment(xi);	    

	lsi = (lsi + 1) % (2*width+1);

	if (si1D.isEndValid())
	{
	    si1D.getEndValues(&si, &ei);
	    si2D[lsi] = _StencilIterator1D(width, parent->yIndex, parent->aa2D, parent->va2D, si, ei);
	    minY = si2D[lsi].getStartIndex();
	}
	else
	{
	    // the _StencilIterator1D at this "stencil-column" should be invalidated
	    si2D[lsi].invalidate();
	    minY = std::numeric_limits<Index>::max();
	}

	for (k = (lsi+1)%(2*width+1), m=0; m<(2*width); m++, k=(k+1)%(2*width+1))
	{
	    si2D[k].reset();
	    if (si2D[k].isValid() && si2D[k].getStartIndex() < minY)
	    {
		minY = si2D[k].getStartIndex();
	    }
	}

	s = minY;
	e = minY+2*width;

	for (k=0; k<(2*width+1); k++)
	{
	    si2D[k].setStartAndEndIndex(s, e);
	}

    }

}



template<class Traits>
void DTGrid<Traits>::_StencilIterator2D::increment()
{
    // first check if we should increment in y-direction
    Index xi, yi;

    yi = std::numeric_limits<Index>::max();

    for (Index k=0; k<2*width+1; k++)
    {
	if (si2D[k].isValid() && !si2D[k].isAtEnd())
	{
	    yi = std::min(static_cast<Index>(si2D[k].getIncrement()), yi);
	}
    }


    if (yi == std::numeric_limits<Index>::max())
    {
	// we should increment only in the x-direction
	xi = si1D.getIncrement();
    }
    else
    {
	xi = 0;
    }




    if (xi == 0)
    {
	// only increment in the y direction

	Index k;

	e += yi;
	s += yi;

	for (k=0; k<(2*width+1); k++)
	{
	    si2D[k].increment(yi);
	}
    }
    else
    {
	// increment in x direction	
	// and initialize in y direction

	Index minY;
	UInt si, ei;
	//UInt k, m;
	Index k, m;

	si1D.increment(xi);	    

	lsi = (lsi + 1) % (2*width+1);

	if (si1D.isEndValid())
	{
	    si1D.getEndValues(&si, &ei);
	    si2D[lsi] = _StencilIterator1D(width, parent->yIndex, parent->aa2D, parent->va2D, si, ei);
	    minY = si2D[lsi].getStartIndex();
	}
	else
	{
	    // the _StencilIterator1D at this "stencil-column" should be invalidated
	    si2D[lsi].invalidate();
	    minY = std::numeric_limits<Index>::max();
	}

	if (si1D.isStartValid())
	{
	    si1D.getStartValues(&si, &ei);
	}

	for (k = (lsi+1)%(2*width+1), m=0; m<(2*width); m++, k=(k+1)%(2*width+1))
	{
	    si2D[k].reset();
	    if (si2D[k].isValid() && si2D[k].getStartIndex() < minY)
	    {
		minY = si2D[k].getStartIndex();
	    }
	}

	s = minY;
	e = minY+2*width;

	for (k=0; k<(2*width+1); k++)
	{
	    si2D[k].setStartAndEndIndex(s, e);
	}

    }

}


///////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF TUBETOPOLOGYITERATOR
///////////////////////////////////////////////////////////////////////////////////

template<class Traits>
DTGrid<Traits>::TubeTopologyIterator::TubeTopologyIterator(const DTGrid *parent, bool begin)
: parent(parent), zIndex(parent==NULL?NULL:parent->zIndex), numVa3D(parent==NULL?0:parent->numVa3D)
{
    if (begin)
    {
	// we have to skip this test for the topology iterator because it may be the case
	// that no values have been pushed yet, so we check on topology instead
	if (numVa3D > 0)  
	{
	    ic1D = ic2D = ic3D = 1;
	    iv1D = iv2D = 1;
	    iv3D = 0;
	    x = parent->xIndex[0];
	    y = parent->yIndex[0];
	    z = zIndex[0];
	    valid = true;
	    inCorrectPColumn = true;
	    zn = zIndex[ic3D];
	}
	else
	{
	    iv3D = 0;
	    valid = false;
	}
    }
    else
    {
	iv3D = numVa3D;
	valid = false;
    }
}


template<class Traits>
bool DTGrid<Traits>::TubeTopologyIterator::operator==(const TubeTopologyIterator& iter) const
{
    return iter.iv3D == iv3D;
}

template<class Traits>
bool DTGrid<Traits>::TubeTopologyIterator::operator!=(const TubeTopologyIterator& iter) const
{
    return iter.iv3D != iv3D;
}

template<class Traits>
typename DTGrid<Traits>::TubeTopologyIterator &DTGrid<Traits>::TubeTopologyIterator::operator++()
{
    // we assume that iv3D < numVa3D when this method is entered

    iv3D++;

    // reached end of connected component in z direction?
    if (z == zIndex[ic3D])
    {
	ic3D++;

	// moved to a new (x,y) column?
	if (parent->va2D[iv2D] == ic3D)
	{
	    iv2D++;

	    // reached end of connected component in y direction?
	    if (y == parent->yIndex[ic2D])                         
	    {
		ic2D++;

		// moved to a new (x) column?
		if (parent->va1D[iv1D] == ic2D)
		{
		    iv1D++;

		    if (x == parent->xIndex[ic1D])
		    {
			ic1D++;
			x = parent->xIndex[ic1D];
			ic1D++;
		    }
		    else
		    {
			x++;
		    }
		}

		y = parent->yIndex[ic2D];
		ic2D++;

	    }
	    else
	    {
		y++;
	    }

	}

	z = zIndex[ic3D];
	ic3D++;
    }
    else
    {
	z++;
    }
    return *this;
}



template<class Traits>
void DTGrid<Traits>::TubeTopologyIterator::incrementFast()
{
    iv3D++;
    z++;
}

template<class Traits>
void DTGrid<Traits>::TubeTopologyIterator::incrementFastZ(Index zn)
{
    // incrementFastZ must only increment the iterator one point forward (therefore, if "z == zIndex[ic3D]" as below and we are exiting a connected componen, 
    // the value should be set to +/- gamma since the iterator will under no circumstances end up in the correct grid point
    // note that it must be possible to continue iteration with incrementUntil following a sequence of incrementFastZ operations

    // first check if this iterator is valid in this (x,y)'th p-column of grid points
    // we have to make this check since a point in the stencil may pass beyond its last point before
    // the center stencil point does! (note that incrementFastZ is only called when the center stencil grid point moves one grid point forward, but
    // we are on the border of the narrow band.

    if ( inCorrectPColumn )   // valid, ie in correct p-column?
    {
	// yes, in correct p-column

	if (z < zn)   
	{
	    // yes, we need to increment

	    if (z == zIndex[ic3D])
	    {

		// we jump to next connected component, ie
		// this stencil iterator jumps more than one grid point, ie.
		// the corresponding grid point lies outside of the narrow band

		valid = false;

		ic3D++;  // next connected component
		z = zIndex[ic3D];

		// we have to make sure that the iterator stays in the same (x,y)'th p-column
		if (parent->va2D[iv2D] == ic3D)
		{
		    // we move to next p-column
		    // adjust 1D and 2D pointers

		    iv2D++;

		    // reached end of connected component in y direction?
		    if (y == parent->yIndex[ic2D])
		    {
			ic2D++;

			// moved to a new (x) column?
			if (parent->va1D[iv1D] == ic2D)
			{
			    iv1D++;

			    if (x == parent->xIndex[ic1D])
			    {
				ic1D++;
				x = parent->xIndex[ic1D];
				ic1D++;
			    }
			    else
			    {
				x++;  
			    }
			}

			y = parent->yIndex[ic2D];
			ic2D++;

		    }
		    else
		    {
			y++;
		    }

		    ic3D++;
		    iv3D++;
		    inCorrectPColumn = false;
		}
		else
		{
		    ic3D++;
		    iv3D++;
		}
	    }
	    else
	    {
		// we stay within the same connected component
		z++;
		iv3D++;
		// no need to set valid to true here, it will already be true
	    }
	}
	else if (z==zn)
	{
	    valid = true;
	}
    }
}


template<class Traits>
void DTGrid<Traits>::TubeTopologyIterator::incrementUntil(Index xn, Index yn, Index zn)
{
    if (x < xn)
    {
	do
	{
	    iv1D++;

	    if (x == parent->xIndex[ic1D])
	    {
		ic1D++;
		x = parent->xIndex[ic1D];
		ic1D++;
	    }
	    else
	    {
		x++;
	    }		
	}
	while (x < xn);

	ic2D = parent->va1D[iv1D-1];  // we know iv1D >= 1
	y = parent->yIndex[ic2D];
	iv2D = parent->aa2D[(ic2D)>>1];	    
	ic2D++;

	ic3D = parent->va2D[iv2D];
	iv2D++;
	z = zIndex[ic3D];
	iv3D = parent->aa3D[(ic3D)>>1];
	ic3D++;
    }
    if (xn==x)
    {

	if (y < yn)  
	{
	    do
	    {
		iv2D++;

		if (y == parent->yIndex[ic2D])
		{
		    ic2D++;

		    // moved to a new (x) column?
		    if (parent->va1D[iv1D] == ic2D)
		    {
			iv1D++;

			if (x == parent->xIndex[ic1D])
			{
			    ic1D++;
			    x = parent->xIndex[ic1D];
			    ic1D++;
			}
			else
			{
			    x++;
			}
		    }

		    y = parent->yIndex[ic2D];
		    ic2D++;
		}
		else
		{
		    y++;
		}		
	    }
	    while (x==xn && y<yn);

	    ic3D = parent->va2D[iv2D-1];
	    z = zIndex[ic3D];
	    iv3D = parent->aa3D[(ic3D)>>1];
	    ic3D++;
	}

	if (xn==x && yn==y)
	{

	    while (xn==x && yn==y && z < zn)  
	    {
		iv3D++;
		if (z == zIndex[ic3D])
		{
		    ic3D++;

		    // moved to a new (x,y) column?
		    if (parent->va2D[iv2D] == ic3D)
		    {
			iv2D++;

			// reached end of connected component in y direction?
			if (y == parent->yIndex[ic2D])
			{
			    ic2D++;

			    // moved to a new (x) column?
			    if (parent->va1D[iv1D] == ic2D)
			    {
				iv1D++;

				if (x == parent->xIndex[ic1D])
				{
				    ic1D++;
				    x = parent->xIndex[ic1D];
				    ic1D++;
				}
				else
				{
				    x++;
				}
			    }

			    y = parent->yIndex[ic2D];
			    ic2D++;

			}
			else
			{
			    y++;
			}

		    }

		    z = zIndex[ic3D];
		    ic3D++;
		}
		else
		{
		    z++;  
		}
	    }

	    if (xn==x && yn==y)
	    {
		inCorrectPColumn = true;

		if (zn==z)
		{
		    valid = true;
		}
		else
		{
		    valid = false;
		}
	    }
	    else
	    {
		inCorrectPColumn = false;
		valid = false;
	    }
	}
	else
	{
	    inCorrectPColumn = false;
	    valid = false;
	}
    }
    else
    {
	inCorrectPColumn = false;
	valid = false;
    }
}



template<class Traits>
bool DTGrid<Traits>::TubeTopologyIterator::hasNext() const
{
    return iv3D < parent->numVa3D;
}




template<class Traits>
void DTGrid<Traits>::TubeTopologyIterator::getIndex(Index *x, Index *y, Index *z) const
{
    *x = this->x;
    *y = this->y;
    *z = this->z;
}


}


#undef FIND_INDEX
#undef FIND_INDEX_EXISTS
#undef INCREMENT_FAST_Z
#undef INCREMENT_UNTIL
