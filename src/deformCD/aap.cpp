/*************************************************************************\

  Copyright 2007 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
   fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             GAMMA Research Group at UNC
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:              {geom, tangm}@cs.unc.edu


\**************************************************************************/


#include <deformCD/aap.h>
// #include "vec4d.h"
#include <assert.h>

#define MAX(a,b)	((a) > (b) ? (a) : (b))
#define MIN(a,b)	((a) < (b) ? (a) : (b))

aap::aap(char xyz, float p)
{
	assert(xyz >= 0 && xyz <= 2);
	this->xyz = xyz;
	this->p = p;
}

aap::aap(const aabb &total)
{
	VEC3F center = total.center();
	char xyz = 2;

	if (total.width() >= total.height() && total.width() >= total.depth()) {
		xyz = 0;
	} else
	if (total.height() >= total.width() && total.height() >= total.depth()) {
		xyz = 1;
	}

	this->xyz = xyz;
	this->p = center[xyz];
}

// float middle_xyz(char xyz, const VEC3F &p1, const VEC3F &p2, const VEC3F &p3)
// {
// 	float t0, t1;

// 	t0 = MIN(p1[xyz], p2[xyz]);
// 	t0 = MIN(t0, p3[xyz]);
// 	t1 = MAX(p1[xyz], p2[xyz]);
// 	t1 = MAX(t1, p3[xyz]);
// 	return (t0+t1)*0.5;
// }

float middle_xyz(char xyz, const VEC3F &p1, const VEC3F &p2, const VEC3F &p3)
{
	float t0, t1;

	t0 = MIN(p1[xyz], p2[xyz]);
	t0 = MIN(t0, p3[xyz]);
	t1 = MAX(p1[xyz], p2[xyz]);
	t1 = MAX(t1, p3[xyz]);
	return (t0+t1)*0.5;
}

bool
aap::inside(const VEC3F &mid) const
{
	return mid[xyz]>p;
}

bool
aap::inside(const VEC3F &p1, const VEC3F &p2, const VEC3F &p3) const
{
	return middle_xyz(xyz, p1, p2, p3)>p;
}

