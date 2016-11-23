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


#include <stdlib.h>
#include <assert.h>

#include "bvh.h"
#include "VEC3.h"
#include "model.h"
#include "aap.h"
// extern int NoDivTriTriIsect(VEC3F& V0,VEC3F& V1,VEC3F& V2, VEC3F& U0,VEC3F& U1,VEC3F& U2);
extern int tri_contact (double *P1, double *P2, double *P3, double *Q1, double *Q2, double *Q3);

static Model *s_mdl1, *s_mdl2;
static bool s_self = false;
#define IsLeaf(node) (node->_left == NULL)

#include <vector>
static std::vector<unsigned int> s_rets;

void
bvh_tree::collide(bvh_tree *other)
{
	s_mdl1 = this->_mdl;
	s_mdl2 = other->_mdl;
	s_self = false;

	s_rets.erase(s_rets.begin(), s_rets.end());
	_root->collide(other->_root);
}

void
bvh_tree::self_collide()
{
	s_mdl1 = this->_mdl;
	s_mdl2 = this->_mdl;
	s_self = true;

	s_rets.erase(s_rets.begin(), s_rets.end());
	_root->self_collide();
}

#define set_vtx_color(mdl, id) {\
			mdl->_colors[id*3] = 0;\
			mdl->_colors[id*3+1] = 255;\
			mdl->_colors[id*3+2] = 0;\
}

bool
bvh_tree::get_contact(int i, unsigned int &id1, unsigned int &id2)
{
	id1 = -1;
	id2 = -1;

	int num = s_rets.size()/2;
	
	if (i<0 || i>=num)
		return false;

	id1 = s_rets[i*2];
	id2 = s_rets[i*2+1];
	return true;

}

void
bvh_tree::color()
{
	int idx=0;
	unsigned int vtx_id;

	for (std::vector<unsigned int>::iterator it=s_rets.begin(); it!=s_rets.end(); it++) {
		unsigned int tri_id = *it;
		if (idx%2) {// odd
			vtx_id = s_mdl2->_tris[tri_id*3];
			set_vtx_color(s_mdl2, vtx_id);
			vtx_id = s_mdl2->_tris[tri_id*3+1];
			set_vtx_color(s_mdl2, vtx_id);
			vtx_id = s_mdl2->_tris[tri_id*3+2];
			set_vtx_color(s_mdl2, vtx_id);
		} else { //even
			vtx_id = s_mdl1->_tris[tri_id*3];
			set_vtx_color(s_mdl1, vtx_id);
			vtx_id = s_mdl1->_tris[tri_id*3+1];
			set_vtx_color(s_mdl1, vtx_id);
			vtx_id = s_mdl1->_tris[tri_id*3+2];
			set_vtx_color(s_mdl1, vtx_id);
		}
		idx++;
	}
}

int
bvh_node::test_triangle(unsigned int id1, unsigned int id2)
{
	VEC3F &p1 = s_mdl1->_vtxs[s_mdl1->_tris[id1*3]];
	VEC3F &p2 = s_mdl1->_vtxs[s_mdl1->_tris[id1*3+1]];
	VEC3F &p3 = s_mdl1->_vtxs[s_mdl1->_tris[id1*3+2]];

	VEC3F &q1 = s_mdl2->_vtxs[s_mdl2->_tris[id2*3]];
	VEC3F &q2 = s_mdl2->_vtxs[s_mdl2->_tris[id2*3+1]];
	VEC3F &q3 = s_mdl2->_vtxs[s_mdl2->_tris[id2*3+2]];

	return tri_contact((double *)&p1, (double *)&p2, (double *)&p3, (double *)&q1, (double *)&q2, (double *)&q3);
}

inline VEC3F interp(VEC3F &p1, VEC3F &p2, double t)
{
	return p1*(1-t)+p2*t;
}

static VEC3F P[3], Q[3];
static float eps_interp = float(10e-5);

inline int tri_contact (VEC3F *p[3], VEC3F *q[3], unsigned int cov)
{
	if (cov == 3)
		return true;

	if (cov == 0) {
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				P[i][j] = (*p[i])[j];
				Q[i][j] = (*q[i])[j];
			}

		return tri_contact(P[0], P[1], P[2], Q[0], Q[1], Q[2]);
	}

	if (cov == 1) {
		for (int i=1; i<3; i++)
			for (int j=0; j<3; j++) {
				P[i][j] = (*p[i])[j];
				Q[i][j] = (*q[i])[j];
			}

		VEC3F mid = interp(*p[1], *p[2], 0.5);
		VEC3F p0 = interp(*p[0], mid, eps_interp);

		mid = interp(*q[1], *q[2], 0.5);
		VEC3F q0 = interp(*q[0], mid, eps_interp);

		for (int j=0; j<3; j++) {
			P[0][j] = p0[j];
			Q[0][j] = q0[j];
		}

		return tri_contact(P[0], P[1], P[2], Q[0], Q[1], Q[2]);
	}

	if (cov == 2) {
		for (int j=0; j<3; j++) {
			P[0][j] = (*p[0])[j];
			Q[0][j] = (*q[0])[j];
		}

		VEC3F p1 = interp(*p[1], *p[0], eps_interp);
		VEC3F p2 = interp(*p[2], *p[0], eps_interp);
		for (int j=0; j<3; j++) {
			P[1][j] = p1[j];
			P[2][j] = p2[j];
		}

		VEC3F q1 = interp(*q[1], *q[0], eps_interp);
		VEC3F q2 = interp(*q[2], *q[0], eps_interp);
		for (int j=0; j<3; j++) {
			Q[1][j] = q1[j];
			Q[2][j] = q2[j];
		}

		return tri_contact(P[0], P[1], P[2], Q[0], Q[1], Q[2]);
	}

	assert(0);
	return false;
}

int
bvh_node::test_triangle(unsigned int id1, unsigned int id2,
						unsigned int cov, unsigned int st1, unsigned int st2)
{
	VEC3F *p[3];
	for (int i=0; i<3; i++)
		p[i] = s_mdl1->_vtxs+(s_mdl1->_tris[id1*3+(st1+i)%3]);

	VEC3F *q[3];
	for (int i=0; i<3; i++)
		q[i] = s_mdl2->_vtxs+(s_mdl2->_tris[id2*3+(st2+i)%3]);

	return tri_contact(p, q, cov);
}

void
bvh_node::self_collide()
{
	if (IsLeaf(this))
		return;

	_left->self_collide();
	_right->self_collide();
	_left->collide(_right);
}

void
bvh_node::collide(bvh_node *other)
{
	assert(other);

	s_mdl1->_num_box_tests++;
	if (!_box.overlaps(other->_box)) {
		return;
	}

	if (IsLeaf(this) && IsLeaf(other)){
		unsigned int st1 = 0, st2 = 0;
		unsigned int cov = 0;
		
		if (s_self)
			cov = s_mdl1->Covertex(_id, other->_id, st1, st2);

		s_mdl1->_num_tri_tests++;
		if (test_triangle(_id, other->_id, cov, st1, st2)) {
			s_mdl1->_num_contacts++;
			s_rets.push_back(_id);
			s_rets.push_back(other->_id);
		}
		return;
	}

	if (IsLeaf(this)) {
		collide(other->_left);
		collide(other->_right);
	} else {
		_left->collide(other);
		_right->collide(other);
	}
}
