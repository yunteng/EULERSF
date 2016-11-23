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


// #include <windows.h>
// #include <gl/gl.h>
//#define USING_OSX 1
#include "SETTINGS.h"
#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "model.h"
// #include "vec4d.h"
#include "bvh.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Model::Model()
{
	_vtxs = NULL;
	_nrms = NULL;
	_colors = NULL;
	_vIDs = NULL;

	_num_tri = 0;
	_tris = NULL;

	_tree = NULL;
	_tri_boxes = NULL;
	_tri_centers = NULL;
	_tri_nrms = NULL;

	_num_box_tests = 0;
	_num_tri_tests = 0;
	_num_contacts = 0;
	_contacts = NULL;
}

Model::~Model()
{
	if(_vtxs != NULL){

		delete[] _vtxs;
		delete[] _vIDs;
		delete[] _colors;
		delete[] _nrms;
		delete[] _tri_nrms;
		delete[] _tris;
		if(_tri_boxes)
			delete[] _tri_boxes;
		if(_tri_centers)
			delete[] _tri_centers;
		if(_contacts)
			delete[] _contacts;

		delete _tree;
	}
}

void
Model::DisplayColor()
{
	glShadeModel(GL_SMOOTH);
	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_NORMAL_ARRAY );

	if (_colors) {
		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		glEnableClientState( GL_COLOR_ARRAY );

		glColorPointer(3, GL_UNSIGNED_BYTE, 0, _colors);
	}

	glVertexPointer( 3, GL_DOUBLE, 0, _vtxs);
	glNormalPointer(GL_DOUBLE, 0, _nrms);

	glDrawElements( GL_TRIANGLES, _num_tri*3, GL_UNSIGNED_INT, _tris);

	// glDisableClientState( GL_VERTEX_ARRAY );
	// glDisableClientState( GL_NORMAL_ARRAY );
	// glDisableClientState( GL_COLOR_ARRAY );
	
	// glDisable(GL_COLOR_MATERIAL);
}

void Model::DisplayBVH(int level)
{
	glDisable(GL_LIGHTING);
	if (_tree)
		_tree->visulization(level);
	glEnable(GL_LIGHTING);
}

void Model::BuildBVH()
{
	_tree = new bvh_tree(this);
	RefitBVH();
}

void Model::DeleteBVH()
{
	delete _tree;
	_tree = NULL;
}

void Model::ResetResult()
{
	// reset results
	_num_box_tests = 0;
	_num_tri_tests = 0;
	_num_contacts = 0;
	if (_contacts) delete [] _contacts;
	_contacts = NULL;
}

void Model::SelfCollide()
{
	ResetResult();

	if (_tree == NULL)
		return;

	_tree->self_collide();
}

void Model::Collide(Model *mdl)
{
	ResetResult();

	if (_tree == NULL || mdl->_tree == NULL)
		return;

	_tree->collide(mdl->_tree);
}

void Model::ResetColor(unsigned char r, unsigned char g, unsigned char b)
{
	for (unsigned int i=0; i<_num_vtx; i++)
	{
		_colors[i*3] = r;
		_colors[i*3+1] = g;
		_colors[i*3+2] = b;
	}

}

bool Model::GetContact(int i, unsigned int &id1, unsigned int &id2)
{
	return _tree->get_contact(i, id1, id2);
}

void Model::ColorCollide()
{
	_tree->color();
}

int Model::RefitBVH()
{
	if(_num_vtx > 0)
		return _tree->refit();
	return 0;
}


inline VEC3F update(VEC3F &v1, VEC3F &v2, VEC3F &v3)
{
	VEC3F s = (v2-v1);
	return cross(s, v3-v1);
}

void
Model::UpdateNorm()
{
	for (unsigned int i=0; i<_num_vtx; i++)
		_nrms[i] = VEC3F(0, 0, 0);
	for (unsigned int i=0; i<_num_tri; i++) {
		unsigned int id0 = _tris[i*3];
		unsigned int id1 = _tris[i*3+1];
		unsigned int id2 = _tris[i*3+2];

		VEC3F v = update(_vtxs[id0], _vtxs[id1], _vtxs[id2]);
		_tri_nrms[i] = v;//VEC3F(v.x, v.y, v.z);

		_nrms[id0] += _tri_nrms[i];
		_nrms[id1] += _tri_nrms[i];
		_nrms[id2] += _tri_nrms[i];
	}

	for (unsigned int i=0; i<_num_vtx; i++)
		_nrms[i].normalize();
}

unsigned int
Model::Covertex(unsigned int id1, unsigned int id2, unsigned int &st1, unsigned int &st2)
{
	unsigned int keeps[4];
	unsigned int num = 0;
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++) {
			if (_tris[id1*3+i] == _tris[id2*3+j]) {
				if (num < 2) {
					keeps[num*2] = i;
					keeps[num*2+1] = j;
				}
				num++;
			}
		}

	assert(num <= 3);
	if (num == 1) {
		st1 = keeps[0];
		st2 = keeps[1];
	} else
	if (num == 2) {
		for (int i=0; i<3; i++) {
			if (i != keeps[0] && i != keeps[2]) {
				st1 = i;
			}
			if (i != keeps[1] && i != keeps[3]) {
				st2 = i;
			}
		}
	} else {
		st1 = st2 = 0;
	}

	return num;
}
