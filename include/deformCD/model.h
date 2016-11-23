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

#pragma once
// class VEC3F;
// class VEC3F;
#include "VEC3.h"

class bvh_tree;
class bvh_node;
class aabb;

class Model {
 protected:
  unsigned int _num_vtx;
  VEC3F *_vtxs;
  VEC3F *_nrms;
  unsigned char *_colors;

  int *_vIDs;

  unsigned int _num_tri;
  unsigned int *_tris;

  // for building BVH
  bvh_tree *_tree;
  VEC3F *_tri_centers;
  aabb *_tri_boxes;
  VEC3F *_tri_nrms;

  // for collide
  unsigned int _num_box_tests;
  unsigned int _num_tri_tests;
  unsigned int _num_contacts;
  unsigned int *_contacts;

  void ResetResult();

 public:
  Model();
  virtual ~Model();

  void UpdateNorm();
  void DisplayColor();

  void BuildBVH();
  void DeleteBVH();
  int RefitBVH();
  void DisplayBVH(int);

  void Collide(Model *);
  void SelfCollide();

  void ColorCollide();
  void ResetColor(unsigned char r, unsigned char g, unsigned char b);

  unsigned int Covertex(unsigned int id1, unsigned int id2, unsigned int &st1,
                        unsigned int &st2);

  int NumBoxTest() { return _num_box_tests; }
  int NumTriTest() { return _num_tri_tests; }
  int NumContact() { return _num_contacts; }
  bool GetContact(int i, unsigned int &id1, unsigned int &id2);

  friend class bvh_tree;
  friend class bvh_node;
};
