/*
This file is part of ESFC (Eulerian Solid-FLuid Coupling).

ESFC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESFC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ESFC.  If not, see <http://www.gnu.org/licenses/>.
*/
// SETTINGS.h: Project-wide options set in one place
//
//////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H
#define SETTINGS_H

#include <cassert>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#define USING_OSX __APPLE__

#define USING_GLVU 1

// select single or double precision
#define DOUBLE_PRECISION 1

#define MATERIAL_DEBUG 0

typedef Eigen::Vector2i VEC2I;
typedef Eigen::Vector3i VEC3I;

#define USING_OPENMP 0
#define HYBRID_DIV_CONSTRAINT 0

#if DOUBLE_PRECISION
typedef Eigen::Vector2d VEC2F;
typedef Eigen::Vector3d VEC3F;
typedef Eigen::Vector4d VEC4F;
typedef Eigen::VectorXd VECTOR;
typedef Eigen::MatrixXd MATRIX;
typedef Eigen::Matrix2d MATRIX2;
typedef Eigen::Matrix3d MATRIX3;
typedef double Real;
#else
typedef Eigen::Vector2f VEC2F;
typedef Eigen::Vector3f VEC3F;
typedef Eigen::Vector4f VEC4F;
typedef Eigen::VectorXf VECTOR;
typedef Eigen::MatrixXf MATRIX;
typedef Eigen::Matrix2f MATRIX2;
typedef Eigen::Matrix3f MATRIX3;
typedef float Real;

#endif
typedef Eigen::Triplet<Real> TRIPLET;
typedef Eigen::SparseMatrix<Real, Eigen::RowMajor> SpMat;

#endif
