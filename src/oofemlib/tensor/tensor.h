/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "tensor/FTensor.hpp"

#pragma once

using namespace FTensor;

static FTensor::Index<'i', 3> i_3;
static FTensor::Index<'j', 3> j_3;
static FTensor::Index<'k', 3> k_3;
static FTensor::Index<'l', 3> l_3;
static FTensor::Index<'m', 3> m_3;
static FTensor::Index<'n', 3> n_3;
/*static FTensor::Index<'o', 3> o_3;
static FTensor::Index<'p', 3> p_3;
static FTensor::Index<'q', 3> q_3;
static FTensor::Index<'r', 3> r_3;
static FTensor::Index<'s', 3> s_3;
static FTensor::Index<'t', 3> t_3;
static FTensor::Index<'u', 3> u_3;
static FTensor::Index<'v', 3> v_3;
static FTensor::Index<'w', 3> w_3;
static FTensor::Index<'x', 3> x_3;
static FTensor::Index<'y', 3> y_3;
static FTensor::Index<'z', 3> z_3;
*/


/**
 * Functions providing basic operations with tensors 
 * based on FTensor library
 * It provides transformation from the second and fourth order tensors in 3d to FloatArrayF<9> and FloatArray<9,9> 
 * It provides calculation of determinant, cofactor, inverse, and tensor cross product
 * More to be added gradually
 * @author Martin Horak
 **/

namespace oofem {


  class Tensor2sym_3d : public Tensor2_symmetric<double, 3>
  {
  public:
    using Tensor2_symmetric<double, 3> :: Tensor2_symmetric;
    Tensor2sym_3d(const oofem::FloatArrayF<6> &array){
      this->data[0] = array.at(1);
      this->data[1] = array.at(6);
      this->data[2] = array.at(5);
      this->data[3] = array.at(6);
      this->data[4] = array.at(2);
      this->data[5] = array.at(4);
      this->data[6] = array.at(5);
      this->data[7] = array.at(4);
      this->data[8] = array.at(3);
    }

    const inline FloatArrayF<6> to_voigt_form() 
    {
      return {
	this->operator()(0, 0),
	  this->operator()(1, 1),
	  this->operator()(2, 2),
	  this->operator()(1, 2),
	  this->operator()(0, 2),
	  this->operator()(0, 1)
	  };
    }

    /*
     * Solves the eigenvalues and eigenvectors of real
     * symmetric matrix by jacobi method.
     *  Written by bp. Inspired by ED WILSON jaco_ procedure.
     *
     * Parameters (input):
     * nf - number of significant figures
     *
     * Output params:
     * eval - eigen values (not sorted)
     * v    - eigenvectors (stored columvise)
     */

    /*  
    std::pair<std::vector<double>, std::vector<Tensor1_3d> > eigs(Tensor2sym_3d, int nf)
      {
	int neq = 3;
	double c_b2 = .10;
	// Function Body
	eval.resize(neq);
	v.resize(neq, neq);

	for ( int i = 1; i <= neq; i++ ) {
	  eval.at(i) = this->at(i, i);
	}

	double tol = pow(c_b2, nf);
	double sum = 0.0;
	for ( int i = 1; i <= neq; ++i ) {
	  for ( int j = 1; j <= neq; ++j ) {
            sum += fabs( this->at(i, j) );
            v.at(i, j) = 0.0;
	  }
	  v.at(i, i) = 1.0;
	}
	
	if ( sum <= 0.0 ) {
	  return 0;
	}


	// ---- REDUCE MATRIX TO DIAGONAL ---------------- 
	int ite = 0;
	double ssum;
	do {
	  ssum = 0.0;
	  for ( int j = 2; j <= neq; ++j ) {
            int ih = j - 1;
            for ( int i = 1; i <= ih; ++i ) {
	      if ( ( fabs( this->at(i, j) ) / sum ) > tol ) {
		ssum += fabs( this->at(i, j) );
		// ---- CALCULATE ROTATION ANGLE ----------------- 
		double aa = atan2( this->at(i, j) * 2.0, eval.at(i) - eval.at(j) ) /  2.0;
		double si = sin(aa);
		double co = cos(aa);
		// ---- MODIFY "I" AND "J" COLUMNS OF "A" AND "V"
		for ( int k = 1; k < i; ++k ) {
		  double tt = this->at(k, i);
		  this->at(k, i) = co * tt + si * this->at(k, j);
		  this->at(k, j) = -si * tt + co * this->at(k, j);
		  tt = v.at(k, i);
		  v.at(k, i) = co * tt + si *v.at(k, j);
		  v.at(k, j) = -si * tt + co *v.at(k, j);
		}
		// diagonal term (i,i)
		double tt = eval.at(i);
		eval.at(i) = co * tt + si * this->at(i, j);
		double aij = -si * tt + co *this->at(i, j);
		tt = v.at(i, i);
		v.at(i, i) = co * tt + si *v.at(i, j);
		v.at(i, j) = -si * tt + co *v.at(i, j);
		for ( int k = i + 1; k < j; ++k ) {
		  double tt = this->at(i, k);
		  this->at(i, k) = co * tt + si * this->at(k, j);
		  this->at(k, j) = -si * tt + co * this->at(k, j);
		  tt = v.at(k, i);
		  v.at(k, i) = co * tt + si *v.at(k, j);
		  v.at(k, j) = -si * tt + co *v.at(k, j);
		}
		// diagonal term (j,j)
		tt = this->at(i, j);
		double aji = co * tt + si *eval.at(j);
		eval.at(j) = -si * tt + co *eval.at(j);
		
		tt = v.at(j, i);
		v.at(j, i) = co * tt + si *v.at(j, j);
		v.at(j, j) = -si * tt + co *v.at(j, j);
		//
		for ( int k = j + 1; k <= neq; ++k ) {
		  double tt = this->at(i, k);
		  this->at(i, k) = co * tt + si *this->at(j, k);
		  this->at(j, k) = -si * tt + co *this->at(j, k);
		  tt = v.at(k, i);
		  v.at(k, i) = co * tt + si *v.at(k, j);
		  v.at(k, j) = -si * tt + co *v.at(k, j);
		}
		
		// ---- MODIFY DIAGONAL TERMS --------------------
		eval.at(i) = co * eval.at(i) + si * aji;
		eval.at(j) = -si * aij + co *eval.at(j);
		this->at(i, j) = 0.0;
	      } else {
		// ---- A(I,J) MADE ZERO BY ROTATION ------------- 
		;
	      }
            }
	  }
	  
	  //
	  if ( ++ite > 50 ) {
            OOFEM_ERROR("too many iterations");
	  }
	} while ( fabs(ssum) / sum > tol );
	
	// restore original matrix
	for ( int i = 1; i <= neq; i++ ) {
	  for ( int j = i; j <= neq; j++ ) {
            this->at(i, j) = this->at(j, i);
	  }
	}
	

      }

*/


    
  };


    class Tensor2_3d : public Tensor2<double, 3, 3>
  {
  public:
    using Tensor2<double, 3,3> :: Tensor2;
    Tensor2_3d(const oofem::FloatArrayF<9> &array){
      this->data[0][0] = array.at(1);
      this->data[0][1] = array.at(6);
      this->data[0][2] = array.at(5);
      this->data[1][0] = array.at(9);
      this->data[1][1] = array.at(2);
      this->data[1][2] = array.at(4);
      this->data[2][0] = array.at(8);
      this->data[2][1] = array.at(7);
      this->data[2][2] = array.at(3);
    }

    const inline FloatArrayF<9> to_voigt_form() 
    {
      return {
	this->operator()(0, 0),
	  this->operator()(1, 1),
	  this->operator()(2, 2),
	  this->operator()(1, 2),
	  this->operator()(0, 2),
	  this->operator()(0, 1),
	  this->operator()(2, 1), 
	  this->operator()(2, 0),
	  this->operator()(1, 0)
	  };
    }


    inline Tensor2_3d compute_cofactor() const
    {

      Tensor2_3d cofF;
      cofF(i_3, j_3) = 0.5 * this->compute_tensor_cross_product(this)(i_3, j_3);
      return cofF;
    }


    inline double compute_determinant() const
    {
      return ( 1./6. *  this->compute_tensor_cross_product(this)(m_3,n_3) * this->operator()(m_3,n_3) );
    }

 


    inline std::pair<double, Tensor2_3d>  compute_determinant_and_cofactor() const
      {

	Tensor2_3d cofF;
	cofF(i_3,j_3)=  0.5 * this->compute_tensor_cross_product(this)(i_3, j_3);
	auto detF =  1./3. * cofF(i_3,j_3) * this->operator()(i_3,j_3);
	return {detF, cofF};
      }
 

    inline Tensor2_3d compute_inverse() const
    {
      Tensor2_3d iF;
      auto [J, cofF] = this->compute_determinant_and_cofactor();
      iF(i_3,j_3) = 1. / J * cofF(j_3,i_3); 
      return iF;
    }


    inline Tensor2_3d compute_tensor_cross_product(const Tensor2_3d &B) const
    {
      Tensor2_3d C (this->operator()(1,1) * B(2,2) - this->operator()(1,2) * B(2,1) - this->operator()(2,1) * B(1,2) + this->operator()(2,2) * B(1,1), this->operator()(1,2)*B(2,0) - this->operator()(1,0)*B(2,2) + this->operator()(2,0)*B(1,2) - this->operator()(2,2)*B(1,0), this->operator()(1,0)*B(2,1) - this->operator()(1,1)*B(2,0) - this->operator()(2,0)*B(1,1) + this->operator()(2,1)*B(1,0), this->operator()(0,2)*B(2,1) - this->operator()(0,1)*B(2,2) + this->operator()(2,1)*B(0,2) - this->operator()(2,2)*B(0,1), this->operator()(0,0)*B(2,2) - this->operator()(0,2)*B(2,0) - this->operator()(2,0)*B(0,2) + this->operator()(2,2)*B(0,0),this->operator()(0,1)*B(2,0) - this->operator()(0,0)*B(2,1) + this->operator()(2,0)*B(0,1) - this->operator()(2,1)*B(0,0), this->operator()(0,1)*B(1,2) - this->operator()(0,2)*B(1,1) - this->operator()(1,1)*B(0,2) + this->operator()(1,2)*B(0,1), this->operator()(0,2)*B(1,0) - this->operator()(0,0)*B(1,2) + this->operator()(1,0)*B(0,2) - this->operator()(1,2)*B(0,0), this->operator()(0,0)*B(1,1) - this->operator()(0,1)*B(1,0) - this->operator()(1,0)*B(0,1) + this->operator()(1,1)*B(0,0));
      return C;    
    }


  
    

    
  };


 class Tensor4_3d : public Tensor4<double, 3, 3, 3, 3>
  {
  public:
    using Tensor4 <double, 3, 3, 3, 3> :: Tensor4;
    Tensor4_3d(const oofem::FloatMatrixF<9,9> &mat){
      // the first column
      this->data[0][0][0][0] = mat.at(1,1);
      this->data[0][1][0][0] = mat.at(6,1);
      this->data[0][2][0][0] = mat.at(5,1);
      this->data[1][0][0][0] = mat.at(9,1);
      this->data[1][1][0][0] = mat.at(2,1);
      this->data[1][2][0][0] = mat.at(4,1);
      this->data[2][0][0][0] = mat.at(8,1);
      this->data[2][1][0][0] = mat.at(7,1);
      this->data[2][2][0][0] = mat.at(3,1);
      // the second column
      this->data[0][0][1][1] = mat.at(1,2);
      this->data[0][1][1][1] = mat.at(6,2);
      this->data[0][2][1][1] = mat.at(5,2);
      this->data[1][0][1][1] = mat.at(9,2);
      this->data[1][1][1][1] = mat.at(2,2);
      this->data[1][2][1][1] = mat.at(4,2);
      this->data[2][0][1][1] = mat.at(8,2);
      this->data[2][1][1][1] = mat.at(7,2);
      this->data[2][2][1][1] = mat.at(3,2);
      // the third column
      this->data[0][0][2][2] = mat.at(1,3);
      this->data[0][1][2][2] = mat.at(6,3);
      this->data[0][2][2][2] = mat.at(5,3);
      this->data[1][0][2][2] = mat.at(9,3);
      this->data[1][1][2][2] = mat.at(2,3);
      this->data[1][2][2][2] = mat.at(4,3);
      this->data[2][0][2][2] = mat.at(8,3);
      this->data[2][1][2][2] = mat.at(7,3);
      this->data[2][2][2][2] = mat.at(3,3);
      // the fourth column
      this->data[0][0][1][2] = mat.at(1,4);
      this->data[0][1][1][2] = mat.at(6,4);
      this->data[0][2][1][2] = mat.at(5,4);
      this->data[1][0][1][2] = mat.at(9,4);
      this->data[1][1][1][2] = mat.at(2,4);
      this->data[1][2][1][2] = mat.at(4,4);
      this->data[2][0][1][2] = mat.at(8,4);
      this->data[2][1][1][2] = mat.at(7,4);
      this->data[2][2][1][2] = mat.at(3,4);
      //the fifth column
      this->data[0][0][0][2] = mat.at(1,5);
      this->data[0][1][0][2] = mat.at(6,5);
      this->data[0][2][0][2] = mat.at(5,5);
      this->data[1][0][0][2] = mat.at(9,5);
      this->data[1][1][0][2] = mat.at(2,5);
      this->data[1][2][0][2] = mat.at(4,5);
      this->data[2][0][0][2] = mat.at(8,5);
      this->data[2][1][0][2] = mat.at(7,5);
      this->data[2][2][0][2] = mat.at(3,5);
      //the sixth column
      this->data[0][0][0][1] = mat.at(1,6);
      this->data[0][1][0][1] = mat.at(6,6);
      this->data[0][2][0][1] = mat.at(5,6);
      this->data[1][0][0][1] = mat.at(9,6);
      this->data[1][1][0][1] = mat.at(2,6);
      this->data[1][2][0][1] = mat.at(4,6);
      this->data[2][0][0][1] = mat.at(8,6);
      this->data[2][1][0][1] = mat.at(7,6);
      this->data[2][2][0][1] = mat.at(3,6);
      //the seventh column
      this->data[0][0][2][1] = mat.at(1,7);
      this->data[0][1][2][1] = mat.at(6,7);
      this->data[0][2][2][1] = mat.at(5,7);
      this->data[1][0][2][1] = mat.at(9,7);
      this->data[1][1][2][1] = mat.at(2,7);
      this->data[1][2][2][1] = mat.at(4,7);
      this->data[2][0][2][1] = mat.at(8,7);
      this->data[2][1][2][1] = mat.at(7,7);
      this->data[2][2][2][1] = mat.at(3,7);
      //the eight column
      this->data[0][0][2][0] = mat.at(1,8);
      this->data[0][1][2][0] = mat.at(6,8);
      this->data[0][2][2][0] = mat.at(5,8);
      this->data[1][0][2][0] = mat.at(9,8);
      this->data[1][1][2][0] = mat.at(2,8);
      this->data[1][2][2][0] = mat.at(4,8);
      this->data[2][0][2][0] = mat.at(8,8);
      this->data[2][1][2][0] = mat.at(7,8);
      this->data[2][2][2][0] = mat.at(3,8);
      //the nineght column
      this->data[0][0][1][0] = mat.at(1,9);
      this->data[0][1][1][0] = mat.at(6,9);
      this->data[0][2][1][0] = mat.at(5,9);
      this->data[1][0][1][0] = mat.at(9,9);
      this->data[1][1][1][0] = mat.at(2,9);
      this->data[1][2][1][0] = mat.at(4,9);
      this->data[2][0][1][0] = mat.at(8,9);
      this->data[2][1][1][0] = mat.at(7,9);
      this->data[2][2][1][0] = mat.at(3,9);
    }

    inline FloatMatrixF<9,9> to_voigt_form()
      {
	return {this->operator()(0,0,0,0),this->operator()(1,1,0,0),this->operator()(2,2,0,0),this->operator()(1,2,0,0),this->operator()(0,2,0,0),this->operator()(0,1,0,0),this->operator()(2,1,0,0),this->operator()(2,0,0,0),this->operator()(1,0,0,0),this->operator()(0,0,1,1),this->operator()(1,1,1,1),this->operator()(2,2,1,1),this->operator()(1,2,1,1),this->operator()(0,2,1,1),this->operator()(0,1,1,1),this->operator()(2,1,1,1),this->operator()(2,0,1,1),this->operator()(1,0,1,1),this->operator()(0,0,2,2),this->operator()(1,1,2,2),this->operator()(2,2,2,2),this->operator()(1,2,2,2),this->operator()(0,2,2,2),this->operator()(0,1,2,2),this->operator()(2,1,2,2),this->operator()(2,0,2,2),this->operator()(1,0,2,2),this->operator()(0,0,1,2),this->operator()(1,1,1,2),this->operator()(2,2,1,2),this->operator()(1,2,1,2),this->operator()(0,2,1,2),this->operator()(0,1,1,2),this->operator()(2,1,1,2),this->operator()(2,0,1,2),this->operator()(1,0,1,2),this->operator()(0,0,0,2),this->operator()(1,1,0,2),this->operator()(2,2,0,2),this->operator()(1,2,0,2),this->operator()(0,2,0,2),this->operator()(0,1,0,2),this->operator()(2,1,0,2),this->operator()(2,0,0,2),this->operator()(1,0,0,2),this->operator()(0,0,0,1),this->operator()(1,1,0,1),this->operator()(2,2,0,1),this->operator()(1,2,0,1),this->operator()(0,2,0,1),this->operator()(0,1,0,1),this->operator()(2,1,0,1),this->operator()(2,0,0,1),this->operator()(1,0,0,1),this->operator()(0,0,2,1),this->operator()(1,1,2,1),this->operator()(2,2,2,1),this->operator()(1,2,2,1),this->operator()(0,2,2,1),this->operator()(0,1,2,1),this->operator()(2,1,2,1),this->operator()(2,0,2,1),this->operator()(1,0,2,1),this->operator()(0,0,2,0),this->operator()(1,1,2,0),this->operator()(2,2,2,0),this->operator()(1,2,2,0),this->operator()(0,2,2,0),this->operator()(0,1,2,0),this->operator()(2,1,2,0),this->operator()(2,0,2,0),this->operator()(1,0,2,0),this->operator()(0,0,1,0),this->operator()(1,1,1,0),this->operator()(2,2,1,0),this->operator()(1,2,1,0),this->operator()(0,2,1,0),this->operator()(0,1,1,0),this->operator()(2,1,1,0),this->operator()(2,0,1,0),this->operator()(1,0,1,0)};
      }

    

  inline Tensor4_3d compute_tensor_cross_product() const
  {
      Tensor4_3d Ax(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
      /*      Ax(1,1,0,0) =  this->operator()(2,2);
      Ax(1,2,0,0) = -this->operator()(2,1);
      Ax(2,1,0,0) = -this->operator()(1,2);
      Ax(2,2,0,0) =  this->operator()(1,1);
    
      Ax(0,1,1,0) = -this->operator()(2,2);
      Ax(0,2,1,0) =  this->operator()(2,1); 
      Ax(2,1,1,0) =  this->operator()(0,2);
      Ax(2,2,1,0) = -this->operator()(0,1);

      Ax(0,1,2,0) =  this->operator()(1,2);
      Ax(0,2,2,0) = -this->operator()(1,1);
      Ax(1,1,2,0) = -this->operator()(0,2);
      Ax(1,2,2,0) =  this->operator()(0,1);

      Ax(1,0,0,1) = -this->operator()(2,2);
      Ax(1,2,0,1) =  this->operator()(2,0);
      Ax(2,0,0,1) =  this->operator()(1,2); 
      Ax(2,2,0,1) = -this->operator()(1,0);
 
      Ax(0,0,1,1) =  this->operator()(2,2);
      Ax(0,2,1,1) = -this->operator()(2,0);
      Ax(2,0,1,1) = -this->operator()(0,2);
      Ax(2,2,1,1) =  this->operator()(0,0);
  
      Ax(0,0,2,1) = -this->operator()(1,2);
      Ax(0,2,2,1) =  this->operator()(1,0);
      Ax(1,0,2,1) =  this->operator()(0,2);
      Ax(1,2,2,1) = -this->operator()(0,0);

      Ax(1,0,0,2) =  this->operator()(2,1);
      Ax(1,1,0,2) = -this->operator()(2,0);
      Ax(2,0,0,2) = -this->operator()(1,1);
      Ax(2,1,0,2) =  this->operator()(1,0);

      Ax(0,0,1,2) = -this->operator()(2,1);   
      Ax(0,1,1,2) =  this->operator()(2,0);
      Ax(2,0,1,2) =  this->operator()(0,1);
      Ax(2,1,1,2) = -this->operator()(0,0);

      Ax(0,0,2,2) =  this->operator()(1,1);
      Ax(0,1,2,2) = -this->operator()(1,0);
      Ax(1,0,2,2) = -this->operator()(0,1);
      Ax(1,1,2,2) =  this->operator()(0,0);
      */
      return Ax;    
    }
 

    
    void __attribute__((noinline)) printYourself() const 
    // Prints the receiver on screen.
    {
      printf("4th order tensor");
      for ( int i = 1; i <= 3; i++ ) {
	for ( int j = 1; j <= 3; j++ ) {
	  for ( int k = 1; k <= 3; k++ ) {
	    for ( int l = 1; l <= 3; l++ ) {
	      printf( "Component %d %d %d %d is: %10.3e \n", i,j,k,l,this->data[i][j][k][l] );
	    }
	  }
	}
      }  
    }
    
    void __attribute__((noinline)) printComponent(int i, int j, int k, int l) const 
    // Prints the receiver on screen.
    {
      printf( "Component %d %d %d %d is: %10.3e\n  ", i,j,k,l,this->data[i][j][k][l] );
    }
    

  };
    
} // end namespace oofem


