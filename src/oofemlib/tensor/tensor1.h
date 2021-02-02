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

#include "floatarray.h"
#include "mathfem.h"
#include "index.h"
#pragma once

using namespace FTensor;



/**
 * Functions providing basic operations with tensors 
 * based on FTensor library
 * It provides transformation from the second and fourth order tensors in 3d to FloatArrayF<9> and FloatArray<9,9> 
 * It provides calculation of determinant, cofactor, inverse, and tensor cross product
 * More to be added gradually
 * @author Martin Horak
 **/

namespace oofem {


  class Tensor1_3d : public Tensor1<double, 3>
  {
  public:
    using Tensor1<double, 3> :: Tensor1;
    Tensor1_3d(const oofem::FloatArrayF<3> &array){
      this->data[0] = array.at(1);
      this->data[1] = array.at(2);
      this->data[2] = array.at(3);
    }

    const inline FloatArrayF<3> to_voigt_form() 
    {
      return {
	this->operator()(0),
	  this->operator()(1),
	  this->operator()(2),
	  };
    }

  };
  
    
} // end namespace oofem


