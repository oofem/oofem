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

#include "thermalcontactelement.h"
#include "feinterpol.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "crosssection.h"


namespace oofem {
REGISTER_Element(ThermalContactElement_QuadLin);
REGISTER_Element(ThermalContactElement_TrLin);


FEI3dTrLin ThermalContactElement_TrLin :: interpolation;
FEI3dQuadLin ThermalContactElement_QuadLin :: interpolation;

  

  ThermalContactElement :: ThermalContactElement(int n, Domain *aDomain) : ContactElement(n, aDomain)
{
}


void
ThermalContactElement :: computeNmatrixAt(const FloatArray &lcoord, FloatMatrix &answer)
{
    FloatArray n;
    this->giveInterpolation()->evalN( n, lcoord,  FEIElementGeometryWrapper(this));
    answer.beNMatrixOf(n, 1);
}     



FloatArray
ThermalContactElement :: computeNormalVectorAt(const FloatArray &lCoords)
{
  FloatArray normal;
  auto norm = this->giveInterpolation()->boundarySurfaceEvalNormal(normal, 0, lCoords, FEIElementGeometryWrapper(this));
  return normal*norm;
}


void ThermalContactElement :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.XF
{
  if ( integrationRulesArray.size() == 0 ) {
    integrationRulesArray.resize(1);
    //integrationRulesArray [ 0 ] = std::make_unique< GaussIntegrationRule >(1, this);
    integrationRulesArray [ 0 ] = std::make_unique<LobattoIntegrationRule>(1, this);
    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
  }
}


  

  

 FEInterpolation *
 ThermalContactElement_QuadLin :: giveInterpolation() const
 {
   return & interpolation;
 }


  FEInterpolation *
  ThermalContactElement_TrLin :: giveInterpolation() const
  {
    return & interpolation;
  }


  
} // end namespace oofem
