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

#ifndef thermalcontactelement_h
#define thermalcontactelement_h

#include "intarray.h"
#include "Contact/contactelement.h"
#include "fei3dquadlin.h"
#include "fei3dtrlin.h"
#include "lobattoir.h"


#define _IFT_ThermalContactElement_QuadLin_Name "thermalcontactelement_quadlin"
#define _IFT_ThermalContactElement_TrLin_Name "thermalcontactelement_trlin"


namespace oofem {
/**
 * @class ThermalContactElement
 * @brief Base class for thermal contact interface elements.
 *
 * This is an abstract contact interface element used in thermal contact
 * formulations. It provides:
 * - Evaluation of the shape function matrix N at local coordinates
 * - Evaluation of the surface normal vector at local coordinates
 * - Definition and construction of Gauss integration points on the interface
 *
 * Derived classes specify the concrete geometry (triangle/quad) and
 * interpolation rule.
 */
/**
 * General thermal contact element.
 */
class ThermalContactElement : public ContactElement
{
protected:
 

public:
  ThermalContactElement(int n, Domain * d);
  void computeNmatrixAt(const FloatArray &lcoord, FloatMatrix &answer) override;
  const char *giveClassName() const override { return "ThermalContactElement"; }
  const char *giveInputRecordName() const override {return "ThermalContactElement";}
  virtual FloatArray computeNormalVectorAt(const FloatArray &lCoords) override;
 
protected:
  void computeGaussPoints() override;


  
};


class FEI3dQuadLin;

class ThermalContactElement_QuadLin : public ThermalContactElement
{
protected:
  static FEI3dQuadLin interpolation;
  //add integration rule typ

public:
  ThermalContactElement_QuadLin(int n, Domain * d): ThermalContactElement(n, d) { numberOfDofMans = 4; numberOfGaussPoints = 4;}
  const char *giveClassName() const override { return "ThermalContactElement_QuadLin"; }
  const char *giveInputRecordName() const override {return _IFT_ThermalContactElement_QuadLin_Name;}
  virtual Element_Geometry_Type giveGeometryType() const override {
    return EGT_quad_1 ;
  }
  FEInterpolation * giveInterpolation() const override;

};


class FEI3dTrLin;

class ThermalContactElement_TrLin : public ThermalContactElement
{
protected:
  static FEI3dTrLin interpolation;
public:
  ThermalContactElement_TrLin(int n, Domain * d): ThermalContactElement(n, d){  numberOfDofMans = 3;numberOfGaussPoints = 3;}
  const char *giveClassName() const override { return "ThermalContactElement_TrLin"; }
  const char *giveInputRecordName() const override {return _IFT_ThermalContactElement_TrLin_Name;}
  virtual Element_Geometry_Type giveGeometryType() const override {
    return EGT_triangle_1 ;
  }
  FEInterpolation * giveInterpolation() const override;
};


  
  

} // end namespace oofem
#endif // thermalcontactelement_quadlin_h
