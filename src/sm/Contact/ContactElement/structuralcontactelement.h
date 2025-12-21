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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef structuralcontactelement_h
#define structuralcontactelement_h

#include "Contact/contactelement.h"
#include "fei3dquadlin.h"
#include "fei3dtrlin.h"
#include "fei2dlinelin.h"


#define _IFT_StructuralContactElement_QuadLin_Name "structuralcontactelement_quadlin"
#define _IFT_StructuralContactElement_LineLin_Name "structuralcontactelement_linelin"
#define _IFT_StructuralContactElement_TrLin_Name "structuralcontactelement_trlin"

namespace oofem {

/**
 * Structural contact elements.
 * 
 */
class StructuralContactElement : public ContactElement
{
protected:
  int nsd = 0;

public:
  StructuralContactElement(int n, Domain * d);
  void computeNmatrixAt(const FloatArray &lcoord, FloatMatrix &answer) override;
  const char *giveClassName() const override { return "StructuralContactElement"; }
  //  const char *giveInputRecordName() const override {return  _IFT_StructuralContactElement_QuadLin_Name;}
protected:
  void computeGaussPoints() override;

};

class FEI3dQuadLin;

class StructuralContactElement_QuadLin : public StructuralContactElement
{
protected:
  static FEI3dQuadLin interpolation;
  int NCP = 0;
  //add integration rule typ
public:
  StructuralContactElement_QuadLin(int n, Domain * d);
  const char *giveClassName() const override { return "StructuralContactElement_QuadLin"; }
  const char *giveInputRecordName() const override {return _IFT_StructuralContactElement_QuadLin_Name;}
  FEInterpolation *giveInterpolation() const override;
  FloatArray computeNormalVectorAt(const FloatArray &lCoords) override;
  void giveDofManDofIDMask(int inode, IntArray &answer) const override;
  Element_Geometry_Type giveGeometryType() const override {return EGT_line_2;}

};


class FEI3dTrLin;

class StructuralContactElement_TrLin : public StructuralContactElement
{
protected:
  static FEI3dTrLin interpolation;
  int NCP = 0;
  //add integration rule typ
public:
  StructuralContactElement_TrLin(int n, Domain * d);
  const char *giveClassName() const override { return "StructuralContactElement_TrLin"; }
  const char *giveInputRecordName() const override {return _IFT_StructuralContactElement_TrLin_Name;}
  FEInterpolation *giveInterpolation() const override;
  FloatArray computeNormalVectorAt(const FloatArray &lCoords) override;
  void giveDofManDofIDMask(int inode, IntArray &answer) const override;
  Element_Geometry_Type giveGeometryType() const override {return EGT_line_2;}

};



  

class FEI2dLineLin;

class StructuralContactElement_LineLin : public StructuralContactElement
{
protected:
  static FEI2dLineLin interpolation;
public:
  StructuralContactElement_LineLin(int n, Domain * d);
  const char *giveClassName() const override { return "StructuralContactElement_LineLin"; }
  const char *giveInputRecordName() const override {return _IFT_StructuralContactElement_LineLin_Name;}
  FEInterpolation *giveInterpolation() const override;
  FloatArray computeNormalVectorAt(const FloatArray &lCoords) override;
  void giveDofManDofIDMask(int inode, IntArray &answer) const override;
  Element_Geometry_Type giveGeometryType() const override {return EGT_line_1;}
  
};

  


} // end namespace oofem
#endif // structuralcontactelement_h
