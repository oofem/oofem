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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#include "winklermodel.h"
#include "../sm/Materials/structuralms.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(WinklerMaterial);

WinklerMaterial:: WinklerMaterial (int n, Domain* d): StructuralMaterial(n, d) 
{ }

WinklerMaterial::~WinklerMaterial()
{ }

IRResultType
WinklerMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, c1, _IFT_WinklerMaterial_C1);
    globalFromulation = true;
    int var=1;
    IR_GIVE_OPTIONAL_FIELD(ir, var,  _IFT_WinklerMaterial_globalFlag);
    globalFromulation=var;

    return StructuralMaterial :: initializeFrom(ir);
}


void
WinklerMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);

    input.setField(this->c1, _IFT_WinklerMaterial_C1);
}


void 
WinklerMaterial::giveRealStressVector_2dPlateSubSoil(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
{
    FloatMatrix tangent;
    this->give2dPlateSubSoilStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, reducedE);

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    status->letTempStrainVectorBe(reducedE);
    status->letTempStressVectorBe(answer);
}

void 
WinklerMaterial::giveRealStressVector_3dBeamSubSoil(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
{
    FloatMatrix tangent;
    this->give3dBeamSubSoilStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, reducedE);

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    status->letTempStrainVectorBe(reducedE);
    status->letTempStressVectorBe(answer);
}


void
WinklerMaterial :: give2dPlateSubSoilStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(3, 3);
    answer.zero();

    answer.at(1, 1) = c1.at(1);
    //answer.at(2, 2) = c2;
    //answer.at(3, 3) = c2;
}


void
WinklerMaterial::give3dBeamSubSoilStiffMtrx(FloatMatrix &answer, MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep)
{

  if (this->c1.giveSize() == 6) {
    answer.beDiagonal(this->c1);

    if (globalFromulation) {
      Beam3dSubsoilMaterialInterface *ei = (Beam3dSubsoilMaterialInterface*) gp->giveElement()->giveInterface(Beam3dSubsoilMaterialInterfaceType);
      FloatMatrix T;
      if (ei) {
	ei->B3SSMI_getUnknownsGtoLRotationMatrix (T);
	answer.rotatedWith(T, 't');
      } else {
	OOFEM_ERROR("Beam3dSubsoilMaterialInterface required from element");
      }
    }
    
  } else {
    OOFEM_ERROR ("C1 attribute size error (shouldequal to 6 for 3dBeamSubsoil mode)");
  }
}


MaterialStatus *
WinklerMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}


int
WinklerMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
  return ((mode == _2dPlateSubSoil) || (mode == _3dBeamSubSoil));
}



  
} // end namespace oofem
