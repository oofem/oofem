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

#include "rankinepm.h"
#include "Materials/isolinearelasticmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(RankinePlasticMaterial);

RankinePlasticMaterial :: RankinePlasticMaterial(int n, Domain *d) : MPlasticMaterial(n, d)
{
    //
    // constructor
    //
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    this->nsurf = 3;
    this->rmType = mpm_CuttingPlane;
}

RankinePlasticMaterial :: ~RankinePlasticMaterial()
{ }


IRResultType
RankinePlasticMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = MPlasticMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;
    result = linearElasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, k, _IFT_RankinePlasticMaterial_ry);
    return IRRT_OK;
}


double
RankinePlasticMaterial :: computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                                              const FloatArray &stressSpaceHardeningVars)
{
    FloatArray princStress(3);
    this->computePrincipalValues(princStress, stressVector, principal_stress);

    return princStress.at(isurf) - this->k;
}

void
RankinePlasticMaterial :: computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                                      const FloatArray &stressSpaceHardeningVars)
{
    FloatArray princStress(3);
    FloatMatrix t(3, 3);

    // compute principal stresses and their directions
    this->computePrincipalValDir(princStress, t, stressVector, principal_stress);

    //derivation through stress transformation. The transformation matrix is stored in t.
    answer.resize(6);
    answer.at(1) = t.at(1, isurf) * t.at(1, isurf); //xx = 11
    answer.at(2) = t.at(2, isurf) * t.at(2, isurf); //yy = 22
    answer.at(3) = t.at(3, isurf) * t.at(3, isurf); //zz = 33
    answer.at(4) = t.at(2, isurf) * t.at(3, isurf); //yz = 23
    answer.at(5) = t.at(1, isurf) * t.at(3, isurf); //xz = 13
    answer.at(6) = t.at(1, isurf) * t.at(2, isurf); //xy = 12

    //crossSection->giveReducedCharacteristicVector(answer, gp, fullAnswer);
}

void
RankinePlasticMaterial :: computeHardeningReducedModuli(FloatMatrix &answer,
                                                        GaussPoint *gp,
                                                        const FloatArray &strainSpaceHardeningVariables,
                                                        TimeStep *tStep)
{
    answer.clear();
}

void
RankinePlasticMaterial :: computeStressSpaceHardeningVarsReducedGradient(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp,
                                                                         const FloatArray &stressVector,
                                                                         const FloatArray &stressSpaceHardeningVars)
{
    answer.clear();
}


void
RankinePlasticMaterial :: computeReducedGradientMatrix(FloatMatrix &answer, int isurf,
                                                       GaussPoint *gp,
                                                       const FloatArray &stressVector,
                                                       const FloatArray &stressSpaceHardeningVars)
{
    answer.clear();
}


void
RankinePlasticMaterial :: computeStressSpaceHardeningVars(FloatArray &answer, GaussPoint *gp,
                                                          const FloatArray &strainSpaceHardeningVariables)
{
    answer.clear();
}


MaterialStatus *
RankinePlasticMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new MPlasticMaterialStatus(1, this->giveDomain(), gp, this->giveSizeOfReducedHardeningVarsVector(gp));
}
} // end namespace oofem
