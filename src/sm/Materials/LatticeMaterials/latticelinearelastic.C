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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "latticelinearelastic.h"
#include "latticematstatus.h"
#include "latticestructuralmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LatticeLinearElastic);

// constructor which creates a dummy material without a status and without random extension interface
LatticeLinearElastic :: LatticeLinearElastic(int n, Domain *d, double e, double a1, double a2) :
    LatticeStructuralMaterial(n, d),
    eNormalMean(e),
    alphaOne(a1),
    alphaTwo(a2)
{}


bool
LatticeLinearElastic :: hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _3dLattice );
}


void
LatticeLinearElastic :: initializeFrom(InputRecord &ir)
{
    LatticeStructuralMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);

    //Young's modulus of the material that the network element is made of
    IR_GIVE_FIELD(ir, this->eNormalMean, _IFT_LatticeLinearElastic_e); // Macro

    //Parameter which relates the shear stiffness to the normal stiffness. Default is 1
    alphaOne = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaOne, _IFT_LatticeLinearElastic_a1); // Macro

    //Parameter which is used for the definition of bending stiffness. Default is 0.
    alphaTwo = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaTwo, _IFT_LatticeLinearElastic_a2); // Macro

    localRandomType = 0; //Default: No local random field
    IR_GIVE_OPTIONAL_FIELD(ir, localRandomType, _IFT_LatticeLinearElastic_localrandomtype); // Macro
    if ( localRandomType == 1 ) { //Gaussian random generator
        coefficientOfVariation = 0.;
        IR_GIVE_FIELD(ir, coefficientOfVariation, _IFT_LatticeLinearElastic_cov); // Macro
    }


    this->cAlpha = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, cAlpha, _IFT_LatticeLinearElastic_calpha);
}

std::unique_ptr<MaterialStatus> 
LatticeLinearElastic :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<LatticeMaterialStatus>(gp);
}

MaterialStatus *
LatticeLinearElastic :: giveStatus(GaussPoint *gp) const
{
    if (!gp->hasMaterialStatus()) {
        // create a new one
        gp->setMaterialStatus(this->CreateStatus(gp));
        this->_generateStatusVariables(gp);
    }

    return static_cast<MaterialStatus*> (gp->giveMaterialStatus());
}


FloatArrayF< 6 >
LatticeLinearElastic :: giveLatticeStress3d(const FloatArrayF< 6 > &strain,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{
    auto status = static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    // subtract stress independent part
    auto reducedStrain = strain;
    FloatArray indepStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( indepStrain.giveSize() > 0 ) {
        reducedStrain -= FloatArrayF< 6 >(indepStrain);
    }

    auto stiffnessMatrix = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);
    auto stress = dot(stiffnessMatrix, reducedStrain);

    //Read in fluid pressures from structural element if this is not a slave problem
    FloatArray pressures;
    if ( !domain->giveEngngModel()->giveMasterEngngModel() ) {
        static_cast< LatticeStructuralElement * >( gp->giveElement() )->givePressures(pressures);
    }

    double waterPressure = 0.;
    for ( int i = 0; i < pressures.giveSize(); i++ ) {
        waterPressure += 1. / pressures.giveSize() * pressures [ i ];
    }

    stress.at(1) += waterPressure;

    //Set all temp values
    status->letTempLatticeStrainBe(strain);
    status->letTempLatticeStressBe(stress);

    return stress;
}


void LatticeLinearElastic :: giveRandomParameters(FloatArray &param)
{
    param.resize(3);
    param.zero();
    param.at(1) = localRandomType;

    if ( localRandomType == 1 ) { //Gaussian
        param.at(2) = coefficientOfVariation;
    } else {
        OOFEM_ERROR("Error: Unknown local random type:\n randomtype 1 = Gaussian\n");
    }
}


Interface *
LatticeLinearElastic :: giveInterface(InterfaceType type)
{
    return nullptr;
}


FloatMatrixF< 6, 6 >
LatticeLinearElastic :: give3dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime) const
{
    //Needed to make sure that status exists before random values are requested for elastic stiffness. Problem is that gp->giveMaterialStatus does not check if status exist already
  
  //static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

    FloatArrayF< 6 >d = {
        1.,
        this->alphaOne, // shear
        this->alphaOne, // shear
        this->alphaTwo, // torsion
        this->alphaTwo, // torsion
        this->alphaTwo, // torsion
    };

    return diag(d * this->give(eNormal_ID, gp) * this->eNormalMean);
}


FloatMatrixF< 3, 3 >
LatticeLinearElastic :: give2dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime) const
{
  //Needed to make sure that status exists before random values are requested for elastic stiffness. Problem is that gp->giveMaterialStatus does not check if status exist already

  //static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

  FloatArrayF< 3 >d = {
        1.,
        this->alphaOne, // shear
        this->alphaTwo, // torsion
    };

    return diag(d * this->give(eNormal_ID, gp) * this->eNormalMean);
}


FloatArrayF< 6 >
LatticeLinearElastic :: giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const
//
// returns a FloatArray(6) of initial strain vector
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    double alpha = this->give(tAlpha, gp);

    //Option to add a eigendisplacement instead of strain
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    alpha += this->cAlpha / length;

    return {
               alpha, 0., 0., 0., 0., 0.
    };
}


double
LatticeLinearElastic :: give(int aProperty, GaussPoint *gp) const
{
   this->giveStatus(gp);
  
    double answer;
    if ( RandomMaterialExtensionInterface :: give(aProperty, gp, answer) ) {
        if ( answer < 0.1 ) { //Introduce cut off to avoid numerical problems
            answer = 0.1;
        } else if ( answer > 10 ) {
            answer = 10;
        }
        return answer;
    } else if ( aProperty == eNormal_ID ) {
        return 1.;
    } else if ( aProperty == 'E' ) {
        return this->eNormalMean;
    } else {
        return LatticeStructuralMaterial :: give(aProperty, gp);
    }
}
}
