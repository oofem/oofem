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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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
LatticeLinearElastic :: LatticeLinearElastic(int n, Domain *d, double e0, double a1, double a2) : LatticeStructuralMaterial(n, d)

{
    eNormalMean = e0;
    alphaOne = a1;
    alphaTwo = a2;
    nu = 0.;
}



LatticeLinearElastic :: ~LatticeLinearElastic()
//
// destructor
//
{}

int
LatticeLinearElastic :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _2dLattice ) || ( mode == _3dLattice ) || ( mode == _1dLattice ) ) {
        return 1;
    }

    return 0;
}


void
LatticeLinearElastic :: initializeFrom(InputRecord &ir)
{

    LatticeStructuralMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);

    //Young's modulus of the material that the network element is made of
    IR_GIVE_FIELD(ir, this->eNormalMean, _IFT_LatticeLinearElastic_eNormal); // Macro

    //Poisson's ratio (Only needed for 3D Bernoulli beam)
    this->nu = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->nu, _IFT_LatticeLinearElastic_nu); // Macro


    //Parameter which relates the shear stiffness to the normal stiffness. Default is 1
    alphaOne = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaOne, _IFT_LatticeLinearElastic_alphaOne); // Macro

    //Parameter which is used for the definition of bending stiffness. Default is 1.
    alphaTwo = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaTwo, _IFT_LatticeLinearElastic_alphaTwo); // Macro

    localRandomType = 0; //Default: No local random field
    IR_GIVE_OPTIONAL_FIELD(ir, localRandomType, _IFT_LatticeLinearElastic_localrandomtype); // Macro
    if ( localRandomType == 1 ) { //Gaussian random generator
        coefficientOfVariation = 0.;
        IR_GIVE_FIELD(ir, coefficientOfVariation, _IFT_LatticeLinearElastic_coefficientOfVariation); // Macro
    }


    double value = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_LatticeLinearElastic_talpha);
    if ( !propertyDictionary.includes(tAlpha) ) {
        //    if (alpha > 0.0 && !propertyDictionary.includes(tAlpha)) {
        // put isotropic thermal expansion coeff into dictionary, if provided
        // and not previosly defined
        propertyDictionary.add(tAlpha, value);
    }

    this->cAlpha = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, cAlpha, _IFT_LatticeLinearElastic_calpha);

}

MaterialStatus *
LatticeLinearElastic :: CreateStatus(GaussPoint *gp) const
{
    return new LatticeLinearElasticMaterialStatus(gp);
}

MaterialStatus *
LatticeLinearElastic :: giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status != NULL ) {
            gp->setMaterialStatus(status);
            this->_generateStatusVariables(gp);
        }
    }

    return status;
}


FloatArrayF< 3 >
LatticeLinearElastic :: giveLatticeStress2d(const FloatArrayF< 3 > &strain,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{
    FloatArray answer;
    answer.resize(3);
    answer.zero();

    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

    FloatArray reducedStrain;

    this->initTempStatus(gp);

    // subtract stress independent part
    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, strain, tStep, VM_Total);

    FloatMatrix stiffnessMatrix;
    stiffnessMatrix.resize(3, 3);
    stiffnessMatrix.zero();
    stiffnessMatrix = this->give2dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    for ( int i = 1; i <= 3; i++ ) { // only diagonal terms matter
        answer.at(i) = stiffnessMatrix.at(i, i) * reducedStrain.at(i);
    }

    //Read in fluid pressures from structural element if this is not a slave problem
    FloatArray pressures;
    if ( !domain->giveEngngModel()->giveMasterEngngModel() ) {
        static_cast< LatticeStructuralElement * >( gp->giveElement() )->givePressures(pressures);
    }

    double waterPressure = 0.;
    for ( int i = 0; i < pressures.giveSize(); i++ ) {
        waterPressure += 1. / pressures.giveSize() * pressures.at(i + 1);
    }

    answer.at(1) += waterPressure;

    //Set all temp values
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);

    return answer;
}


FloatArrayF< 6 >
LatticeLinearElastic :: giveLatticeStress3d(const FloatArrayF< 6 > &strain,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{
    FloatArray answer;
    answer.resize(6);
    answer.zero();

    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

    FloatArray reducedStrain;

    this->initTempStatus(gp);

    // subtract stress independent part
    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, strain, tStep, VM_Total);

    FloatMatrix stiffnessMatrix;
    stiffnessMatrix.resize(6, 6);
    stiffnessMatrix.zero();
    stiffnessMatrix = this->give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    for ( int i = 1; i <= 6; i++ ) { // only diagonal terms matter
        answer.at(i) = stiffnessMatrix.at(i, i) * reducedStrain.at(i);
    }

    //Read in fluid pressures from structural element if this is not a slave problem
    FloatArray pressures;
    if ( !domain->giveEngngModel()->giveMasterEngngModel() ) {
        static_cast< LatticeStructuralElement * >( gp->giveElement() )->givePressures(pressures);
    }

    double waterPressure = 0.;
    for ( int i = 0; i < pressures.giveSize(); i++ ) {
        waterPressure += 1. / pressures.giveSize() * pressures.at(i + 1);
    }

    answer.at(1) += waterPressure;

    //Set all temp values
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);

    return answer;
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

    return;
}


Interface *
LatticeLinearElastic :: giveInterface(InterfaceType type)
{
    return NULL;
}

FloatMatrixF< 1, 1 >
LatticeLinearElastic :: give1dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    /* Returns elastic moduli in reduced stress-strain space*/
    FloatMatrix answer;
    answer.resize(1, 1);
    answer.zero();

    //  answer.at(1, 1) = this->give(eNormal_ID, gp) * this->eNormalMean;
    //  answer.at(1, 1) = this->give(eNormal_ID, gp) * this->eNormalMean;

    return answer;
}

FloatMatrixF< 3, 3 >
LatticeLinearElastic :: give2dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    /* Returns elastic moduli in reduced stress-strain space*/
    FloatMatrix answer;
    answer.resize(3, 3);
    answer.zero();

    answer.at(1, 1) = 1.;
    answer.at(2, 2) = this->alphaOne; // shear
    answer.at(3, 3) = this->alphaTwo; // torsion

    answer.times(this->give(eNormal_ID, gp) * this->eNormalMean);

    return answer;
}


FloatMatrixF< 6, 6 >
LatticeLinearElastic :: give3dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime) const
{
    /* Returns elastic moduli in reduced stress-strain space*/
    FloatMatrix answer;
    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) = 1.;
    answer.at(2, 2) = this->alphaOne; // shear
    answer.at(3, 3) = this->alphaOne; // shear
    answer.at(4, 4) = this->alphaTwo; // torsion
    answer.at(5, 5) = this->alphaTwo; // torsion
    answer.at(6, 6) = this->alphaTwo; // torsion

    answer.times(this->give(eNormal_ID, gp) * this->eNormalMean);

    return answer;
}



void
LatticeLinearElastic :: giveThermalDilatationVector(FloatArray &answer,
                                                    GaussPoint *gp,  TimeStep *tStep)
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

    answer.resize(6);
    answer.zero();

    answer.at(1) = alpha;
}


double
LatticeLinearElastic :: give(int aProperty, GaussPoint *gp) const
{
    double answer;
    if ( static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) )->_giveProperty(aProperty, answer) ) {
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
    } else if ( aProperty == 'n' ) {
        return this->nu;
    } else {
        return LatticeStructuralMaterial :: give(aProperty, gp);
    }
}


int
LatticeLinearElastic :: giveIPValue(FloatArray &answer,
                                    GaussPoint *gp,
                                    InternalStateType type,
                                    TimeStep *atTime)
{
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

    if ( type == IST_CharacteristicLength ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
        return 1;
    } else if ( type == IST_DeltaDissWork ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDeltaDissipation();
        return 1;
    } else if ( type == IST_CrackWidth ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveCrackWidth();
        return 1;
    } else if ( type == IST_DissWork ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDissipation();
        return 1;
    }  else if ( type == IST_CrackStatuses ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveCrackFlag();
        return 1;
    } else {
        return LatticeStructuralMaterial :: giveIPValue(answer, gp, type, atTime);
    }
}


//Status

LatticeLinearElasticMaterialStatus :: LatticeLinearElasticMaterialStatus(GaussPoint *g) : LatticeMaterialStatus(g)
{}



void
LatticeLinearElasticMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    FloatArray forces;

    if ( strcmp(this->gp->giveElement()->giveClassName(), "LatticeBeam3d") == 0 ) {
        static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveInternalForcesVector(forces, tStep, 0);
        fprintf(file, "LatticeBeam forces = %e %e %e %e %e %e.\n", forces.at(7), forces.at(8), forces.at(9), forces.at(10), forces.at(11), forces.at(12) );
    } else if ( strcmp(this->gp->giveElement()->giveClassName(), "LatticeBeam3dBoundary") == 0 ) {
        static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveInternalForcesVector(forces, tStep, 0);
        fprintf(file, "LatticeBeam3dBounday forces = %e %e %e %e %e %e.\n", forces.at(7), forces.at(8), forces.at(9), forces.at(10), forces.at(11), forces.at(12) );
    } else {
        LatticeMaterialStatus :: printOutputAt(file, tStep);
    }
    return;
}
}
