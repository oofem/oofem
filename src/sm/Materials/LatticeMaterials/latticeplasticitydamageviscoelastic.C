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

#include "latticeplasticitydamageviscoelastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(LatticePlasticityDamageViscoelastic);

LatticePlasticityDamageViscoelastic :: LatticePlasticityDamageViscoelastic(int n, Domain *d) : LatticePlasticityDamage(n, d)
{}


bool
LatticePlasticityDamageViscoelastic :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _3dLattice;
}


void
LatticePlasticityDamageViscoelastic :: initializeFrom(InputRecord &ir)
{
    LatticePlasticityDamage :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, slaveMat, _IFT_LatticePlasticityDamageViscoelastic_slaveMat); // number of slave material

    RheoChainMaterial *rChM = giveViscoelasticMaterial();

    // override elastic modulus by the one given by the compliance function
    // provided values of the strain at peak stress must be changed

    double E28 = 1. / rChM->computeCreepFunction(28.01, 28., NULL, NULL); // modulus of elasticity evaluated at 28 days, duration of loading 15 min
    this->eNormalMean = E28; // swap elastic modulus/stiffness


    // must match with the viscoelastic material
    if ( rChM->giveAlphaOne() != this->alphaOne ) {
        OOFEM_ERROR("a1 must be set to the same value in both master and viscoelastic slave materials");
    }

    // must match with the viscoelastic material
    if ( rChM->giveAlphaTwo() != this->alphaTwo ) {
        OOFEM_ERROR("a2 must be set to the same value in both master and viscoelastic slave materials");
    }
}


MaterialStatus *
LatticePlasticityDamageViscoelastic :: CreateStatus(GaussPoint *gp) const
{
    return new LatticePlasticityDamageViscoelasticStatus(1, LatticePlasticityDamageViscoelastic :: domain, gp, slaveMat);
}

FloatArrayF<6>
LatticePlasticityDamageViscoelastic :: giveReducedStrain(GaussPoint *gp,
                                                         TimeStep *tStep) const
{
    auto status = static_cast< LatticePlasticityDamageStatus * >( this->giveStatus(gp) );

    auto elasticStiffnessMatrix = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    double omega = status->giveDamage();
    auto strain = status->giveLatticeStress() * (1. / ( 1. - omega ) );


    for ( int i = 1; i <= 6; i++ ) { // only diagonal terms matter
        strain.at(i) /= elasticStiffnessMatrix.at(i, i);
    }

    const auto &oldPlasticStrain = status->givePlasticLatticeStrain();
    for ( int i = 1; i <= 3; i++ ) {
        strain.at(i) += oldPlasticStrain.at(i);
    }
    return strain;
}

void
LatticePlasticityDamageViscoelastic :: giveRealStressVector(FloatArray &answer,
                                                            GaussPoint *gp,
                                                            const FloatArray &totalStrain,
                                                            TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    double tol = 1.e-12;

    auto status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );

    //    this->initGpForNewStep(gp);

    // total deformation - temperature - plastic strain
    FloatArray reducedStrainForViscoMat;

    // VE stress / stiffness
    FloatArray quasiElasticStrain;

    // VE stress / stiffness + plastic strain
    FloatArray reducedStrain;

    RheoChainMaterial *rChM = giveViscoelasticMaterial();

    MaterialMode mMode = gp->giveMaterialMode();
    StressVector tempEffectiveStress(mMode);

    GaussPoint *rChGP = status->giveViscoelasticGaussPoint();

    FloatArrayF<6> tempPlasticStrain;

    auto elasticStiffnessMatrix = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    FloatArrayF<3> strain;

    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector(gp->giveMaterialMode() );

    int itercount = 0;
    double prevIterKappaP = 0.;

    FloatArrayF<3> stress3;
    while ( ( fabs(status->giveTempKappaP() - prevIterKappaP) >= tol ) || ( itercount == 0 ) ) {
        itercount++;

        if ( itercount > 100 ) {
            OOFEM_ERROR("Algorithm not converging");
        }

        prevIterKappaP = status->giveTempKappaP();

        // subtract stress independent part = temperature
        // slave (viscoelastic) material uses incremental formulation
        // must subtract not only temperature but also increment of plastic deformation
        this->giveStressDependentPartOfStrainVector(reducedStrainForViscoMat, gp, totalStrain, tStep, VM_Total);

        tempPlasticStrain = status->giveTempPlasticLatticeStrain();

        for ( int i = 1; i <= 3; i++ ) {
            reducedStrainForViscoMat.at(i) -= tempPlasticStrain.at(i);
        }

        rChM->giveRealStressVector(tempEffectiveStress, rChGP, reducedStrainForViscoMat, tStep);

        // transform effective stress computed by the viscoelastic model into strain by
        // dividing it by a stiffness matrix (which is diagonal)
        quasiElasticStrain.resize(rsize);
        for ( int i = 1; i <= rsize; i++ ) { // only diagonal terms matter
            quasiElasticStrain.at(i) = tempEffectiveStress.at(i) / elasticStiffnessMatrix.at(i, i);
        }

        // total deformation is needed = visco stress / stiffness + plastic strain
        for ( int i = 1; i <= 3; i++ ) {
            strain.at(i) = quasiElasticStrain.at(i) + tempPlasticStrain.at(i);
        }

        // this (parent) method needs a total, temperature-free strain (with plastic strain)
        // from now on new = updated plastic strain
        stress3 = this->performPlasticityReturn(gp, strain, tStep);
    } // end while loop

    reducedStrain = quasiElasticStrain;
    for ( int i = 1; i <= 3; i++ ) {
        reducedStrain.at(i) += tempPlasticStrain.at(i);
    }

    // performDamageEvolution input "strain" is used only to compute crack width. Therefore
    // it should correspond to stress-dependent strain deformation = total - shrinkage - temperature
    // but with plastic deformation

    FloatArray shrinkageStrain;
    rChM->giveShrinkageStrainVector(shrinkageStrain, rChGP, tStep, VM_Total);

    strain = reducedStrainForViscoMat - FloatArrayF<3>(shrinkageStrain) + status->giveTempPlasticLatticeStrain()[{0, 1, 2}];

    // strain is passed only to compute crack width, nothing else
    this->performDamageEvaluation(gp, strain);

    double omega = status->giveTempDamage();

    FloatArrayF<6> stress;
    stress.at(1) = ( 1. - omega ) * stress3.at(1);
    stress.at(2) = ( 1. - omega ) * stress3.at(2);
    stress.at(3) = ( 1. - omega ) * stress3.at(3);
    stress.at(4) = ( 1. - omega ) * quasiElasticStrain.at(4) * this->alphaTwo * this->eNormalMean;
    stress.at(5) = ( 1. - omega ) * quasiElasticStrain.at(5) * this->alphaTwo * this->eNormalMean;
    stress.at(6) = ( 1. - omega ) * quasiElasticStrain.at(6) * this->alphaTwo * this->eNormalMean;

    status->letTempLatticeStrainBe(totalStrain);
    status->letTempReducedLatticeStrainBe(reducedStrain);
    status->letTempLatticeStressBe(stress);
}

RheoChainMaterial *
LatticePlasticityDamageViscoelastic :: giveViscoelasticMaterial()
{
    auto mat = domain->giveMaterial(slaveMat);
    auto rChMat = dynamic_cast< RheoChainMaterial * >( mat );
    return rChMat;
}

FloatMatrixF< 3, 3 >
LatticePlasticityDamageViscoelastic :: give2dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const
{
    /* Returns elastic moduli in reduced stress-strain space*/
    OOFEM_ERROR("2d mode not implemented");
}

FloatMatrixF< 6, 6 >
LatticePlasticityDamageViscoelastic :: give3dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const
{
    FloatMatrixF<6,6> tangent;

    //@todo: this has to be rewritten for the 3d case
    // LatticePlasticityDamageViscoelasticStatus *status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );
    // GaussPoint *slaveGp;
    // RheoChainMaterial *rChMat;
    // double Eincr;

    // slaveGp = status->giveViscoelasticGaussPoint();
    // //    rChMat = giveViscoelasticMaterial();

    // answer = LatticePlasticityDamage :: give3dLatticeStiffnessMatrix(rmode, gp, tStep);

    // Eincr = rChMat->giveEModulus(slaveGp, tStep);

    // answer.times(Eincr / this->eNormalMean);

    return tangent;
}

int
LatticePlasticityDamageViscoelastic :: giveIPValue(FloatArray &answer,
                                                   GaussPoint *gp,
                                                   InternalStateType type,
                                                   TimeStep *tStep)
{
    return LatticePlasticityDamage :: giveIPValue(answer, gp, type, tStep);
}


LatticePlasticityDamageViscoelasticStatus :: LatticePlasticityDamageViscoelasticStatus(int n, Domain *d, GaussPoint *g, int s) :
    LatticePlasticityDamageStatus(n, d, g), slaveMat(s)
{
    viscoelasticGP = new GaussPoint(g->giveIntegrationRule(), g->giveNumber(), g->giveNaturalCoordinates(), g->giveWeight(), g->giveMaterialMode() );
}

void
LatticePlasticityDamageViscoelasticStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticePlasticityDamageStatus :: initTempStatus();
}

void
LatticePlasticityDamageViscoelasticStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    MaterialStatus *mS = this->giveViscoelasticMatStatus();

    LatticePlasticityDamageStatus :: printOutputAt(file, tStep);
    fprintf(file, "\nViscoelastic material:");

    mS->printOutputAt(file, tStep);

    fprintf(file, "\n");
}

MaterialStatus *
LatticePlasticityDamageViscoelasticStatus :: giveViscoelasticMatStatus() const
{
    //    Material *mat;
    //    RheoChainMaterial *rChMat;
    //    GaussPoint *rChGP;

    // @todo: Don't now how to fix this!
    //    mat = domain->Givematerial(slaveMat);
    //    rChMat = dynamic_cast< RheoChainMaterial * >(mat);

    //    rChGP = this->giveViscoelasticGaussPoint();

    //  MaterialStatus *mS = rChMat->giveStatus(rChGP);

    return nullptr;
}


void
LatticePlasticityDamageViscoelasticStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    MaterialStatus *mS = this->giveViscoelasticMatStatus();
    mS->updateYourself(tStep);

    LatticePlasticityDamageStatus :: updateYourself(tStep);
}

void
LatticePlasticityDamageViscoelasticStatus :: saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    LatticePlasticityDamageStatus :: saveContext(stream, mode);

    this->giveViscoelasticMatStatus()->saveContext(stream, mode);
}

void
LatticePlasticityDamageViscoelasticStatus :: restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    LatticePlasticityDamageStatus :: restoreContext(stream, mode);

    this->giveViscoelasticMatStatus()->restoreContext(stream, mode);
}
}     // end namespace oofem
