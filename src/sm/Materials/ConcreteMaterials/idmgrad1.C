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

#include "idmgrad1.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "Materials/isolinearelasticmaterial.h"
#include "error.h"
#include "nonlocalmaterialext.h"
#include "datastream.h"
#include "contextioerr.h"
#include "stressvector.h"
#include "strainvector.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(IDGMaterial);

IDGMaterial :: IDGMaterial(int n, Domain *d) : IsotropicDamageMaterial1(n, d), GradDpMaterialExtensionInterface(d)
    //
    // constructor
    //
{
    internalLength = 0;
    averType = 0;
}


IDGMaterial :: ~IDGMaterial()
//
// destructor
//
{ }

IRResultType
IDGMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    result = IsotropicDamageMaterial1 :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    result = GradDpMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    return this->mapper.initializeFrom(ir);
}

/////////////////////////////////////////////////////////////////////////////



int
IDGMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _1dMat || mode == _PlaneStress || mode == _PlaneStrain;
}

void
IDGMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
                                   MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}


void
IDGMaterial :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,  TimeStep *tStep)
{
    IsotropicDamageMaterialStatus *status = static_cast< IsotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    double om;
    om = status->giveTempDamage();
    om = min(om, maxOmega);
    answer.resize(1, 1);
    answer.at(1, 1) = lmat->give('E', gp);
    answer.times(1.0 - om);
}

void
IDGMaterial :: give1dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );

    double kappa = status->giveKappa();
    double tempKappa = status->giveTempKappa();
    FloatArray totalStrain =  status->giveTempStrainVector();

    answer.resize(1, 1);
    if ( tempKappa > kappa ) {
        FloatArray eta;
        this->computeEta(eta, totalStrain, gp, tStep);
        answer.at(1, 1) = eta.at(1);
    }
}

void
IDGMaterial :: give1dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();

    double damage = status->giveDamage();
    double tempDamage = status->giveTempDamage();
    double E = lmat->give('E', gp);

    answer.resize(1, 1);
    if ( tempDamage > damage ) {
        double nlKappa =  status->giveNonlocalCumulatedStrain();
        answer.at(1, 1) = E * status->giveTempStrainVector().at(1);
        double gPrime = damageFunctionPrime(nlKappa, gp);
        answer.times(gPrime);
    } else {
        answer.zero();
    }
}



void
IDGMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
    double tempDamage;
    if ( mode == ElasticStiffness ) {
        tempDamage = 0.0;
    } else {
        tempDamage = status->giveTempDamage();
        if ( tempDamage > 0.0 ) {
            tempDamage = min(tempDamage, maxOmega);
        }
    }

    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    answer.times(1.0 - tempDamage);
}


void
IDGMaterial :: givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1, 3);
    answer.zero();
    if ( mode == TangentStiffness ) {
        IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
        FloatArray eta;
        FloatArray totalStrain = status->giveTempStrainVector();
        this->computeEta(eta, totalStrain, gp, tStep);
        answer.at(1, 1) = eta.at(1);
        answer.at(1, 2) = eta.at(2);
        answer.at(1, 3) = eta.at(3);
    }
}

void
IDGMaterial :: givePlaneStressGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(3, 1);
    if ( mode == TangentStiffness ) {
        IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
        double tempDamage = status->giveTempDamage();
        tempDamage = min(tempDamage, maxOmega);
        double nlKappa =  status->giveNonlocalCumulatedStrain();
        double gPrime =  this->damageFunctionPrime(nlKappa, gp);
        FloatArray stress =  status->giveTempStressVector();

        answer.at(1, 1) = stress.at(1) / ( 1 - tempDamage );
        answer.at(2, 1) = stress.at(2) / ( 1 - tempDamage );
        answer.at(3, 1) = stress.at(3) / ( 1 - tempDamage );
        answer.times(gPrime);
    }
}


void
IDGMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
    double tempDamage;
    if ( mode == ElasticStiffness ) {
        tempDamage = 0.0;
    } else {
        tempDamage = status->giveTempDamage();
        if ( tempDamage > 0.0 ) {
            tempDamage = min(tempDamage, maxOmega);
        }
    }

    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
    answer.times(1.0 - tempDamage);
}


void
IDGMaterial :: givePlaneStrainKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
    FloatArray totalStrain =  status->giveTempStrainVector();
    FloatArray eta;

    answer.resize(1, 4);
    this->computeEta(eta, totalStrain, gp, tStep);
    answer.at(1, 1) = eta.at(1);
    answer.at(1, 2) = eta.at(2);
    answer.at(1, 3) = eta.at(3);
    answer.at(1, 4) = eta.at(4);
}

void
IDGMaterial :: givePlaneStrainGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
    double damage = status->giveDamage();
    double tempDamage = status->giveTempDamage();
    tempDamage = min(tempDamage, maxOmega);
    answer.resize(4, 1);
    if ( tempDamage > damage ) {
        double nlEquivStrain =  status->giveNonlocalCumulatedStrain();
        double gPrime =  this->damageFunctionPrime(nlEquivStrain, gp);
        FloatArray stress =  status->giveTempStressVector();
        answer.at(1, 1) = stress.at(1) / ( 1 - tempDamage );
        answer.at(2, 1) = stress.at(2) / ( 1 - tempDamage );
        answer.at(3, 1) = stress.at(3) / ( 1 - tempDamage );
        answer.at(4, 1) = stress.at(4) / ( 1 - tempDamage );
        answer.times(gPrime);
    }
}


void
IDGMaterial :: giveInternalLength(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    if ( averType == 0 ) {
        answer.resize(1, 1);
        answer.at(1, 1) = cl;
    } else if ( averType == 1 ) {
        answer.resize(1, 1);
        FloatArray gpCoords;
        if ( gp->giveElement()->computeGlobalCoordinates( gpCoords, gp->giveNaturalCoordinates() ) == 0 ) {
            OOFEM_ERROR("computeGlobalCoordinates of GP failed");
        }

        this->giveDistanceBasedCharacteristicLength(gpCoords);
        answer.at(1, 1) = cl;
    } else if ( averType == 2 ) {
        MaterialMode mMode = gp->giveMaterialMode();
        if ( mMode == _PlaneStress ) {
            answer.resize(2, 2);
            answer.zero();
            IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
            FloatArray stress = status->giveTempStressVector();
            StressVector fullStress(stress, _PlaneStress);
            FloatArray sigPrinc;
            FloatMatrix nPrinc;
            // get principal stresses (ordered) and principal stress directions
            fullStress.computePrincipalValDir(sigPrinc, nPrinc);
            if ( nPrinc.giveNumberOfRows() == 0 ) {
                nPrinc.resize(2, 2);
                nPrinc.at(1, 1) = 1;
                nPrinc.at(2, 2) = 1;
            }

            // principal internal lengths
            double l1 = cl0;
            double l2 = cl0;
            double gamma;
            if ( sigPrinc.at(1) > 0 ) {
                if ( sigPrinc.at(2) > 0 ) {
                    gamma = beta + ( 1 - beta ) * sigPrinc.at(2) * sigPrinc.at(2) / ( sigPrinc.at(1) * sigPrinc.at(1) );
                } else {
                    gamma = beta;
                }
            } else {
                gamma = 1;
            }

            l2 = l2 * gamma;

            // compose the internal Length matrix in global coordinates
            //   the first subscript refers to coordinate
            //   the second subscript refers to eigenvalue
            double n11 = nPrinc.at(1, 1);
            double n12 = nPrinc.at(1, 2);
            double n21 = nPrinc.at(2, 1);
            double n22 = nPrinc.at(2, 2);

            answer.at(1, 1) = l1 * l1 * n11 * n11 + l2 * l2 * n12 * n12;
            answer.at(1, 2) = l1 * l1 * n11 * n21 + l2 * l2 * n12 * n22;
            answer.at(2, 1) = l1 * l1 * n11 * n21 + l2 * l2 * n12 * n22;
            answer.at(2, 2) = l1 * l1 * n21 * n21 + l2 * l2 * n22 * n22;
        } else {
            OOFEM_ERROR("Unknown material mode.");
        }
    }
}


void
IDGMaterial :: giveInternalLengthDerivative(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    if ( mMode == _PlaneStress ) {
        answer.resize(5, 5);
        answer.zero();
        if ( mode == TangentStiffness ) {
            IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
            FloatArray stress = status->giveTempStressVector();
            stress.resize(3);
            StressVector fullStress(stress, _PlaneStress);
            FloatArray sigPrinc;
            FloatMatrix nPrinc;
            // get principal trial stresses (ordered) and principal stress directions
            fullStress.computePrincipalValDir(sigPrinc, nPrinc);
            if ( sigPrinc.at(1) > 0 ) {
                if ( sigPrinc.at(2) > 0 && sigPrinc.at(1) != sigPrinc.at(2) ) {
                    double gamma = beta + ( 1 - beta ) * sigPrinc.at(2) * sigPrinc.at(2) / ( sigPrinc.at(1) * sigPrinc.at(1) );
                    double dL2dS1 =  4. * gamma * cl * cl * ( beta - 1. ) * sigPrinc.at(2) * sigPrinc.at(2) / ( sigPrinc.at(1) * sigPrinc.at(1) * sigPrinc.at(1) );
                    double dL2dS2 =  4. * gamma * cl * cl * ( 1. - beta ) * sigPrinc.at(2) / ( sigPrinc.at(1) * sigPrinc.at(1) );
                    double dLdN = ( gamma * gamma - 1. ) * cl * cl / ( sigPrinc.at(2) - sigPrinc.at(1) );
                    answer.at(1, 1) = nPrinc.at(1, 1);
                    answer.at(1, 2) = nPrinc.at(1, 2);
                    answer.at(2, 1) = nPrinc.at(2, 1);
                    answer.at(2, 2) = nPrinc.at(2, 2);
                    answer.at(3, 3) = dL2dS1;
                    answer.at(4, 4) = dL2dS2;
                    answer.at(5, 5) = dLdN;
                }
            }
        }
    } else {
        OOFEM_ERROR("Unknown material mode.");
    }
}






void
IDGMaterial :: giveRealStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    IDGMaterialStatus *status = static_cast< IDGMaterialStatus * >( this->giveStatus(gp) );
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray reducedTotalStrainVector, strain;

    FloatMatrix de;
    double f, equivStrain, tempKappa = 0.0, omega = 0.0;

    this->initTempStatus(gp);
    
    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, tStep, VM_Total);


    // compute equivalent strain
    this->computeEquivalentStrain(equivStrain, reducedTotalStrainVector, gp, tStep);

    if ( llcriteria == idm_strainLevelCR ) {
        // compute value of loading function if strainLevel crit apply
        f = nonlocalCumulatedStrain - status->giveKappa();

        if ( f <= 0.0 ) {
            // damage does not grow
            tempKappa = status->giveKappa();
            omega     = status->giveDamage();
        } else {
            // damage grow
            this->initDamaged(tempKappa, reducedTotalStrainVector, gp);
            // evaluate damage parameter
            this->computeDamageParam(omega, nonlocalCumulatedStrain, reducedTotalStrainVector, gp);
        }
    } else if ( llcriteria == idm_damageLevelCR ) {
        // evaluate damage parameter first
        tempKappa = nonlocalCumulatedStrain;
        this->initDamaged(tempKappa, strain, gp);
        this->computeDamageParam(omega, tempKappa, reducedTotalStrainVector, gp);
        if ( omega < status->giveDamage() ) {
            // unloading takes place
            omega = status->giveDamage();
            //printf (".");
        }
    } else {
        OOFEM_ERROR("unsupported loading/uloading criteria");
    }


    lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
    de.times(1.0 - omega);
    answer1.beProductOf(de, strain);
    answer2 = equivStrain;

    // update gp

    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer1);
    status->setNonlocalCumulatedStrain(nonlocalCumulatedStrain);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
}




MaterialStatus *
IDGMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new IDGMaterialStatus(1, IDGMaterial :: domain, gp);
}


IDGMaterialStatus :: IDGMaterialStatus(int n, Domain *d, GaussPoint *g) : IsotropicDamageMaterial1Status(n, d, g)
{ }


IDGMaterialStatus :: ~IDGMaterialStatus()
{ }




void
IDGMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    IsotropicDamageMaterial1Status :: initTempStatus();
}



void
IDGMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    IsotropicDamageMaterial1Status :: updateYourself(tStep);
}



contextIOResultType
IDGMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = IsotropicDamageMaterial1Status :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream.write(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
IDGMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = IsotropicDamageMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream.read(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

void
IDGMaterial :: givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMat:
        give1dStressStiffMtrx(answer, mode, gp, tStep);
        break;
    case _PlaneStress:
        givePlaneStressStiffMtrx(answer, mode, gp, tStep);
        break;
    case _PlaneStrain:
        givePlaneStrainStiffMtrx(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
IDGMaterial :: givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMat:
        give1dKappaMatrix(answer, mode, gp, tStep);
        break;
    case _PlaneStress:
        givePlaneStressKappaMatrix(answer, mode, gp, tStep);
        break;
    case _PlaneStrain:
        givePlaneStrainKappaMatrix(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
IDGMaterial :: givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMat:
        give1dGprime(answer, mode, gp, tStep);
        break;
    case _PlaneStress:
        givePlaneStressGprime(answer, mode, gp, tStep);
        break;
    case _PlaneStrain:
        givePlaneStrainGprime(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
IDGMaterial :: givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMat:
        giveInternalLength(answer, mode, gp, tStep);
        break;
    case _PlaneStress:
        giveInternalLength(answer, mode, gp, tStep);
        break;
    case _PlaneStrain:
        giveInternalLength(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

void
IDGMaterial :: givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _PlaneStress:
        giveInternalLengthDerivative(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}
}     // end namespace oofem
