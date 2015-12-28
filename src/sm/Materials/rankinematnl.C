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

#include "rankinematnl.h"
#include "../sm/Elements/structuralelement.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "sparsemtrx.h"
#include "error.h"
#include "nonlocalmaterialext.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "datastream.h"

namespace oofem {
REGISTER_Material(RankineMatNl);

RankineMatNl :: RankineMatNl(int n, Domain *d) : RankineMat(n, d), StructuralNonlocalMaterialExtensionInterface(d), NonlocalMaterialStiffnessInterface()
    //
    // constructor
    //
{ }

void
RankineMatNl :: giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp,
                                                 const FloatArray &totalStrain, TimeStep *tStep)
{
    RankineMatNlStatus *nlStatus = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );
    
    double tempDam;
    //mj performPlasticityReturn(gp, totalStrain, mode);
    // nonlocal method "computeDamage" performs the plastic return
    // for all Gauss points when it is called for the first time
    // in the iteration
    tempDam = this->computeDamage(gp, tStep);
    answer.beScaled(1.0 - tempDam, nlStatus->giveTempEffectiveStress());
    nlStatus->setTempDamage(tempDam);
    nlStatus->letTempStrainVectorBe(totalStrain);
    nlStatus->letTempStressVectorBe(answer);
#ifdef keep_track_of_dissipated_energy
    double gf = sig0 * sig0 / E; // only estimated, but OK for this purpose
    nlStatus->computeWork_PlaneStress(gp, gf);
#endif
}

void
RankineMatNl :: giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp,
                                        const FloatArray &totalStrain, TimeStep *tStep)
{
    RankineMatNlStatus *nlStatus = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );

    double tempDam;
    //mj performPlasticityReturn(gp, totalStrain, mode);
    // nonlocal method "computeDamage" performs the plastic return
    // for all Gauss points when it is called for the first time
    // in the iteration
    tempDam = this->computeDamage(gp, tStep);
    answer.beScaled(1.0 - tempDam, nlStatus->giveTempEffectiveStress());
    nlStatus->setTempDamage(tempDam);
    nlStatus->letTempStrainVectorBe(totalStrain);
    nlStatus->letTempStressVectorBe(answer);
#ifdef keep_track_of_dissipated_energy
    double gf = sig0 * sig0 / E; // only estimated, but OK for this purpose
    nlStatus->computeWork_1d(gp, gf);
#endif
}


void
RankineMatNl :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
        return;
    }

    RankineMatNlStatus *status = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );

    if ( mode == SecantStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
        double damage = status->giveTempDamage();
        answer.times(1. - damage);
        return;
    }

    if ( mode == TangentStiffness ) {
        double tempDamage = status->giveTempDamage();
        double damage = status->giveDamage();
        double gprime;
        if ( tempDamage <= damage ) { // unloading
            gprime = 0.;
        } else { // loading
            double kappa;
            computeCumPlasticStrain(kappa, gp, tStep);
            gprime = computeDamageParamPrime(kappa);
            gprime *= ( 1. - mm );
        }

        evaluatePlaneStressStiffMtrx(answer, mode, gp, tStep, gprime);
        return;
    }

    OOFEM_ERROR("unknown type of stiffness");
}

void
RankineMatNl :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep)
{
    /* Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     */

    double cumPlasticStrain;
    RankineMatNlStatus *nlstatus = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);
    this->performPlasticityReturn(gp, strainVector);
    this->computeLocalCumPlasticStrain(cumPlasticStrain, gp, tStep);
    // standard formulation based on averaging of equivalent strain
    nlstatus->setLocalCumPlasticStrainForAverage(cumPlasticStrain);
}

// returns in "kappa" the value of kappa_hat = m*kappa_nonlocal + (1-m)*kappa_local
void
RankineMatNl :: computeCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *tStep)
{
    double nonlocalContribution, nonlocalCumPlasticStrain = 0.0;
    RankineMatNlStatus *nonlocStatus, *status = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(tStep);
    double localCumPlasticStrain = status->giveLocalCumPlasticStrainForAverage();
    // compute nonlocal cumulative plastic strain
    std :: list< localIntegrationRecord > *list = this->giveIPIntegrationList(gp);

    for ( auto &lir: *list ) {
        nonlocStatus = static_cast< RankineMatNlStatus * >( this->giveStatus(lir.nearGp) );
        nonlocalContribution = nonlocStatus->giveLocalCumPlasticStrainForAverage();
        if ( nonlocalContribution > 0 ) {
            nonlocalContribution *= lir.weight;
        }

        nonlocalCumPlasticStrain += nonlocalContribution;
    }

    double scale = status->giveIntegrationScale();
    if ( scaling == ST_Standard ) { // standard rescaling
        nonlocalCumPlasticStrain *= 1. / scale;
    } else if ( scaling == ST_Borino ) { // Borino modification
        if ( scale > 1. ) {
            nonlocalCumPlasticStrain *= 1. / scale;
        } else {
            nonlocalCumPlasticStrain += ( 1. - scale ) * status->giveLocalCumPlasticStrainForAverage();
        }
    }

    kappa = mm * nonlocalCumPlasticStrain + ( 1. - mm ) * localCumPlasticStrain;
    status->setKappa_nl(nonlocalCumPlasticStrain);
    status->setKappa_hat(kappa);
}

Interface *
RankineMatNl :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialExtensionInterface * >(this);
    } else if ( type == NonlocalMaterialStiffnessInterfaceType ) {
        return static_cast< NonlocalMaterialStiffnessInterface * >(this);
    } else {
        return NULL;
    }
}


IRResultType
RankineMatNl :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    return RankineMat :: initializeFrom(ir);
}


void
RankineMatNl :: giveInputRecord(DynamicInputRecord &input)
{
    RankineMat :: giveInputRecord(input);
    StructuralNonlocalMaterialExtensionInterface :: giveInputRecord(input);
}



double
RankineMatNl :: computeDamage(GaussPoint *gp, TimeStep *tStep)
{
    RankineMatNlStatus *nlStatus = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );
    double nlKappa;
    this->computeCumPlasticStrain(nlKappa, gp, tStep);
    double dam, tempDam;
    dam = nlStatus->giveDamage();
    tempDam = this->computeDamageParam(nlKappa);
    if ( tempDam < dam ) {
        tempDam = dam;
    }

    return tempDam;
}

// this method is running over all neighboring Gauss points
// and computes the contribution of the nonlocal interaction to tangent stiffness
// it uses methods
//   giveLocalNonlocalStiffnessContribution
//   giveRemoteNonlocalStiffnessContribution
// the first one computes m*gprime*Btransposed*sigmaeff for the present point
// the second one computes Btransposed*eta for the neighboring point
// (where eta is the derivative of cum. plastic strain wrt final strain)
// THIS METHOD CAN BE USED IN THE SAME FORM FOR ALL NONLOCAL DAMAGE-PLASTIC MODELS
void
RankineMatNl :: NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                                     GaussPoint *gp, TimeStep *tStep)
{
    double coeff;
    RankineMatNlStatus *status = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );
    std :: list< localIntegrationRecord > *list = status->giveIntegrationDomainList();
    RankineMatNl *rmat;
    FloatArray rcontrib, lcontrib;
    IntArray loc, rloc;

    FloatMatrix contrib;

    if ( this->giveLocalNonlocalStiffnessContribution(gp, loc, s, lcontrib, tStep) == 0 ) {
        return;
    }

    for ( auto &lir: *list ) {
        rmat = dynamic_cast< RankineMatNl * >( lir.nearGp->giveMaterial() );
        if ( rmat ) {
            rmat->giveRemoteNonlocalStiffnessContribution(lir.nearGp, rloc, s, rcontrib, tStep);
            coeff = gp->giveElement()->computeVolumeAround(gp) * lir.weight / status->giveIntegrationScale();

            contrib.clear();
            contrib.plusDyadUnsym(lcontrib, rcontrib, - 1.0 * coeff);
            dest.assemble(loc, rloc, contrib);
        }
    }
}

std :: list< localIntegrationRecord > *
RankineMatNl :: NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp)
{
    RankineMatNlStatus *status = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );
    this->buildNonlocalPointTable(gp);
    return status->giveIntegrationDomainList();
}


// computes m*gprime*Btransposed*sigmaeff for the given Gauss point
// and returns 0 of damage is not growing, 1 if it is growing
// (if damage is not growing, the contribution is not considered at all)
int
RankineMatNl :: giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                                       FloatArray &lcontrib, TimeStep *tStep)
{
    int nrows, nsize;
    double sum, nlKappa, damage, tempDamage;
    RankineMatNlStatus *status = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );
    StructuralElement *elem = static_cast< StructuralElement * >( gp->giveElement() );
    FloatMatrix b;

    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    if ( tempDamage <= damage ) {
        return 0; // no contribution if damage is not growing
    }

    elem->giveLocationArray(loc, s);
    const FloatArray &stress = status->giveTempEffectiveStress();
    elem->computeBmatrixAt(gp, b);
    this->computeCumPlasticStrain(nlKappa, gp, tStep);
    double factor = computeDamageParamPrime(nlKappa);
    factor *= mm; // this factor is m*gprime
    nrows = b.giveNumberOfColumns();
    nsize = stress.giveSize();
    lcontrib.resize(nrows);
    // compute the product Btransposed*stress and multiply by factor
    for ( int i = 1; i <= nrows; i++ ) {
        sum = 0.0;
        for ( int j = 1; j <= nsize; j++ ) {
            sum += b.at(j, i) * stress.at(j);
        }

        lcontrib.at(i) = sum * factor;
    }

    return 1; // contribution will be considered
}

// computes Btransposed*eta for the given Gauss point
// (where eta is the derivative of cum. plastic strain wrt final strain)
void
RankineMatNl :: giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                        FloatArray &rcontrib, TimeStep *tStep)
{
    RankineMatNlStatus *status = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );
    StructuralElement *elem = static_cast< StructuralElement * >( gp->giveElement() );
    elem->giveLocationArray(rloc, s);
    FloatMatrix b;
    elem->computeBmatrixAt(gp, b);
    int ncols = b.giveNumberOfColumns();
    rcontrib.resize(ncols);

    double kappa = status->giveCumulativePlasticStrain();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    if ( tempKappa <= kappa ) {
        rcontrib.zero();
        return;
    }

    double sum;
    int nsize = 3;
    FloatArray eta(3);
    computeEta(eta, status);
    for ( int i = 1; i <= ncols; i++ ) {
        sum = 0.;
        for ( int j = 1; j <= nsize; j++ ) {
            sum += eta.at(j) * b.at(j, i);
        }

        rcontrib.at(i) = sum;
    }
}


int
RankineMatNl :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_CumPlasticStrain_2 ) {
        answer.resize(1);
        double dummy;
        // this method also stores the nonlocal kappa in status ... kappa_nl
        computeCumPlasticStrain(dummy, gp, tStep);
        RankineMatNlStatus *status = static_cast< RankineMatNlStatus * >( this->giveStatus(gp) );
        answer.at(1) = status->giveKappa_nl();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        computeCumPlasticStrain(answer.at(1), gp, tStep);
        return 1;
    } else {
        return RankineMat :: giveIPValue(answer, gp, type, tStep);
    }
}


//*******************************
//*************status************
//*******************************

RankineMatNlStatus :: RankineMatNlStatus(int n, Domain *d, GaussPoint *g) :
    RankineMatStatus(n, d, g), StructuralNonlocalMaterialStatusExtensionInterface()
{
    localCumPlasticStrainForAverage = 0.0;
}


RankineMatNlStatus :: ~RankineMatNlStatus()
{ }


void
RankineMatNlStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "damage %g, kappa %g, kappa_nl %g, kappa_hat %g", damage, kappa, kappa_nl, kappa_hat);
#ifdef keep_track_of_dissipated_energy
    fprintf(file, ", dissW %g, freeE %g, stressW %g", this->dissWork, ( this->stressWork ) - ( this->dissWork ), this->stressWork);
#endif
    fprintf(file, " }\n");
}

void
RankineMatNlStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    RankineMatStatus :: initTempStatus();
}


void
RankineMatNlStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    RankineMatStatus :: updateYourself(tStep);
}


contextIOResultType
RankineMatNlStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = RankineMatStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //if (!stream.write(&localEquivalentStrainForAverage,1)) THROW_CIOERR(CIO_IOERR);
    return CIO_OK;
}

contextIOResultType
RankineMatNlStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = RankineMatStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    //if (!stream.read (&localEquivalentStrainForAverage,1)) THROW_CIOERR(CIO_IOERR);

    return CIO_OK;
}

Interface *
RankineMatNlStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialStatusExtensionInterface * >(this);
    } else {
        return RankineMatStatus :: giveInterface(type);
    }
}


int
RankineMatNl :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip)
{
    RankineMatNlStatus *nlStatus = static_cast< RankineMatNlStatus * >( this->giveStatus(ip) );

    this->buildNonlocalPointTable(ip);
    this->updateDomainBeforeNonlocAverage(tStep);

    return buff.write( nlStatus->giveLocalCumPlasticStrainForAverage() );
}

int
RankineMatNl :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip)
{
    int result;
    RankineMatNlStatus *nlStatus = static_cast< RankineMatNlStatus * >( this->giveStatus(ip) );
    double localCumPlasticStrainForAverage;

    result = buff.read(localCumPlasticStrainForAverage);
    nlStatus->setLocalCumPlasticStrainForAverage(localCumPlasticStrainForAverage);
    return result;
}

int
RankineMatNl :: estimatePackSize(DataStream &buff, GaussPoint *ip)
{
    // Note: nlStatus localStrainVectorForAverage memeber must be properly sized!
    // IDNLMaterialStatus *nlStatus = (IDNLMaterialStatus*) this -> giveStatus (ip);
    return buff.givePackSizeOfDouble(1);
}

} // end namespace oofem
