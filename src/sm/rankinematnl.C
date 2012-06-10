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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "rankinematnl.h"
#include "structuralelement.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "sparsemtrx.h"
#include "dynalist.h"
#include "error.h"
#include "nonlocalmaterialext.h"
#include "contextioerr.h"

#ifdef __PARALLEL_MODE
 #include "combuff.h"
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
RankineMatNl :: RankineMatNl(int n, Domain *d) : RankineMat(n, d), StructuralNonlocalMaterialExtensionInterface(d), NonlocalMaterialStiffnessInterface()
//
// constructor
//
{}

void
RankineMatNl :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                     const FloatArray &totalStrain, TimeStep *atTime)
{
    RankineMatNlStatus *nlStatus = ( RankineMatNlStatus * ) this->giveStatus(gp);
    //mj this->initGpForNewStep(gp);

    double tempDam;
    FloatArray tempEffStress, totalStress;
    MaterialMode mode = gp->giveMaterialMode();
    //mj performPlasticityReturn(gp, totalStrain, mode);
    // nonlocal method "computeDamage" performs the plastic return
    // for all Gauss points when it is called for the first time
    // in the iteration
    tempDam = this->computeDamage(gp, atTime);
    nlStatus->giveTempEffectiveStress(tempEffStress);
    answer.beScaled( 1.0 - tempDam, tempEffStress);
    nlStatus->setTempDamage(tempDam);
    nlStatus->letTempStrainVectorBe(totalStrain);
    nlStatus->letTempStressVectorBe(answer);
#ifdef keep_track_of_dissipated_energy
    double gf = sig0 * sig0 / E; // only estimated, but OK for this purpose
    nlStatus->computeWork(gp, mode, gf);
#endif
}

void
RankineMatNl :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
        return;
    }

    RankineMatNlStatus *status = ( RankineMatNlStatus * ) this->giveStatus(gp);

    if ( mode == SecantStiffness ) {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
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
            computeCumPlasticStrain(kappa, gp, atTime);
            gprime = computeDamageParamPrime(kappa);
            gprime *= ( 1. - mm );
        }

        evaluatePlaneStressStiffMtrx(answer, form, mode, gp, atTime, gprime);
        return;
    }

    _error("RankineMatNl :: givePlaneStressStiffMtrx ... unknown type of stiffness\n");
}

void
RankineMatNl :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime)
{
    /* Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     */

    double cumPlasticStrain;
    RankineMatNlStatus *nlstatus = ( RankineMatNlStatus * ) this->giveStatus(gp);

    this->initTempStatus(gp);
    //mj this->initGpForNewStep(gp);
    MaterialMode mode = gp->giveMaterialMode();
    this->performPlasticityReturn(gp, strainVector, mode);
    this->computeLocalCumPlasticStrain(cumPlasticStrain, gp, atTime);
    // standard formulation based on averaging of equivalent strain
    nlstatus->setLocalCumPlasticStrainForAverage(cumPlasticStrain);
}

// returns in "kappa" the value of kappa_hat = m*kappa_nonlocal + (1-m)*kappa_local
void
RankineMatNl :: computeCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *atTime)
{
    double nonlocalContribution, nonlocalCumPlasticStrain = 0.0;
    RankineMatNlStatus *nonlocStatus, *status = ( RankineMatNlStatus * ) this->giveStatus(gp);

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(atTime);
    double localCumPlasticStrain = status->giveLocalCumPlasticStrainForAverage();
    // compute nonlocal cumulative plastic strain
    dynaList< localIntegrationRecord > *list = this->giveIPIntegrationList(gp);
    dynaList< localIntegrationRecord > :: iterator pos;

    for ( pos = list->begin(); pos != list->end(); ++pos ) {
        nonlocStatus = ( RankineMatNlStatus * ) this->giveStatus( ( * pos ).nearGp );
        nonlocalContribution = nonlocStatus->giveLocalCumPlasticStrainForAverage();
        if ( nonlocalContribution > 0 ) {
            nonlocalContribution *= ( * pos ).weight;
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
        return ( StructuralNonlocalMaterialExtensionInterface * ) this;
    } else if ( type == NonlocalMaterialStiffnessInterfaceType ) {
        return ( NonlocalMaterialStiffnessInterface * ) this;
    } else {
        return NULL;
    }
}


IRResultType
RankineMatNl :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                // Required by IR_GIVE_FIELD macro

    RankineMat :: initializeFrom(ir);
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);

    return IRRT_OK;
}


int
RankineMatNl :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    RankineMat :: giveInputRecordString(str, keyword);
    StructuralNonlocalMaterialExtensionInterface :: giveInputRecordString(str, false);
    sprintf(buff, " r %e", this->cl);
    str += buff;

    return 1;
}


double
RankineMatNl :: computeDamage(GaussPoint *gp, TimeStep *atTime)
{
    RankineMatNlStatus *nlStatus = ( RankineMatNlStatus * ) this->giveStatus(gp);
    double nlKappa;
    this->computeCumPlasticStrain(nlKappa, gp, atTime);
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
                                                                     GaussPoint *gp, TimeStep *atTime)
{
    double coeff;
    RankineMatNlStatus *status = ( RankineMatNlStatus * ) this->giveStatus(gp);
    dynaList< localIntegrationRecord > *list = status->giveIntegrationDomainList();
    dynaList< localIntegrationRecord > :: iterator pos;
    RankineMatNl *rmat;
    FloatArray rcontrib, lcontrib;
    IntArray loc, rloc;

    FloatMatrix contrib;

    if ( this->giveLocalNonlocalStiffnessContribution(gp, loc, s, lcontrib, atTime) == 0 ) {
        return;
    }

    for ( pos = list->begin(); pos != list->end(); ++pos ) {
        rmat = ( RankineMatNl * )( ( * pos ).nearGp )->giveMaterial();
        if ( rmat->giveClassID() == this->giveClassID() ) {
            rmat->giveRemoteNonlocalStiffnessContribution( ( * pos ).nearGp, rloc, s, rcontrib, atTime );
            coeff = gp->giveElement()->computeVolumeAround(gp) * ( * pos ).weight / status->giveIntegrationScale();

            int i, j, dim1 = loc.giveSize(), dim2 = rloc.giveSize();
            contrib.resize(dim1, dim2);
            for ( i = 1; i <= dim1; i++ ) {
                for ( j = 1; j <= dim2; j++ ) {
                    contrib.at(i, j) = -1.0 * lcontrib.at(i) * rcontrib.at(j) * coeff;
                }
            }

            dest.assemble(loc, rloc, contrib);
        }
    }
}

dynaList< localIntegrationRecord > *
RankineMatNl :: NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp)
{
    RankineMatNlStatus *status = ( RankineMatNlStatus * ) this->giveStatus(gp);
    this->buildNonlocalPointTable(gp);
    return status->giveIntegrationDomainList();
}


// computes m*gprime*Btransposed*sigmaeff for the given Gauss point
// and returns 0 of damage is not growing, 1 if it is growing
// (if damage is not growing, the contribution is not considered at all)
int
RankineMatNl :: giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                                       FloatArray &lcontrib, TimeStep *atTime)
{
    int nrows, nsize, i, j;
    double sum, nlKappa, damage, tempDamage;
    RankineMatNlStatus *status = ( RankineMatNlStatus * ) this->giveStatus(gp);
    StructuralElement *elem = ( StructuralElement * )( gp->giveElement() );
    FloatMatrix b;
    FloatArray stress;

    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    if ( tempDamage <= damage ) {
        return 0; // no contribution if damage is not growing
    }

    elem->giveLocationArray(loc, EID_MomentumBalance, s);
    status->giveTempEffectiveStress(stress);
    elem->computeBmatrixAt(gp, b);
    this->computeCumPlasticStrain(nlKappa, gp, atTime);
    double factor = computeDamageParamPrime(nlKappa);
    factor *= mm; // this factor is m*gprime
    nrows = b.giveNumberOfColumns();
    nsize = stress.giveSize();
    lcontrib.resize(nrows);
    // compute the product Btransposed*stress and multiply by factor
    for ( i = 1; i <= nrows; i++ ) {
        sum = 0.0;
        for ( j = 1; j <= nsize; j++ ) {
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
                                                        FloatArray &rcontrib, TimeStep *atTime)
{
    RankineMatNlStatus *status = ( RankineMatNlStatus * ) this->giveStatus(gp);
    StructuralElement *elem = ( StructuralElement * )( gp->giveElement() );
    elem->giveLocationArray(rloc, EID_MomentumBalance, s);
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

    int i, j;
    double sum;
    int nsize = 3;
    FloatArray eta(3);
    computeEta(eta, status);
    for ( i = 1; i <= ncols; i++ ) {
        sum = 0.;
        for ( j = 1; j <= nsize; j++ ) {
            sum += eta.at(j) * b.at(j, i);
        }

        rcontrib.at(i) = sum;
    }
}



int
RankineMatNl :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    if ( type == IST_CumPlasticStrain_2 ) {
        answer.resize(1);
        double dummy;
        // this method also stores the nonlocal kappa in status ... kappa_nl
        computeCumPlasticStrain(dummy, aGaussPoint, atTime);
        RankineMatNlStatus *status = ( RankineMatNlStatus * ) this->giveStatus(aGaussPoint);
        answer.at(1) = status->giveKappa_nl();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        computeCumPlasticStrain(answer.at(1), aGaussPoint, atTime);
        return 1;
    } else {
        return RankineMat :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

InternalStateValueType
RankineMatNl :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_CumPlasticStrain_2 || type == IST_MaxEquivalentStrainLevel ) {
        return ISVT_SCALAR;
    } else {
        return RankineMat :: giveIPValueType(type);
    }
}

int
RankineMatNl :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_CumPlasticStrain_2 || type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return RankineMat :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}

int
RankineMatNl :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if ( type == IST_CumPlasticStrain_2 || type == IST_MaxEquivalentStrainLevel ) {
        return 1;
    } else {
        return RankineMat :: giveIPValueSize(type, gp);
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
RankineMatNlStatus :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    RankineMatStatus :: updateYourself(atTime);
}


contextIOResultType
RankineMatNlStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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

    //if (!stream->write(&localEquivalentStrainForAverage,1)) THROW_CIOERR(CIO_IOERR);
    return CIO_OK;
}

contextIOResultType
RankineMatNlStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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
    //if (!stream->read (&localEquivalentStrainForAverage,1)) THROW_CIOERR(CIO_IOERR);

    return CIO_OK;
}

Interface *
RankineMatNlStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return ( StructuralNonlocalMaterialStatusExtensionInterface * ) this;
    } else {
        return RankineMatStatus :: giveInterface(type);
    }
}


#ifdef __PARALLEL_MODE
int
RankineMatNl :: packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
{
    RankineMatNlStatus *nlStatus = ( RankineMatNlStatus * ) this->giveStatus(ip);

    this->buildNonlocalPointTable(ip);
    this->updateDomainBeforeNonlocAverage(stepN);

    return buff.packDouble( nlStatus->giveLocalCumPlasticStrainForAverage() );
}

int
RankineMatNl :: unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
{
    int result;
    RankineMatNlStatus *nlStatus = ( RankineMatNlStatus * ) this->giveStatus(ip);
    double localCumPlasticStrainForAverage;

    result = buff.unpackDouble(localCumPlasticStrainForAverage);
    nlStatus->setLocalCumPlasticStrainForAverage(localCumPlasticStrainForAverage);
    return result;
}

int
RankineMatNl :: estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip)
{
    // Note: nlStatus localStrainVectorForAverage memeber must be properly sized!
    // IDNLMaterialStatus *nlStatus = (IDNLMaterialStatus*) this -> giveStatus (ip);
    return buff.givePackSize(MPI_DOUBLE, 1);
}
#endif
} // end namespace oofem
