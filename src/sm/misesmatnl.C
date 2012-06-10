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

#include "misesmatnl.h"
#include "structuralelement.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "dynalist.h"
#include "nonlocalmaterialext.h"
#include "contextioerr.h"

#ifdef __PARALLEL_MODE
 #include "combuff.h"
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
double sign(double number)
{
    if ( number > 0 ) {
        return 1;
    } else if ( number < 0 ) {
        return -1;
    } else {
        return 0;
    }
}


MisesMatNl :: MisesMatNl(int n, Domain *d) : MisesMat(n, d), StructuralNonlocalMaterialExtensionInterface(d), NonlocalMaterialStiffnessInterface()
//
// constructor
//
{
    Rf = 0.;
    exponent = 1.;
    averType = 0;
}


MisesMatNl :: ~MisesMatNl()
//
// destructor
//
{ }


void
MisesMatNl :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                   const FloatArray &totalStrain, TimeStep *atTime)
{
    MisesMatNlStatus *nlStatus = ( MisesMatNlStatus * ) this->giveStatus(gp);
    this->initGpForNewStep(gp);

    double tempDam;
    FloatArray tempEffStress, totalStress;
    MaterialMode mode = gp->giveMaterialMode();
    performPlasticityReturn(gp, totalStrain, mode);
    tempDam = this->computeDamage(gp, atTime);
    nlStatus->giveTempEffectiveStress(tempEffStress);
    answer.beScaled( 1.0 - tempDam, tempEffStress);
    nlStatus->setTempDamage(tempDam);
    nlStatus->letTempStrainVectorBe(totalStrain);
    nlStatus->letTempStressVectorBe(answer);
}


void
MisesMatNl :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    answer.resize(1, 1);
    answer.zero();
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    double E = lmat->give('E', gp);
    double nlKappa;
    MisesMatNlStatus *status = ( MisesMatNlStatus * ) this->giveStatus(gp);
    double kappa = status->giveCumulativePlasticStrain();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double tempDamage = status->giveTempDamage();
    double damage = status->giveDamage();
    FloatArray stressVector;
    answer.at(1, 1) = ( 1 - tempDamage ) * E;
    if ( mode != TangentStiffness ) {
        return;
    }

    if ( tempKappa <= kappa ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return;
    }

    // === plastic loading ===
    status->giveTempEffectiveStress(stressVector);
    double stress = stressVector.at(1);
    answer.at(1, 1) = ( 1. - tempDamage ) * E * H / ( E + H );
    if ( tempDamage > damage ) {
        this->computeCumPlasticStrain(nlKappa, gp, atTime);
        answer.at(1, 1) = answer.at(1, 1) - ( 1 - mm ) * computeDamageParamPrime(nlKappa) * E / ( E + H ) * stress * sign(stress);
    }
}


void
MisesMatNl :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime)
{
    /* Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     */

    double cumPlasticStrain;
    MisesMatNlStatus *nlstatus = ( MisesMatNlStatus * ) this->giveStatus(gp);

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);
    MaterialMode mode = gp->giveMaterialMode();
    this->performPlasticityReturn(gp, strainVector, mode);
    this->computeLocalCumPlasticStrain(cumPlasticStrain, gp, atTime);
    // standard formulation based on averaging of equivalent strain
    nlstatus->setLocalCumPlasticStrainForAverage(cumPlasticStrain);

    // influence of damage on weight function
    if ( averType >= 2 && averType <= 5 ) {
        this->modifyNonlocalWeightFunctionAround(gp);
    }
}


void
MisesMatNl :: modifyNonlocalWeightFunctionAround(GaussPoint *gp)
{
    MisesMatNlStatus *nonlocStatus, *status = ( MisesMatNlStatus * ) this->giveStatus(gp);
    dynaList< localIntegrationRecord > *list = this->giveIPIntegrationList(gp);
    dynaList< localIntegrationRecord > :: iterator pos, postarget;

    // find the current Gauss point (target) in the list of it neighbors
    for ( pos = list->begin(); pos != list->end(); ++pos ) {
        if ( ( * pos ).nearGp == gp ) {
            postarget = pos;
        }
    }

    Element *elem = gp->giveElement();
    FloatArray coords;
    elem->computeGlobalCoordinates( coords, * ( gp->giveCoordinates() ) );
    double xtarget = coords.at(1);

    double w, wsum = 0., x, xprev, damage, damageprev;
    int n;
    Element *nearElem;

    // process the list from the target to the end
    double distance = 0.; // distance modified by damage
    xprev = xtarget;
    for ( pos = postarget; pos != list->end(); ++pos ) {
        nearElem = ( ( * pos ).nearGp )->giveElement();
        nearElem->computeGlobalCoordinates( coords, * ( ( ( * pos ).nearGp )->giveCoordinates() ) );
        x = coords.at(1);
        nonlocStatus = ( MisesMatNlStatus * ) this->giveStatus( ( * pos ).nearGp );
        damage = nonlocStatus->giveTempDamage();
        if ( pos != postarget ) {
            distance += ( x - xprev ) * 0.5 * ( computeDistanceModifier(damage) + computeDistanceModifier(damageprev) );
        }

        w = computeWeightFunction(distance) * nearElem->computeVolumeAround( ( * pos ).nearGp );
        ( * pos ).weight = w;
        wsum += w;
        xprev = x;
        damageprev = damage;
    }

    // process the list from the target to the beginning
    distance = 0.;
    for ( pos = postarget; pos != list->begin(); --pos ) {
        nearElem = ( ( * pos ).nearGp )->giveElement();
        nearElem->computeGlobalCoordinates( coords, * ( ( ( * pos ).nearGp )->giveCoordinates() ) );
        x = coords.at(1);
        nonlocStatus = ( MisesMatNlStatus * ) this->giveStatus( ( * pos ).nearGp );
        damage = nonlocStatus->giveTempDamage();
        if ( pos != postarget ) {
            distance += ( xprev - x ) * 0.5 * ( computeDistanceModifier(damage) + computeDistanceModifier(damageprev) );
            w = computeWeightFunction(distance) * nearElem->computeVolumeAround( ( * pos ).nearGp );
            ( * pos ).weight = w;
            wsum += w;
        }

        xprev = x;
        damageprev = damage;
    }

    // the beginning must be treated separately
    pos = list->begin();
    if ( pos != postarget ) {
        nearElem = ( ( * pos ).nearGp )->giveElement();
        nearElem->computeGlobalCoordinates( coords, * ( ( ( * pos ).nearGp )->giveCoordinates() ) );
        x = coords.at(1);
        nonlocStatus = ( MisesMatNlStatus * ) this->giveStatus( ( * pos ).nearGp );
        damage = nonlocStatus->giveTempDamage();
        n = ( ( * pos ).nearGp )->giveElement()->giveNumber();
        distance += ( xprev - x ) * 0.5 * ( computeDistanceModifier(damage) + computeDistanceModifier(damageprev) );
        w = computeWeightFunction(distance) * nearElem->computeVolumeAround( ( * pos ).nearGp );
        ( * pos ).weight = w;
        wsum += w;
    }

    status->setIntegrationScale(wsum);
}

double
MisesMatNl :: computeDistanceModifier(double damage)
{
    switch ( averType ) {
    case 2: return 1. / ( Rf / cl + ( 1. - Rf / cl ) * pow(1. - damage, exponent) );

    case 3: if ( damage == 0. ) {
            return 1.;
    } else {
            return 1. / ( 1. - ( 1. - Rf / cl ) * pow(damage, exponent) );
    }

    case 4: return 1. / pow(Rf / cl, damage);

    case 5: return ( 2. * cl ) / ( cl + Rf + ( cl - Rf ) * cos(3.1415926 * damage) );

    default: return 1.;
    }
}

void
MisesMatNl :: computeCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *atTime)
{
    double nonlocalContribution, nonlocalCumPlasticStrain = 0.0;
    MisesMatNlStatus *nonlocStatus, *status = ( MisesMatNlStatus * ) this->giveStatus(gp);

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(atTime);
    double localCumPlasticStrain = status->giveLocalCumPlasticStrainForAverage();
    // compute nonlocal cumulative plastic strain
    dynaList< localIntegrationRecord > *list = this->giveIPIntegrationList(gp);
    dynaList< localIntegrationRecord > :: iterator pos;

    for ( pos = list->begin(); pos != list->end(); ++pos ) {
        nonlocStatus = ( MisesMatNlStatus * ) this->giveStatus( ( * pos ).nearGp );
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
}

Interface *
MisesMatNl :: giveInterface(InterfaceType type)
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
MisesMatNl :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    MisesMat :: initializeFrom(ir);
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);

    averType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, averType, IFT_MisesMatNl_averagingtype, "averagingtype");
    if ( averType == 2 ) {
        exponent = 0.5; // default value for averaging type 2
    }

    if ( averType == 3 ) {
        exponent = 1.; // default value for averaging type 3
    }

    if ( averType == 2 || averType == 3 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, exponent, IFT_MisesMatNl_averagingtype, "exp");
    }

    if ( averType >= 2 && averType <= 5 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, Rf, IFT_MisesMatNl_averagingtype, "rf");
    }

    return IRRT_OK;
}


int
MisesMatNl :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    MisesMat :: giveInputRecordString(str, keyword);
    StructuralNonlocalMaterialExtensionInterface :: giveInputRecordString(str, false);
    sprintf(buff, " r %e", this->cl);
    str += buff;

    return 1;
}


double
MisesMatNl :: computeDamage(GaussPoint *gp, TimeStep *atTime)
{
    MisesMatNlStatus *nlStatus = ( MisesMatNlStatus * ) this->giveStatus(gp);
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


void
MisesMatNl :: NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                                   GaussPoint *gp, TimeStep *atTime)
{
    double coeff;
    MisesMatNlStatus *status = ( MisesMatNlStatus * ) this->giveStatus(gp);
    dynaList< localIntegrationRecord > *list = status->giveIntegrationDomainList();
    dynaList< localIntegrationRecord > :: iterator pos;
    MisesMatNl *rmat;
    FloatArray rcontrib, lcontrib;
    IntArray loc, rloc;

    FloatMatrix contrib;

    if ( this->giveLocalNonlocalStiffnessContribution(gp, loc, s, lcontrib, atTime) == 0 ) {
        return;
    }

    for ( pos = list->begin(); pos != list->end(); ++pos ) {
        rmat = ( MisesMatNl * )( ( * pos ).nearGp )->giveMaterial();
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
MisesMatNl :: NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp)
{
    MisesMatNlStatus *status = ( MisesMatNlStatus * ) this->giveStatus(gp);
    this->buildNonlocalPointTable(gp);
    return status->giveIntegrationDomainList();
}


int
MisesMatNl :: giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                                     FloatArray &lcontrib, TimeStep *atTime)
{
    int nrows, nsize, i, j;
    double sum, nlKappa, damage, tempDamage, dDamF;
    MisesMatNlStatus *status = ( MisesMatNlStatus * ) this->giveStatus(gp);
    StructuralElement *elem = ( StructuralElement * )( gp->giveElement() );
    FloatMatrix b;
    FloatArray stress;

    this->computeCumPlasticStrain(nlKappa, gp, atTime);
    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    if ( ( tempDamage - damage ) > 0 ) {
        elem->giveLocationArray(loc, EID_MomentumBalance, s);
        status->giveTempEffectiveStress(stress);
        elem->computeBmatrixAt(gp, b);
        dDamF = computeDamageParamPrime(nlKappa);
        nrows = b.giveNumberOfColumns();
        nsize = stress.giveSize();
        lcontrib.resize(nrows);

        for ( i = 1; i <= nrows; i++ ) {
            sum = 0.0;
            for ( j = 1; j <= nsize; j++ ) {
                sum += b.at(j, i) * stress.at(j);
            }

            lcontrib.at(i) = sum * mm * dDamF;
        }
    }

    return 1;
}


void
MisesMatNl :: giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                      FloatArray &rcontrib, TimeStep *atTime)
{
    int ncols, nsize, i, j;
    double sum, kappa, tempKappa;
    MisesMatNlStatus *status = ( MisesMatNlStatus * ) this->giveStatus(gp);
    StructuralElement *elem = ( StructuralElement * )( gp->giveElement() );
    FloatMatrix b;
    FloatArray stress;
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    double E = lmat->give('E', gp);

    elem->giveLocationArray(rloc, EID_MomentumBalance, s);
    elem->computeBmatrixAt(gp, b);
    ncols = b.giveNumberOfColumns();
    rcontrib.resize(ncols);
    kappa = status->giveCumulativePlasticStrain();
    tempKappa = status->giveTempCumulativePlasticStrain();

    if ( ( tempKappa - kappa ) > 0 ) {
        status->giveTempEffectiveStress(stress);
        if ( gp->giveMaterialMode() == _1dMat ) {
            nsize = stress.giveSize();
            double coeff = sign( stress.at(1) ) * E / ( E + H );
            for ( i = 1; i <= ncols; i++ ) {
                sum = 0.;
                for ( j = 1; j <= nsize; j++ ) {
                    sum += stress.at(j) * coeff * b.at(j, i);
                }

                rcontrib.at(i) = sum;
            }
        }
    } else {
        rcontrib.zero();
    }
}


/*********************************************status**************************************************************/

MisesMatNlStatus :: MisesMatNlStatus(int n, Domain *d, GaussPoint *g) :
    MisesMatStatus(n, d, g), StructuralNonlocalMaterialStatusExtensionInterface()
{
    localCumPlasticStrainForAverage = 0.0;
}


MisesMatNlStatus :: ~MisesMatNlStatus()
{ }


void
MisesMatNlStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
    fprintf(file, "}\n");
}


void
MisesMatNlStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    MisesMatStatus :: initTempStatus();
}


void
MisesMatNlStatus :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    MisesMatStatus :: updateYourself(atTime);
}


contextIOResultType
MisesMatNlStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = MisesMatStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //if (!stream->write(&localEquivalentStrainForAverage,1)) THROW_CIOERR(CIO_IOERR);
    return CIO_OK;
}


contextIOResultType
MisesMatNlStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = MisesMatStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    //if (!stream->read (&localEquivalentStrainForAverage,1)) THROW_CIOERR(CIO_IOERR);

    return CIO_OK;
}


Interface *
MisesMatNlStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return ( StructuralNonlocalMaterialStatusExtensionInterface * ) this;
    } else {
        return MisesMatStatus :: giveInterface(type);
    }
}


#ifdef __PARALLEL_MODE
int
MisesMatNl :: packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
{
    MisesMatNlStatus *nlStatus = ( MisesMatNlStatus * ) this->giveStatus(ip);

    this->buildNonlocalPointTable(ip);
    this->updateDomainBeforeNonlocAverage(stepN);

    return buff.packDouble( nlStatus->giveLocalCumPlasticStrainForAverage() );
}


int
MisesMatNl :: unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
{
    int result;
    MisesMatNlStatus *nlStatus = ( MisesMatNlStatus * ) this->giveStatus(ip);
    double localCumPlasticStrainForAverage;

    result = buff.unpackDouble(localCumPlasticStrainForAverage);
    nlStatus->setLocalCumPlasticStrainForAverage(localCumPlasticStrainForAverage);
    return result;
}


int
MisesMatNl :: estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip)
{
    // Note: nlStatus localStrainVectorForAverage memeber must be properly sized!
    // IDNLMaterialStatus *nlStatus = (IDNLMaterialStatus*) this -> giveStatus (ip);
    return buff.givePackSize(MPI_DOUBLE, 1);
}
#endif
} // end namespace oofem
