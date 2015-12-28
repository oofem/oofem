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

#include "misesmatnl.h"
#include "../sm/Elements/structuralelement.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "nonlocalmaterialext.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "datastream.h"

namespace oofem {
REGISTER_Material(MisesMatNl);

MisesMatNl :: MisesMatNl(int n, Domain *d) : MisesMat(n, d), StructuralNonlocalMaterialExtensionInterface(d), NonlocalMaterialStiffnessInterface()
{
    Rf = 0.;
    exponent = 1.;
    averType = 0;
}


MisesMatNl :: ~MisesMatNl()
{ }


void
MisesMatNl :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                   const FloatArray &totalStrain, TimeStep *tStep)
{
    OOFEM_ERROR("3D mode not supported");
}


void
MisesMatNl :: giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp,
                                   const FloatArray &totalStrain, TimeStep *tStep)
{
    MisesMatNlStatus *nlStatus = static_cast< MisesMatNlStatus * >( this->giveStatus(gp) );
    
    double tempDam;
    performPlasticityReturn(gp, totalStrain);
    tempDam = this->computeDamage(gp, tStep);
    answer.beScaled(1.0 - tempDam, nlStatus->giveTempEffectiveStress());
    nlStatus->setTempDamage(tempDam);
    nlStatus->letTempStrainVectorBe(totalStrain);
    nlStatus->letTempStressVectorBe(answer);
}


void
MisesMatNl :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1, 1);
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    double E = lmat->give('E', gp);
    MisesMatNlStatus *status = static_cast< MisesMatNlStatus * >( this->giveStatus(gp) );
    double kappa = status->giveCumulativePlasticStrain();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double tempDamage = status->giveTempDamage();
    double damage = status->giveDamage();
    answer.at(1, 1) = ( 1 - tempDamage ) * E;
    if ( mode != TangentStiffness ) {
        return;
    }

    if ( tempKappa <= kappa ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return;
    }

    // === plastic loading ===
    const FloatArray &stressVector = status->giveTempEffectiveStress();
    double stress = stressVector.at(1);
    answer.at(1, 1) = ( 1. - tempDamage ) * E * H / ( E + H );
    if ( tempDamage > damage ) {
        double nlKappa;
        this->computeCumPlasticStrain(nlKappa, gp, tStep);
        answer.at(1, 1) -= ( 1 - mm ) * computeDamageParamPrime(nlKappa) * E / ( E + H ) * stress * sgn(stress);
    }
}


void
MisesMatNl :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep)
{
    /* Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     */

    double cumPlasticStrain;
    MisesMatNlStatus *nlstatus = static_cast< MisesMatNlStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);
    this->performPlasticityReturn(gp, strainVector);
    this->computeLocalCumPlasticStrain(cumPlasticStrain, gp, tStep);
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
    MisesMatNlStatus *nonlocStatus, *status = static_cast< MisesMatNlStatus * >( this->giveStatus(gp) );
    std :: list< localIntegrationRecord > *list = this->giveIPIntegrationList(gp);
    std :: list< localIntegrationRecord > :: iterator pos, postarget;

    // find the current Gauss point (target) in the list of it neighbors
    for ( pos = list->begin(); pos != list->end(); ++pos ) {
        if ( pos->nearGp == gp ) {
            postarget = pos;
        }
    }

    Element *elem = gp->giveElement();
    FloatArray coords;
    elem->computeGlobalCoordinates( coords, gp->giveNaturalCoordinates() );
    double xtarget = coords.at(1);

    double w, wsum = 0., x, xprev, damage, damageprev = 0.0;
    Element *nearElem;

    // process the list from the target to the end
    double distance = 0.; // distance modified by damage
    xprev = xtarget;
    for ( pos = postarget; pos != list->end(); ++pos ) {
        nearElem = ( pos->nearGp )->giveElement();
        nearElem->computeGlobalCoordinates( coords, pos->nearGp->giveNaturalCoordinates() );
        x = coords.at(1);
        nonlocStatus = static_cast< MisesMatNlStatus * >( this->giveStatus(pos->nearGp) );
        damage = nonlocStatus->giveTempDamage();
        if ( pos != postarget ) {
            distance += ( x - xprev ) * 0.5 * ( computeDistanceModifier(damage) + computeDistanceModifier(damageprev) );
        }

        w = computeWeightFunction(distance) * nearElem->computeVolumeAround(pos->nearGp);
        pos->weight = w;
        wsum += w;
        xprev = x;
        damageprev = damage;
    }

    // process the list from the target to the beginning
    distance = 0.;
    for ( pos = postarget; pos != list->begin(); --pos ) {
        nearElem = ( pos->nearGp )->giveElement();
        nearElem->computeGlobalCoordinates( coords, pos->nearGp->giveNaturalCoordinates() );
        x = coords.at(1);
        nonlocStatus = static_cast< MisesMatNlStatus * >( this->giveStatus(pos->nearGp) );
        damage = nonlocStatus->giveTempDamage();
        if ( pos != postarget ) {
            distance += ( xprev - x ) * 0.5 * ( computeDistanceModifier(damage) + computeDistanceModifier(damageprev) );
            w = computeWeightFunction(distance) * nearElem->computeVolumeAround(pos->nearGp);
            pos->weight = w;
            wsum += w;
        }

        xprev = x;
        damageprev = damage;
    }

    // the beginning must be treated separately
    pos = list->begin();
    if ( pos != postarget ) {
        nearElem = ( pos->nearGp )->giveElement();
        nearElem->computeGlobalCoordinates( coords, pos->nearGp->giveNaturalCoordinates() );
        x = coords.at(1);
        nonlocStatus = static_cast< MisesMatNlStatus * >( this->giveStatus(pos->nearGp) );
        damage = nonlocStatus->giveTempDamage();
        distance += ( xprev - x ) * 0.5 * ( computeDistanceModifier(damage) + computeDistanceModifier(damageprev) );
        w = computeWeightFunction(distance) * nearElem->computeVolumeAround(pos->nearGp);
        pos->weight = w;
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

    case 5: return ( 2. * cl ) / ( cl + Rf + ( cl - Rf ) * cos(M_PI * damage) );

    default: return 1.;
    }
}

void
MisesMatNl :: computeCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *tStep)
{
    double nonlocalContribution, nonlocalCumPlasticStrain = 0.0;
    MisesMatNlStatus *nonlocStatus, *status = static_cast< MisesMatNlStatus * >( this->giveStatus(gp) );

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(tStep);
    double localCumPlasticStrain = status->giveLocalCumPlasticStrainForAverage();
    // compute nonlocal cumulative plastic strain
    std :: list< localIntegrationRecord > *list = this->giveIPIntegrationList(gp);

    for ( auto &lir: *list ) {
        nonlocStatus = static_cast< MisesMatNlStatus * >( this->giveStatus(lir.nearGp) );
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
}

Interface *
MisesMatNl :: giveInterface(InterfaceType type)
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
MisesMatNl :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = MisesMat :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;
    result = StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    averType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, averType, _IFT_MisesMatNl_averagingtype);
    if ( averType == 2 ) {
        exponent = 0.5; // default value for averaging type 2
    }

    if ( averType == 3 ) {
        exponent = 1.; // default value for averaging type 3
    }

    if ( averType == 2 || averType == 3 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, exponent, _IFT_MisesMatNl_exp);
    }

    if ( averType >= 2 && averType <= 5 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, Rf, _IFT_MisesMatNl_rf);
    }

    return IRRT_OK;
}


void
MisesMatNl :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    StructuralNonlocalMaterialExtensionInterface :: giveInputRecord(input);

    input.setField(averType, _IFT_MisesMatNl_averagingtype);

    if ( averType == 2 || averType == 3 ) {
        input.setField(exponent, _IFT_MisesMatNl_exp);
    }

    if ( averType >= 2 && averType <= 5 ) {
        input.setField(Rf, _IFT_MisesMatNl_rf);
    }
}


double
MisesMatNl :: computeDamage(GaussPoint *gp, TimeStep *tStep)
{
    MisesMatNlStatus *nlStatus = static_cast< MisesMatNlStatus * >( this->giveStatus(gp) );
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


void
MisesMatNl :: NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                                   GaussPoint *gp, TimeStep *tStep)
{
    double coeff;
    MisesMatNlStatus *status = static_cast< MisesMatNlStatus * >( this->giveStatus(gp) );
    std :: list< localIntegrationRecord > *list = status->giveIntegrationDomainList();
    MisesMatNl *rmat;
    FloatArray rcontrib, lcontrib;
    IntArray loc, rloc;

    FloatMatrix contrib;

    if ( this->giveLocalNonlocalStiffnessContribution(gp, loc, s, lcontrib, tStep) == 0 ) {
        return;
    }

    for ( auto &lir: *list ) {
        rmat = dynamic_cast< MisesMatNl * >( lir.nearGp->giveMaterial() );
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
MisesMatNl :: NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp)
{
    MisesMatNlStatus *status = static_cast< MisesMatNlStatus * >( this->giveStatus(gp) );
    this->buildNonlocalPointTable(gp);
    return status->giveIntegrationDomainList();
}


int
MisesMatNl :: giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                                     FloatArray &lcontrib, TimeStep *tStep)
{
    double nlKappa, damage, tempDamage, dDamF;
    MisesMatNlStatus *status = static_cast< MisesMatNlStatus * >( this->giveStatus(gp) );
    StructuralElement *elem = static_cast< StructuralElement * >( gp->giveElement() );
    FloatMatrix b;

    this->computeCumPlasticStrain(nlKappa, gp, tStep);
    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    if ( ( tempDamage - damage ) > 0 ) {
        const FloatArray &stress = status->giveTempEffectiveStress();
        elem->giveLocationArray(loc, s);
        elem->computeBmatrixAt(gp, b);
        dDamF = computeDamageParamPrime(nlKappa);
        lcontrib.clear();
        lcontrib.plusProduct(b, stress, mm * dDamF);
    }

    return 1;
}


void
MisesMatNl :: giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                      FloatArray &rcontrib, TimeStep *tStep)
{
    double kappa, tempKappa;
    MisesMatNlStatus *status = static_cast< MisesMatNlStatus * >( this->giveStatus(gp) );
    StructuralElement *elem = static_cast< StructuralElement * >( gp->giveElement() );
    FloatMatrix b;

    elem->giveLocationArray(rloc, s);
    elem->computeBmatrixAt(gp, b);
    kappa = status->giveCumulativePlasticStrain();
    tempKappa = status->giveTempCumulativePlasticStrain();

    rcontrib.clear();
    if ( ( tempKappa - kappa ) > 0 ) {
        LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
        double E = lmat->give('E', gp);
        const FloatArray &stress = status->giveTempEffectiveStress();
        if ( gp->giveMaterialMode() == _1dMat ) {
            double coeff = sgn( stress.at(1) ) * E / ( E + H );
            rcontrib.plusProduct(b, stress, coeff);
            return;
        }
    }
    rcontrib.resize(b.giveNumberOfColumns());
    rcontrib.zero();
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
MisesMatNlStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    MisesMatStatus :: updateYourself(tStep);
}


contextIOResultType
MisesMatNlStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
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

    //if (!stream.write(&localEquivalentStrainForAverage,1)) THROW_CIOERR(CIO_IOERR);
    return CIO_OK;
}


contextIOResultType
MisesMatNlStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
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
    //if (!stream.read (&localEquivalentStrainForAverage,1)) THROW_CIOERR(CIO_IOERR);

    return CIO_OK;
}


Interface *
MisesMatNlStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialStatusExtensionInterface * >(this);
    } else {
        return MisesMatStatus :: giveInterface(type);
    }
}


int
MisesMatNl :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip)
{
    MisesMatNlStatus *nlStatus = static_cast< MisesMatNlStatus * >( this->giveStatus(ip) );

    this->buildNonlocalPointTable(ip);
    this->updateDomainBeforeNonlocAverage(tStep);

    return buff.write( nlStatus->giveLocalCumPlasticStrainForAverage() );
}


int
MisesMatNl :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip)
{
    int result;
    MisesMatNlStatus *nlStatus = static_cast< MisesMatNlStatus * >( this->giveStatus(ip) );
    double localCumPlasticStrainForAverage;

    result = buff.read(localCumPlasticStrainForAverage);
    nlStatus->setLocalCumPlasticStrainForAverage(localCumPlasticStrainForAverage);
    return result;
}


int
MisesMatNl :: estimatePackSize(DataStream &buff, GaussPoint *ip)
{
    // Note: nlStatus localStrainVectorForAverage memeber must be properly sized!
    // IDNLMaterialStatus *nlStatus = (IDNLMaterialStatus*) this -> giveStatus (ip);
    return buff.givePackSizeOfDouble(1);
}

} // end namespace oofem
