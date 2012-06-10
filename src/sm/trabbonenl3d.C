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

#include "trabbonenl3d.h"
#include "structuralelement.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "dynalist.h"
#include "nonlocalmaterialext.h"

#ifdef __PARALLEL_MODE
 #include "idmnl1.h"
 #include "combuff.h"
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {

TrabBoneNL3D :: TrabBoneNL3D(int n, Domain *d) : TrabBone3D(n, d), StructuralNonlocalMaterialExtensionInterface(d), NonlocalMaterialStiffnessInterface()
{
    R = 0.;
}


TrabBoneNL3D :: ~TrabBoneNL3D()
{}


void
TrabBoneNL3D :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime)
{
    FloatArray SDstrainVector, fullSDStrainVector;
    double cumPlastStrain;
    TrabBoneNL3DStatus *nlStatus = ( TrabBoneNL3DStatus * ) this->giveStatus(gp);

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector, atTime, VM_Total);

    nlStatus->letTempStrainVectorBe(strainVector);

    StrainVector strain( strainVector, gp->giveMaterialMode() );

    this->performPlasticityReturn(gp, strain);
    this->computeLocalCumPlastStrain(cumPlastStrain, strain, gp, atTime);
    nlStatus->setLocalCumPlastStrainForAverage(cumPlastStrain);
}


void
TrabBoneNL3D :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseForm form, MatResponseMode mode, GaussPoint *gp,
                                              TimeStep *atTime)
{
    TrabBoneNL3DStatus *nlStatus = ( TrabBoneNL3DStatus * ) this->giveStatus(gp);

    double tempDam, beta, nlKappa;
    FloatArray tempEffectiveStress, tempTensor2, prodTensor, plasFlowDirec;
    FloatMatrix elasticity, compliance, SSaTensor, secondTerm, thirdTerm, tangentMatrix;

    if ( mode == ElasticStiffness ) {
        this->constructAnisoComplTensor(compliance);
        elasticity.beInverseOf(compliance);

        answer = elasticity;
    } else if ( mode == SecantStiffness )     {
        this->constructAnisoComplTensor(compliance);
        elasticity.beInverseOf(compliance);
        tempDam = nlStatus->giveTempDam();

        answer = elasticity;
        answer.times(1.0 - tempDam);
    } else if ( mode == TangentStiffness )     {
        double kappa = nlStatus->giveKappa();
        double tempKappa = nlStatus->giveTempKappa();
        double dKappa = tempKappa - kappa;
        if ( dKappa < 10.e-9 ) {
            dKappa = 0;
        }

        if ( dKappa > 0.0 ) {
            // Imports
            tempEffectiveStress = * nlStatus->giveTempEffectiveStress();
            this->computeCumPlastStrain(nlKappa, gp, atTime);
            tempDam = nlStatus->giveTempDam();
            double dam = nlStatus->giveDam();
            plasFlowDirec = * nlStatus->givePlasFlowDirec();
            SSaTensor = * nlStatus->giveSSaTensor();
            beta = nlStatus->giveBeta();
            // Construction of the dyadic product tensor
            prodTensor.beTProductOf(SSaTensor, plasFlowDirec);
            // Construction of the tangent stiffness third term
            if ( tempDam - dam > 0 ) {
                thirdTerm.beDyadicProductOf(tempEffectiveStress, prodTensor);
                thirdTerm.times(-expDam * critDam * exp(-expDam * nlKappa) * ( 1.0 - mParam ) / beta);
            } else   {
                thirdTerm.resize(6, 6);
            }

            // Construction of the tangent stiffness second term
            tempTensor2.beProductOf(SSaTensor, plasFlowDirec);
            secondTerm.beDyadicProductOf(tempTensor2, prodTensor);
            secondTerm.times(-( 1.0 - tempDam ) / beta);
            // Construction of the tangent stiffness
            tangentMatrix = SSaTensor;
            tangentMatrix.times(1.0 - tempDam);
            tangentMatrix.add(secondTerm);
            tangentMatrix.add(thirdTerm);

            answer = tangentMatrix;
        } else   {
            // Import of state variables
            tempDam = nlStatus->giveTempDam();
            // Construction of the tangent stiffness
            this->constructAnisoComplTensor(compliance);
            elasticity.beInverseOf(compliance);
            answer = elasticity;
            answer.times(1.0 - tempDam);
        }
    }

    nlStatus->setSmtrx(answer);
}


void
TrabBoneNL3D :: NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s, GaussPoint *gp, TimeStep *atTime)
{
    TrabBoneNL3DStatus *nlStatus = ( TrabBoneNL3DStatus * ) this->giveStatus(gp);
    dynaList< localIntegrationRecord > *list = nlStatus->giveIntegrationDomainList();
    dynaList< localIntegrationRecord > :: iterator pos;
    TrabBoneNL3D *rmat;

    double coeff;
    FloatArray rcontrib, lcontrib;
    IntArray loc, rloc;

    FloatMatrix contrib;

    if ( this->giveLocalNonlocalStiffnessContribution(gp, loc, s, lcontrib, atTime) == 0 ) {
        return;
    }

    for ( pos = list->begin(); pos != list->end(); ++pos ) {
        rmat = ( TrabBoneNL3D * ) ( ( * pos ).nearGp )->giveMaterial();
        if ( rmat->giveClassID() == this->giveClassID() ) {
            rmat->giveRemoteNonlocalStiffnessContribution( ( * pos ).nearGp, rloc, s, rcontrib, atTime );
            coeff = gp->giveElement()->computeVolumeAround(gp) * ( * pos ).weight / nlStatus->giveIntegrationScale();

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
TrabBoneNL3D :: NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp)
{
    TrabBoneNL3DStatus *nlStatus = ( TrabBoneNL3DStatus * ) this->giveStatus(gp);
    this->buildNonlocalPointTable(gp);
    return nlStatus->giveIntegrationDomainList();
}


int
TrabBoneNL3D :: giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                                       FloatArray &lcontrib, TimeStep *atTime)
{
    TrabBoneNL3DStatus *nlStatus = ( TrabBoneNL3DStatus * ) this->giveStatus(gp);
    StructuralElement *elem = ( StructuralElement * ) ( gp->giveElement() );

    int nrows, nsize, i, j;
    double sum, nlKappa, dDamFunc, dam, tempDam;
    FloatArray localNu;
    FloatMatrix b;

    this->computeCumPlastStrain(nlKappa, gp, atTime);
    dam = nlStatus->giveDam();
    tempDam = nlStatus->giveTempDam();

    if ( ( tempDam - dam ) > 0.0 ) {
        elem->giveLocationArray(loc, EID_MomentumBalance, s);
        localNu = * nlStatus->giveTempEffectiveStress();

        elem->giveLocationArray( loc, EID_MomentumBalance, EModelDefaultEquationNumbering() );
        elem->computeBmatrixAt(gp, b);
        dDamFunc = expDam * critDam * exp(-expDam * nlKappa);

        nrows = b.giveNumberOfColumns();
        nsize = localNu.giveSize();
        lcontrib.resize(nrows);
        for ( i = 1; i <= nrows; i++ ) {
            sum = 0.0;
            for ( j = 1; j <= nsize; j++ ) {
                sum += b.at(j, i) * localNu.at(j);
            }

            lcontrib.at(i) = mParam * dDamFunc * sum;
        }

        return 1;
    } else {
        loc.resize(0);
        return 0;
    }

    return 0;
}


void
TrabBoneNL3D :: giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                        FloatArray &rcontrib, TimeStep *atTime)
{
    TrabBoneNL3DStatus *nlStatus = ( TrabBoneNL3DStatus * ) this->giveStatus(gp);
    StructuralElement *elem = ( StructuralElement * ) ( gp->giveElement() );

    int ncols, nsize, i, j;
    double sum, beta;
    FloatArray remoteNu, plasFlowDirec, prodTensor;
    FloatMatrix b, SSaTensor;

    elem->giveLocationArray(rloc, EID_MomentumBalance, s);
    elem->computeBmatrixAt(gp, b);

    ncols = b.giveNumberOfColumns();
    rcontrib.resize(ncols);


    double kappa = nlStatus->giveKappa();
    double tempKappa = nlStatus->giveTempKappa();
    double dKappa = tempKappa - kappa;
    if ( dKappa < 10.e-9 ) {
        dKappa = 0;
    }

    if ( dKappa > 0.0 ) {
        plasFlowDirec = * nlStatus->givePlasFlowDirec();
        SSaTensor = * nlStatus->giveSSaTensor();
        beta = nlStatus->giveBeta();

        prodTensor.beTProductOf(SSaTensor, plasFlowDirec);
        remoteNu.beScaled( 1.0 / beta, prodTensor);
        nsize = remoteNu.giveSize();

        for ( i = 1; i <= ncols; i++ ) {
            sum = 0.0;
            for ( j = 1; j <= nsize; j++ ) {
                sum += remoteNu.at(j) * b.at(j, i);
            }

            rcontrib.at(i) = sum;
        }
    } else   {
        for ( i = 1; i <= ncols; i++ ) {
            rcontrib.at(i) = 0.;
        }
    }
}


void
TrabBoneNL3D :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                     const FloatArray &totalStrain, TimeStep *atTime)
{
    TrabBoneNL3DStatus *nlStatus = ( TrabBoneNL3DStatus * ) this->giveStatus(gp);
    this->initGpForNewStep(gp);

    int i;
    double tempDam;
    FloatArray effStress, totalStress, densStress;

    performPlasticityReturn(gp, totalStrain);
    tempDam = computeDamage(gp, atTime);
    effStress = * nlStatus->giveTempEffectiveStress();

    totalStress.beScaled( 1.0 - tempDam, effStress);

    for ( i = 1; i <= 6; i++ ) {
        if ( sqrt( totalStress.at(i) * totalStress.at(i) ) < pow(10.0, -8.0) ) {
            totalStress.at(i) = 0.;
        }
    }

    computePlasStrainEnerDensity(gp, totalStrain, totalStress);

    if ( JCrit != 0. ) {
        computeDensificationStress(densStress, gp, totalStrain, atTime);
    } else   {
        densStress.resize(6);
    }

    answer = totalStress;
    nlStatus->setTempDam(tempDam);
    nlStatus->letTempStrainVectorBe(totalStrain);
    nlStatus->letTempStressVectorBe(answer);
}


void
TrabBoneNL3D :: computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime)
{
    double nonlocalContribution, nonlocalCumPlastStrain = 0.0, coeff = 0.0;
    TrabBoneNL3DStatus *nonlocStatus, *nlStatus = ( TrabBoneNL3DStatus * ) this->giveStatus(gp);

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(atTime);

    dynaList< localIntegrationRecord > *list = nlStatus->giveIntegrationDomainList();
    dynaList< localIntegrationRecord > :: iterator pos;

    for ( pos = list->begin(); pos != list->end(); ++pos ) {
        coeff = gp->giveElement()->computeVolumeAround(gp) * ( * pos ).weight / nlStatus->giveIntegrationScale();

        nonlocStatus = ( TrabBoneNL3DStatus * ) this->giveStatus( ( * pos ).nearGp );
        nonlocalContribution = nonlocStatus->giveLocalCumPlastStrainForAverage();
        nonlocalContribution *= ( * pos ).weight;
        nonlocalCumPlastStrain += nonlocalContribution;
    }

    nonlocalCumPlastStrain *= 1. / nlStatus->giveIntegrationScale();

    double localCumPlastStrain = nlStatus->giveLocalCumPlastStrainForAverage();
    kappa = mParam * nonlocalCumPlastStrain + ( 1 - mParam ) * localCumPlastStrain;
}


Interface *
TrabBoneNL3D :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return ( StructuralNonlocalMaterialExtensionInterface * ) this;
    } else if ( type == NonlocalMaterialStiffnessInterfaceType ) {
        return ( NonlocalMaterialStiffnessInterface * ) this;
    }
    //
    else {
        return NULL;
    }
}


IRResultType
TrabBoneNL3D :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    TrabBone3D :: initializeFrom(ir);
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, R, IFT_TrabBoneNL3D_r, "r"); // Macro
    if ( R < 0.0 ) {
        R = 0.0;
    }

    mParam = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, IFT_TrabBoneNL3D_m, "m"); // Macro

    return IRRT_OK;
}


int
TrabBoneNL3D :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    TrabBone3D :: giveInputRecordString(str, keyword);
    StructuralNonlocalMaterialExtensionInterface :: giveInputRecordString(str, false);
    sprintf(buff, " r %e", this->R);
    str += buff;

    return 1;
}


double
TrabBoneNL3D :: computeWeightFunction(const FloatArray &src, const FloatArray &coord)
{
    double dist = src.distance(coord);

    if ( ( dist >= 0. ) && ( dist <= this->R ) ) {
        double help = ( 1. - dist * dist / ( R * R ) );
        return help * help;
    }

    return 0.0;
}


TrabBoneNL3DStatus :: TrabBoneNL3DStatus(int n, Domain *d, GaussPoint *g) :
    TrabBone3DStatus(n, d, g), StructuralNonlocalMaterialStatusExtensionInterface()
{
    localCumPlastStrainForAverage = 0.0;
}


TrabBoneNL3DStatus :: ~TrabBoneNL3DStatus()
{}


void
TrabBoneNL3DStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status {");
    fprintf(file, " plastrains ");
    int n = plasDef.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", plasDef.at(i) );
    }

    fprintf(file, " ,");
    fprintf(file, " kappa % .4e ,", kappa);
    fprintf(file, " dam  % .4e ,", tempDam);
    fprintf(file, " esed  % .4e ,", this->tempTSED - this->tempPSED);
    fprintf(file, " psed  % .4e ,", this->tempPSED);
    fprintf(file, " tsed  % .4e", this->tempTSED);
    fprintf(file, "}\n");
}


void
TrabBoneNL3DStatus :: initTempStatus()
{
    TrabBone3DStatus :: initTempStatus();
}


void
TrabBoneNL3DStatus :: updateYourself(TimeStep *atTime)
{
    TrabBone3DStatus :: updateYourself(atTime);
}


Interface *
TrabBoneNL3DStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return this;
    } else {
        return NULL;
    }
}


#ifdef __PARALLEL_MODE
int
TrabBoneNL3D :: packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
{
    abort();
 #if 0
    IDNLMaterialStatus *nlStatus = ( IDNLMaterialStatus * ) this->giveStatus(ip);

    this->buildNonlocalPointTable(ip);
    this->updateDomainBeforeNonlocAverage(stepN);

    return buff.packDouble( nlStatus->giveLocalEquivalentStrainForAverage() );

 #endif
}


int
TrabBoneNL3D :: unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
{
    abort();
 #if 0
    int result;
    IDNLMaterialStatus *nlStatus = ( IDNLMaterialStatus * ) this->giveStatus(ip);
    double localEquivalentStrainForAverage;

    result = buff.unpackDouble(localEquivalentStrainForAverage);
    nlStatus->setLocalEquivalentStrainForAverage(localEquivalentStrainForAverage);
    return result;

 #endif
}


int
TrabBoneNL3D :: estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip)
{
    abort();
 #if 0
    // Note: nlStatus localStrainVectorForAverage memeber must be properly sized!
    // IDNLMaterialStatus *nlStatus = (IDNLMaterialStatus*) this -> giveStatus (ip);
    return buff.givePackSize(MPI_DOUBLE, 1);

 #endif
}
#endif

} // end namespace oofem
