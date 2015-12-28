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
#include "trabbonenl3d.h"
#include "../sm/Elements/structuralelement.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "nonlocalmaterialext.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "unknownnumberingscheme.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

#include <cstdlib>

namespace oofem {
REGISTER_Material(TrabBoneNL3D);

TrabBoneNL3D :: TrabBoneNL3D(int n, Domain *d) : TrabBone3D(n, d), StructuralNonlocalMaterialExtensionInterface(d), NonlocalMaterialStiffnessInterface()
{
    R = 0.;
}


TrabBoneNL3D :: ~TrabBoneNL3D()
{ }


void
TrabBoneNL3D :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray SDstrainVector;
    double cumPlastStrain;
    TrabBoneNL3DStatus *nlStatus = static_cast< TrabBoneNL3DStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector, tStep, VM_Total);

    nlStatus->letTempStrainVectorBe(strainVector);

    this->performPlasticityReturn(gp, strainVector, tStep);
    this->computeLocalCumPlastStrain(cumPlastStrain, strainVector, gp, tStep);
    nlStatus->setLocalCumPlastStrainForAverage(cumPlastStrain);
}


void
TrabBoneNL3D :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseMode mode, GaussPoint *gp,
                                              TimeStep *tStep)
{
    TrabBoneNL3DStatus *nlStatus = static_cast< TrabBoneNL3DStatus * >( this->giveStatus(gp) );

    double tempDam, beta, nlKappa;
    FloatArray tempEffectiveStress, tempTensor2, prodTensor, plasFlowDirec;
    FloatMatrix elasticity, compliance, SSaTensor, secondTerm, thirdTerm, tangentMatrix;

    if ( mode == ElasticStiffness ) {
        this->constructAnisoComplTensor(compliance);
        elasticity.beInverseOf(compliance);

        answer = elasticity;
    } else if ( mode == SecantStiffness ) {
        this->constructAnisoComplTensor(compliance);
        elasticity.beInverseOf(compliance);
        tempDam = nlStatus->giveTempDam();

        answer = elasticity;
        answer.times(1.0 - tempDam);
    } else if ( mode == TangentStiffness ) {
        double kappa = nlStatus->giveKappa();
        double tempKappa = nlStatus->giveTempKappa();
        double dKappa = tempKappa - kappa;
        if ( dKappa < 10.e-9 ) {
            dKappa = 0;
        }

        if ( dKappa > 0.0 ) {
            // Imports
            tempEffectiveStress = nlStatus->giveTempEffectiveStress();
            this->computeCumPlastStrain(nlKappa, gp, tStep);
            tempDam = nlStatus->giveTempDam();
            double dam = nlStatus->giveDam();
            plasFlowDirec = nlStatus->givePlasFlowDirec();
            SSaTensor = nlStatus->giveSSaTensor();
            beta = nlStatus->giveBeta();
            // Construction of the dyadic product tensor
            prodTensor.beTProductOf(SSaTensor, plasFlowDirec);
            // Construction of the tangent stiffness third term
            if ( tempDam - dam > 0 ) {
                thirdTerm.beDyadicProductOf(tempEffectiveStress, prodTensor);
                thirdTerm.times(-expDam * critDam * exp(-expDam * nlKappa) * ( 1.0 - mParam ) / beta);
            } else {
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
        } else {
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
TrabBoneNL3D :: NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s, GaussPoint *gp, TimeStep *tStep)
{
    TrabBoneNL3DStatus *nlStatus = static_cast< TrabBoneNL3DStatus * >( this->giveStatus(gp) );
    std :: list< localIntegrationRecord > *list = nlStatus->giveIntegrationDomainList();
    TrabBoneNL3D *rmat;

    double coeff;
    FloatArray rcontrib, lcontrib;
    IntArray loc, rloc;

    FloatMatrix contrib;

    if ( this->giveLocalNonlocalStiffnessContribution(gp, loc, s, lcontrib, tStep) == 0 ) {
        return;
    }

    for ( auto &lir: *list ) {
        rmat = dynamic_cast< TrabBoneNL3D * >( lir.nearGp->giveMaterial() );
        if ( rmat ) {
            rmat->giveRemoteNonlocalStiffnessContribution(lir.nearGp, rloc, s, rcontrib, tStep);
            coeff = gp->giveElement()->computeVolumeAround(gp) * lir.weight / nlStatus->giveIntegrationScale();

            contrib.clear();
            contrib.plusDyadUnsym(lcontrib, rcontrib, - 1.0 * coeff);
            dest.assemble(loc, rloc, contrib);
        }
    }
}

std :: list< localIntegrationRecord > *
TrabBoneNL3D :: NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp)
{
    TrabBoneNL3DStatus *nlStatus = static_cast< TrabBoneNL3DStatus * >( this->giveStatus(gp) );
    this->buildNonlocalPointTable(gp);
    return nlStatus->giveIntegrationDomainList();
}

int
TrabBoneNL3D :: giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                                       FloatArray &lcontrib, TimeStep *tStep)
{
    TrabBoneNL3DStatus *nlStatus = static_cast< TrabBoneNL3DStatus * >( this->giveStatus(gp) );
    StructuralElement *elem = static_cast< StructuralElement * >( gp->giveElement() );

    int nrows, nsize;
    double sum, nlKappa, dDamFunc, dam, tempDam;
    FloatArray localNu;
    FloatMatrix b;

    this->computeCumPlastStrain(nlKappa, gp, tStep);
    dam = nlStatus->giveDam();
    tempDam = nlStatus->giveTempDam();

    if ( ( tempDam - dam ) > 0.0 ) {
        elem->giveLocationArray(loc, s);
        localNu = nlStatus->giveTempEffectiveStress();

        elem->giveLocationArray(loc, EModelDefaultEquationNumbering() );
        elem->computeBmatrixAt(gp, b);
        dDamFunc = expDam * critDam * exp(-expDam * nlKappa);

        nrows = b.giveNumberOfColumns();
        nsize = localNu.giveSize();
        lcontrib.resize(nrows);
        for ( int i = 1; i <= nrows; i++ ) {
            sum = 0.0;
            for ( int j = 1; j <= nsize; j++ ) {
                sum += b.at(j, i) * localNu.at(j);
            }

            lcontrib.at(i) = mParam * dDamFunc * sum;
        }

        return 1;
    } else {
        loc.clear();
        return 0;
    }
}

void
TrabBoneNL3D :: giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                        FloatArray &rcontrib, TimeStep *tStep)
{
    TrabBoneNL3DStatus *nlStatus = static_cast< TrabBoneNL3DStatus * >( this->giveStatus(gp) );
    StructuralElement *elem = static_cast< StructuralElement * >( gp->giveElement() );

    FloatMatrix b;

    elem->giveLocationArray(rloc, s);
    elem->computeBmatrixAt(gp, b);

    double kappa = nlStatus->giveKappa();
    double tempKappa = nlStatus->giveTempKappa();
    double dKappa = tempKappa - kappa;
    if ( dKappa < 10.e-9 ) {
        dKappa = 0;
    }

    if ( dKappa > 0.0 ) {
        FloatArray remoteNu, prodTensor;
        const FloatArray &plasFlowDirec = nlStatus->givePlasFlowDirec();
        const FloatMatrix &SSaTensor = nlStatus->giveSSaTensor();
        double beta = nlStatus->giveBeta();

        prodTensor.beTProductOf(SSaTensor, plasFlowDirec);
        remoteNu = 1 / beta * prodTensor;
        rcontrib.beTProductOf(b, remoteNu);
    } else {
        rcontrib.resize(b.giveNumberOfColumns());
        rcontrib.zero();
    }
}


void
TrabBoneNL3D :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                        const FloatArray &totalStrain, TimeStep *tStep)
{
    TrabBoneNL3DStatus *nlStatus = static_cast< TrabBoneNL3DStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    double tempDam;
    FloatArray effStress, totalStress, densStress;

    performPlasticityReturn(gp, totalStrain, tStep);
    tempDam = computeDamage(gp, tStep);
    effStress = nlStatus->giveTempEffectiveStress();

    totalStress = ( 1 - tempDam ) * effStress;

    for ( int i = 1; i <= 6; i++ ) {
        if ( sqrt( totalStress.at(i) * totalStress.at(i) ) < 1e-8 ) {
            totalStress.at(i) = 0.;
        }
    }

    computePlasStrainEnerDensity(gp, totalStrain, totalStress);

    if ( densCrit != 0. ) {
        computeDensificationStress(densStress, gp, totalStrain, tStep);
        answer.add(densStress);
    }

    answer = totalStress;
    nlStatus->setTempDam(tempDam);
    nlStatus->letTempStrainVectorBe(totalStrain);
    nlStatus->letTempStressVectorBe(answer);
}


void
TrabBoneNL3D :: computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep)
{
    double nonlocalContribution, nonlocalCumPlastStrain = 0.0;
    TrabBoneNL3DStatus *nonlocStatus, *nlStatus = static_cast< TrabBoneNL3DStatus * >( this->giveStatus(gp) );

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(tStep);

    std :: list< localIntegrationRecord > *list = nlStatus->giveIntegrationDomainList();

    for ( auto &lir: *list ) {
        nonlocStatus = static_cast< TrabBoneNL3DStatus * >( this->giveStatus(lir.nearGp) );
        nonlocalContribution = nonlocStatus->giveLocalCumPlastStrainForAverage();
        nonlocalContribution *= lir.weight;
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
        return static_cast< StructuralNonlocalMaterialExtensionInterface * >(this);
    } else if ( type == NonlocalMaterialStiffnessInterfaceType ) {
        return static_cast< NonlocalMaterialStiffnessInterface * >(this);
    } else {
        return NULL;
    }
}


IRResultType
TrabBoneNL3D :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    result = TrabBone3D :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    result = StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    IR_GIVE_FIELD(ir, R, _IFT_TrabBoneNL3D_r);
    if ( R < 0.0 ) {
        R = 0.0;
    }

    mParam = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, _IFT_TrabBoneNL3D_m);

    return IRRT_OK;
}


void
TrabBoneNL3D :: giveInputRecord(DynamicInputRecord &input)
{
    TrabBone3D :: giveInputRecord(input);
    input.setField(this->R, _IFT_TrabBoneNL3D_r);
    input.setField(this->mParam, _IFT_TrabBoneNL3D_m);
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


/////////////////////////////////////////////////////////////////
/////////TRABECULAR BONE NONLOCAL MATERIAL STATUS////////////////
/////////////////////////////////////////////////////////////////

TrabBoneNL3DStatus :: TrabBoneNL3DStatus(int n, Domain *d, GaussPoint *g) :
    TrabBone3DStatus(n, d, g), StructuralNonlocalMaterialStatusExtensionInterface()
{
    localCumPlastStrainForAverage = 0.0;
}


TrabBoneNL3DStatus :: ~TrabBoneNL3DStatus()
{ }


void
TrabBoneNL3DStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status {");
    fprintf(file, " plastrains ");
    for ( auto &val : plasDef ) {
        fprintf( file, " %.4e", val );
    }

    fprintf(file, " ,");
    fprintf(file, " kappa %.4e ,", kappa);
    fprintf(file, " dam  %.4e ,", tempDam);
    fprintf(file, " esed  %.4e ,", this->tempTSED - this->tempPSED);
    fprintf(file, " psed  %.4e ,", this->tempPSED);
    fprintf(file, " tsed  %.4e", this->tempTSED);
    fprintf(file, "}\n");
}


void
TrabBoneNL3DStatus :: initTempStatus()
{
    TrabBone3DStatus :: initTempStatus();
}


void
TrabBoneNL3DStatus :: updateYourself(TimeStep *tStep)
{
    TrabBone3DStatus :: updateYourself(tStep);
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


int
TrabBoneNL3D :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip)
{
    abort();
    return 0;
 #if 0
    IDNLMaterialStatus *nlStatus = static_cast< IDNLMaterialStatus * >( this->giveStatus(ip) );

    this->buildNonlocalPointTable(ip);
    this->updateDomainBeforeNonlocAverage(tStep);

    return buff.write( nlStatus->giveLocalEquivalentStrainForAverage() );

 #endif
}

int
TrabBoneNL3D :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip)
{
    abort();
    return 0;
 #if 0
    int result;
    IDNLMaterialStatus *nlStatus = static_cast< IDNLMaterialStatus * >( this->giveStatus(ip) );
    double localEquivalentStrainForAverage;

    result = buff.read(localEquivalentStrainForAverage);
    nlStatus->setLocalEquivalentStrainForAverage(localEquivalentStrainForAverage);
    return result;

 #endif
}

int
TrabBoneNL3D :: estimatePackSize(DataStream &buff, GaussPoint *ip)
{
    abort();
    return 0;
 #if 0
    // Note: nlStatus localStrainVectorForAverage memeber must be properly sized!
    // IDNLMaterialStatus *nlStatus = (IDNLMaterialStatus*) this -> giveStatus (ip);
    return buff.givePackSizeOfDouble(1);

 #endif
}

//
// END: PARALLEL MODE OPTION
/////////////////////////////////////////////////////////////////
} // end namespace oofem
