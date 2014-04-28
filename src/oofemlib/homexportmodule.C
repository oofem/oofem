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

#include "homexportmodule.h"
#include "timestep.h"
#include "element.h"
#include "femcmpnn.h"
#include "materialinterface.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "crosssection.h"

#ifdef __SM_MODULE
 #include "../sm/structuralelement.h"
 #include "../sm/structuralmaterial.h"
 #include "../sm/truss3d.h"
 #include "../sm/stressvector.h"
 #include "../sm/strainvector.h"
#endif

#ifdef __TM_MODULE
 #include "../tm/transportmaterial.h"
#endif


namespace oofem {
REGISTER_ExportModule(HOMExportModule)

//inherit LinearElasticMaterial for accessing stress/strain transformation functions
HOMExportModule :: HOMExportModule(int n, EngngModel *e) : ExportModule(n, e)
{
    this->matnum.clear();
    this->internalSourceEnergy.clear();
    this->capacityEnergy.clear();
}


HOMExportModule :: ~HOMExportModule()
{ }

IRResultType
HOMExportModule :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    IRResultType val;
    this->scale = 1.;
    ExportModule :: initializeFrom(ir);
    val = IR_GIVE_OPTIONAL_FIELD(ir, this->scale, _IFT_HOMExportModule_scale);
    if ( val == IRRT_NOTFOUND ) {
        this->scale = 1.;
    }

    val = IR_GIVE_OPTIONAL_FIELD(ir, this->matnum, _IFT_HOMExportModule_matnum);
    internalSourceEnergy.resize( emodel->giveDomain(1)->giveNumberOfCrossSectionModels() );
    return IRRT_OK;
}


void
HOMExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    Element *elem;
    TransportElement *transpElem;
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    Domain *d  = emodel->giveDomain(1);
    int nelem = d->giveNumberOfElements();
    double dV, VolTot = 0.;
    double sumState = 0.;
    FloatArray vecState, vecFlow, sumFlow(3);
    FloatArray internalSource, capacity;
    FloatArray answerArr, answerArr1;
    FloatMatrix answerMtrx;
    GaussPoint *gp;
    IntegrationRule *iRule;
    IntArray Mask;
    FloatMatrix baseGCS;

    sumFlow.zero();
    domainType domType = d->giveDomainType();

    if ( domType == _HeatTransferMode || domType == _HeatMass1Mode ) {
        for ( int ielem = 1; ielem <= nelem; ielem++ ) {
            elem = d->giveElement(ielem);
            if ( this->matnum.giveSize() == 0 || this->matnum.contains( elem->giveMaterial()->giveNumber() ) ) {
                iRule = elem->giveDefaultIntegrationRulePtr();
                for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
                    gp  = iRule->getIntegrationPoint(i);
                    dV  = elem->computeVolumeAround(gp);
                    VolTot += dV;
                    elem->giveIPValue(vecState, gp, IST_Temperature, tStep);
                    elem->giveIPValue(vecFlow, gp, IST_TemperatureFlow, tStep);
                    sumState += vecState.at(1) * dV;
                    vecFlow.resize(3);
                    vecFlow.times(dV);
                    sumFlow += vecFlow;
                }
            }
        }

        sumState *= ( 1. / VolTot * this->scale );
        fprintf(this->stream, "%e  % e       ", tStep->giveTargetTime(), sumState);
        sumFlow.times(1. / VolTot * this->scale);
        fprintf( this->stream, "% e  % e  % e       ", sumFlow.at(1), sumFlow.at(2), sumFlow.at(3) );
        //add total heat for each material - accumulated energy due to material capacity (J for heat) and internal source (J for heat)
        internalSource.resize( d->giveNumberOfCrossSectionModels() );
        capacity.resize( d->giveNumberOfCrossSectionModels() );
        for ( int ielem = 1; ielem <= nelem; ielem++ ) {
            transpElem = static_cast< TransportElement* >( d->giveElement(ielem) );
            transpElem->computeInternalSourceRhsVectorAt(answerArr, tStep, VM_Total);
            internalSource.at( transpElem->giveCrossSection()->giveNumber() ) += answerArr.sum();
            
            transpElem->giveCharacteristicMatrix(answerMtrx,CapacityMatrix,tStep);
            transpElem->computeVectorOf(EID_ConservationEquation, VM_Incremental, tStep, answerArr);
            answerArr1.beProductOf(answerMtrx,answerArr);
            capacity.at( transpElem->giveCrossSection()->giveNumber() ) -= answerArr1.sum();
        }
        internalSource.times(tStep->giveTimeIncrement());
        internalSourceEnergy.add(internalSource);
        capacityEnergy.add(capacity);
        for ( int i = 1; i <= internalSourceEnergy.giveSize(); i++ ) {
            fprintf( this->stream, "% e ", internalSourceEnergy.at(i) );
        }
        fprintf( this->stream, "     " );
        for ( int i = 1; i <= capacityEnergy.giveSize(); i++ ) {
            fprintf( this->stream, "% e ", capacityEnergy.at(i) );
        }
        fprintf( this->stream, "\n" );
    } else { //structural analysis
#ifdef __SM_MODULE
        StructuralElement *structElem;
        //int nnodes = d->giveNumberOfDofManagers();
        FloatArray VecStrain, VecStress, VecEigStrain, VecEigStrainReduced, tempStrain(6), tempStress(6), tempEigStrain(6), SumStrain(6), SumStress(6), SumEigStrain(6), tempFloatAr, damage;
        double sumDamage = 0.;
        //stress and strain vectors are always in global c.s.
        SumStrain.zero(); //xx, yy, zz, yz, zx, xy
        SumStress.zero();
        SumEigStrain.zero();

        for ( int ielem = 1; ielem <= nelem; ielem++ ) {
            tempStrain.zero();
            tempStress.zero();
            tempEigStrain.zero();

            elem = d->giveElement(ielem);
            if ( this->matnum.giveSize() == 0 || this->matnum.contains( elem->giveMaterial()->giveNumber() ) ) {
                iRule = elem->giveDefaultIntegrationRulePtr();
                for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
                    gp  = iRule->getIntegrationPoint(i);
                    structElem = static_cast< StructuralElement * >(elem);
                    // structElem->computeResultingIPEigenstrainAt(VecEigStrain, tStep, gp, VM_Incremental);
                    structElem->computeResultingIPEigenstrainAt(VecEigStrainReduced, tStep, gp, VM_Total);
		    if( VecEigStrainReduced.giveSize() == 0 ){
                        VecEigStrain.resize(0);
                    } else {
                        ((StructuralMaterial*)structElem->giveMaterial())->giveFullSymVectorForm(VecEigStrain, VecEigStrainReduced,  gp->giveMaterialMode());
                    }
                    dV  = elem->computeVolumeAround(gp);
                    elem->giveIPValue(VecStrain, gp, IST_StrainTensor, tStep);
                    elem->giveIPValue(VecStress, gp, IST_StressTensor, tStep);
                    elem->giveIPValue(damage, gp, IST_DamageTensor, tStep);

                    //truss element has strains and stresses in the first array so transform them to global coordinates
                    ///@todo Should this be the job of giveIPValue to ensure that vectors are given in global c.s.?
                    if ( dynamic_cast< Truss3d * >(elem) ) {
                        MaterialMode mmode = _3dMat;
                        tempStress.at(1) = VecStress.at(1);
                        tempStrain.at(1) = VecStrain.at(1);
                        if ( VecEigStrain.giveSize() ) {
                            tempEigStrain.at(1) = VecEigStrain.at(1);
                        }
                        StressVector stressVector1(tempStress, mmode); //convert from array
                        StressVector stressVector2(mmode);
                        StrainVector strainVector1(tempStrain, mmode); //convert from array
                        StrainVector strainVector2(mmode);
                        StrainVector strainEigVector1(tempEigStrain, mmode); //convert from array
                        StrainVector strainEigVector2(mmode);
                        elem->giveLocalCoordinateSystem(baseGCS);

                        stressVector1.transformTo(stressVector2, baseGCS, 0);
                        strainVector1.transformTo(strainVector2, baseGCS, 0);
                        strainEigVector1.transformTo(strainEigVector2, baseGCS, 0);

                        stressVector2.convertToFullForm(VecStress);
                        strainVector2.convertToFullForm(VecStrain);
                        strainEigVector2.convertToFullForm(VecEigStrain);
                    }

                    VolTot += dV;
                    VecStrain.times(dV);
                    VecStress.times(dV);
                    VecEigStrain.times(dV);
                    if ( damage.giveSize() ) { //only for models with damage
                        sumDamage += damage.at(1) * dV;
                    }

                    for ( int j = 1; j <= 6; j++ ) {
                        SumStrain.at(j) += VecStrain.at(j);
                        if ( VecEigStrain.giveSize() ) {
                            SumEigStrain.at(j) += VecEigStrain.at(j);
                        }

                        SumStress.at(j) += VecStress.at(j);
                    }

                    //VecStrain.printYourself();
                    //VecEigStrain.printYourself();
                    //SumStrain.printYourself();
                }
            }
        }

        //averaging
        SumStrain.times(1. / VolTot * this->scale);
        SumEigStrain.times(1. / VolTot * this->scale);
        SumStress.times(1. / VolTot * this->scale);
        sumDamage *= ( 1. / VolTot );
        //SumStrain.printYourself();
        //SumEigStrain.printYourself();
        //SumStress.printYourself();
        fprintf( this->stream, "%f   ", tStep->giveTargetTime() );
        for ( int j = 0; j < 6; j++ ) { //strain
            fprintf( this->stream, "%06.5e ", SumStrain(j) );
        }

        fprintf(this->stream, "    ");
        for ( int j = 0; j < 6; j++ ) { //stress
            fprintf( this->stream, "%06.5e ", SumStress(j) );
        }

        fprintf(this->stream, "    ");
        for ( int j = 0; j < 6; j++ ) { //eigenstrain
            fprintf( this->stream, "%06.5e ", SumEigStrain(j) );
        }

        fprintf(this->stream, "Vol %06.5e ", VolTot);
        fprintf(this->stream, "AvrDamage %06.5e ", sumDamage);
        fprintf(this->stream, "\n");

#endif
    }
    fflush(this->stream);
}

void
HOMExportModule :: initialize()
{
    Domain *d  = emodel->giveDomain(1);
    domainType domType = d->giveDomainType();

    std :: string fileName = emodel->giveOutputBaseFileName() + ".hom";
    if ( ( this->stream = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str());
    }

    if ( domType == _HeatTransferMode || domType == _HeatMass1Mode ) {
        fprintf(this->stream, "#Time          AvrState           AvrFlow (xx,yy,zz)                               internalSourceEnergy(over_each_crosssection)   capacityEnergy(over_each_crosssection)\n");
    } else {
        fprintf(this->stream, "#Time   AvrStrain (xx, yy, zz, yz, zx, xy)                                           AvrStress (xx, yy, zz, yz, zx, xy)                                             AvrEigenstrain (xx, yy, zz, yz, zx, xy)\n");
    }
    fflush(this->stream);
}


void
HOMExportModule :: terminate()
{
    fclose(this->stream);
}
} // end namespace oofem
