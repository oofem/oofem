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

#include "homexportmodule.h"
#include "timestep.h"
#include "structuralelement.h"
#include "materialinterface.h"
#include "gausspnt.h"
#ifndef __MAKEDEPEND
 #include <vector>
#endif

namespace oofem {
//inherit LinearElasticMaterial for accessing stress/strain transformation functions
HOMExportModule :: HOMExportModule(int n, EngngModel *e) : ExportModule(n, e), LinearElasticMaterial(0, NULL)
{
    this->matnum.resize(0);
}


HOMExportModule :: ~HOMExportModule()
{ }

IRResultType
HOMExportModule :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    IRResultType val;
    this->scale = 1.;
    ExportModule :: initializeFrom(ir);
    val = IR_GIVE_OPTIONAL_FIELD(ir, this->scale, IFT_HOMExportModule_scale, "scale"); // Macro
    if ( val == IRRT_NOTFOUND ) {
        this->scale = 1.;
    }

    val = IR_GIVE_OPTIONAL_FIELD(ir, this->matnum, IFT_HOMExportModule_matnum, "matnum"); // Macro
    return IRRT_OK;
}


void
HOMExportModule :: doOutput(TimeStep *tStep)
{
    int i, j, index, ielem;
    Element *elem;
    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    Domain *d  = emodel->giveDomain(1);
    int nelem = d->giveNumberOfElements();
    double dV, VolTot = 0.;
    double sumState = 0.;
    FloatArray vecState, vecFlow, sumFlow(3);
    GaussPoint *gp;
    IntegrationRule *iRule;
    IntArray Mask;
    FloatMatrix baseGCS;

    sumFlow.zero();
    domainType domType = d->giveDomainType();

    if ( domType == _HeatTransferMode || domType == _HeatMass1Mode ) {
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            elem = d->giveElement(ielem);
            if ( this->matnum.giveSize() == 0 || this->matnum.contains( elem->giveMaterial()->giveNumber() ) ) {
                iRule = elem->giveDefaultIntegrationRulePtr();
                for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
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
        fprintf( this->stream, "% e  % e  % e\n", sumFlow.at(1), sumFlow.at(2), sumFlow.at(3) );
    } else {
        StructuralElement *structElem;
        //int nnodes = d->giveNumberOfDofManagers();
        FloatArray VecStrain, VecStress, VecEigStrain, SumStrain(6), SumStress(6), SumEigStrain(6), tempFloatAr;
        //stress and strain vectors are always in global c.s.
        SumStrain.zero(); //xx, yy, zz, yz, zx, xy
        SumStress.zero();
        SumEigStrain.zero();

        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            elem = d->giveElement(ielem);
            if ( this->matnum.giveSize() == 0 || this->matnum.contains( elem->giveMaterial()->giveNumber() ) ) {
                iRule = elem->giveDefaultIntegrationRulePtr();
                for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
                    gp  = iRule->getIntegrationPoint(i);
                    structElem = ( StructuralElement * ) gp->giveElement();
                    structElem->computeResultingIPEigenstrainAt(VecEigStrain, tStep, gp, VM_Incremental);
                    //gp->giveCoordinate(1);
                    dV  = elem->computeVolumeAround(gp);
                    //OOFEM_LOG_INFO("Element %d GP %d Vol %f\n", elem->giveNumber(), gp->giveNumber(), dV);
                    //fprintf(this->stream, "Element %d GP %d stress %f\n", elem->giveNumber(), gp->giveNumber(), 0.0);
                    //((StructuralCrossSection*) gp->giveCrossSection())->giveFullCharacteristicVector(helpVec, gp, strainVector);
                    elem->giveIPValue(VecStrain, gp, IST_StrainTensor, tStep);
                    elem->giveIPValue(VecStress, gp, IST_StressTensor, tStep);
                    elem->giveIntVarCompFullIndx(Mask, IST_StrainTensor);

                    //truss element has strains and stresses in the first array so transform them to global coordinates
                    if ( elem->giveClassID() == Truss3dClass ) {
                        elem->giveLocalCoordinateSystem(baseGCS);
                        //this->giveStressVectorTranformationMtrx(transfMatrix,baseGCS,1);//from structuralMaterial
                        //baseGCS.printYourself();
                        VecStress.resize(6);
                        VecStrain.resize(6);
                        VecEigStrain.resize(6);
                        this->transformStressVectorTo(tempFloatAr, baseGCS, VecStress, 0);
                        VecStress = tempFloatAr;
                        this->transformStrainVectorTo(tempFloatAr, baseGCS, VecStrain, 0);
                        VecStrain = tempFloatAr;
                        tempFloatAr = VecEigStrain;
                        VecEigStrain = tempFloatAr;

                        for ( j = 1; j <= 6; j++ ) {
                            Mask.at(j) = j;
                        }

                        //VecStrain.printYourself();
                        //VecEigStrain.printYourself();
                        //VecStress.printYourself();
                    }

                    VolTot += dV;
                    VecStrain.times(dV);
                    VecStress.times(dV);
                    VecEigStrain.times(dV);

                    for ( j = 0; j < 6; j++ ) {
                        index = Mask(j); //indexes in Mask from 1
                        if ( index ) {
                            SumStrain(j) += VecStrain(index - 1);
                            if ( VecEigStrain.giveSize() ) {
                                SumEigStrain(j) += VecEigStrain(index - 1);
                            }

                            SumStress(j) += VecStress(index - 1);
                        }
                    }

                    //VecStrain.printYourself();
                    //VecEigStrain.printYourself();
                    //SumStrain.printYourself();
                }
            }
        }

        //average
        SumStrain.times(1. / VolTot * this->scale);
        SumEigStrain.times(1. / VolTot * this->scale);
        SumStress.times(1. / VolTot * this->scale);
        //SumStrain.printYourself();
        //SumEigStrain.printYourself();
        //SumStress.printYourself();
        fprintf( this->stream, "%f   ", tStep->giveTargetTime() );
        for ( j = 0; j < 6; j++ ) { //strain
            fprintf( this->stream, "%06.5e ", SumStrain(j) );
        }

        fprintf(this->stream, "    ");
        for ( j = 0; j < 6; j++ ) { //stress
            fprintf( this->stream, "%06.5e ", SumStress(j) );
        }

        fprintf(this->stream, "    ");
        for ( j = 0; j < 6; j++ ) { //eigenstrain
            fprintf( this->stream, "%06.5e ", SumEigStrain(j) );
        }

        fprintf(this->stream, "Vol %06.5e ", VolTot);
        fprintf(this->stream, "\n");
    }
    fflush(this->stream);
}

void
HOMExportModule :: initialize()
{
    Domain *d  = emodel->giveDomain(1);
    domainType domType = d->giveDomainType();

    std::string fileName = emodel->giveOutputBaseFileName() + ".hom";
    if ( ( this->stream = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR2("HOMExportModule::giveOutputStream: failed to open file %s", fileName.c_str());
    }

    if ( domType == _HeatTransferMode || domType == _HeatMass1Mode ) {
        fprintf(this->stream, "#Time          AvrState           AvrFlow (xx,yy,zz)\n");
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
