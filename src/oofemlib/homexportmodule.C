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
#include "gausspoint.h"
#include "engngm.h"
#include "material.h"
#include "classfactory.h"
#include "tm/EngineeringModels/transienttransportproblem.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Materials/structuralmaterial.h"

namespace oofem {
REGISTER_ExportModule(HOMExportModule)

HOMExportModule :: HOMExportModule(int n, EngngModel *e) : ExportModule(n, e) { }

HOMExportModule :: ~HOMExportModule() { }

void
HOMExportModule :: initializeFrom(InputRecord &ir)
{
    ExportModule :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, this->ists, _IFT_HOMExportModule_ISTs);
    this->reactions = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->reactions, _IFT_HOMExportModule_reactions);
    this->scale = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->scale, _IFT_HOMExportModule_scale);
    this->strainEnergy = ir.hasField(_IFT_HOMExportModule_strain_energy);
    this->strainEnergySum = 0;
}

void
HOMExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }
    bool outputVol = false;
    double volTot;
    this->stream << std::scientific << tStep->giveTargetTime()*this->timeScale << "   ";
    
    //assemble list of eligible elements. Elements can be present more times in a list but averaging goes just once over each element.
    elements.resize(0);
    for ( int ireg = 1; ireg <= this->giveNumberOfRegions(); ireg++ ) {
        elements.followedBy(this->giveRegionSet(ireg)->giveElementList());
    }
    //elements.printYourself();
    
    if (!ists.isEmpty()) {
        for ( int ist: ists ) {
            FloatArray answer;
            average(answer, volTot, ist, tStep);
            
            if(!outputVol){
                this->stream << std::scientific << volTot << "    ";
                outputVol = true;
            }
            
            if(answer.giveSize() == 0) { //value not defined
                this->stream << "-";
            } else if (answer.giveSize() == 1) {
                this->stream << answer.at(1);
            } else {
                this->stream << answer.giveSize() << " ";
                for ( auto s: answer ) {
                    this->stream << std::scientific << s << " ";
                }
            }
            this->stream << "    ";
        }
    }
    if (reactions) {
       ///@todo need reaction forces from engngm model, not separately from sm/tm modules as it is implemented now
      /*
       TransientTransportProblem *TTP = dynamic_cast< TransientTransportProblem * >(emodel);
       FieldPtr fld = TTP->giveField(FT_HumidityConcentration, tStep);
      */
    
    }
    
    if (strainEnergy) {
        FloatArray averageStress, averageStrain, avSig, dEps;
        double strainEnergyIncr = 0.;
        
        average(averageStress, volTot, IST_StressTensor, tStep);
        average(averageStrain, volTot, IST_StrainTensor, tStep);
        
        if(lastAverageStress.isEmpty()){
            lastAverageStress.resize(averageStress.giveSize());
            lastAverageStress.zero();
            lastAverageStrain.resize(averageStrain.giveSize());
            lastAverageStrain.zero();
        }
        
        //Average stress for integration
        avSig = 0.5*lastAverageStress + 0.5*averageStress;
        dEps.beDifferenceOf(averageStrain, lastAverageStrain);
        strainEnergyIncr = volTot*avSig.dotProduct(dEps);
        
        lastAverageStress = averageStress;
        lastAverageStrain = averageStrain;
        
        strainEnergySum += strainEnergyIncr;
        this->stream << strainEnergyIncr << " " << strainEnergySum << " ";
    }
    
    this->stream << "\n" << std::flush;
}


void 
HOMExportModule :: average(FloatArray &answer, double &volTot, int ist, TimeStep *tStep)
{
       
    FloatArray ipState;
    volTot = 0.;
    Domain *d = emodel->giveDomain(1);
    for ( auto &elem : d->giveElements() ) {
        //printf("%d ", elem -> giveNumber());
        if ( elements.contains(elem -> giveNumber()) ){
            for ( GaussPoint *gp: *elem->giveDefaultIntegrationRulePtr() ) {
                double dV = elem->computeVolumeAround(gp);
                volTot += dV;
                elem->giveGlobalIPValue(ipState, gp, (InternalStateType)ist, tStep);
                answer.add(dV, ipState);
            }
        }
    }

    answer.times( 1. / volTot * this->scale );
}


void
HOMExportModule :: initialize()
{
    char numStr[32];
    sprintf(numStr, "%02d", this->number);
    std::string fileName = emodel->giveOutputBaseFileName() + "." + numStr + ".hom";
    stream = std::ofstream(fileName);
    if ( !stream.good() ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }
    
    this->stream << "#Time          Volume          ";
    for ( int var: this->ists ) {
        this->stream << __InternalStateTypeToString( ( InternalStateType ) var) << "        ";
    }
    
    if(this->strainEnergy){
        this->stream << "StrainEnergyIncr StrainEnergySum      ";
    }
    
    this->stream << "\n" << std::flush;

    initializeElementSet();

    ExportModule :: initialize();
}

void
HOMExportModule :: terminate()
{
    
    if (this->stream){
        this->stream.close();
    }
}
} // end namespace oofem
