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

#include "crackexportmodule.h"
#include "gausspoint.h"
#include "material.h"
#include "element.h"
#include "integrationrule.h"
#include "timestep.h"
#include "engngm.h"
#include "classfactory.h"
#include "Materials/structuralmaterial.h"
#include "crosssection.h"
#include "floatarray.h"

#include <fstream>
#include <sstream>
#include <boost/concept_check.hpp>

namespace oofem {

REGISTER_ExportModule( CrackExportModule )


CrackExportModule :: CrackExportModule(int n, EngngModel *e) : ExportModule(n, e)
{
}


CrackExportModule :: ~CrackExportModule()
{
}


IRResultType
CrackExportModule :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro


    IR_GIVE_FIELD(ir, crossSections, _IFT_CrackExportModule_cs);
    this->threshold = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, threshold, _IFT_CrackExportModule_threshold);

    return ExportModule :: initializeFrom(ir);
}


void
CrackExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    double crackLength_projection;
    double crackLength_sqrt;
    double crackWidth;
    double crackAngle;
    double elementLength;
   
    double damage;
    FloatArray crackVector;
    FloatMatrix princDir;
    FloatMatrix rotMatrix;
    FloatArray strainVector, princStrain;

    double weight, totWeight;

    Domain *d  = emodel->giveDomain(1);

    std::vector< FloatArray > pointsVector;

    for ( auto &elem : d->giveElements() ) {

        int csNumber = elem->giveCrossSection()->giveNumber();
        if ( this->crossSections.containsSorted( csNumber ) ) {

            totWeight = 0.;
            crackWidth = 0.;
            crackLength_projection = 0.;
            crackLength_sqrt = 0.;
            crackAngle = 0.;

            for ( int i = 0; i < elem->giveNumberOfIntegrationRules(); i++ ) {
                IntegrationRule *iRule = elem->giveIntegrationRule(i);

                for ( auto &gp: *iRule ) {

                    FloatArray tmp;
                    elem->giveIPValue(tmp, gp, IST_DamageScalar, tStep);
                    damage = tmp.at(1);
                    elem->giveIPValue(strainVector, gp, IST_StrainTensor, tStep);

                    if ( damage > 0. ) {

                        weight = elem->computeVolumeAround(gp);
                        totWeight += weight;

                        elem->giveIPValue(tmp, gp, IST_CharacteristicLength, tStep);
                        elementLength = tmp.at(1);
                        elem->giveIPValue(crackVector, gp, IST_CrackVector, tStep);
                        crackVector.times(1./damage);

                        princDir.resize(2,2);

                        princDir.at(1,1) = crackVector.at(2);
                        princDir.at(2,1) = -crackVector.at(1);

                        princDir.at(1,2) = crackVector.at(1);
                        princDir.at(2,2) = crackVector.at(2);

                        // modify shear strain in order to allow transformation with the stress transformation matrix
                        strainVector = {strainVector(0), strainVector(1), strainVector(6)*0.5};

                        StructuralMaterial ::  givePlaneStressVectorTranformationMtrx(rotMatrix, princDir, false);
                        princStrain.beProductOf(rotMatrix, strainVector);

                        crackWidth += elementLength * princStrain.at(1) * damage * weight;
                        crackLength_projection += elem->giveCharacteristicSize(gp, crackVector, ECSM_Projection) * weight;
                        crackLength_sqrt += elem->giveCharacteristicSize(gp, crackVector, ECSM_SquareRootOfArea) * weight;

                        if ( crackVector.at(1) != 0. ) {
                            double contrib = atan( crackVector.at(2) / crackVector.at(1) ) * 180./3.1415926;
                            if ( contrib < 0. ) {
                                contrib += 180.;
                            } else if ( contrib > 180. ) {
                                contrib -= 180.;
                            }
                            crackAngle += contrib  * weight;
                        } else {
                            crackAngle += 90. * weight;
                        }

#if 0
                        double length1 = elem->giveCharacteristicSize(gp, crackVector, ECSM_SquareRootOfArea);
                        double length2 = elem->giveCharacteristicSize(gp, crackVector, ECSM_ProjectionCentered);
                        double length3 = elem->giveCharacteristicSize(gp, crackVector, ECSM_Oliver1);
                        double length4 = elem->giveCharacteristicSize(gp, crackVector, ECSM_Oliver1modified);
                        double length5 = elem->giveCharacteristicSize(gp, crackVector, ECSM_Projection);
#endif

                    } else {
                        crackWidth += 0.;
                        crackLength_projection += 0.;
                        crackLength_sqrt += 0.;
                        crackAngle += 0.;
                    }
                }

                if ( totWeight > 0. ) {
                    crackWidth /= totWeight;
                    crackLength_projection /= totWeight;
                    crackLength_sqrt /= totWeight;
                    crackAngle /= totWeight;
                }

                if ( crackWidth >= threshold ) {

                    pointsVector.emplace_back({
                            elem->giveNumber(),
                            crackWidth,
                            crackAngle,
                            crackLength_projection,
                            crackLength_sqrt
                    });
                }
            }
        }
    }
    // vyblejt do outputu
    std :: stringstream strCracks;
    strCracks  << ".dat";
    std :: string nameCracks = this->giveOutputBaseFileName(tStep) + strCracks.str();
    writeToOutputFile(nameCracks, pointsVector);
}


void CrackExportModule :: writeToOutputFile(const std :: string &iName, const std :: vector< FloatArray > &iPoints)
{
    std :: ofstream file;
    file.open( iName.data() );

    // Set some output options
    file << std :: scientific;

    file << "#elem\twidth\tangle\tlength_project\tlength_sqrt\n";

    for ( auto posVec: iPoints ) {
        for ( auto &val : posVec ) {
            file << val << "\t";
        }
        file << "\n";
    }

    file.close();
}


void
CrackExportModule :: initialize()
{
    ExportModule :: initialize();
}


void
CrackExportModule :: terminate()
{}
  
}
