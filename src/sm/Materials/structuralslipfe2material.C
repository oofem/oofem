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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#include "structuralslipfe2material.h"
#include "gausspoint.h"
#include "engngm.h"
#include "oofemtxtdatareader.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "util.h"
#include "contextioerr.h"
#include "generalboundarycondition.h"
#include "prescribedgradienthomogenization.h"
#include "exportmodulemanager.h"
#include "vtkxmlexportmodule.h"
#include "nummet.h"
#include "sm/EngineeringModels/xfemsolverinterface.h"
#include "sm/EngineeringModels/staticstructural.h"
#include "unknownnumberingscheme.h"
#include "xfem/xfemstructuremanager.h"
#include "mathfem.h"

#include "dynamicdatareader.h"

#include <sstream>

namespace oofem {
REGISTER_Material(StructuralSlipFE2Material);

StructuralSlipFE2Material :: StructuralSlipFE2Material(int n, Domain *d) : StructuralFE2Material(n, d)
{}


void
StructuralSlipFE2Material :: initializeFrom(InputRecord &ir)
{
    StructuralFE2Material :: initializeFrom(ir);

    useExtStiff = ir.hasField(_IFT_StructuralSlipFE2Material_useExternalStiffness);
    allGPRes = ir.hasField(_IFT_StructuralSlipFE2Material_allGPResults);
    IR_GIVE_OPTIONAL_FIELD(ir, outputSelected, _IFT_StructuralSlipFE2Material_outputSelectedResults);
    IR_GIVE_OPTIONAL_FIELD(ir, givendStressdEpsTangent, _IFT_StructuralSlipFE2Material_dStressdEps);
}


void
StructuralSlipFE2Material :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralFE2Material :: giveInputRecord(input);

    if ( useExtStiff ) {
        input.setField(_IFT_StructuralSlipFE2Material_useExternalStiffness);
        input.setField(givendStressdEpsTangent, _IFT_StructuralSlipFE2Material_dStressdEps);
    }

    input.setField(allGPRes, _IFT_StructuralSlipFE2Material_allGPResults);
    input.setField(outputSelected, _IFT_StructuralSlipFE2Material_outputSelectedResults);
}


MaterialStatus *
StructuralSlipFE2Material :: CreateStatus(GaussPoint *gp) const
{
    int rank = -1;
    auto emodel = this->domain->giveEngngModel();
    if ( emodel->isParallel() && emodel->giveNumberOfProcesses() > 1 ) {
        rank = emodel->giveRank();
    }

    if ( !(outputSelected.giveSize()==0) ) {
        std::vector<GaussPoint*> gpArray;
        gpArray.resize(outputSelected.giveSize()/2);
        for (int i = 1; i <= outputSelected.giveSize()/2; i++) {
            GaussPoint* gpOut = domain->giveElement(outputSelected.at(i*2-1))->giveIntegrationRule(0)->getIntegrationPoint(outputSelected.at(i*2) - 1);
            gpArray[i - 1] = gpOut;
        }

        if ( std::find(gpArray.begin(), gpArray.end(), gp) != gpArray.end() ) {
            return new StructuralSlipFE2MaterialStatus(rank, gp, this->inputfile);
        } else {
            return new StructuralFE2MaterialStatus(rank, gp, this->inputfile);
        }
    } else if ( allGPRes ) {
        return new StructuralSlipFE2MaterialStatus(rank, gp, this->inputfile);
    } else {
        return new StructuralFE2MaterialStatus(rank, gp, this->inputfile);
    }

}

FloatMatrixF<3, 3>
StructuralSlipFE2Material::givePlaneStressStiffMtrx( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const
{
    if (useExtStiff) {
        auto status = static_cast<StructuralSlipFE2MaterialStatus*>(this->giveStatus(gp));
        FloatMatrixF<3,3> answer;
        answer = this->givendStressdEpsTangent;
        status->setTangent(answer);
        return answer;
    } else {
        return StructuralFE2Material::givePlaneStressStiffMtrx(rMode, gp, tStep);
    }
}

//=============================================================================


StructuralSlipFE2MaterialStatus :: StructuralSlipFE2MaterialStatus( int rank, GaussPoint * g,  const std :: string & inputfile ) :
    StructuralFE2MaterialStatus(rank, g, inputfile)
{
    nel = gp->giveElement()->giveGlobalNumber();
    gpn = gp->giveNumber();
    std :: ostringstream name;
    std :: string oldName = this->rve->giveOutputBaseFileName();
    if ( !oldName.empty() ) {
        name << oldName.substr(0, oldName.size()-4) << "-el" << nel << "-gp" << gpn;
        this->rve->letOutputBaseFileNameBe( name.str() );
    }
}

} // end namespace oofem
