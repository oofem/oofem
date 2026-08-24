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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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

#include "sm/Materials/LatticeMaterials/latticestructuralmaterial.h"
#include "sm/Materials/LatticeMaterials/latticematstatus.h"
#include "domain.h"
#include "verbose.h"
#include "sm/Materials/structuralms.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Elements/nlstructuralelement.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "mathfem.h"
#include "engngm.h"
#include "fieldmanager.h"
#include "dynamicinputrecord.h"

namespace oofem {
LatticeStructuralMaterial :: LatticeStructuralMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }


bool
LatticeStructuralMaterial :: hasMaterialModeCapability(MaterialMode mode) const
//
// returns whether receiver supports given mode
//
{
    return mode == _3dLattice || mode == _2dLattice || mode == _1dLattice;
}

void
LatticeStructuralMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
                                                 MatResponseMode rMode,
                                                 GaussPoint *gp, TimeStep *tStep) const
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    answer = this->give3dLatticeStiffnessMatrix(rMode, gp, tStep);
}


int
LatticeStructuralMaterial :: giveIPValue(FloatArray &answer,
                                         GaussPoint *gp,
                                         InternalStateType type,
                                         TimeStep *atTime)
{
    auto status = static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

        if ( type == IST_LatticeForce ) {
            auto help = status->giveLatticeStress();
            answer.resize(3);
            answer.at(1) = help.at(1);
            answer.at(2) = help.at(2);
            answer.at(3) = help.at(3);
            return 1;
        } else if ( type == IST_LatticeMoment )   {
            auto help = status->giveLatticeStress();
            answer.resize(3);
            answer.at(1) = help.at(4);
            answer.at(2) = help.at(5);
            answer.at(3) = help.at(6);
        return 1;
    } else if  ( type == IST_LatticeStrain ) {
            auto help = status->giveLatticeStrain();
            answer.resize(3);
            answer.at(1) = help.at(1);
            answer.at(2) = help.at(2);
            answer.at(3) = help.at(3);
            return 1;
        } else if ( type == IST_LatticeCurvature )   {
            auto help = status->giveLatticeStrain();
            answer.resize(3);
            answer.at(1) = help.at(4);
            answer.at(2) = help.at(5);
            answer.at(3) = help.at(6);
            return 1;
        } else if ( type == IST_PlasticLatticeStrain )   {
            auto help = status->givePlasticLatticeStrain();
            answer.resize(3);
            answer.at(1) = help.at(1);
            answer.at(2) = help.at(2);
            answer.at(3) = help.at(3);
            return 1;
        } else if ( type == IST_PlasticLatticeCurvature )   {
            auto help = status->givePlasticLatticeStrain();
            answer.resize(3);
            answer.at(1) = help.at(4);
            answer.at(2) = help.at(5);
            answer.at(3) = help.at(6);
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, atTime);
    }
}


void
LatticeStructuralMaterial::initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    StructuralMaterial::initializeFrom(ir);

    //temperature threshold used to compute reduction
    this->tCrit = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->tCrit, _IFT_LatticeStructuralMaterial_tcrit); // Macro


}



double
LatticeStructuralMaterial :: giveLatticeStress1d(double strain, GaussPoint *gp, TimeStep *tStep)
{
    FloatArrayF< 6 >tempStrain;
    tempStrain [ 0 ] = strain;
    auto answer = giveLatticeStress3d(tempStrain, gp, tStep);
        return answer [ { 0 } ];
}

FloatArrayF< 3 >
LatticeStructuralMaterial :: giveLatticeStress2d(const FloatArrayF< 3 > &strain, GaussPoint *gp, TimeStep *tStep)
{
    auto answer = giveLatticeStress3d(assemble< 6 >(strain, { 0, 1, 5 }), gp, tStep);
    return answer [ { 0, 1, 5 } ];
}

FloatArrayF< 6 >
LatticeStructuralMaterial :: giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("3dLattice mode not supported");
}

    FloatArrayF < 6 >
    LatticeStructuralMaterial::giveFrameForces3d(const FloatArrayF < 6 > & strain, GaussPoint * gp, TimeStep * tStep)
    {
        OOFEM_ERROR("3dFrame mode not supported");
    }



FloatMatrixF< 1, 1 >
LatticeStructuralMaterial :: give1dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
//
// return material stiffness matrix for 1dlattice
//
{
    OOFEM_ERROR("No general implementation provided");
}

FloatMatrixF< 3, 3 >
LatticeStructuralMaterial :: give2dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
//
// return material stiffness matrix for 2dlattice
//
{
    OOFEM_ERROR("No general implementation provided");
}

FloatMatrixF< 6, 6 >
LatticeStructuralMaterial :: give3dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
//
// return material stiffness matrix for 2dlattice
//
{
    OOFEM_ERROR("No general implementation provided");
}

    FloatMatrixF < 6, 6 >
    LatticeStructuralMaterial::give3dFrameStiffnessMatrix(MatResponseMode mode, GaussPoint * gp, TimeStep * tStep) const
    //
    // return material stiffness matrix for 2dlattice
    //
    {
        OOFEM_ERROR("No general implementation provided");
    }


 double LatticeStructuralMaterial::computeTemperatureReductionFactor(GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
 {
   double reductionFactor = 1.;

    FloatArray et;

    if ( gp->giveIntegrationRule() == NULL ) {
        ///@todo Hack for loose gausspoints. We shouldn't ask for "gp->giveElement()". FIXME
        return reductionFactor;
    }

    Element *elem = gp->giveElement();
    StructuralElement *selem = dynamic_cast< StructuralElement * >( gp->giveElement() );

    if ( tStep->giveIntrinsicTime() < this->castingTime ) {
        return reductionFactor;
    }

    //sum up all prescribed temperatures over an element
    if ( selem ) {
        selem->computeResultingIPTemperatureAt(et, tStep, gp, mode);
    }

    /* add external source, if provided */
    FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    FieldPtr tf = fm->giveField(FT_Temperature);
    if ( tf ) {
        // temperature field registered
        Coordinates gcoords;
        FloatArray et2;
        elem->computeGlobalCoordinates(gcoords, gp->giveNaturalCoordinates() );
        int err;
        if ( ( err = tf->evaluateAt(et2, gcoords, mode, tStep) ) ) {
            OOFEM_ERROR("tf->evaluateAt failed, element %d, error code %d", elem->giveNumber(), err);
        }

        if ( et2.isNotEmpty() ) {
            if ( et.isEmpty() ) {
                et = et2;
            } else {
                et.at(1) += et2.at(1);
            }
        }
    }

    //Compute reductionFactor
    if(et.isNotEmpty()) {
        if ( et.at( 1 ) > this->referenceTemperature && et.at( 1 ) > 0. ) {
            reductionFactor = exp( -pow( ( et.at( 1 ) - this->referenceTemperature ) / ( tCrit - this->referenceTemperature ), 2. ) );
        }
    }
    return reductionFactor;
 }

} // end namespace oofem
