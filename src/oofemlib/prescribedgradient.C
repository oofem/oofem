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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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

#include "prescribedgradient.h"
#include "dofiditem.h"
#include "dofmanager.h"
#include "dof.h"
#include "valuemodetype.h"
#include "classtype.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {

double PrescribedGradient :: give(Dof *dof, ValueModeType mode, TimeStep *tStep)
{
    DofIDItem id = dof->giveDofID();
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();

    if ( coords == NULL || coords->giveSize() != this->centerCoord.giveSize() ) {
        OOFEM_ERROR("give: Size of coordinatesystem different from center coordinate in b.c.");
    }

    // Reminder: u_i = d_ij . (x_j - xb_j) = d_ij . dx_j
    FloatArray dx(*coords);
    dx -= this->centerCoord;

    FloatArray u;
    u.beProductOf(gradient, dx);

    switch ( id ) {
    case D_u:
    case V_u:
    case P_f:
    case T_f:
        return u.at(1);

    case D_v:
    case V_v:
        if ( u.giveSize() >= 2 ) {
            return u.at(2);
        } else {
            OOFEM_ERROR("give: Prescribed tensor dimensions to small for D_v or V_v.");
        }

    case D_w:
    case V_w:
        if ( u.giveSize() >= 3 ) {
            return u.at(3);
        } else {
            OOFEM_ERROR("give: Prescribed tensor dimensions to small for D_w or V_w.");
        }

    default:
        return 0.0;
    }
}

void PrescribedGradient :: setPrescribedGradientVoigt(const FloatArray &t)
{
    int n = t.giveSize();
    if ( n == 3 ) { // Then 2D
        this->gradient.resize(2, 2, 0);
        this->gradient.at(1, 1) = t.at(1);
        this->gradient.at(2, 2) = t.at(2);
        this->gradient.at(1, 2) = this->gradient.at(2, 1) = t.at(3);
    } else if ( n == 6 )     { // Then 3D
        OOFEM_ERROR("setPrescribedTensorVoigt: Check the order of voigt vectors in OOFEM");
        this->gradient.resize(2, 2, 0);
        this->gradient.at(1, 1) = t.at(1);
        this->gradient.at(2, 2) = t.at(2);
        this->gradient.at(2, 2) = t.at(3);
        // In voigt form, assuming the use of gamma_12 instead of eps_12
        this->gradient.at(1, 2) = this->gradient.at(2, 1) = t.at(6) / 2;
        this->gradient.at(1, 3) = this->gradient.at(3, 1) = t.at(5) / 2;
        this->gradient.at(2, 3) = this->gradient.at(3, 2) = t.at(4) / 2;
    } else   {
        OOFEM_ERROR("setPrescribedTensorVoigt: Tensor is in strange voigt format. Should be 3 or 6. Use setPrescribedTensor directly if needed.");
    }
}

IRResultType PrescribedGradient :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    GeneralBoundaryCondition :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->gradient, IFT_PrescribedTensor_gradient, "gradient");
    IRResultType rt = IR_GIVE_OPTIONAL_FIELD(ir, this->centerCoord, IFT_PrescribedTensor_centercoords, "ccoord")
    if (rt == IRRT_OK) {
        this->centerCoord.resize(this->gradient.giveNumberOfColumns());
        this->centerCoord.zero();
    }

    return IRRT_OK;
}

int PrescribedGradient :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    GeneralBoundaryCondition :: giveInputRecordString(str, keyword);

    sprintf( buff, " gradient %d %d ", this->gradient.giveNumberOfRows(), this->gradient.giveNumberOfColumns() );
    for ( int i = 1; i <= this->gradient.giveNumberOfColumns(); i++ ) {
        for ( int j = 1; j <= this->gradient.giveNumberOfRows(); j++ ) {
            sprintf( buff, " %e", this->gradient.at(i, j) );
            str += buff;
        }
    }

    sprintf( buff, " ccord %d", this->centerCoord.giveSize() );
    for ( int i = 1; i <= this->centerCoord.giveSize(); i++ ) {
        sprintf( buff, " %e", this->centerCoord.at(i) );
        str += buff;
    }

    return 1;
}
} // end namespace oofem

