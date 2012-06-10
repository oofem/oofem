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

#include "lumpedmasselement.h"
#include "dofmanager.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
LumpedMassElement :: LumpedMassElement(int n, Domain *aDomain) : StructuralElement(n, aDomain)
// Constructor.
{
    numberOfDofMans = 1;
}


void
LumpedMassElement :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    int _i, _ndofs = this->computeNumberOfDofs(EID_MomentumBalance);
    answer.resize(_ndofs, _ndofs);
    answer.zero();

    for ( _i = 1; _i <= _ndofs; _i++ ) {
        answer.at(_i, _i) = this->components.at(_i);
    }
}


IRResultType
LumpedMassElement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, components, IFT_LumpedMassElement_components, "components"); // Macro

    return IRRT_OK;
}


int
LumpedMassElement :: checkConsistency()
//
// check internal consistency
//
{
    int _result = StructuralElement :: checkConsistency();
    int _ndofs = this->computeNumberOfDofs(EID_MomentumBalance);
    if ( _ndofs != this->components.giveSize() ) {
        _warning("checkConsistency : component array size mismatch");
        _result = 0;
    }

    return _result;
}


int
LumpedMassElement :: computeNumberOfDofs(EquationID ut)
{
    DofManager *dman = this->giveDofManager(1);
    int _i, _ndof = dman->giveNumberOfDofs();
    int answer = 0;
    DofIDItem _dofid;

    // simply count all "structural" dofs of element node
    for ( _i = 1; _i <= _ndof; _i++ ) {
        _dofid = dman->giveDof(_i)->giveDofID();
        if ( ( _dofid == D_u ) || ( _dofid == D_v ) || ( _dofid == D_w ) ||
            ( _dofid == R_u ) || ( _dofid == R_v ) || ( _dofid == R_w ) ) {
            answer++;
        }
    }

    return answer;
}


void
LumpedMassElement :: giveDofManDofIDMask(int inode, EquationID eid, IntArray &answer) const
{
    answer.resize(0, 6);
    DofManager *dman = this->giveDofManager(inode);
    int _i, _ndof = dman->giveNumberOfDofs();
    DofIDItem _dofid;

    // simply collect all "structural" dofs of element node
    for ( _i = 1; _i <= _ndof; _i++ ) {
        _dofid = dman->giveDof(_i)->giveDofID();
        if ( ( _dofid == D_u ) || ( _dofid == D_v ) || ( _dofid == D_w ) ||
            ( _dofid == R_u ) || ( _dofid == R_v ) || ( _dofid == R_w ) ) {
            answer.followedBy(_dofid);
        }
    }
}

#ifdef __OOFEG
void LumpedMassElement :: drawRawGeometry(oofegGraphicContext &gc)
{
    GraphicObj *go;
    WCRec p [ 1 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetMType(SQUARE_MARKER);
    EASValsSetMSize(8);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    go = CreateMarker3D(p);
    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | MTYPE_MASK | MSIZE_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void LumpedMassElement :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();
    WCRec p [ 1 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetMType(SQUARE_MARKER);
    EASValsSetMSize(8);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
    go = CreateMarker3D(p);
    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | MTYPE_MASK | MSIZE_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void LumpedMassElement :: drawScalar(oofegGraphicContext &context)
{}

#endif
} // end namespace oofem
