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

#include "../sm/Elements/lumpedmasselement.h"
#include "dofmanager.h"
#include "dof.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "node.h"
#endif

namespace oofem {
REGISTER_Element(LumpedMassElement);

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
    int ndofs = this->computeNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    answer.beDiagonal(this->components);
}


IRResultType
LumpedMassElement :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, components, _IFT_LumpedMassElement_components);

    return StructuralElement :: initializeFrom(ir);
}


int
LumpedMassElement :: checkConsistency()
//
// check internal consistency
//
{
    int _result = StructuralElement :: checkConsistency();
    int _ndofs = this->computeNumberOfDofs();
    if ( _ndofs != this->components.giveSize() ) {
        OOFEM_WARNING("component array size mismatch");
        _result = 0;
    }

    return _result;
}


int
LumpedMassElement :: computeNumberOfDofs()
{
    DofManager *dman = this->giveDofManager(1);
    int answer = 0;

    // simply count all "structural" dofs of element node
    for ( Dof *dof: *dman ) {
        DofIDItem _dofid = dof->giveDofID();
        if ( ( _dofid == D_u ) || ( _dofid == D_v ) || ( _dofid == D_w ) ||
            ( _dofid == R_u ) || ( _dofid == R_v ) || ( _dofid == R_w ) ) {
            answer++;
        }
    }

    return answer;
}


void
LumpedMassElement :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer.resize(6);
    answer.clear();
    DofManager *dman = this->giveDofManager(inode);

    // simply collect all "structural" dofs of element node
    for ( Dof *dof: *dman ) {
        DofIDItem _dofid = dof->giveDofID();
        if ( ( _dofid == D_u ) || ( _dofid == D_v ) || ( _dofid == D_w ) ||
            ( _dofid == R_u ) || ( _dofid == R_v ) || ( _dofid == R_w ) ) {
            answer.followedBy(_dofid);
        }
    }
}

#ifdef __OOFEG
void LumpedMassElement :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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


void LumpedMassElement :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;
    double defScale = gc.getDefScale();
    WCRec p [ 1 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetMType(SQUARE_MARKER);
    EASValsSetMSize(8);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);
    go = CreateMarker3D(p);
    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | MTYPE_MASK | MSIZE_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void LumpedMassElement :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{ }

#endif
} // end namespace oofem
