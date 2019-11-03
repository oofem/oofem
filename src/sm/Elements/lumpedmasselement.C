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

#include "sm/Elements/lumpedmasselement.h"
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
{
    numberOfDofMans = 1;
}


void
LumpedMassElement :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    int ndofs = this->computeNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    answer.beDiagonal(this->components);
}


void
LumpedMassElement :: initializeFrom(InputRecord &ir)
{
    StructuralElement :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, dofs, _IFT_LumpedMassElement_dofs);
    IR_GIVE_FIELD(ir, components, _IFT_LumpedMassElement_components);
}


int
LumpedMassElement :: checkConsistency()
{
    int _result = StructuralElement :: checkConsistency();
    if ( this->dofs.giveSize() != this->components.giveSize() ) {
        OOFEM_WARNING("component/dof array size mismatch");
        _result = 0;
    }

    return _result;
}


int
LumpedMassElement :: computeNumberOfDofs()
{
    return dofs.giveSize();
}


void
LumpedMassElement :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = dofs;
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
