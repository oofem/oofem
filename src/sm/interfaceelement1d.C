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

#include "interfaceelement1d.h"
#include "domain.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "lobattoir.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"

#include "engngm.h"
#ifndef __MAKEDEPEND
 #include <stdlib.h>
 #include <math.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #ifndef __MAKEDEPEND
  #include "Emarkwd3d.h"
 #endif
#endif

namespace oofem {
InterfaceElem1d :: InterfaceElem1d(int n, Domain *aDomain) :
    StructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    referenceNode = 0;
    normal.resize(3);
    normal.zero();
}


void
InterfaceElem1d :: setCoordMode()
{
    switch ( domain->giveDomainType() ) {
    case _2dPlaneStressMode:
        this->mode = ie1d_2d;
        break;
    case _PlaneStrainMode:
        this->mode = ie1d_2d;
        break;
    case _3dMode:
        this->mode = ie1d_3d;
        break;
    case _3dAxisymmMode:
        this->mode = ie1d_3d;
        break;
    case _2dTrussMode:
        this->mode = ie1d_2d;
        break;
    case _1dTrussMode:
        this->mode = ie1d_1d;
        break;
    case _2dBeamMode:
        this->mode = ie1d_2d;
        break;
    default:
        _error("setCoordMode: Unsupported domain type")
    }
}


MaterialMode
InterfaceElem1d :: giveMaterialMode()
{
    setCoordMode();
    switch ( mode ) {
    case ie1d_1d: return _1dInterface;

    case ie1d_2d: return _2dInterface;

    case ie1d_3d: return _3dInterface;

    default: _error("giveMaterialMode: Unsupported coord mode");
    }
    return _1dInterface; // to make the compiler happy
}


void
InterfaceElem1d :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix (which is zero matrix) of the receiver. This expression is
// valid in both local and global axes.
{
    int ndofs = this->computeNumberOfDofs(EID_MomentumBalance);
    answer.resize(ndofs, ndofs);
    answer.zero();
}


void
InterfaceElem1d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    setCoordMode();
    double area = this->giveCrossSection()->give(CS_Area);
    this->computeLocalSlipDir(normal);
    FloatMatrix bLoc, lcs;
    evaluateLocalCoordinateSystem(lcs);
    switch ( mode ) {
    case ie1d_1d:
        bLoc.resize(1, 2);
        bLoc.at(1, 1) = -1.0;
        bLoc.at(1, 2) = +1.0;
        break;
    case ie1d_2d:
        bLoc.resize(2, 4);
        bLoc.zero();
        bLoc.at(1, 1) = -1.0;
        bLoc.at(1, 3) = +1.0;
        bLoc.at(2, 2) = -1.0;
        bLoc.at(2, 4) = +1.0;
        break;
    case ie1d_3d:
        bLoc.resize(3, 6);
        bLoc.zero();
        bLoc.at(1, 1) = -1.0;
        bLoc.at(1, 4) = +1.0;
        bLoc.at(2, 2) = -1.0;
        bLoc.at(2, 5) = +1.0;
        bLoc.at(3, 3) = -1.0;
        bLoc.at(3, 6) = +1.0;
        break;
    default:
        _error("giveDofManDofIDMask: unsupported mode");
    }

    bLoc.times(area);
    answer.beProductOf(lcs, bLoc);
}


void
InterfaceElem1d :: evaluateLocalCoordinateSystem(FloatMatrix &lcs)
//
// Computes unit vectors of local coordinate system, stored by rows.
//
{
    setCoordMode();
    switch ( mode ) {
    case ie1d_1d:
        lcs.resize(1, 1);
        lcs.at(1, 1) = 1.;
        return;

    case ie1d_2d:
        lcs.resize(2, 2);
        lcs.at(1, 1) =  normal.at(1);
        lcs.at(1, 2) =  normal.at(2);
        lcs.at(2, 1) = -normal.at(2);
        lcs.at(2, 2) =  normal.at(1);
        return;

    case ie1d_3d:
    {
        FloatArray ly(3), lz(3);
        normal.normalize();
        ly.zero();
        if ( abs( normal.at(1) ) > abs( normal.at(2) ) ) {
            ly.at(2) = 1.;
        } else {
            ly.at(1) = 1.;
        }

        lz.beVectorProductOf(normal, ly);
        lz.normalize();
        ly.beVectorProductOf(lz, normal);
        ly.normalize();

        lcs.resize(3, 3);
        int i;
        for ( i = 1; i <= 3; i++ ) {
            lcs.at(1, i) = normal.at(i);
            lcs.at(2, i) = ly.at(i);
            lcs.at(3, i) = lz.at(i);
        }

        return;
    }

    default:
        _error("giveDofManDofIDMask: unsupported mode");
    }
}


void
InterfaceElem1d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        //integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 1, this->giveMaterialMode());
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 1, _1dInterface);
    }
}


int
InterfaceElem1d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    answer.resize(3);
    answer.at(1) = this->giveNode(1)->giveCoordinate(1);
    answer.at(2) = this->giveNode(1)->giveCoordinate(2);
    answer.at(3) = this->giveNode(1)->giveCoordinate(3);

    return 1;
}


int
InterfaceElem1d :: computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords)
{
    _error("Not implemented");
    return 0;
}


double
InterfaceElem1d :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    return 1.0;
}


IRResultType
InterfaceElem1d :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom(ir);
    IR_GIVE_OPTIONAL_FIELD(ir, referenceNode, IFT_Beam3d_refnode, "refnode"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, normal, IFT_Node_coords, "normal");
    if ( referenceNode == 0 && normal.at(1) == 0 && normal.at(2) == 0 && normal.at(1) == 0 && normal.at(3) == 0 ) {
        _error("instanciateFrom: wrong reference node or normal specified");
    }

    this->computeGaussPoints();
    return IRRT_OK;
}


int
InterfaceElem1d :: computeNumberOfDofs(EquationID)
{
    setCoordMode();
    switch (mode) {
        case ie1d_1d:
            return 2;
        case ie1d_2d:
            return 4;
        case ie1d_3d:
            return 6;
        default:
            _error("giveDofManDofIDMask: unsupported mode");
    }

    return 0; // to suppress compiler warning
}


void
InterfaceElem1d :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    /*
     * cmode mode = giveCoordMode();
     * if ( mode == ie1d_1d ) {
     * answer.resize(1);
     * answer.at(1) = D_u;
     * } else if ( mode == ie1d_2d ) {
     * answer.resize(2);
     * answer.at(1) = D_u;
     * answer.at(2) = D_v;
     * } else if ( mode == ie1d_3d ) {
     * answer.resize(3);
     * answer.at(1) = D_u;
     * answer.at(2) = D_v;
     * answer.at(3) = D_w;
     * } else {
     */
    switch ( domain->giveDomainType() ) {
    case _2dPlaneStressMode:
    case _PlaneStrainMode:
    case _3dAxisymmMode:
        answer.setValues(2, D_u, D_v);
        break;
    case _3dMode:
        answer.setValues(3, D_u, D_v, D_w);
        break;
    case _2dTrussMode:
    case _2dBeamMode:
        answer.setValues(2, D_u, D_w);
        break;
    case _1dTrussMode:
        answer.setValues(1, D_u);
        break;
    default:
        _error("giveDofManDofIDMask: unsupported mode");
    }
}


void
InterfaceElem1d :: computeLocalSlipDir(FloatArray &normal)
{
    normal.resize(3);
    if ( this->referenceNode ) {
        // tangent
        normal.at(1) = domain->giveNode(this->referenceNode)->giveCoordinate(1) - this->giveNode(1)->giveCoordinate(1);
        normal.at(2) = domain->giveNode(this->referenceNode)->giveCoordinate(2) - this->giveNode(1)->giveCoordinate(2);
        normal.at(3) = domain->giveNode(this->referenceNode)->giveCoordinate(3) - this->giveNode(1)->giveCoordinate(3);
    } else   {
        if ( normal.at(1) == 0 && normal.at(2) == 0 && normal.at(3) == 0 ) {
            _error("computeLocalSlipDir: normal is not defined (referenceNode=0,normal=(0,0,0))");
        }
    }

    normal.normalize();
}


#ifdef __OOFEG
void InterfaceElem1d :: drawRawGeometry(oofegGraphicContext &gc)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 1 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    p [ 0 ].x = ( FPNum )( this->giveNode(1)->giveCoordinate(1) );
    p [ 0 ].y = ( FPNum )( this->giveNode(1)->giveCoordinate(2) );
    p [ 0 ].z = ( FPNum )( this->giveNode(1)->giveCoordinate(3) );

    EASValsSetMType(CIRCLE_MARKER);
    go = CreateMarker3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void InterfaceElem1d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 1 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale) +
                                 this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale) );
    p [ 0 ].y = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale) +
                                 this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale) );
    p [ 0 ].z = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale) +
                                 this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale) );

    EASValsSetMType(CIRCLE_MARKER);
    go = CreateMarker3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void InterfaceElem1d :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray gcoord(3), v1;
    WCRec p [ 1 ];
    IntArray map;
    GraphicObj *go;
    double val [ 1 ];

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( context.getInternalVarsDefGeoFlag() ) {
        double defScale = context.getDefScale();
        p [ 0 ].x = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale) +
                                     this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale) );
        p [ 0 ].y = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale) +
                                     this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale) );
        p [ 0 ].z = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale) +
                                     this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale) );
    } else {
        p [ 0 ].x = ( FPNum )( this->giveNode(1)->giveCoordinate(1) );
        p [ 0 ].y = ( FPNum )( this->giveNode(1)->giveCoordinate(2) );
        p [ 0 ].z = ( FPNum )( this->giveNode(1)->giveCoordinate(3) );
    }

    result += giveIPValue(v1, iRule->getIntegrationPoint(0), context.giveIntVarType(), tStep);


    for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        result = 0;
        gp  = iRule->getIntegrationPoint(i);
        result += giveIPValue(v1, gp, context.giveIntVarType(), tStep);
        result += this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
        if ( result != 2 ) {
            continue;
        }

        if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
            return;
        }

        val [ 0 ] = v1.at(indx);
        context.updateFringeTableMinMax(val, 1);

        EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
        EASValsSetMType(FILLED_CIRCLE_MARKER);
        go = CreateMarkerWD3D(p, val [ 0 ]);
        EGWithMaskChangeAttributes(LAYER_MASK | FILL_MASK | MTYPE_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
        //}
    }
}

#endif
} // end namespace oofem
