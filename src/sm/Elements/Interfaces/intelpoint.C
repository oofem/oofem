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

#include "../sm/Elements/Interfaces/intelpoint.h"
#include "../sm/CrossSections/structuralinterfacecrosssection.h"
#include "domain.h"
#include "node.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"

 #include <Emarkwd3d.h>
#endif

namespace oofem {
REGISTER_Element(IntElPoint);

IntElPoint :: IntElPoint(int n, Domain *aDomain) :
    StructuralInterfaceElement(n, aDomain)
{
    numberOfDofMans = 2;
    referenceNode = 0;
    normal.resize(3);
    normal.zero();
    area = 1.0;
}


void
IntElPoint :: setCoordMode()
{
    switch ( domain->giveNumberOfSpatialDimensions() ) {
    case 1:
        this->mode = ie1d_1d;
        break;
    case 2:
        this->mode = ie1d_2d;
        break;
    case 3:
        this->mode = ie1d_3d;
        break;
    default:
        OOFEM_ERROR("Unsupported domain type")
    }
}


MaterialMode
IntElPoint :: giveMaterialMode()
{
    setCoordMode();
    switch ( mode ) {
    case ie1d_1d: return _1dInterface;

    case ie1d_2d: return _2dInterface;

    case ie1d_3d: return _3dInterface;

    default: OOFEM_ERROR("Unsupported coord mode");
    }
    return _1dInterface; // to make the compiler happy
}


void
IntElPoint :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    setCoordMode();
    

    switch ( mode ) {
    case ie1d_1d:
        answer.resize(1, 2);
        answer.at(1, 1) = -1.0;
        answer.at(1, 2) = +1.0;
        break;
    case ie1d_2d:
        answer.resize(2, 4);
        answer.zero();
        answer.at(1, 1) = -1.0;
        answer.at(1, 3) = +1.0;
        answer.at(2, 2) = -1.0;
        answer.at(2, 4) = +1.0;
        break;
    case ie1d_3d:
        answer.resize(3, 6);
        answer.zero();
        answer.at(1, 1) = -1.0;
        answer.at(1, 4) = +1.0;
        answer.at(2, 2) = -1.0;
        answer.at(2, 5) = +1.0;
        answer.at(3, 3) = -1.0;
        answer.at(3, 6) = +1.0;
        break;
    default:
        OOFEM_ERROR("Unsupported mode");
    }

}


void
IntElPoint :: computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Computes transformation matrix to local coordinate system.
    setCoordMode();
    switch ( mode ) {
    case ie1d_1d:
        answer.resize(1, 1);
        answer.at(1, 1) = 1.;
        return;

    case ie1d_2d:
        answer.resize(2, 2);
        answer.at(1, 1) =  normal.at(1);
        answer.at(1, 2) =  normal.at(2);
        answer.at(2, 1) = -normal.at(2);
        answer.at(2, 2) =  normal.at(1);
        return;

    case ie1d_3d:
    {
        //FloatMatrix test;
        //test.beLocalCoordSys(normal.normalize());

        FloatArray ly(3), lz(3);
        normal.normalize();
        ly.zero();
        if ( fabs( normal.at(1) ) > fabs( normal.at(2) ) ) {
            ly.at(2) = 1.;
        } else {
            ly.at(1) = 1.;
        }

        lz.beVectorProductOf(normal, ly);
        lz.normalize();
        ly.beVectorProductOf(lz, normal);
        ly.normalize();

        answer.resize(3, 3);
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(1, i) = normal.at(i);
            answer.at(2, i) = ly.at(i);
            answer.at(3, i) = lz.at(i);
        }

        return;
    }

    default:
        OOFEM_ERROR("Unsupported mode");
    }
}


void
IntElPoint :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        integrationRulesArray [ 0 ]->setUpIntegrationPoints( _Line, 1, this->giveMaterialMode() );
    }
}


int
IntElPoint :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    answer = *this->giveNode(1)->giveCoordinates();
    answer.add(*this->giveNode(2)->giveCoordinates());
    answer.times(0.5);
    return 1;
}



double
IntElPoint :: computeAreaAround(GaussPoint *gp)
{
    // The modeled area/extension around the connected nodes. 
    // Compare with the cs area of a bar. ///@todo replace with cs-property? /JB
    return this->area;
}


IRResultType
IntElPoint :: initializeFrom(InputRecord *ir)
{
    IRResultType result;  // Required by IR_GIVE_FIELD macro

    result = StructuralInterfaceElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    if ( ir->hasField(_IFT_IntElPoint_refnode) &&  ir->hasField(_IFT_IntElPoint_normal) ) {
        OOFEM_WARNING("Ambiguous input: 'refnode' and 'normal' cannot both be specified");
        return IRRT_BAD_FORMAT;
    }
    IR_GIVE_OPTIONAL_FIELD(ir, referenceNode, _IFT_IntElPoint_refnode);
    IR_GIVE_OPTIONAL_FIELD(ir, normal, _IFT_IntElPoint_normal);
    
    this->area = 1.0; // Default area ///@todo Make non-optional? /JB
    IR_GIVE_OPTIONAL_FIELD(ir, this->area, _IFT_IntElPoint_area);
    
    this->computeLocalSlipDir(normal); 
    return IRRT_OK;
}


int
IntElPoint :: computeNumberOfDofs()
{
    setCoordMode();
    switch ( mode ) {
    case ie1d_1d:
        return 2;

    case ie1d_2d:
        return 4;

    case ie1d_3d:
        return 6;

    default:
        OOFEM_ERROR("Unsupported mode");
    }

    return 0; // to suppress compiler warning
}


void
IntElPoint :: giveDofManDofIDMask(int inode, IntArray &answer) const
{

    switch ( domain->giveNumberOfSpatialDimensions() ) {
    case 1:
        answer = IntArray{ D_u };
        break;
    case 2:
        answer = { D_u, D_v };
        break;
    case 3:
        answer = { D_u, D_v, D_w };
        break;
    default:
        OOFEM_ERROR("Unsupported mode");
    }

}


void
IntElPoint :: computeLocalSlipDir(FloatArray &normal)
{
    normal.resizeWithValues(3);
    if ( this->referenceNode ) {
        // normal
        normal.beDifferenceOf(*domain->giveNode(this->referenceNode)->giveCoordinates(), *this->giveNode(1)->giveCoordinates());
    } else {
        if ( normal.at(1) == 0 && normal.at(2) == 0 && normal.at(3) == 0 ) {
            OOFEM_ERROR("Normal is not defined (referenceNode=0,normal=(0,0,0))");
        }
    }

    normal.normalize();
}



#ifdef __OOFEG
void IntElPoint :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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
    p [ 0 ].x = ( FPNum ) ( this->giveNode(1)->giveCoordinate(1) );
    p [ 0 ].y = ( FPNum ) ( this->giveNode(1)->giveCoordinate(2) );
    p [ 0 ].z = ( FPNum ) ( this->giveNode(1)->giveCoordinate(3) );

    EASValsSetMType(CIRCLE_MARKER);
    go = CreateMarker3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void IntElPoint :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 1 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    double defScale = gc.getDefScale();

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale) +
                                 this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale) );
    p [ 0 ].y = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale) +
                                 this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale) );
    p [ 0 ].z = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale) +
                                 this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale) );

    EASValsSetMType(CIRCLE_MARKER);
    go = CreateMarker3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void IntElPoint :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int indx, result = 0;
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    FloatArray gcoord(3), v1;
    WCRec p [ 1 ];
    GraphicObj *go;
    double val [ 1 ];

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.getInternalVarsDefGeoFlag() ) {
        double defScale = gc.getDefScale();
        p [ 0 ].x = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale) +
                                     this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale) );
        p [ 0 ].y = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale) +
                                     this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale) );
        p [ 0 ].z = ( FPNum ) 0.5 * ( this->giveNode(1)->giveUpdatedCoordinate(3, tStep,  defScale) +
                                     this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale) );
    } else {
        p [ 0 ].x = ( FPNum ) ( this->giveNode(1)->giveCoordinate(1) );
        p [ 0 ].y = ( FPNum ) ( this->giveNode(1)->giveCoordinate(2) );
        p [ 0 ].z = ( FPNum ) ( this->giveNode(1)->giveCoordinate(3) );
    }

    result += giveIPValue(v1, iRule->getIntegrationPoint(0), gc.giveIntVarType(), tStep);


    for ( GaussPoint *gp: *iRule ) {
        result = 0;
        result += giveIPValue(v1, gp, gc.giveIntVarType(), tStep);
        if ( result != 1 ) {
            continue;
        }

        indx = gc.giveIntVarIndx();

        val [ 0 ] = v1.at(indx);
        gc.updateFringeTableMinMax(val, 1);

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
