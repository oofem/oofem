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

#include "Elements/Axisymmetry/q4axisymm.h"
#include "fei2dquadquad.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#include "mathfem.h"
#include "crosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Element(Q4Axisymm);

FEI2dQuadQuad Q4Axisymm :: interp(1, 2);

Q4Axisymm :: Q4Axisymm(int n, Domain *aDomain) :
    AxisymElement(n, aDomain), ZZNodalRecoveryModelInterface(this)
{
    numberOfDofMans = 8;
    numberOfGaussPoints = 4;
    numberOfFiAndShGaussPoints = 1;
}


Q4Axisymm :: ~Q4Axisymm()
{ }


FEInterpolation *
Q4Axisymm :: giveInterpolation() const
{
    return & interp;
}



IRResultType
Q4Axisymm :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    numberOfGaussPoints = 4;
    result = StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    numberOfFiAndShGaussPoints = 1;
    ///@todo only works for 1 //JB
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfFiAndShGaussPoints, _IFT_Q4Axisymm_nipfish);

    return result;
}


void
Q4Axisymm :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    // Returns the [ 6 x (nno*2) ] strain-displacement matrix {B} of the receiver,
    // evaluated at gp. Uses reduced integration
    // (epsilon_x,epsilon_y,...,Gamma_xy) = B . r
    // r = ( u1,v1,u2,v2,u3,v3,u4,v4)

    if ( numberOfFiAndShGaussPoints == 1 ) { // Reduced integration
        FEInterpolation *interp = this->giveInterpolation();
        
        FloatArray N, NRed, redCoord = {0.0, 0.0}; // eval in centroid
        interp->evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        interp->evalN( NRed, redCoord, FEIElementGeometryWrapper(this) );
        
        // Evaluate radius at center
        double r = 0.0;
        for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
            double x = this->giveNode(i)->giveCoordinate(1);
            r += x * NRed.at(i);
        } 
        
        FloatMatrix dNdx, dNdxRed;
        interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        interp->evaldNdx( dNdxRed, redCoord, FEIElementGeometryWrapper(this) );
        answer.resize(6, dNdx.giveNumberOfRows() * 2);
        answer.zero();

        for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
            answer.at(1, i * 2 - 1) = dNdx.at(i, 1);
            answer.at(2, i * 2 - 0) = dNdx.at(i, 2);
            answer.at(3, i * 2 - 1) = NRed.at(i) / r;
            answer.at(6, 2 * i - 1) = dNdxRed.at(i, 2);
            answer.at(6, 2 * i - 0) = dNdxRed.at(i, 1);
        }
    }
}





Interface *
Q4Axisymm :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    }

    return NULL;
}


#ifdef __OOFEG

void Q4Axisymm :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 4 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);

    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Q4Axisymm :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 4 ];
    GraphicObj *go;
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}

#endif
} // end namespace oofem
