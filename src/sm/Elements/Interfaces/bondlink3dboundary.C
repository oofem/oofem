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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "domain.h"
#include "../sm/Elements/Interfaces/bondlink3dboundary.h"
#include "../sm/Materials/structuralms.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "../sm/Elements/structuralelement.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "../sm/Materials/structuralmaterial.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "parametermanager.h"
#include "paramkey.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(BondLink3dBoundary);
ParamKey BondLink3dBoundary::IPK_BondLink3dBoundary_location("location");

BondLink3dBoundary :: BondLink3dBoundary(int n, Domain *aDomain) : BondLink3d(n, aDomain), location(2)
{
    numberOfDofMans     = 3;
    geometryFlag = 0;
}

BondLink3dBoundary :: ~BondLink3dBoundary()
{}

void
BondLink3dBoundary :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                             TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    FloatMatrix d, b, bt, db;
    FloatArray u, slip;

    this->computeVectorOf(VM_Total, tStep, u);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    answer.clear();

    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    this->computeBmatrixAt(gp, b);
    bt.beTranspositionOf(b);

    if ( !this->isActivated(tStep) ) {
        slip.resize(StructuralMaterial :: giveSizeOfVoigtSymVector(gp->giveMaterialMode() ) );
        slip.zero();
    }
    slip.beProductOf(b, u);

    answer.resize(9, 9);
    answer.zero();

    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    db.beProductOf(d, b);
    answer.beProductOf(bt, db);

    //Introduce integration of bond strength
    double area = this->computeVolumeAround(gp) / this->giveLength();
    answer.times(area);


    return;
}

void
BondLink3dBoundary :: computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    FloatArray unitCellSize;
    unitCellSize.resize(3);
    unitCellSize.at(1) = this->giveNode(3)->giveCoordinate(1);
    unitCellSize.at(2) = this->giveNode(3)->giveCoordinate(2);
    unitCellSize.at(3) = this->giveNode(3)->giveCoordinate(3);

    IntArray switches1, switches2;
    this->giveSwitches(switches1, this->location.at(1) );
    this->giveSwitches(switches2, this->location.at(2) );

    FloatMatrix k1, k2;
    k1.resize(6, 9);
    k2.resize(6, 9);

    for ( int i = 1; i <= 3; i++ ) {
        k1.at(i, 3 * i - 2) = unitCellSize.at(1) * switches1.at(1);
        k1.at(i, 3 * i - 1) = unitCellSize.at(2) * switches1.at(2);
        k1.at(i, 3 * i) = unitCellSize.at(3) * switches1.at(3);
    }

    for ( int i = 1; i <= 3; i++ ) {
        k2.at(i, 3 * i - 2) = unitCellSize.at(1) * switches2.at(1);
        k2.at(i, 3 * i - 1) = unitCellSize.at(2) * switches2.at(2);
        k2.at(i, 3 * i) = unitCellSize.at(3) * switches2.at(3);
    }

    answer.resize(12, 12);
    answer.beUnitMatrix();
    answer.resizeWithData(12, 21);

    answer.assemble(k1, { 1, 2, 3, 4, 5, 6 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21 });
    answer.assemble(k2, { 7, 8, 9, 10, 11, 12 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21 });
}


bool
BondLink3dBoundary :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    int i, j;

    answer.resize(9, 9);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
            answer.at(i + 6, j + 6) = lcs.at(i, j);
            //            answer.at(i + 9, j + 9) = lcs.at(i, j);
        }
    }

    return 1;
}

int
BondLink3dBoundary :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->localCoordinateSystem;

    return 1;
}



void
BondLink3dBoundary ::   giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode == 1 ) {
        answer = { D_u, D_v, D_w, R_u, R_v, R_w };
    } else if ( inode == 2 )        {
        answer = { D_u, D_v, D_w };
    } else if ( inode == 3 )        {
        answer = { E_xx, E_xy, E_xz, E_yx, E_yy, E_yz, E_zx, E_zy, E_zz };
    }
}

void
BondLink3dBoundary :: initializeFrom(InputRecord &ir, int priority)
{
    // first call parent
    StructuralElement :: initializeFrom(ir, priority);
    ParameterManager &ppm = giveDomain()->elementPPM;
    PM_UPDATE_PARAMETER(location, ppm, ir, this->number, IPK_BondLink3dBoundary_location, priority);
}



void
BondLink3dBoundary :: computeGeometryProperties()
{
    //coordinates of the two nodes
    Node *nodeA, *nodeB;
    FloatArray coordsA(3), coordsB(3);

    //Order of nodes. Important, because continuum node does not have rotational DOFs.
    //Beam node
    nodeA  = this->giveNode(1);
    //Continuum node
    nodeB  = this->giveNode(2);


    //Calculate components of distance from continuum node to lattice node.
    for ( int i = 0; i < 3; i++ ) {
        coordsA.at(i + 1) =  nodeA->giveCoordinate(i + 1);
        coordsB.at(i + 1) =  nodeB->giveCoordinate(i + 1);
    }

    FloatArray rigidGlobal(3);

    //Calculate normal vector
    for ( int i = 0; i < 3; i++ ) {
        rigidGlobal.at(i + 1) = coordsA.at(i + 1) - coordsB.at(i + 1);
    }

    //Construct an initial temporary local coordinate system
    FloatArray normal(3), s(3), t(3);

    //Calculate normal vector
    normal = this->directionVector;
    normal.normalize();

    //Construct two perpendicular axis so that n is normal to the plane which they create
    //Check, if one of the components of the normal-direction is zero
    if ( normal.at(1) == 0 ) {
        s.at(1) = 0.;
        s.at(2) = normal.at(3);
        s.at(3) = -normal.at(2);
    } else if ( normal.at(2) == 0 ) {
        s.at(1) = normal.at(3);
        s.at(2) = 0.;
        s.at(3) = -normal.at(1);
    } else {
        s.at(1) = normal.at(2);
        s.at(2) = -normal.at(1);
        s.at(3) = 0.;
    }

    s.normalize();

    t.beVectorProductOf(normal, s);
    t.normalize();

    //Set up rotation matrix
    FloatMatrix lcs(3, 3);

    this->localCoordinateSystem.resize(3, 3);
    this->localCoordinateSystem.zero();

    for ( int i = 1; i <= 3; i++ ) {
        this->localCoordinateSystem.at(1, i) = normal.at(i);
        this->localCoordinateSystem.at(2, i) = s.at(i);
        this->localCoordinateSystem.at(3, i) = t.at(i);
    }

    // Rotate rigidarm vector into local coordinate system

    this->rigid.beProductOf(localCoordinateSystem, rigidGlobal);


    this->globalCentroid.resize(3);
    for ( int i = 1; i <= 3; i++ ) {
        this->globalCentroid.at(i) = nodeB->giveCoordinate(i);
        ;
    }

    this->geometryFlag = 1;

    return;
}

void
BondLink3dBoundary :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralElement :: saveContext(stream, mode);
}


void
BondLink3dBoundary :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralElement :: restoreContext(stream, mode);
}


void
BondLink3dBoundary :: giveInternalForcesVector(FloatArray &answer,
                                               TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix b, bt;
    FloatArray u, stress(6), slip(6);

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    for ( GaussPoint *gp: * this->giveDefaultIntegrationRulePtr() ) {
        StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
        this->computeBmatrixAt(gp, b);
        bt.beTranspositionOf(b);

        if ( useUpdatedGpRecord == 1 ) {
            stress = matStat->giveStressVector();
        } else {
            if ( !this->isActivated(tStep) ) {
                slip.resize(6);
                slip.zero();
            }
            slip.beProductOf(b, u);
            this->computeStressVector(stress, slip, gp, tStep);
        }

        answer.beProductOf(bt, stress);

        //Introduce integration of bond strength
        double area = this->computeVolumeAround(gp) / this->giveLength();
        answer.times(area);
    }

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}



void
BondLink3dBoundary :: giveSwitches(IntArray &answer, int location) {
    int counter = 1;
    for ( int x = -1; x <  2; x++ ) {
        for ( int y = -1; y <  2; y++ ) {
            for ( int z = -1; z <  2; z++ ) {
                if ( !( z == 0 && y == 0 && x == 0 ) ) {
                    if ( counter == location ) {
                        answer(0) = x;
                        answer(1) = y;
                        answer(2) = z;
                    }
                    counter++;
                }
            }
        }
    }
    return;
}


void
BondLink3dBoundary :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveCharMaterialStiffnessMatrix(answer, rMode, gp, tStep);
}

void
BondLink3dBoundary :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveRealStresses(strain, gp, tStep);
}



#ifdef __OOFEG

void
BondLink3dBoundary :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc, tStep);
    } else if ( mode == OGC_deformedGeometry ) {
        this->drawDeformedGeometry(gc, tStep, DisplacementVector);
    } else if ( mode == OGC_eigenVectorGeometry ) {
        this->drawDeformedGeometry(gc, tStep, EigenVector);
    } else if ( mode == OGC_scalarPlot ) {
        this->drawScalar(gc, tStep);
    } else if ( mode == OGC_elemSpecial ) {
        this->drawSpecial(gc, tStep);
    } else {
        OOFEM_ERROR("unsupported mode");
    }
}



void BondLink3dBoundary :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    WCRec p [ 2 ]; /* points */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void BondLink3dBoundary :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    double defScale = gc.getDefScale();

    WCRec p [ 2 ]; /* points */

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}

#endif
} // end namespace oofem
