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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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
#include "../sm/Elements/LatticeElements/latticebeam3dboundary.h"
#include "../sm/Materials/LatticeMaterials/latticematstatus.h"
#include "../sm/Materials/LatticeMaterials/latticelinearelastic.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "../sm/Elements/LatticeElements/latticestructuralelement.h"
#include "classfactory.h"
#include "../sm/Materials/structuralmaterial.h"
#include "contextioerr.h"
#include "datastream.h"
#include "sm/CrossSections/latticecrosssection.h"
#include "latticebeam3d.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "dof.h"

#endif

namespace oofem {
REGISTER_Element(LatticeBeam3dBoundary);

LatticeBeam3dBoundary :: LatticeBeam3dBoundary(int n, Domain *aDomain) : LatticeBeam3d(n, aDomain)
{
    numberOfDofMans = 3;
    geometryFlag = 0;
}

LatticeBeam3dBoundary :: ~LatticeBeam3dBoundary()
{}


void
LatticeBeam3dBoundary :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    //    double dV;
    FloatMatrix d, bi, bj, dbj, dij, bjt;
    FloatMatrix t(12, 18), tt;
    FloatMatrix answerTemp(12, 12), ttk(18, 12);
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answerTemp.zero();
    t.zero();

    //Call the parent and compute stiffness
    LatticeBeam3d :: computeStiffnessMatrix(answerTemp, rMode, tStep);

    answer.resize(computeNumberOfDofs(), computeNumberOfDofs() );
    answer.zero();

    for ( int m = 1; m <= 12; m++ ) {
        for ( int k = 1; k <= 12; k++ ) {
            answer.at(m, k) = answerTemp.at(m, k);
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }

    // Rotate to global system add the projections and rotate back to local system.

    FloatMatrix R;
    if ( this->giveRotationMatrix(R) ) {
        answer.rotatedWith(R);
    }

    for ( int m = 1; m <= 12; m++ ) {
        for ( int k = 1; k <= 12; k++ ) {
            answerTemp.at(m, k) = answer.at(m, k);
        }
    }

    IntArray projectionComponentNodeOne(3);
    projectionComponentNodeOne.zero();
    if ( location.at(1) != 0 ) {
        giveSwitches(projectionComponentNodeOne, location.at(1) );
    }

    IntArray projectionComponentNodeTwo(3);
    projectionComponentNodeTwo.zero();
    if ( location.at(2) != 0 ) {
        giveSwitches(projectionComponentNodeTwo, location.at(2) );
    }

    for ( int k = 1; k <= 12; k++ ) {
        t.at(k, k) = 1.;
    }
    t.at(1, 13) = projectionComponentNodeOne.at(1);
    t.at(2, 14) = projectionComponentNodeOne.at(2);
    t.at(3, 15) = projectionComponentNodeOne.at(3);


    t.at(2, 16) = projectionComponentNodeOne.at(3);
    t.at(1, 17) = projectionComponentNodeOne.at(3);
    t.at(1, 18) = projectionComponentNodeOne.at(2);


    t.at(7, 13) = projectionComponentNodeTwo.at(1);
    t.at(8, 14) = projectionComponentNodeTwo.at(2);
    t.at(9, 15) = projectionComponentNodeTwo.at(3);

    t.at(8, 16) = projectionComponentNodeTwo.at(3);
    t.at(7, 17) = projectionComponentNodeTwo.at(3);
    t.at(7, 18) = projectionComponentNodeTwo.at(2);


    tt.beTranspositionOf(t);

    ttk.beProductOf(tt, answerTemp);
    answer.beProductOf(ttk, t);

    //Rotate back to local system
    FloatMatrix Rtranspose;
    Rtranspose.beTranspositionOf(R);
    answer.rotatedWith(Rtranspose);

    return;
}



double
LatticeBeam3dBoundary :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->area * this->length;
}

void
LatticeBeam3dBoundary :: giveVTKCoordinates(int nodeNumber, FloatArray &coords) {
    coords.resize(3);
    coords.zero();
    Node *node;

    FloatArray specimenDimension(3);
    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);
    specimenDimension.at(3) =  this->giveNode(3)->giveCoordinate(3);

    IntArray projectionComponent(3);
    projectionComponent.zero();

    if ( nodeNumber == 1 ) {
        node  = this->giveNode(1);
        if ( location.at(1) != 0 ) {
            giveSwitches(projectionComponent, location.at(1) );
        }
    } else if ( nodeNumber == 2 ) {
        node  = this->giveNode(2);
        if ( location.at(2) != 0 ) {
            giveSwitches(projectionComponent, location.at(2) );
        }
    } else {
        OOFEM_ERROR("wrong element used in the vtk module");
    }

    for ( int i = 0; i < 3; i++ ) {
        coords.at(i + 1) =  node->giveCoordinate(i + 1) + projectionComponent.at(i + 1) * specimenDimension.at(i + 1);
    }
    return;
}

bool
LatticeBeam3dBoundary :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    int i, j;

    answer.resize(18, 18);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
            answer.at(i + 6, j + 6) = lcs.at(i, j);
            answer.at(i + 9, j + 9) = lcs.at(i, j);
        }
    }
    answer.at(13, 13) = 1.;
    answer.at(14, 14) = 1.;
    answer.at(15, 15) = 1.;
    answer.at(16, 16) = 1.;
    answer.at(17, 17) = 1.;
    answer.at(18, 18) = 1.;

    return 1;
}

int
LatticeBeam3dBoundary :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->localCoordinateSystem;

    return 1;
}



void
LatticeBeam3dBoundary ::   giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

void
LatticeBeam3dBoundary :: initializeFrom(InputRecord &ir)
{
    LatticeBeam3d :: initializeFrom(ir);

    location.resize(2);
    IR_GIVE_FIELD(ir, location, _IFT_LatticeBeam3dBoundary_location);
}



void
LatticeBeam3dBoundary :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix stiffness;
    FloatArray u;

    this->computeVectorOf(VM_Total, tStep, u);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // //Compute projection vector
    // IntArray projectionComponentNodeOne(3);
    // projectionComponentNodeOne.zero();
    // if ( location.at(1) != 0 ) {
    //     giveSwitches( projectionComponentNodeOne, location.at(1) );
    // }

    // IntArray projectionComponentNodeTwo(3);
    // projectionComponentNodeTwo.zero();

    // if ( location.at(2) != 0 ) {
    //     giveSwitches( projectionComponentNodeTwo, location.at(2) );
    // }

    // //Rotate displacements into global coordinate system
    // FloatMatrix rotationMatrix;
    // if ( this->computeGtoLRotationMatrix(rotationMatrix) ) {
    //     u.rotatedWith(rotationMatrix, 't');
    // }

    // //General expressions for the corrected displacements. Rotations at nodes are not influenced. Only translations
    // u.at(1) = u.at(1) + projectionComponentNodeOne.at(1) * u.at(13) + projectionComponentNodeOne.at(2) * u.at(18) + projectionComponentNodeOne.at(3) * u.at(17);
    // u.at(2) = u.at(2) + projectionComponentNodeOne.at(2) * u.at(14) + projectionComponentNodeOne.at(3) * u.at(16);
    // u.at(3) = u.at(3) + projectionComponentNodeOne.at(3) * u.at(15);

    // u.at(7) = u.at(7) + projectionComponentNodeTwo.at(1) * u.at(13) + projectionComponentNodeTwo.at(2) * u.at(18) + projectionComponentNodeTwo.at(3) * u.at(17);
    // u.at(8) = u.at(8) + projectionComponentNodeTwo.at(2) * u.at(14) + projectionComponentNodeTwo.at(3) * u.at(16);
    // u.at(9) = u.at(9) + projectionComponentNodeTwo.at(3) * u.at(15);


    // if ( this->computeGtoLRotationMatrix(rotationMatrix) ) {
    //     u.rotatedWith(rotationMatrix, 'n');
    // }

    this->computeStiffnessMatrix(stiffness, ElasticStiffness, tStep);

    // zero answer will resize accordingly when adding first contribution
    answer.clear();
    answer.beProductOf(stiffness, u);

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
    return;
}


void
LatticeBeam3dBoundary :: giveSwitches(IntArray &answer, int location) {
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
LatticeBeam3dBoundary :: computeGeometryProperties()
{
    Node *nodeA, *nodeB;

    FloatArray coordsA(3);
    FloatArray coordsB(3);

    nodeA  = this->giveNode(1);
    nodeB  = this->giveNode(2);

    FloatArray specimenDimension(3);
    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);
    specimenDimension.at(3) =  this->giveNode(3)->giveCoordinate(3);

    IntArray projectionComponentNodeOne(3);
    projectionComponentNodeOne.zero();
    if ( location.at(1) != 0 ) {
        giveSwitches(projectionComponentNodeOne, location.at(1) );
    }

    IntArray projectionComponentNodeTwo(3);
    projectionComponentNodeTwo.zero();
    if ( location.at(2) != 0 ) {
        giveSwitches(projectionComponentNodeTwo, location.at(2) );
    }

    for ( int i = 0; i < 3; i++ ) {
        coordsA.at(i + 1) =  nodeA->giveCoordinate(i + 1) + projectionComponentNodeOne.at(i + 1) * specimenDimension.at(i + 1);
        coordsB.at(i + 1) =  nodeB->giveCoordinate(i + 1) + projectionComponentNodeTwo.at(i + 1) * specimenDimension.at(i + 1);
    }

    //Calculate normal vector
    FloatArray s(3), t(3);
    this->midPoint.resize(3);

    //Calculate normal vector
    this->normal.resize(3);
    for ( int i = 0; i < 3; i++ ) {
        this->normal.at(i + 1) = coordsB.at(i + 1) - coordsA.at(i + 1);
    }

    this->length  = sqrt(pow(this->normal.at(1), 2.) + pow(this->normal.at(2), 2.) + pow(this->normal.at(3), 2.) );

    for ( int i = 0; i < 3; i++ ) {
        this->normal.at(i + 1) /= length;
    }

    // Compute midpoint
    this->midPoint.resize(3);
    for ( int i = 0; i < 3; i++ ) {
        this->midPoint.at(i + 1) = 0.5 * ( coordsB.at(i + 1) + coordsA.at(i + 1) );
    }

    this->globalCentroid = this->midPoint;

    computeCrossSectionProperties();

    //Set geometry flag to 1 so that this is done only once
    this->geometryFlag = 1;

    return;
}


void
LatticeBeam3dBoundary :: saveContext(DataStream &stream, ContextMode mode)
{
    LatticeBeam3d :: saveContext(stream, mode);

    contextIOResultType iores;

    if ( ( mode & CM_Definition ) ) {
        if ( ( iores = location.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }
}

void
LatticeBeam3dBoundary :: restoreContext(DataStream &stream, ContextMode mode)
{
    LatticeBeam3d :: restoreContext(stream, mode);

    contextIOResultType iores;

    if ( mode & CM_Definition ) {
        if ( ( iores = this->location.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }
}


int
LatticeBeam3dBoundary :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer.resize(3);
    answer = this->globalCentroid;

    return 1;
}



#ifdef __OOFEG

void
LatticeBeam3dBoundary :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc, tStep);
    } else if ( mode == OGC_deformedGeometry ) {
        this->drawDeformedGeometry(gc, tStep, DisplacementVector);
    } else if ( mode == OGC_elemSpecial ) {
        this->drawSpecial(gc, tStep);
    } else {
        OOFEM_ERROR("drawYourself : unsupported mode");
    }
}


void LatticeBeam3dBoundary :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    WCRec p [ 2 ]; /* points */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getActiveCrackColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);


    FloatArray specimenDimension(3);
    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);
    specimenDimension.at(3) =  this->giveNode(3)->giveCoordinate(3);

    //Compute projection vector
    IntArray projectionComponentNodeOne(3);
    projectionComponentNodeOne.zero();
    if ( location.at(1) != 0 ) {
        giveSwitches(projectionComponentNodeOne, location.at(1) );
    }

    IntArray projectionComponentNodeTwo(3);
    projectionComponentNodeTwo.zero();

    if ( location.at(2) != 0 ) {
        giveSwitches(projectionComponentNodeTwo, location.at(2) );
    }

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1) + projectionComponentNodeOne.at(1) * specimenDimension.at(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2) + projectionComponentNodeOne.at(2) * specimenDimension.at(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3) + projectionComponentNodeOne.at(3) * specimenDimension.at(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1) + projectionComponentNodeTwo.at(1) * specimenDimension.at(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2) + projectionComponentNodeTwo.at(2) * specimenDimension.at(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3) + projectionComponentNodeTwo.at(3) * specimenDimension.at(3);

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void LatticeBeam3dBoundary :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
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

    FloatArray specimenDimension(3);
    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);
    specimenDimension.at(3) =  this->giveNode(3)->giveCoordinate(3);

    FloatArray dispOne(6);
    dispOne.at(1) = this->giveNode(1)->giveDofWithID(D_u)->giveUnknown(VM_Total, tStep);
    dispOne.at(2) = this->giveNode(1)->giveDofWithID(D_v)->giveUnknown(VM_Total, tStep);
    dispOne.at(3) = this->giveNode(1)->giveDofWithID(D_w)->giveUnknown(VM_Total, tStep);
    dispOne.at(4) = this->giveNode(1)->giveDofWithID(R_u)->giveUnknown(VM_Total, tStep);
    dispOne.at(5) = this->giveNode(1)->giveDofWithID(R_v)->giveUnknown(VM_Total, tStep);
    dispOne.at(6) = this->giveNode(1)->giveDofWithID(R_w)->giveUnknown(VM_Total, tStep);

    FloatArray dispTwo(6);
    dispTwo.at(1) = this->giveNode(2)->giveDofWithID(D_u)->giveUnknown(VM_Total, tStep);
    dispTwo.at(2) = this->giveNode(2)->giveDofWithID(D_v)->giveUnknown(VM_Total, tStep);
    dispTwo.at(3) = this->giveNode(2)->giveDofWithID(D_w)->giveUnknown(VM_Total, tStep);
    dispTwo.at(4) = this->giveNode(2)->giveDofWithID(R_u)->giveUnknown(VM_Total, tStep);
    dispTwo.at(5) = this->giveNode(2)->giveDofWithID(R_v)->giveUnknown(VM_Total, tStep);
    dispTwo.at(6) = this->giveNode(2)->giveDofWithID(R_w)->giveUnknown(VM_Total, tStep);

    FloatArray dispThree(6);
    dispThree.at(1) = this->giveNode(3)->giveDofWithID(D_u)->giveUnknown(VM_Total, tStep);
    dispThree.at(2) = this->giveNode(3)->giveDofWithID(D_v)->giveUnknown(VM_Total, tStep);
    dispThree.at(3) = this->giveNode(3)->giveDofWithID(D_w)->giveUnknown(VM_Total, tStep);
    dispThree.at(4) = this->giveNode(3)->giveDofWithID(R_u)->giveUnknown(VM_Total, tStep);
    dispThree.at(5) = this->giveNode(3)->giveDofWithID(R_v)->giveUnknown(VM_Total, tStep);
    dispThree.at(6) = this->giveNode(3)->giveDofWithID(R_w)->giveUnknown(VM_Total, tStep);

    IntArray projectionComponentNodeOne(3);
    projectionComponentNodeOne.zero();
    if ( location.at(1) != 0 ) {
        giveSwitches(projectionComponentNodeOne, location.at(1) );
    }

    IntArray projectionComponentNodeTwo(3);
    projectionComponentNodeTwo.zero();
    if ( location.at(2) != 0 ) {
        giveSwitches(projectionComponentNodeTwo, location.at(2) );
    }


    //Modify dispOne and dispTwo
    dispOne.at(1) = dispOne.at(1) + projectionComponentNodeOne.at(1) * dispThree.at(1) + projectionComponentNodeOne.at(2) * dispThree.at(4) + projectionComponentNodeOne.at(3) * dispThree.at(5);
    dispOne.at(2) = dispOne.at(2) + projectionComponentNodeOne.at(2) * dispThree.at(2) + projectionComponentNodeOne.at(3) * dispThree.at(6);
    dispOne.at(3) = dispOne.at(3) + projectionComponentNodeOne.at(3) * dispThree.at(3);


    dispTwo.at(1) = dispTwo.at(1) + projectionComponentNodeTwo.at(1) * dispThree.at(1) + projectionComponentNodeTwo.at(2) * dispThree.at(4) + projectionComponentNodeTwo.at(3) * dispThree.at(5);
    dispTwo.at(2) = dispTwo.at(2) + projectionComponentNodeTwo.at(2) * dispThree.at(2) + projectionComponentNodeTwo.at(3) * dispThree.at(6);
    dispTwo.at(3) = dispTwo.at(3) + projectionComponentNodeTwo.at(3) * dispThree.at(3);

    double x1, y1, z1, x2, y2, z2;
    x1 = this->giveNode(1)->giveCoordinate(1) + projectionComponentNodeOne.at(1) * specimenDimension.at(1);
    y1 = this->giveNode(1)->giveCoordinate(2) + projectionComponentNodeOne.at(2) * specimenDimension.at(2);
    z1 = this->giveNode(1)->giveCoordinate(3) + projectionComponentNodeOne.at(3) * specimenDimension.at(3);

    x2 = this->giveNode(2)->giveCoordinate(1) + projectionComponentNodeTwo.at(1) * specimenDimension.at(1);
    y2 = this->giveNode(2)->giveCoordinate(2) + projectionComponentNodeTwo.at(2) * specimenDimension.at(2);
    z2 = this->giveNode(2)->giveCoordinate(3) + projectionComponentNodeTwo.at(3) * specimenDimension.at(3);

    p [ 0 ].x = ( FPNum ) x1 + defScale * dispOne.at(1);
    p [ 0 ].y = ( FPNum ) y1 + defScale * dispOne.at(2);
    p [ 0 ].z = ( FPNum ) z1 + defScale * dispOne.at(3);

    p [ 1 ].x = ( FPNum ) x2 + defScale * dispTwo.at(1);
    p [ 1 ].y = ( FPNum ) y2 + defScale * dispTwo.at(2);
    p [ 1 ].z = ( FPNum ) z2 + defScale * dispTwo.at(3);

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}

#endif
} // end namespace oofem
