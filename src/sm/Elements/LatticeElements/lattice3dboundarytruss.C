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
#include "lattice3dboundarytruss.h"
#include "../sm/Materials/LatticeMaterials/latticematstatus.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "latticestructuralelement.h"
#include "classfactory.h"
#include "../sm/Materials/structuralmaterial.h"
#include "contextioerr.h"
#include "datastream.h"
#include "crosssection.h"
#include "dof.h"
#include "parametermanager.h"
#include "paramkey.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Lattice3dBoundaryTruss);
ParamKey Lattice3dBoundaryTruss::IPK_Lattice3dBoundaryTruss_location("location");

Lattice3dBoundaryTruss :: Lattice3dBoundaryTruss(int n, Domain *aDomain) : Lattice3dBoundary(n, aDomain), location (2)
{
    numberOfDofMans = 3;
    geometryFlag = 0;
}

Lattice3dBoundaryTruss :: ~Lattice3dBoundaryTruss()
{}


void
Lattice3dBoundaryTruss :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    //Assemble Bmatrix (used to compute strains and rotations}
    answer.resize(6, 12);
    answer.zero();

    //Normal displacement jump in x-direction
    //First node
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 0.;
    answer.at(1, 3) = 0.;
    answer.at(1, 4) = 0.;
    answer.at(1, 5) = -this->eccT;
    answer.at(1, 6) = this->eccS;
    //Second node
    answer.at(1, 7) = 1.;
    answer.at(1, 8) = 0.;
    answer.at(1, 9) = 0.;
    answer.at(1, 10) = 0.;
    answer.at(1, 11) = this->eccT;
    answer.at(1, 12) = -this->eccS;

    //Shear displacement jump in y-plane
    //first node
    answer.at(2, 1) = 0.;
    answer.at(2, 2) = -1.;
    answer.at(2, 3) =  0.;
    answer.at(2, 4) = this->eccT;
    answer.at(2, 5) = 0;
    answer.at(2, 6) = -this->length / 2.;
    //Second node
    answer.at(2, 7) = 0.;
    answer.at(2, 8) = 1.;
    answer.at(2, 9) =  0.;
    answer.at(2, 10) = -this->eccT;
    answer.at(2, 11) = 0;
    answer.at(2, 12) = -this->length / 2.;

    //Shear displacement jump in z-plane
    //first node
    answer.at(3, 1) = 0.;
    answer.at(3, 2) = 0.;
    answer.at(3, 3) = -1.;
    answer.at(3, 4) = -this->eccS;
    answer.at(3, 5) = this->length / 2.;
    answer.at(3, 6) = 0.;
    //Second node
    answer.at(3, 7) = 0.;
    answer.at(3, 8) = 0.;
    answer.at(3, 9) =  1.;
    answer.at(3, 10) = this->eccS;
    answer.at(3, 11) = this->length / 2.;
    answer.at(3, 12) = 0.;

    //Rotation around x-axis
    //First node
    answer.at(4, 1) = 0.;
    answer.at(4, 2) = 0;
    answer.at(4, 3) = 0.;
    answer.at(4, 4) = -sqrt(Ip / this->area);
    answer.at(4, 5) = 0.;
    answer.at(4, 6) = 0.;
    //Second node
    answer.at(4, 7) = 0.;
    answer.at(4, 8) = 0.;
    answer.at(4, 9) = 0.;
    answer.at(4, 10) = sqrt(Ip / this->area);
    answer.at(4, 11) = 0.;
    answer.at(4, 12) = 0.;

    //Rotation around y-axis
    //First node
    answer.at(5, 1) = 0.;
    answer.at(5, 2) = 0.;
    answer.at(5, 3) = 0.;
    answer.at(5, 4) = 0.;
    answer.at(5, 5) = -sqrt(I1 / this->area);
    answer.at(5, 6) = 0.;
    //Second node
    answer.at(5, 7) = 0.;
    answer.at(5, 8) = 0.;
    answer.at(5, 9) =  0.;
    answer.at(5, 10) = 0.;
    answer.at(5, 11) = sqrt(I1 / this->area);
    answer.at(5, 12) = 0.;

    //Rotation around z-axis
    //First node
    answer.at(6, 1) = 0.;
    answer.at(6, 2) = 0.;
    answer.at(6, 3) = 0.;
    answer.at(6, 4) = 0.;
    answer.at(6, 5) = 0.;
    answer.at(6, 6) = -sqrt(I2 / this->area);
    //Second node
    answer.at(6, 7) = 0.;
    answer.at(6, 8) = 0.;
    answer.at(6, 9) =  0.;
    answer.at(6, 10) = 0.;
    answer.at(6, 11) = 0.;
    answer.at(6, 12) = sqrt(I2 / this->area);

    answer.times(1. / this->length);

    return;
}

void
Lattice3dBoundaryTruss :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                 TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    //    double dV;
    FloatMatrix d, bi, bj, dbj, dij, bjt;
    FloatMatrix t(12, 13), tt;
    FloatMatrix answerTemp(12, 12), answerHelp, ttk(13, 12);
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answerTemp.zero();
    answerHelp.zero();
    t.zero();

    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    double volume = this->computeVolumeAround(integrationRulesArray [ 0 ]->getIntegrationPoint(0) );

    this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bj);
    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    for ( int i = 1; i <= 6; i++ ) {
        d.at(i, i) *= volume;
    }

    dbj.beProductOf(d, bj);
    bjt.beTranspositionOf(bj);
    answerTemp.beProductOf(bjt, dbj);

    answer.resize(computeNumberOfDofs(), computeNumberOfDofs() );
    answer.zero();

    answerHelp.resize(computeNumberOfDofs(), computeNumberOfDofs() );
    answerHelp.zero();

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
    t.at(7, 13) = projectionComponentNodeTwo.at(1);

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
Lattice3dBoundaryTruss :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->area * this->length;
}

void
Lattice3dBoundaryTruss :: recalculateCoordinates(int nodeNumber, FloatArray &coords) {
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

void
Lattice3dBoundaryTruss :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u;

    //Compute strain vector
    //Get the 13 components of the displacement vector of this element
    this->computeVectorOf(VM_Total, stepN, u);
    this->computeBmatrixAt(gp, b);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

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

    FloatMatrix rotationMatrix;
    if ( this->computeGtoLRotationMatrix(rotationMatrix) ) {
        u.rotatedWith(rotationMatrix, 't');
    }

    //General expressions for the corrected displacements. Rotations at nodes are not influenced. Only translations
    u.at(1) = u.at(1) + projectionComponentNodeOne.at(1) * u.at(13);

    u.at(7) = u.at(7) + projectionComponentNodeTwo.at(1) * u.at(13);

    if ( this->computeGtoLRotationMatrix(rotationMatrix) ) {
        u.rotatedWith(rotationMatrix, 'n');
    }
    //define a temp displacement vector
    FloatArray uTemp(12);

    for ( int i = 1; i <= 12; i++ ) {
        uTemp.at(i) = u.at(i);
    }

    answer.beProductOf(b, uTemp);
}

bool
Lattice3dBoundaryTruss :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    int i, j;

    answer.resize(13, 13);
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

    return 1;
}

int
Lattice3dBoundaryTruss :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->localCoordinateSystem;

    return 1;
}



void
Lattice3dBoundaryTruss ::   giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode == 3 ) {
        answer = { E_xx };
    } else {
        answer = { D_u, D_v, D_w, R_u, R_v, R_w };
    }
}

void
Lattice3dBoundaryTruss :: initializeFrom(InputRecord &ir, int priority)
{
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    Lattice3d :: initializeFrom(ir, priority);
    PM_UPDATE_PARAMETER(location, ppm, ir, this->number, IPK_Lattice3dBoundaryTruss_location, priority) ;
}

void
Lattice3dBoundaryTruss :: postInitialize()
{
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    // first call parent
    Lattice3d :: postInitialize();
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_Lattice3dBoundaryTruss_location) ;
}


void
Lattice3dBoundaryTruss :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    Material *mat = this->giveMaterial();

    FloatMatrix b, bt, A, R, GNT;
    FloatArray bs, TotalStressVector, u, strain;
    double dV;

    this->computeVectorOf(VM_Total, tStep, u);

    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    answer.resize(13);
    answer.zero();

    this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), b);

    bt.beTranspositionOf(b);
    if ( useUpdatedGpRecord == 1 ) {
        TotalStressVector = ( ( LatticeMaterialStatus * ) mat->giveStatus(integrationRulesArray [ 0 ]->getIntegrationPoint(0) ) )
                            ->giveLatticeStress();
    } else
    if ( !this->isActivated(tStep) ) {
        strain.resize(StructuralMaterial :: giveSizeOfVoigtSymVector(integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialMode() ) );
        strain.zero();
    }
    this->computeStrainVector(strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    this->computeStressVector(TotalStressVector, strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    dV  = this->computeVolumeAround(integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
    bs.beProductOf(bt, TotalStressVector);
    bs.times(dV);

    for ( int m = 1; m <= 12; m++ ) {
        answer.at(m) = bs.at(m);
    }

    if ( this->giveRotationMatrix(R) ) {
        answer.rotatedWith(R, 't');
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

    //Normal stresses
    answer.at(13) = projectionComponentNodeOne.at(1) * answer.at(1) + projectionComponentNodeTwo.at(1) * answer.at(7);

    //Rotate to local system
    if ( this->giveRotationMatrix(R) ) {
        answer.rotatedWith(R, 'n');
    }

    return;
}


void
Lattice3dBoundaryTruss :: giveSwitches(IntArray &answer, int location) {
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
Lattice3dBoundaryTruss :: computeGeometryProperties()
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

    computeCrossSectionProperties();

    //Set geometry flag to 1 so that this is done only once
    this->geometryFlag = 1;

    return;
}


void
Lattice3dBoundaryTruss :: saveContext(DataStream &stream, ContextMode mode)
{
    Lattice3d :: saveContext(stream, mode);

    contextIOResultType iores;

    if ( ( mode & CM_Definition ) ) {
        if ( ( iores = location.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }
}



void
Lattice3dBoundaryTruss :: restoreContext(DataStream &stream, ContextMode mode)
{
    Lattice3d :: restoreContext(stream, mode);

    contextIOResultType iores;

    if ( mode & CM_Definition ) {
        if ( ( iores = this->location.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }
}



#ifdef __OOFEG

void
Lattice3dBoundaryTruss :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc, tStep);
        this->drawRawCrossSections(gc, tStep);
    } else if ( mode == OGC_deformedGeometry ) {
        this->drawDeformedGeometry(gc, tStep, DisplacementVector);
    } else if ( mode == OGC_elemSpecial ) {
        this->drawSpecial(gc, tStep);
    } else {
        OOFEM_ERROR("drawYourself : unsupported mode");
    }
}




void Lattice3dBoundaryTruss :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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


void
Lattice3dBoundaryTruss :: drawRawCrossSections(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    //  if (!go) { // create new one
    //Create as many points as we have polygon vertices
    numberOfPolygonVertices = this->polygonCoords.giveSize() / 3.;
    WCRec p [ numberOfPolygonVertices ]; /* poin */

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    //  EASValsSetColor(gc.getNodeColor());
    EASValsSetLayer(OOFEG_RAW_CROSSSECTION_LAYER);

    for ( int i = 0; i < numberOfPolygonVertices; i++ ) {
        p [ i ].x = ( FPNum ) polygonCoords(3 * i);
        p [ i ].y = ( FPNum ) polygonCoords(3 * i + 1);
        p [ i ].z = ( FPNum ) polygonCoords(3 * i + 2);
    }


    WCRec pTemp [ 2 ]; /* points */
    for ( int i = 0; i < numberOfPolygonVertices; i++ ) {
        if ( i < numberOfPolygonVertices - 1 ) {
            pTemp [ 0 ] = p [ i ];
            pTemp [ 1 ] = p [ i + 1 ];
        } else {
            pTemp [ 0 ] = p [ i ];
            pTemp [ 1 ] = p [ 0 ];
        }

        go = CreateLine3D(pTemp);
        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
        EGAttachObject(go, ( EObjectP ) this);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}


void Lattice3dBoundaryTruss :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    //That seems to be wrong. The strain field should be ordered exx, eyy, ezz, gyz, gzx, gyx
    //Therefore, the x displacement should include 5th and 6th strain components.
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

    FloatArray dispThree(1);
    dispThree.at(1) = this->giveNode(3)->giveDofWithID(E_xx)->giveUnknown(VM_Total, tStep);

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
    dispOne.at(1) = dispOne.at(1) + projectionComponentNodeOne.at(1) * dispThree.at(1);
    dispOne.at(2) = dispOne.at(2);
    dispOne.at(3) = dispOne.at(3);


    dispTwo.at(1) = dispTwo.at(1) + projectionComponentNodeTwo.at(1) * dispThree.at(1);
    dispTwo.at(2) = dispTwo.at(2);
    dispTwo.at(3) = dispTwo.at(3);

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
