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

#include "../sm/Elements/LatticeElements/lattice2dboundary.h"
#include "../sm/Elements/LatticeElements/lattice2d.h"
#include "../sm/Materials/LatticeMaterials/latticematstatus.h"
#include "domain.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dof.h"
#include "../sm/Materials/structuralmaterial.h"
#include "parametermanager.h"
#include "paramkey.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Lattice2dBoundary);
ParamKey Lattice2dBoundary::IPK_Lattice2dBoundary_location("location");

Lattice2dBoundary :: Lattice2dBoundary(int n, Domain *aDomain) : Lattice2d(n, aDomain)
    // Constructor.
{
    numberOfDofMans     = 3;

    kappa = -1; // set kappa to undef value (should be always > 0.)

    pitch               = 10.;  // a dummy value
}

Lattice2dBoundary :: ~Lattice2dBoundary()
{
    // destructor
}

void
Lattice2dBoundary :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    double length = this->giveLength();
    double x1, x2, y1, y2, xp, yp;

    FloatArray projectionComponent(2);
    giveSwitches(projectionComponent);


    xp = this->gpCoords.at(1);
    yp = this->gpCoords.at(2);


    //specimen dimension
    FloatArray specimenDimension(2);
    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);


    //Coordinates of the nodes
    x1 = this->giveNode(1)->giveCoordinate(1);
    y1 = this->giveNode(1)->giveCoordinate(2);
    x2 = this->giveNode(2)->giveCoordinate(1) + projectionComponent.at(1) * specimenDimension.at(1);
    y2 = this->giveNode(2)->giveCoordinate(2) + projectionComponent.at(2) * specimenDimension.at(2);


    double ecc;

    double areaHelp = 0.5 * ( x1 * y2 + x2 * yp + xp * y1 - ( xp * y2 + yp * x1 + x2 * y1 ) );

    ecc = 2 * areaHelp / length;

    //    Assemble Bmatrix (used to compute strains and rotations
    answer.resize(3, 6);
    answer.zero();
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 0.;
    answer.at(1, 3) = ecc;
    answer.at(1, 4) = 1.;
    answer.at(1, 5) = 0.;
    answer.at(1, 6) = -ecc;

    answer.at(2, 1) = 0.;
    answer.at(2, 2) = -1.;
    answer.at(2, 3) = -length / 2.;
    answer.at(2, 4) = 0.;
    answer.at(2, 5) = 1.;
    answer.at(2, 6) = -length / 2.;

    answer.at(3, 1) = 0.;
    answer.at(3, 2) = 0.;
    answer.at(3, 3) = -this->width / sqrt(12.);

    answer.at(3, 4) = 0.;
    answer.at(3, 5) = 0.;
    answer.at(3, 6) = this->width / sqrt(12.);


    answer.times(1. / length);

    return;
}

const IntArray
Lattice2dBoundary :: giveLocation()
{
    //Map 2D locations to 3D system
    int newLocation = 0;
    if ( this->location == 1 ) {
        newLocation = 22;
    } else if ( this->location == 2 ) {
        newLocation = 25;
    } else if ( this->location == 3 ) {
        newLocation = 16;
    } else if ( this->location == 4 ) {
        newLocation = 8;
    } else if ( this->location == 5 ) {
        newLocation = 5;
    } else if ( this->location == 6 ) {
        newLocation = 2;
    } else if ( this->location == 7 ) {
        newLocation = 11;
    } else if ( this->location == 8 ) {
        newLocation = 19;
    }

    IntArray tempLocation(2);
    tempLocation.at(1) = 0;
    tempLocation.at(2) = newLocation;
    return tempLocation;
}

void
Lattice2dBoundary :: recalculateCoordinates(int nodeNumber, FloatArray &coords) {
    coords.resize(3);
    coords.zero();
    Node *node;

    FloatArray specimenDimension(3);
    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);
    specimenDimension.at(3) =  this->giveNode(3)->giveCoordinate(3);

    FloatArray projectionComponent(3);
    projectionComponent.zero();

    node  = this->giveNode(2);
    if ( location != 0 ) {
        giveSwitches(projectionComponent);
    }

    for ( int i = 0; i < 3; i++ ) {
        coords.at(i + 1) =  node->giveCoordinate(i + 1) + projectionComponent.at(i + 1) * specimenDimension.at(i + 1);
    }
    return;
}


void
Lattice2dBoundary :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,                                      TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    double dV;
    FloatMatrix d, bi, bj, dbj, dij;

    FloatMatrix answerTemp(6, 6);
    answerTemp.zero();


    this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bj);
    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    dV = this->computeVolumeAround(integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
    dbj.beProductOf(d, bj);
    answerTemp.plusProductUnsym(bj, dbj, dV);

    answer.resize(computeNumberOfDofs(), computeNumberOfDofs() );
    answer.zero();

    for ( int m = 1; m <= 6; m++ ) {
        for ( int k = 1; k <= 6; k++ ) {
            answer.at(m, k) = answerTemp.at(m, k);
        }
    }


    //Rotate to global system
    FloatMatrix R;
    if ( this->giveRotationMatrix(R) ) {
        answer.rotatedWith(R);
    }

    //Convert normal stiffness matrix to the modified one

    FloatArray projectionComponent(2);
    giveSwitches(projectionComponent);

    answer.at(1, 7) = projectionComponent.at(1) * answer.at(1, 4);
    answer.at(1, 8) = projectionComponent.at(2) * answer.at(1, 5);
    answer.at(1, 9) = projectionComponent.at(1) * answer.at(1, 5);

    answer.at(2, 7) = projectionComponent.at(1) * answer.at(2, 4);
    answer.at(2, 8) = projectionComponent.at(2) * answer.at(2, 5);
    answer.at(2, 9) = projectionComponent.at(1) * answer.at(2, 5);

    answer.at(3, 7) = projectionComponent.at(1) * answer.at(3, 4);
    answer.at(3, 8) = projectionComponent.at(2) * answer.at(3, 5);
    answer.at(3, 9) = projectionComponent.at(1) * answer.at(3, 5);

    answer.at(4, 7) = projectionComponent.at(1) * answer.at(4, 4);
    answer.at(4, 8) = projectionComponent.at(2) * answer.at(4, 5);
    answer.at(4, 9) = projectionComponent.at(1) * answer.at(4, 5);

    answer.at(5, 7) = projectionComponent.at(1) * answer.at(5, 4);
    answer.at(5, 8) = projectionComponent.at(2) * answer.at(5, 5);
    answer.at(5, 9) = projectionComponent.at(1) * answer.at(5, 5);

    answer.at(6, 7) = projectionComponent.at(1) * answer.at(6, 4);
    answer.at(6, 8) = projectionComponent.at(2) * answer.at(6, 5);
    answer.at(6, 9) = projectionComponent.at(1) * answer.at(6, 5);

    answer.at(7, 1) = projectionComponent.at(1) * answer.at(4, 1);
    answer.at(7, 2) = projectionComponent.at(1) * answer.at(4, 2);
    answer.at(7, 3) = projectionComponent.at(1) * answer.at(4, 3);
    answer.at(7, 4) = projectionComponent.at(1) * answer.at(4, 4);
    answer.at(7, 5) = projectionComponent.at(1) * answer.at(4, 5);
    answer.at(7, 6) = projectionComponent.at(1) * answer.at(4, 6);
    answer.at(7, 7) = projectionComponent.at(1) * projectionComponent.at(1) * answer.at(4, 4);
    answer.at(7, 8) = projectionComponent.at(1) * projectionComponent.at(2) * answer.at(4, 5);
    answer.at(7, 9) = projectionComponent.at(1) * projectionComponent.at(1) * answer.at(4, 5);

    answer.at(8, 1) = projectionComponent.at(2) * answer.at(5, 1);
    answer.at(8, 2) = projectionComponent.at(2) * answer.at(5, 2);
    answer.at(8, 3) = projectionComponent.at(2) * answer.at(5, 3);
    answer.at(8, 4) = projectionComponent.at(2) * answer.at(5, 4);
    answer.at(8, 5) = projectionComponent.at(2) * answer.at(5, 5);
    answer.at(8, 6) = projectionComponent.at(2) * answer.at(5, 6);
    answer.at(8, 7) = projectionComponent.at(1) * projectionComponent.at(2) * answer.at(5, 4);
    answer.at(8, 8) = projectionComponent.at(2) * projectionComponent.at(2) * answer.at(5, 5);
    answer.at(8, 9) = projectionComponent.at(1) * projectionComponent.at(2) * answer.at(5, 5);

    answer.at(9, 1) = projectionComponent.at(1) * answer.at(5, 1);
    answer.at(9, 2) = projectionComponent.at(1) * answer.at(5, 2);
    answer.at(9, 3) = projectionComponent.at(1) * answer.at(5, 3);
    answer.at(9, 4) = projectionComponent.at(1) * answer.at(5, 4);
    answer.at(9, 5) = projectionComponent.at(1) * answer.at(5, 5);
    answer.at(9, 6) = projectionComponent.at(1) * answer.at(5, 6);
    answer.at(9, 7) = projectionComponent.at(1) * projectionComponent.at(1) * answer.at(5, 4);
    answer.at(9, 8) = projectionComponent.at(1) * projectionComponent.at(2) * answer.at(5, 5);
    answer.at(9, 9) = projectionComponent.at(1) * projectionComponent.at(1) * answer.at(5, 5);
    //Rotate back to local system
    FloatMatrix Rtranspose;
    Rtranspose.beTranspositionOf(R);
    answer.rotatedWith(Rtranspose);
}

void
Lattice2dBoundary :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u;

    //Compute strain vector
    //Get the 9 components of the displacement vector of this element
    this->computeVectorOf(VM_Total, stepN, u);
    this->computeBmatrixAt(gp, b);


    //Compute projection vector
    FloatArray projectionComponent(2);
    giveSwitches(projectionComponent);

    FloatMatrix rotationMatrix;
    if ( this->computeGtoLRotationMatrix(rotationMatrix) ) {
        u.rotatedWith(rotationMatrix, 't');
    }

    //General expressions for the corrected displacements
    u.at(4) = u.at(4) + projectionComponent.at(1) * u.at(7);
    u.at(5) = u.at(5) + projectionComponent.at(2) * u.at(8) + projectionComponent.at(1) * u.at(9);


    if ( this->computeGtoLRotationMatrix(rotationMatrix) ) {
        u.rotatedWith(rotationMatrix, 'n');
    }


    //define a temp displacement vector
    FloatArray uTemp(6);

    for ( int i = 1; i <= 6; i++ ) {
        uTemp.at(i) = u.at(i);
    }

    answer.beProductOf(b, uTemp);
}


bool
Lattice2dBoundary :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    double sine, cosine;
    answer.resize(9, 9);
    answer.zero();

    sine           = sin(this->givePitch() );
    cosine         = cos(pitch);
    answer.at(1, 1) =  cosine;
    answer.at(1, 2) =  sine;
    answer.at(2, 1) = -sine;
    answer.at(2, 2) =  cosine;
    answer.at(3, 3) =  1.;
    answer.at(4, 4) =  cosine;
    answer.at(4, 5) =  sine;
    answer.at(5, 4) = -sine;
    answer.at(5, 5) =  cosine;
    answer.at(6, 6) =  1.;
    answer.at(7, 7) = 1.;
    answer.at(8, 8) = 1.;
    answer.at(9, 9) = 1.;

    return 1;
}


double
Lattice2dBoundary :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    double weight  = aGaussPoint->giveWeight();
    return weight * 0.5 * this->giveLength() * this->width * this->thickness;
}


void
Lattice2dBoundary ::   giveDofManDofIDMask(int inode, IntArray &answer) const {
    if ( inode == 3 ) {
        answer = { E_xx, E_yy, G_xy };
    } else {
        answer = { D_u, D_v, R_w };
    }
}



double Lattice2dBoundary :: giveLength()
// Returns the length of the receiver.
{
    //Compute strain vector
    FloatArray projectionComponent(2);
    giveSwitches(projectionComponent);

    //specimen dimension
    FloatArray specimenDimension(2);
    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);

    double dx, dy;
    Node *nodeA, *nodeB;

    nodeA   = this->giveNode(1);
    nodeB   = this->giveNode(2);

    dx      = nodeB->giveCoordinate(1) + projectionComponent.at(1) * specimenDimension.at(1) - nodeA->giveCoordinate(1);
    dy      = nodeB->giveCoordinate(2) + projectionComponent.at(2) * specimenDimension.at(2) - nodeA->giveCoordinate(2);

    return sqrt(dx * dx + dy * dy);
}


double Lattice2dBoundary :: givePitch()
// Returns the pitch of the receiver.
{
    double xA, xB, yA, yB;
    Node *nodeA, *nodeB;

    //specimen dimension
    FloatArray specimenDimension(2);
    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);

    FloatArray projectionComponent(2);
    giveSwitches(projectionComponent);

    if ( pitch == 10. ) {            // 10. : dummy initialization value
        nodeA  = this->giveNode(1);
        nodeB  = this->giveNode(2);
        xA     = nodeA->giveCoordinate(1);
        xB     = nodeB->giveCoordinate(1) + projectionComponent.at(1) * specimenDimension.at(1);
        yA     = nodeA->giveCoordinate(2);
        yB     = nodeB->giveCoordinate(2) + projectionComponent.at(2) * specimenDimension.at(2);
        pitch  = atan2(yB - yA, xB - xA);
    }
    return pitch;
}



void
Lattice2dBoundary :: initializeFrom(InputRecord &ir, int priority)
{
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    // first call parent
    Lattice2d :: initializeFrom(ir, priority);
    PM_UPDATE_PARAMETER(location, ppm, ir, this->number, IPK_Lattice2dBoundary_location, priority);
}

void
Lattice2dBoundary :: postInitialize()
{
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    // first call parent
    Lattice2d :: postInitialize();
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_Lattice2dBoundary_location);
    this->computeGaussPoints();
}



void
Lattice2dBoundary :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    Material *mat = this->giveMaterial();

    FloatMatrix b, bt, A, R, GNT;
    FloatArray bs, TotalStressVector, u, strain;
    double dV;


    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    answer.resize(9);
    answer.zero();

    this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), b);

    bt.beTranspositionOf(b);
    // TotalStressVector = gp->giveStressVector() ;
    if ( useUpdatedGpRecord == 1 ) {
        TotalStressVector = ( ( LatticeMaterialStatus * ) mat->giveStatus(integrationRulesArray [ 0 ]->getIntegrationPoint(0) ) )
                            ->giveLatticeStress();
    } else
    if ( !this->isActivated(tStep) ) {
        strain.resize(StructuralMaterial :: giveSizeOfVoigtSymVector(integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialMode() ) );
        strain.zero();
    }
    this->computeStrainVector(strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    //    strain.beProductOf(b, u);

    this->computeStressVector(TotalStressVector, strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    //
    // compute nodal representation of internal forces using f = B^T*Sigma dV
    //
    dV  = this->computeVolumeAround(integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
    bs.beProductOf(bt, TotalStressVector);
    bs.times(dV);

    for ( int m = 1; m <= 6; m++ ) {
        answer.at(m) = bs.at(m);
    }

    //Rotate to global system

    if ( this->giveRotationMatrix(R) ) {
        answer.rotatedWith(R, 't');
    }
    //Compute strain vector

    FloatArray projectionComponent(2);
    giveSwitches(projectionComponent);

    answer.at(7) = projectionComponent.at(1) * answer.at(4);
    answer.at(8) = projectionComponent.at(2) * answer.at(5);
    answer.at(9) = projectionComponent.at(1) * answer.at(5);

    //Rotate to local system
    if ( this->giveRotationMatrix(R) ) {
        answer.rotatedWith(R, 'n');
    }

    return;
}


void
Lattice2dBoundary :: giveSwitches(FloatArray &answer) {
    if ( this->location == 1 ) {
        answer(0) = 1;
        answer(1) = 0;
    } else if ( this->location == 2 ) {
        answer(0) = 1;
        answer(1) = 1;
    } else if ( this->location == 3 ) {
        answer(0) = 0;
        answer(1) = 1;
    } else if ( this->location == 4 ) {
        answer(0) = -1;
        answer(1) = 1;
    } else if ( this->location == 5 ) {
        answer(0) = -1;
        answer(1) = 0;
    } else if ( this->location == 6 ) {
        answer(0) = -1;
        answer(1) = -1;
    } else if ( this->location == 7 ) {
        answer(0) = 0;
        answer(1) = -1;
    } else if ( this->location == 8 ) {
        answer(0) = 1;
        answer(1) = -1;
    }

    return;
}


void Lattice2dBoundary :: saveContext(DataStream &stream, ContextMode mode)
{
    Lattice2d :: saveContext(stream, mode);


    if ( ( mode & CM_Definition ) ) {
        if ( !stream.write(location) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}


void
Lattice2dBoundary :: restoreContext(DataStream &stream, ContextMode mode)
{
    Lattice2d :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.read(location) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}


#ifdef __OOFEG

void
Lattice2dBoundary :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
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
        OOFEM_ERROR("unsupported mode");
    }
}


void Lattice2dBoundary :: drawRawCrossSections(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getCrossSectionColor() );
    EASValsSetLayer(OOFEG_RAW_CROSSSECTION_LAYER);

    //The cs definition is not needed anymore

    FloatArray gpCoords;
    giveGpCoordinates(gpCoords);

    FloatArray projectionComponent(2);
    giveSwitches(projectionComponent);

    FloatArray specimenDimension(2);

    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);

    // double width = cs->give(WIDTH);

    double x1, y1, x2, y2;
    x1 = this->giveNode(1)->giveCoordinate(1);
    y1 = this->giveNode(1)->giveCoordinate(2);

    x2 = this->giveNode(2)->giveCoordinate(1) + projectionComponent.at(1) * specimenDimension.at(1);
    y2 = this->giveNode(2)->giveCoordinate(2) + projectionComponent.at(2) * specimenDimension.at(2);

    //Compute normal and shear direction
    FloatArray normalDirection;
    FloatArray shearDirection;
    normalDirection.resize(2);
    normalDirection.zero();
    shearDirection.resize(2);
    shearDirection.zero();
    normalDirection.at(1) = x2 - x1;
    normalDirection.at(2) = y2 - y1;
    normalDirection.normalize();
    if ( normalDirection.at(2) == 0. ) {
        shearDirection.at(1) = 0.;
        shearDirection.at(2) = 1.;
    } else {
        shearDirection.at(1) = 1.0;
        shearDirection.at(2) =
            -normalDirection.at(1) / normalDirection.at(2);
    }
    shearDirection.normalize();

    p [ 0 ].x = ( FPNum ) this->gpCoords.at(1) - shearDirection.at(1) * this->width / 2.;
    p [ 0 ].y = ( FPNum ) this->gpCoords.at(2) - shearDirection.at(2) * this->width / 2.;
    p [ 0 ].z = 0.;

    p [ 1 ].x = ( FPNum ) this->gpCoords.at(1) + shearDirection.at(1) * this->width / 2.;
    p [ 1 ].y = ( FPNum ) this->gpCoords.at(2) + shearDirection.at(2) * this->width / 2.;
    p [ 1 ].z = 0.;

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Lattice2dBoundary :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);

    FloatArray projectionComponent(2);
    giveSwitches(projectionComponent);


    FloatArray specimenDimension(2);

    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1) + projectionComponent.at(1) * specimenDimension.at(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2) + projectionComponent.at(2) * specimenDimension.at(2);
    p [ 1 ].z = 0.;

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}



void Lattice2dBoundary :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);


    FloatArray projectionComponent(2);
    giveSwitches(projectionComponent);

    //specimen dimension
    FloatArray specimenDimension(2);
    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);

    double x1, y1, x2, y2;
    x1 = this->giveNode(1)->giveCoordinate(1);
    y1 = this->giveNode(1)->giveCoordinate(2);
    x2 = this->giveNode(2)->giveCoordinate(1) + projectionComponent.at(1) * specimenDimension.at(1);
    y2 = this->giveNode(2)->giveCoordinate(2) + projectionComponent.at(2) * specimenDimension.at(2);

    //Node One
    FloatArray dispOne(3);
    dispOne.at(1) = this->giveNode(1)->giveDofWithID(D_u)->giveUnknown(VM_Total, tStep);
    dispOne.at(2) = this->giveNode(1)->giveDofWithID(D_v)->giveUnknown(VM_Total, tStep);
    dispOne.at(3) = this->giveNode(1)->giveDofWithID(R_w)->giveUnknown(VM_Total, tStep);

    //Node Two
    FloatArray dispTwo(3);
    dispTwo.at(1) = this->giveNode(2)->giveDofWithID(D_u)->giveUnknown(VM_Total, tStep);
    dispTwo.at(2) = this->giveNode(2)->giveDofWithID(D_v)->giveUnknown(VM_Total, tStep);
    dispTwo.at(3) = this->giveNode(2)->giveDofWithID(R_w)->giveUnknown(VM_Total, tStep);

    //Node Three
    FloatArray dispThree(3);
    dispThree.at(1) = this->giveNode(3)->giveDofWithID(E_xx)->giveUnknown(VM_Total, tStep);
    dispThree.at(2) = this->giveNode(3)->giveDofWithID(E_yy)->giveUnknown(VM_Total, tStep);
    dispThree.at(3) = this->giveNode(3)->giveDofWithID(G_xy)->giveUnknown(VM_Total, tStep);

    //Modify dispOne and dispTwo
    dispTwo.at(1) = dispTwo.at(1) + projectionComponent.at(1) * dispThree.at(1);
    dispTwo.at(2) = dispTwo.at(2) + projectionComponent.at(2) * dispThree.at(2) + projectionComponent.at(1) * dispThree.at(3);

    p [ 0 ].x = ( FPNum ) x1 + defScale * dispOne.at(1);
    p [ 0 ].y = ( FPNum ) y1 + defScale * dispOne.at(2);
    p [ 0 ].z = 0.;

    p [ 1 ].x = ( FPNum ) x2 + defScale * dispTwo.at(1);
    p [ 1 ].y = ( FPNum ) y2 + defScale * dispTwo.at(2);
    p [ 1 ].z = 0.;

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}






void
Lattice2dBoundary :: drawSpecial(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec l [ 2 ];
    GraphicObj *tr;
    StructuralMaterial *mat = ( StructuralMaterial * ) this->giveMaterial();
    GaussPoint *gp;
    FloatArray crackStatuses, cf;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }
    if ( gc.giveIntVarType() == IST_CrackState ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        mat->giveIPValue(crackStatuses, gp, IST_CrackStatuses, tStep); //    printf("crackFlag= %e\n",crackStatuses(0));
        if ( crackStatuses(0) == 1. || crackStatuses(0) == 2. ) {
            FloatArray gpCoords;
            this->giveGpCoordinates(gpCoords);

            //specimen dimension
            FloatArray specimenDimension(2);
            specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
            specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);


            FloatArray projectionComponent(2);
            giveSwitches(projectionComponent);


            double x1, y1, x2, y2;
            x1 = this->giveNode(1)->giveCoordinate(1);
            y1 = this->giveNode(1)->giveCoordinate(2);
            x2 = this->giveNode(2)->giveCoordinate(1) + projectionComponent.at(1) * specimenDimension.at(1);
            y2 = this->giveNode(2)->giveCoordinate(2) + projectionComponent.at(2) * specimenDimension.at(2);

            //Compute normal and shear direction
            FloatArray normalDirection;
            FloatArray shearDirection;
            normalDirection.resize(2);
            normalDirection.zero();
            shearDirection.resize(2);
            shearDirection.zero();
            normalDirection.at(1) = x2 - x1;
            normalDirection.at(2) = y2 - y1;
            normalDirection.normalize();
            if ( normalDirection.at(2) == 0. ) {
                shearDirection.at(1) = 0.;
                shearDirection.at(2) = 1.;
            } else {
                shearDirection.at(1) = 1.0;
                shearDirection.at(2) =
                    -normalDirection.at(1) / normalDirection.at(2);
            }
            shearDirection.normalize();

            l [ 0 ].x = ( FPNum ) this->gpCoords.at(1) - shearDirection.at(1) * this->width / 2.;
            l [ 0 ].y = ( FPNum ) this->gpCoords.at(2) - shearDirection.at(2) * this->width / 2.;
            l [ 0 ].z = 0.;
            l [ 1 ].x = ( FPNum ) this->gpCoords.at(1) + shearDirection.at(1) * this->width / 2.;
            ;
            l [ 1 ].y = ( FPNum ) this->gpCoords.at(2) + shearDirection.at(2) * this->width / 2.;
            l [ 1 ].z = 0.;

            EASValsSetLayer(OOFEG_CRACK_PATTERN_LAYER);
            EASValsSetLineWidth(OOFEG_CRACK_PATTERN_WIDTH);
            if ( ( crackStatuses(0) == 1. ) ) {
                EASValsSetColor(gc.getActiveCrackColor() );
            } else if ( crackStatuses(0) == 2. ) {
                EASValsSetColor(gc.getCrackPatternColor() );
            }
            tr = CreateLine3D(l);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}

void
Lattice2dBoundary :: giveCrossSectionCoordinates(FloatArray &coords)
{
    FloatArray gpCoords;
    this->giveGpCoordinates(gpCoords);
    FloatArray projectionComponent(2);
    giveSwitches(projectionComponent);

    FloatArray specimenDimension(2);

    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);

    double x1, y1, x2, y2;
    x1 = this->giveNode(1)->giveCoordinate(1);
    y1 = this->giveNode(1)->giveCoordinate(2);

    x2 = this->giveNode(2)->giveCoordinate(1) + projectionComponent.at(1) * specimenDimension.at(1);
    y2 = this->giveNode(2)->giveCoordinate(2) + projectionComponent.at(2) * specimenDimension.at(2);

    //Compute normal and shear direction
    FloatArray normalDirection;
    FloatArray shearDirection;
    normalDirection.resize(2);
    normalDirection.zero();
    shearDirection.resize(2);
    shearDirection.zero();
    normalDirection.at(1) = x2 - x1;
    normalDirection.at(2) = y2 - y1;
    normalDirection.normalize();
    if ( normalDirection.at(2) == 0. ) {
        shearDirection.at(1) = 0.;
        shearDirection.at(2) = 1.;
    } else {
        shearDirection.at(1) = 1.0;
        shearDirection.at(2) =
            -normalDirection.at(1) / normalDirection.at(2);
    }
    shearDirection.normalize();

    coords.resize(6);
    coords.at(1) = this->gpCoords.at(1) - shearDirection.at(1) * this->width / 2.;
    coords.at(2) = this->gpCoords.at(2) - shearDirection.at(2) * this->width / 2.;
    coords.at(3) = 0.;
    coords.at(4) = this->gpCoords.at(1) + shearDirection.at(1) * this->width / 2.;
    coords.at(5) = this->gpCoords.at(2) + shearDirection.at(2) * this->width / 2.;
    coords.at(6) = 0.;

    return;
}

#endif
} // end namespace oofem
