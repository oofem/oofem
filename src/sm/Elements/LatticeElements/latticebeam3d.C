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
#include "../sm/Elements/LatticeElements/latticebeam3d.h"
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
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "sm/CrossSections/latticecrosssection.h"
#include "parametermanager.h"
#include "paramkey.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "../sm/Materials/structuralmaterial.h"
#endif

namespace oofem {
REGISTER_Element(LatticeBeam3d);
ParamKey LatticeBeam3d::IPK_LatticeBeam3d_diameter("diameter");

LatticeBeam3d :: LatticeBeam3d(int n, Domain *aDomain) : LatticeStructuralElement(n, aDomain)
{
    geometryFlag = 0;
    numberOfDofMans     = 2;
    myPi = 3.14159;
}

LatticeBeam3d :: ~LatticeBeam3d()
{}


void
LatticeBeam3d :: giveGPCoordinates(FloatArray &coords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    coords.resize(3);
    coords = this->globalCentroid;

    return;
}

double
LatticeBeam3d :: giveLength()
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->length;
}

void
LatticeBeam3d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give3dStiffnessMatrix(rMode, gp, tStep);
}


void
LatticeBeam3d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                        TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    FloatMatrix d, bi, bj, bjt, dbj, dij;

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    answer.resize(12, 12);
    answer.zero();

    //Stiffness matrix of Bernoulli beam from McGuire
    double a = pow(this->diameter / 2., 2.) * myPi;
    double iz = pow(this->diameter / 2., 4.) * myPi / 4.;
    double iy = iz;
    double j = iz + iy;
    double l = this->giveLength();

    double nu = ( static_cast< LatticeLinearElastic * >( this->giveMaterial() ) )->give('n', gp);
    double e = ( static_cast< LatticeLinearElastic * >( this->giveMaterial() ) )->give('E', gp);

    answer.at(1, 1) = a / l;
    answer.at(1, 7) = -a / l;

    answer.at(2, 2) = 12. * iz / pow(l, 3.);
    answer.at(2, 6) = 6. * iz / pow(l, 2.);
    answer.at(2, 8) = -12. * iz / pow(l, 3.);
    answer.at(2, 12) = 6. * iz / pow(l, 2.);

    answer.at(3, 3) = 12 * iy / pow(l, 3.);
    answer.at(3, 5) = -6. * iy / pow(l, 2.);
    answer.at(3, 9) = -12. * iy / pow(l, 3.);
    answer.at(3, 11) = -6. * iy / pow(l, 2.);

    answer.at(4, 4) = j / ( 2. * ( 1. + nu ) * l );
    answer.at(4, 10) = -j / ( 2. * ( 1. + nu ) * l );

    answer.at(5, 3) = answer.at(3, 5);
    answer.at(5, 5) = 4 * iy / l;
    answer.at(5, 9) = 6. * iy / pow(l, 2.);
    answer.at(5, 11) = 2. * iy / l;

    answer.at(6, 2) = answer.at(2, 6);
    answer.at(6, 6) = 4. * iz / l;
    answer.at(6, 8) = -6. * iz / pow(l, 2.);
    answer.at(6, 12) = 2. * iz / l;

    answer.at(7, 1) = answer.at(1, 7);
    answer.at(7, 7) = a / l;

    answer.at(8, 2) = answer.at(2, 8);
    answer.at(8, 6) = answer.at(6, 8);
    answer.at(8, 8) = 12. * iz / pow(l, 3.);
    answer.at(8, 12) = -6. * iz / pow(l, 2.);

    answer.at(9, 3) = answer.at(3, 9);
    answer.at(9, 5) = answer.at(5, 9);
    answer.at(9, 9) = 12 * iy / pow(l, 3.);
    answer.at(9, 11) = 6 * iy / pow(l, 2.);

    answer.at(10, 4) = answer.at(4, 10.);
    answer.at(10, 10) = j / ( 2 * ( 1 + nu ) * l );

    answer.at(11, 3) = answer.at(3, 11);
    answer.at(11, 5) = answer.at(5, 11);
    answer.at(11, 9) = answer.at(9, 11);
    answer.at(11, 11) = 4 * iy / l;

    answer.at(12, 2) = answer.at(2, 12);
    answer.at(12, 6) = answer.at(6, 12);
    answer.at(12, 8) = answer.at(8, 12);
    answer.at(12, 12) = 4 * iz / l;

    answer.times(e);

    return;
}

void LatticeBeam3d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ].reset(new GaussIntegrationRule(1, this, 1, 3) );
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _3dLattice);
}



double LatticeBeam3d :: giveArea() {
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->area;
}


bool
LatticeBeam3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    int i, j;

    answer.resize(12, 12);
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

    return 1;
}

int
LatticeBeam3d :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->localCoordinateSystem;

    return 1;
}


double
LatticeBeam3d :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->area * this->length;
}


void
LatticeBeam3d ::   giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

void
LatticeBeam3d :: initializeFrom(InputRecord &ir, int priority)
{
    ParameterManager &ppm =  this->giveDomain()->elementPPM;
    LatticeStructuralElement :: initializeFrom(ir, priority);
    PM_UPDATE_PARAMETER(diameter, ppm, ir, this->number, IPK_LatticeBeam3d_diameter, priority) ;
}

void
LatticeBeam3d :: postInitialize()
{
    ParameterManager &ppm =  this->giveDomain()->elementPPM;
    this->LatticeStructuralElement :: postInitialize();
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_LatticeBeam3d_diameter) ;
}

int
LatticeBeam3d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer.resize(3);
    answer = this->globalCentroid;

    return 1;
}


void
LatticeBeam3d :: computeGeometryProperties()
{
    //coordinates of the two nodes
    Node *nodeA, *nodeB;
    FloatArray coordsA(3), coordsB(3);

    nodeA  = this->giveNode(1);
    nodeB  = this->giveNode(2);

    for ( int i = 0; i < 3; i++ ) {
        coordsA.at(i + 1) =  nodeA->giveCoordinate(i + 1);
        coordsB.at(i + 1) =  nodeB->giveCoordinate(i + 1);
    }

    //Construct an initial temporary local coordinate system
    FloatArray s(3), t(3);

    //Calculate normal vector
    this->normal.resize(3);
    for ( int i = 0; i < 3; i++ ) {
        this->normal.at(i + 1) = coordsB.at(i + 1) - coordsA.at(i + 1);
    }

    this->length  = sqrt(pow(normal.at(1), 2.) + pow(normal.at(2), 2.) + pow(normal.at(3), 2.) );

    // Compute midpoint
    this->midPoint.resize(3);
    for ( int i = 0; i < 3; i++ ) {
        this->midPoint.at(i + 1) = 0.5 * ( coordsB.at(i + 1) + coordsA.at(i + 1) );
    }

    for ( int i = 0; i < 3; i++ ) {
        this->normal.at(i + 1) /= length;
    }

    this->globalCentroid = this->midPoint;

    computeCrossSectionProperties();

    this->geometryFlag = 1;

    return;
}

void
LatticeBeam3d :: computeCrossSectionProperties() {
    //Construct two perpendicular axis so that n is normal to the plane which they create
    //Check, if one of the components of the normal-direction is zero
    FloatArray s(3), t(3);
    if ( this->normal.at(1) == 0 ) {
        s.at(1) = 0.;
        s.at(2) = this->normal.at(3);
        s.at(3) = -this->normal.at(2);
    } else if ( this->normal.at(2) == 0 ) {
        s.at(1) = this->normal.at(3);
        s.at(2) = 0.;
        s.at(3) = -this->normal.at(1);
    } else {
        s.at(1) = this->normal.at(2);
        s.at(2) = -this->normal.at(1);
        s.at(3) = 0.;
    }

    s.normalize();

    t.beVectorProductOf(this->normal, s);
    t.normalize();

    this->localCoordinateSystem.resize(3, 3);
    this->localCoordinateSystem.zero();

    //Set up rotation matrix
    for ( int i = 1; i <= 3; i++ ) {
        this->localCoordinateSystem.at(1, i) = this->normal.at(i);
        this->localCoordinateSystem.at(2, i) = s.at(i);
        this->localCoordinateSystem.at(3, i) = t.at(i);
    }
}

void
LatticeBeam3d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->giveLatticeStress3d(strain, gp, tStep);
}


void
LatticeBeam3d :: giveInternalForcesVector(FloatArray &answer,
                                          TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatArray u;
    FloatMatrix stiffness;

    FloatArray strain(6), stress(6);

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    //    printf("Check what u is.\n");
    //    u.printYourself();

    this->computeStiffnessMatrix(stiffness, ElasticStiffness, tStep);

    // zero answer will resize accordingly when adding first contribution
    answer.clear();
    answer.beProductOf(stiffness, u);
    //    printf("answer.at(1) = %e, answer.at(7) = %e\n", answer.at(1), answer.at(7));

    double area = pow(this->diameter / 2., 2.) * myPi;
    //Apply yield limit to axial component
    strain.zero();
    strain.at(1) = ( u.at(7) - u.at(1) ) / this->giveLength();
    this->computeStressVector(stress, strain, gp, tStep);
    answer.at(1) = -stress.at(1) * area;
    answer.at(7) = -answer.at(1);
    //     printf("answerNew.at(1) = %e, answerNew.at(7) = %e\n", answer.at(1), answer.at(7));

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
    return;
}

#ifdef __OOFEG

void
LatticeBeam3d :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
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



void LatticeBeam3d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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


void LatticeBeam3d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
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
