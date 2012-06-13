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

#include "domain.h"
#include "lattice2d.h"
#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "latticematstatus.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "flotarry.h"
#include "mathfem.h"
#include "latticestructuralelement.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "structuralmaterial.h"
#endif

namespace oofem {
Lattice2d :: Lattice2d(int n, Domain *aDomain) : LatticeStructuralElement(n, aDomain)
    // Constructor.
{
    numberOfDofMans     = 2;

    length              = 0.;
    pitch               = 10.;  // a dummy value
}

Lattice2d :: ~Lattice2d()
{}


int
Lattice2d :: giveCrackFlag()
{
    LatticeMaterialStatus *status;

    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    gp = iRule->getIntegrationPoint(0);
    Material *mat = this->giveMaterial();

    status = ( LatticeMaterialStatus * ) mat->giveStatus(gp);
    int crackFlag = 0;
    crackFlag = status->giveCrackFlag();

    return crackFlag;
}


double
Lattice2d :: giveCrackWidth()
{
    LatticeMaterialStatus *status;

    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    gp = iRule->getIntegrationPoint(0);
    Material *mat = this->giveMaterial();

    status = ( LatticeMaterialStatus * ) mat->giveStatus(gp);
    double crackWidth = 0;
    crackWidth = status->giveCrackWidth();

    return crackWidth;
}

double
Lattice2d :: giveDissipation()
{
    LatticeMaterialStatus *status;

    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    gp = iRule->getIntegrationPoint(0);
    Material *mat = this->giveMaterial();

    status = ( LatticeMaterialStatus * ) mat->giveStatus(gp);
    double dissipation = 0;
    dissipation = status->giveDissipation();

    return dissipation;
}


double
Lattice2d :: giveDeltaDissipation()
{
    LatticeMaterialStatus *status;

    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    gp = iRule->getIntegrationPoint(0);
    Material *mat = this->giveMaterial();

    status = ( LatticeMaterialStatus * ) mat->giveStatus(gp);
    double deltaDissipation = 0;
    deltaDissipation = status->giveDeltaDissipation();

    return deltaDissipation;
}

void
Lattice2d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    double l = this->giveLength();
    double x1, x2, y1, y2, xp, yp;

    //Coordinates of the nodes
    x1 = this->giveNode(1)->giveCoordinate(1);
    y1 = this->giveNode(1)->giveCoordinate(2);
    x2 = this->giveNode(2)->giveCoordinate(1);
    y2 = this->giveNode(2)->giveCoordinate(2);

    xp = this->gpCoords.at(1);
    yp = this->gpCoords.at(2);

    double ecc;

    double areaHelp = 0.5 * ( x1 * y2 + x2 * yp + xp * y1 - ( xp * y2 + yp * x1 + x2 * y1 ) );

    ecc = 2 * areaHelp / l;

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
    answer.at(2, 3) =  -l / 2.;
    answer.at(2, 4) = 0.;
    answer.at(2, 5) = 1.;
    answer.at(2, 6) = -l / 2.;

    answer.at(3, 1) = 0.;
    answer.at(3, 2) = 0.;
    answer.at(3, 3) = -pow(this->width, 2.) / ( sqrt(12.) * l );
    //  answer.at(3,3) = -this->width;
    answer.at(3, 4) = 0.;
    answer.at(3, 5) = 0.;

    answer.at(3, 6) = pow(this->width, 2.) / ( sqrt(12.) * l );
    //  answer.at(3,6) = this->width;

    answer.times(1. / l);
}


void
Lattice2d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                    TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    double dV;
    FloatMatrix d, bi, bj, dbj, dij;
    answer.resize(6, 6);
    answer.zero();
    this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bj);
    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    dV = this->computeVolumeAround( integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
    dbj.beProductOf(d, bj);
    answer.plusProductUnsym(bj, dbj, dV);
}


void Lattice2d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    // the gauss point is used only when methods from crosssection and/or material
    // classes are requested
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ 1 ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
    integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 1, _2dLattice);
}


bool
Lattice2d :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    double sine, cosine;
    answer.resize(6, 6);
    answer.zero();

    sine           = sin( this->givePitch() );
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
    return 1;
}


double
Lattice2d :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    double area = this->width * this->thickness;
    double weight  = aGaussPoint->giveWeight();
    return weight * 0.5 * this->giveLength() * area;
}


void
Lattice2d ::   giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const {
    answer.resize(3);
    answer.at(1) = D_u;
    answer.at(2) = D_v;
    answer.at(3) = R_w;
}



double Lattice2d :: giveLength()
// Returns the length of the receiver.
{
    double dx, dy;
    Node *nodeA, *nodeB;

    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
        length  = sqrt(dx * dx + dy * dy);
    }

    return length;
}


double Lattice2d :: givePitch()
// Returns the pitch of the receiver.
{
    double xA, xB, yA, yB;
    Node *nodeA, *nodeB;

    if ( pitch == 10. ) {            // 10. : dummy initialization value
        nodeA  = this->giveNode(1);
        nodeB  = this->giveNode(2);
        xA     = nodeA->giveCoordinate(1);
        xB     = nodeB->giveCoordinate(1);
        yA     = nodeA->giveCoordinate(2);
        yB     = nodeB->giveCoordinate(2);
        pitch  = atan2(yB - yA, xB - xA);
    }

    return pitch;
}


int
Lattice2d :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    double sine, cosine;

    answer.resize(3, 3);
    answer.zero();

    sine           = sin( this->givePitch() );
    cosine         = cos(pitch);

    answer.at(1, 1) = cosine;
    answer.at(1, 2) = sine;
    answer.at(2, 1) = -sine;
    answer.at(2, 2) = cosine;
    answer.at(3, 3) = 1.0;

    return 1;
}

IRResultType
Lattice2d :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    // first call parent
    LatticeStructuralElement :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, thickness, IFT_Lattice2d_thick, "thick"); // Macro

    IR_GIVE_OPTIONAL_FIELD(ir, width, IFT_Lattice2d_width, "width"); // Macro

    IR_GIVE_OPTIONAL_FIELD(ir, gpCoords, IFT_Lattice2d_gpcoords, "gpcoords"); // Macro

    this->computeGaussPoints();
    return IRRT_OK;
}


int
Lattice2d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    answer.resize(3);
    answer.at(1) = this->gpCoords.at(1);
    answer.at(2) = this->gpCoords.at(2);

    return 1;
}


#ifdef __OOFEG

void
Lattice2d :: drawYourself(oofegGraphicContext &gc)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc);
    } else if ( mode == OGC_deformedGeometry ) {
        this->drawDeformedGeometry(gc, DisplacementVector);
    } else if ( mode == OGC_eigenVectorGeometry ) {
        this->drawDeformedGeometry(gc, EigenVector);
    } else if ( mode == OGC_scalarPlot ) {
        this->drawScalar(gc);
    } else if ( mode == OGC_elemSpecial ) {
        this->drawSpecial(gc);
    } else {
        _error("drawYourself : unsupported mode");
    }
}




void Lattice2d :: drawRawGeometry(oofegGraphicContext &gc)
{
    GraphicObj *go;

    WCRec p [ 2 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void Lattice2d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();

    WCRec p [ 2 ]; /* poin */
    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = 0.;

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].z = 0.;

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}

void
Lattice2d :: drawSpecial(oofegGraphicContext &gc)
{
    WCRec l [ 2 ];
    GraphicObj *tr;
    Material *mat = ( StructuralMaterial * ) this->giveMaterial();
    GaussPoint *gp;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    FloatArray crackStatuses, cf;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarType() == IST_CrackState ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        mat->giveIPValue(crackStatuses, gp, IST_CrackStatuses, tStep);
        if ( crackStatuses(0) == 1. || crackStatuses(0) == 2. || crackStatuses(0) == 3 || crackStatuses(0) == 4 ) {
            double x1, y1, x2, y2;
            x1 = this->giveNode(1)->giveCoordinate(1);
            y1 = this->giveNode(1)->giveCoordinate(2);
            x2 = this->giveNode(2)->giveCoordinate(1);
            y2 = this->giveNode(2)->giveCoordinate(2);

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
                EASValsSetColor( gc.getActiveCrackColor() );
            } else if ( crackStatuses(0) == 2. ) {
                EASValsSetColor( gc.getCrackPatternColor() );
            } else if ( crackStatuses(0) == 3. ) {
                EASValsSetColor( gc.getActiveCrackColor() );
            } else if ( crackStatuses(0) == 4. ) {
                EASValsSetColor( gc.getActiveCrackColor() );
            }


            tr = CreateLine3D(l);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}
#endif
} // end namespace oofem
