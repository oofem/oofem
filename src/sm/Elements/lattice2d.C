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

#include "../sm/Elements/lattice2d.h"
#include "../sm/Materials/latticematstatus.h"
#include "../sm/Elements/latticestructuralelement.h"
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

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "../sm/Materials/structuralmaterial.h"
#endif

namespace oofem {
REGISTER_Element(Lattice2d);

Lattice2d :: Lattice2d(int n, Domain *aDomain) : LatticeStructuralElement(n, aDomain)
{
    numberOfDofMans = 2;

    length = 0.;
    pitch = 10.;  // a dummy value
    couplingNumbers.zero();
}

Lattice2d :: ~Lattice2d()
{ }


int
Lattice2d :: giveCrackFlag()
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    int crackFlag = 0;
    crackFlag = status->giveCrackFlag();

    return crackFlag;
}


double
Lattice2d :: giveCrackWidth()
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    double crackWidth = 0;
    crackWidth = status->giveCrackWidth();

    return crackWidth;
}

double
Lattice2d :: giveOldCrackWidth()
{
    LatticeMaterialStatus *status;

    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    GaussPoint *gp = iRule->getIntegrationPoint(0);
    status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    double crackWidth = 0;
    crackWidth = status->giveOldCrackWidth();

    return crackWidth;
}



double
Lattice2d :: giveDissipation()
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    double dissipation = 0;
    dissipation = status->giveDissipation();

    return dissipation;
}


double
Lattice2d :: giveDeltaDissipation()
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    double deltaDissipation = 0;
    deltaDissipation = status->giveDeltaDissipation();

    return deltaDissipation;
}

void
Lattice2d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
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
    answer.at(3, 3) = -this->width / sqrt(12.);
    answer.at(3, 4) = 0.;
    answer.at(3, 5) = 0.;
    answer.at(3, 6) = this->width / sqrt(12.);

    answer.times(1. / l);
}


void
Lattice2d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                    TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    double dV;
    FloatMatrix d, bj, dbj;
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
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _2dLattice);
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
Lattice2d :: computeVolumeAround(GaussPoint *gp)
{
    double area = this->width * this->thickness;
    double weight  = gp->giveWeight();
    return weight * 0.5 * this->giveLength() * area;
}


void
Lattice2d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, R_w
    };
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

double
Lattice2d :: giveNormalStress()
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    double normalStress = 0;
    normalStress = status->giveNormalStress();

    return normalStress;
}


double
Lattice2d :: giveOldNormalStress()
{
    LatticeMaterialStatus *status;

    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    GaussPoint *gp = iRule->getIntegrationPoint(0);
    status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    double normalStress = 0;
    normalStress = status->giveOldNormalStress();

    return normalStress;
}

int
Lattice2d :: hasBeenUpdated()
{
    LatticeMaterialStatus *status;

    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    GaussPoint *gp = iRule->getIntegrationPoint(0);
    status = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
    int updateFlag = 0;
    updateFlag = status->hasBeenUpdated();

    return updateFlag;
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
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    IR_GIVE_OPTIONAL_FIELD(ir, thickness, _IFT_Lattice2d_thick);

    IR_GIVE_OPTIONAL_FIELD(ir, width, _IFT_Lattice2d_width);

    IR_GIVE_OPTIONAL_FIELD(ir, gpCoords, _IFT_Lattice2d_gpcoords);

    couplingFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, couplingFlag, _IFT_Lattice2d_couplingflag);

    couplingNumbers.resize(1);
    couplingNumbers.zero();
    if ( couplingFlag == 1 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, couplingNumbers.at(1), _IFT_Lattice2d_couplingnumber);
    }

    // first call parent
    return LatticeStructuralElement :: initializeFrom(ir);
}


int
Lattice2d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    answer.resize(3);
    answer.at(1) = this->gpCoords.at(1);
    answer.at(2) = this->gpCoords.at(2);

    return 1;
}


void
Lattice2d :: giveGpCoordinates(FloatArray &answer)
{
    answer.resize(3);
    answer.at(1) = this->gpCoords.at(1);
    answer.at(2) = this->gpCoords.at(2);

    return;
}


contextIOResultType Lattice2d :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores =  LatticeStructuralElement :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( mode & CM_Definition ) ) {

        if ( !stream.write(width) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.write(thickness) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.write(couplingFlag) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = couplingNumbers.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = gpCoords.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}



contextIOResultType Lattice2d :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = LatticeStructuralElement :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
     
        if ( !stream.read(width) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(thickness) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(couplingFlag) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = couplingNumbers.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = gpCoords.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

    }

    return CIO_OK;
}

#ifdef __OOFEG

void
Lattice2d :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc, tStep);
        this->drawRawCrossSections(gc, tStep);
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

void
Lattice2d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    WCRec p [ 2 ]; /* points */
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

void
Lattice2d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    double defScale = gc.getDefScale();

    WCRec p [ 2 ]; /* poin */
    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.;

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.;

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void
Lattice2d :: drawSpecial(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 2 ];
    GraphicObj *tr;
    GaussPoint *gp;
    FloatArray crackStatuses, cf;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarType() == IST_CrackState ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        this->giveIPValue(crackStatuses, gp, IST_CrackStatuses, tStep);
        if ( crackStatuses(0) == 1. || crackStatuses(0) == 2. || crackStatuses(0) == 3 || crackStatuses(0) == 4 ) {
	  FloatArray coords;
	  this->giveCrossSectionCoordinates(coords);

	    p [ 0 ].x = ( FPNum ) coords.at(1);
	    p [ 0 ].y = ( FPNum ) coords.at(2);
	    p [ 0 ].z = ( FPNum ) coords.at(3);
	    p [ 1 ].x = ( FPNum ) coords.at(4);
	    p [ 1 ].y = ( FPNum ) coords.at(5);
	    p [ 1 ].z = ( FPNum ) coords.at(6);


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


            tr = CreateLine3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
	    EGAttachObject(tr, ( EObjectP ) this);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}

void
Lattice2d :: drawRawCrossSections(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getCrossSectionColor() );
    EASValsSetLayer(OOFEG_RAW_CROSSSECTION_LAYER);

    FloatArray coords;
    this->giveCrossSectionCoordinates(coords);

    p [ 0 ].x = ( FPNum ) coords.at(1);
    p [ 0 ].y = ( FPNum ) coords.at(2);
    p [ 0 ].z = ( FPNum ) coords.at(3);
    p [ 1 ].x = ( FPNum ) coords.at(4);
    p [ 1 ].y = ( FPNum ) coords.at(5);
    p [ 1 ].z = ( FPNum ) coords.at(6);

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void
Lattice2d :: giveCrossSectionCoordinates(FloatArray &coords)
{
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
