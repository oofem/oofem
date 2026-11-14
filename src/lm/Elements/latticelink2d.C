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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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
#include "Elements/latticelink2d.h"
#include "Materials/latticematstatus.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "Elements/latticestructuralelement.h"
#include "parametermanager.h"
#include "paramkey.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "sm/Materials/structuralmaterial.h"
#include "lm/CrossSections/latticecrosssection.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(LatticeLink2d);

ParamKey LatticeLink2d::IPK_LatticeLink2d_length("length");
ParamKey LatticeLink2d::IPK_LatticeLink2d_diameter("diameter");
ParamKey LatticeLink2d::IPK_LatticeLink2d_dirvector("dirvector");

LatticeLink2d :: LatticeLink2d(int n, Domain *aDomain) : LatticeStructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
}

LatticeLink2d :: ~LatticeLink2d()
{}

double
LatticeLink2d :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    // Bond area times bond length
    return pow(this->bondLength, 2.) * this->bondDiameter * M_PI;
}


void
LatticeLink2d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    //Assemble Bmatrix based on the rigid arm from node 1 to node 2
    //rigid.at(1) (along the bar), rigid.at(2) (lateral)
    answer.resize(3, 6);
    answer.zero();

    //Slip: jump along the bar
    //First node
    answer.at(1, 1) = -1.;
    answer.at(1, 3) = this->rigid.at(2);
    //Second node
    answer.at(1, 4) = 1.;

    //Lateral jump perpendicular to the bar
    //First node
    answer.at(2, 2) = -1.;
    answer.at(2, 3) = -this->rigid.at(1);
    //Second node
    answer.at(2, 5) = 1.;

    //Rotation around z-axis
    //First node
    answer.at(3, 3) = -1.;
    //Second node
    answer.at(3, 6) = 1.;
}

void
LatticeLink2d :: giveGpCoordinates(FloatArray &coords, GaussPoint *gp)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    coords = this->globalCentroid;
}

void
LatticeLink2d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                        TimeStep *tStep)
{
    FloatMatrix d, b, bt, db;

    answer.clear();

    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    this->computeBmatrixAt(gp, b);
    bt.beTranspositionOf(b);

    this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);

    db.beProductOf(d, b);
    answer.beProductOf(bt, db);

    // Scale by bond area
    double area = this->computeVolumeAround(gp) / this->giveLength();
    answer.times(area);
}

void LatticeLink2d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ].reset(new GaussIntegrationRule(1, this, 1, 3) );
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _2dLattice);
}

bool
LatticeLink2d :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;

    answer.resize(6, 6);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( int i = 1; i <= 2; i++ ) {
        for ( int j = 1; j <= 2; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
        }
    }
    //Rotations about z are the same in the global and local system
    answer.at(3, 3) = 1.;
    answer.at(6, 6) = 1.;

    return 1;
}


int
LatticeLink2d :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->localCoordinateSystem;

    return 1;
}

double
LatticeLink2d :: giveBondLength()
{
    return this->bondLength;
}


double
LatticeLink2d :: giveBondDiameter()
{
    return this->bondDiameter;
}


void
LatticeLink2d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, R_w
    };
}

void
LatticeLink2d :: initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority)
{
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    LatticeStructuralElement :: initializeFrom(ir, priority);

    PM_UPDATE_PARAMETER(bondLength, ppm, ir, this->number, IPK_LatticeLink2d_length, priority) ;
    PM_UPDATE_PARAMETER(bondDiameter, ppm, ir, this->number, IPK_LatticeLink2d_diameter, priority) ;
    PM_UPDATE_PARAMETER(directionVector, ppm, ir, this->number, IPK_LatticeLink2d_dirvector, priority) ;
}

void
LatticeLink2d :: postInitialize()
{
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    LatticeStructuralElement :: postInitialize();
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_LatticeLink2d_length) ;
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_LatticeLink2d_diameter) ;
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_LatticeLink2d_dirvector) ;
}


int
LatticeLink2d :: computeGlobalCoordinates(Coordinates &answer, const FloatArray &lcoords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->globalCentroid;

    return 1;
}


void
LatticeLink2d :: computeGeometryProperties()
{
    //Node 1 carries the rigid arm; node 2 is where the slip is evaluated.
    Node *nodeA = this->giveNode(1);
    Node *nodeB = this->giveNode(2);

    //Rigid arm from node 1 to node 2 (zero for coinciding nodes).
    FloatArray rigidGlobal(2);
    for ( int i = 1; i <= 2; i++ ) {
        rigidGlobal.at(i) = nodeB->giveCoordinate(i) - nodeA->giveCoordinate(i);
    }

    //Direction of the bar in the x-y plane (a third component of dirvector is ignored)
    if ( this->directionVector.giveSize() < 2 ) {
        OOFEM_ERROR("LatticeLink2d %d: dirvector needs at least 2 components.", this->giveGlobalNumber() );
    }
    FloatArray normal(2);
    normal.at(1) = this->directionVector.at(1);
    normal.at(2) = this->directionVector.at(2);
    if ( normal.computeNorm() < 1e-12 ) {
        OOFEM_ERROR("LatticeLink2d %d: dirvector has no component in the x-y plane.", this->giveGlobalNumber() );
    }
    normal.normalize();

    //Local system: first axis along the bar, second axis perpendicular to it in the x-y plane
    this->localCoordinateSystem.resize(2, 2);
    this->localCoordinateSystem.at(1, 1) = normal.at(1);
    this->localCoordinateSystem.at(1, 2) = normal.at(2);
    this->localCoordinateSystem.at(2, 1) = -normal.at(2);
    this->localCoordinateSystem.at(2, 2) = normal.at(1);

    // Rotate rigid arm vector into local coordinate system
    this->rigid.beProductOf(localCoordinateSystem, rigidGlobal);

    this->globalCentroid.resize(3);
    for ( int i = 1; i <= 3; i++ ) {
        this->globalCentroid.at(i) = nodeB->giveCoordinate(i);
    }

    this->geometryFlag = 1;
}

void
LatticeLink2d :: saveContext(DataStream &stream, ContextMode mode)
{
    LatticeStructuralElement :: saveContext(stream, mode);
}


void
LatticeLink2d :: restoreContext(DataStream &stream, ContextMode mode)
{
    LatticeStructuralElement :: restoreContext(stream, mode);
}


void
LatticeLink2d :: giveInternalForcesVector(FloatArray &answer,
                                          TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix b, bt;
    FloatArray u, stress(3), strain(3);

    this->computeVectorOf(VM_Total, tStep, u);
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    answer.clear();

    for ( GaussPoint *gp: * this->giveDefaultIntegrationRulePtr() ) {
        this->computeBmatrixAt(gp, b);
        bt.beTranspositionOf(b);

        if ( useUpdatedGpRecord == 1 ) {
            //The status stores the 3d lattice stress; slip, lateral and rotation are components 1, 2 and 6
            LatticeMaterialStatus *matStat = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
            const auto &stress3d = matStat->giveLatticeStress();
            stress.at(1) = stress3d.at(1);
            stress.at(2) = stress3d.at(2);
            stress.at(3) = stress3d.at(6);
        } else {
            strain.beProductOf(b, u);
            if ( !this->isActivated(tStep) ) {
                strain.zero();
            }
            this->computeStressVector(stress, strain, gp, tStep);
        }

        if ( stress.giveSize() == 0 ) {
            break;
        }

        // f = B^T * stress
        answer.beProductOf(bt, stress);

        // Scale by bond area
        double area = this->computeVolumeAround(gp) / this->giveLength();
        answer.times(area);
    }

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}


double
LatticeLink2d :: giveLength()
{
    // Bond length
    return this->bondLength;
}

void
LatticeLink2d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give2dStiffnessMatrix(rMode, gp, tStep);
}

void
LatticeLink2d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->giveLatticeStress2d(strain, gp, tStep);
}



#ifdef __OOFEG

void
LatticeLink2d :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
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
void LatticeLink2d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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

void LatticeLink2d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
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
