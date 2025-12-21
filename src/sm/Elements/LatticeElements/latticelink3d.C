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
#include "../sm/Elements/LatticeElements/latticelink3d.h"
#include "../sm/Materials/LatticeMaterials/latticematstatus.h"
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
#include "../sm/Materials/structuralmaterial.h"
#include "sm/CrossSections/latticecrosssection.h"
#include "parametermanager.h"
#include "paramkey.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(LatticeLink3d);
ParamKey LatticeLink3d::IPK_LatticeLink3d_length("length");
ParamKey LatticeLink3d::IPK_LatticeLink3d_diameter("diameter");
ParamKey LatticeLink3d::IPK_LatticeLink3d_dirvector("dirvector");
ParamKey LatticeLink3d::IPK_LatticeLink3d_l_end("l_end");

LatticeLink3d :: LatticeLink3d(int n, Domain *aDomain) : LatticeStructuralElement(n, aDomain)
{
    numberOfDofMans     = 2;
    geometryFlag = 0;
}

LatticeLink3d :: ~LatticeLink3d()
{}

double
LatticeLink3d :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    //Returns artifical volume (bond area times bond length) so that general parts of post processing work
    return pow(this->bondLength, 2.) * this->bondDiameter * M_PI;
}


void
LatticeLink3d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    //Assemble Bmatrix based on three rigid arm components
    //rigid.at(1) (tangential), rigid.at(2) (lateral), rigid.at(3) (lateral)
    answer.resize(6, 12);
    answer.zero();

    //Normal displacement jump in x-direction
    //First node
    answer.at(1, 1) = -1.;
    answer.at(1, 5) = -this->rigid.at(3);
    answer.at(1, 6) = this->rigid.at(2);
    //Second node
    answer.at(1, 7) = 1.;

    //Shear displacement jump in y-plane
    //first node
    answer.at(2, 2) = -1.;
    answer.at(2, 4) = this->rigid.at(3);
    answer.at(2, 6) = -this->rigid.at(1);
    //Second node
    answer.at(2, 8) = 1.;

    //Shear displacement jump in z-plane
    //first node
    answer.at(3, 3) = -1.;
    answer.at(3, 4) = -this->rigid.at(2);
    answer.at(3, 5) = this->rigid.at(1);
    //Second node
    answer.at(3, 9) =  1.;

    //Rotation around x-axis
    //First node
    answer.at(4, 4) = -1.;
    //Second node
    answer.at(4, 10) = 1.;

    //Rotation around y-axis
    //First node
    answer.at(5, 5) = -1.;
    //Second node
    answer.at(5, 11) = 1.;

    //Rotation around z-axis
    //First node
    answer.at(6, 6) = -1.;
    //Second node
    answer.at(6, 12) = 1.;

    return;
}

void
LatticeLink3d :: giveGPCoordinates(FloatArray &coords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }
    coords.resize(3);
    coords = this->globalCentroid;
    return;
}


void
LatticeLink3d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                        TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    FloatMatrix d, b, bt, db;
    FloatArray u, strain;

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();


    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    this->computeBmatrixAt(gp, b);
    bt.beTranspositionOf(b);

    if ( !this->isActivated(tStep) ) {
        strain.resize(StructuralMaterial :: giveSizeOfVoigtSymVector(gp->giveMaterialMode() ) );
        strain.zero();
    }
    strain.beProductOf(b, u);


    answer.resize(12, 12);
    answer.zero();

    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);


    db.beProductOf(d, b);
    answer.beProductOf(bt, db);

    //Introduce integration of bond strength
    double area = this->computeVolumeAround(gp) / this->giveLength();
    answer.times(area);

    return;
}

void LatticeLink3d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ].reset(new GaussIntegrationRule(1, this, 1, 3) );
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _3dLattice);
}




bool
LatticeLink3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
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
LatticeLink3d :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->localCoordinateSystem;

    return 1;
}

double
LatticeLink3d :: giveBondLength()
{
    return this->bondLength;
}

double
LatticeLink3d :: giveBondEndLength()
{
    return this->bondEndLength;
}


double
LatticeLink3d :: giveBondDiameter()
{
    return this->bondDiameter;
}




void
LatticeLink3d ::   giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

void
LatticeLink3d :: initializeFrom(InputRecord &ir, int priority)
{
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    // first call parent
    LatticeStructuralElement :: initializeFrom(ir, priority);

    PM_UPDATE_PARAMETER(bondLength, ppm, ir, this->number, IPK_LatticeLink3d_length, priority) ;
    PM_UPDATE_PARAMETER(bondDiameter, ppm, ir, this->number, IPK_LatticeLink3d_diameter, priority) ;
    PM_UPDATE_PARAMETER(directionVector, ppm, ir, this->number, IPK_LatticeLink3d_dirvector, priority) ;
    PM_UPDATE_PARAMETER(bondEndLength, ppm, ir, this->number, IPK_LatticeLink3d_l_end, priority) ;

}

void
LatticeLink3d :: postInitialize()
{
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    LatticeStructuralElement :: postInitialize();
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_LatticeLink3d_length) ;
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_LatticeLink3d_diameter) ;
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_LatticeLink3d_dirvector) ;
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_LatticeLink3d_l_end) ;
}

int
LatticeLink3d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer.resize(3);
    answer = this->globalCentroid;

    return 1;
}


void
LatticeLink3d :: computeGeometryProperties()
{
    //coordinates of the two nodes
    Node *nodeA, *nodeB;
    FloatArray coordsA(3), coordsB(3);

    //Order of nodes. However, might not matter.
    //Reinforcement node
    nodeA  = this->giveNode(1);
    //Lattice node
    nodeB  = this->giveNode(2);

    //Calculate components of distance from reinforcement node to lattice node.

    for ( int i = 0; i < 3; i++ ) {
        coordsA.at(i + 1) =  nodeA->giveCoordinate(i + 1);
        coordsB.at(i + 1) =  nodeB->giveCoordinate(i + 1);
    }

    FloatArray rigidGlobal(3);

    //Calculate normal vector
    for ( int i = 0; i < 3; i++ ) {
        rigidGlobal.at(i + 1) = coordsB.at(i + 1) - coordsA.at(i + 1);
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
        this->globalCentroid.at(i) = nodeA->giveCoordinate(i);
        ;
    }

    this->geometryFlag = 1;

    return;
}
void
LatticeLink3d :: saveContext(DataStream &stream, ContextMode mode)
{
    LatticeStructuralElement :: saveContext(stream, mode);
}


void
LatticeLink3d :: restoreContext(DataStream &stream, ContextMode mode)
{
    LatticeStructuralElement :: restoreContext(stream, mode);
}


void
LatticeLink3d :: giveInternalForcesVector(FloatArray &answer,
                                          TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix b, bt;
    FloatArray u, stress(6), strain(6);

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    for ( GaussPoint *gp: * this->giveDefaultIntegrationRulePtr() ) {
        LatticeMaterialStatus *matStat = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
        this->computeBmatrixAt(gp, b);
        bt.beTranspositionOf(b);

        if ( useUpdatedGpRecord == 1 ) {
            stress = matStat->giveLatticeStress();
        } else {
            if ( !this->isActivated(tStep) ) {
                strain.resize(StructuralMaterial :: giveSizeOfVoigtSymVector(gp->giveMaterialMode() ) );
                strain.zero();
            }
            strain.beProductOf(b, u);
            this->computeStressVector(stress, strain, gp, tStep);
        }

        // updates gp stress and strain record  acording to current
        // increment of displacement
        if ( stress.giveSize() == 0 ) {
            break;
        }

        // now every gauss point has real stress vector
        // compute nodal representation of internal forces using f = B^T*Sigma dV
        //	double dV = this->computeVolumeAround(gp);
        //		double dV = 1.;
        if ( stress.giveSize() == 6 ) {
            // It may happen that e.g. plane strain is computed
            // using the default 3D implementation. If so,
            // the stress needs to be reduced.
            // (Note that no reduction will take place if
            //  the simulation is actually 3D.)
            FloatArray stressTemp;
            StructuralMaterial :: giveReducedSymVectorForm(stressTemp, stress, gp->giveMaterialMode() );

            answer.beProductOf(bt, stressTemp);
        } else {
            answer.beProductOf(bt, stress);
        }

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


double
LatticeLink3d :: giveLength()
{
    //returns the bond length, not the length between the two nodes
    return this->bondLength;
}

void
LatticeLink3d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give3dStiffnessMatrix(rMode, gp, tStep);
}

void
LatticeLink3d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->giveLatticeStress3d(strain, gp, tStep);
}



#ifdef __OOFEG

void
LatticeLink3d :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
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



void LatticeLink3d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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

void LatticeLink3d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
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
