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
#include "../sm/Elements/Interfaces/bondlink2dtruss.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "mathfem.h"
#include "../sm/Elements/structuralelement.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "../sm/Materials/structuralmaterial.h"
#include "sm/CrossSections/structuralinterfacecrosssection.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(BondLink2dTruss);

BondLink2dTruss::BondLink2dTruss(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
}


double
BondLink2dTruss::computeVolumeAround(GaussPoint *aGaussPoint)
{
    // Bond area times bond length
    return pow(this->bondLength, 2.) * this->bondDiameter * M_PI;
}


void
BondLink2dTruss::computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    //Assemble Bmatrix in the local system (displacements are rotated by computeGtoLRotationMatrix)
    answer.resize(2, 4);
    answer.zero();

    //Slip: jump along the bar
    //First node
    answer.at(1, 1) = 1.;
    //Second node
    answer.at(1, 3) = -1.;

    //Lateral jump perpendicular to the bar
    //first node
    answer.at(2, 2) = 1.;
    //Second node
    answer.at(2, 4) = -1.;

    return;
}

void
BondLink2dTruss::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                        TimeStep *tStep)
{
    FloatMatrix d, b, bt, db;
    FloatArray u, slip;

    this->computeVectorOf(VM_Total, tStep, u);
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    answer.clear();

    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    this->computeBmatrixAt(gp, b);
    bt.beTranspositionOf(b);

    if ( !this->isActivated(tStep) ) {
        slip.resize( StructuralMaterial::giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        slip.zero();
    }
    slip.beProductOf(b, u);

    answer.resize(4, 4);
    answer.zero();

    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    db.beProductOf(d, b);
    answer.beProductOf(bt, db);

    // Scale by bond area
    double area = this->computeVolumeAround(gp) / this->giveLength();
    answer.times(area);


    return;
}

void BondLink2dTruss::computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, this->giveMaterialMode() );
}




bool
BondLink2dTruss::computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    int i, j;

    answer.resize(4, 4);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( i = 1; i <= 2; i++ ) {
        for ( j = 1; j <= 2; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 2, j + 2) = lcs.at(i, j);
        }
    }

    return 1;
}


int
BondLink2dTruss::giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->localCoordinateSystem;

    return 1;
}

double
BondLink2dTruss::giveBondLength()
{
    return this->bondLength;
}



double
BondLink2dTruss::giveBondDiameter()
{
    return this->bondDiameter;
}


void
BondLink2dTruss::giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { D_u, D_v };
}

void
BondLink2dTruss::initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority)
{
    StructuralElement::initializeFrom(ir, priority);

    IR_GIVE_FIELD(ir, this->bondLength, _IFT_BondLink2dTruss_length);

    IR_GIVE_FIELD(ir, this->bondDiameter, _IFT_BondLink2dTruss_diameter);

    IR_GIVE_FIELD(ir, this->directionVector, _IFT_BondLink2dTruss_dirvector);
}


int
BondLink2dTruss::computeGlobalCoordinates(Coordinates &answer, const FloatArray &lcoords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer = this->globalCentroid;

    return 1;
}


void
BondLink2dTruss::computeGeometryProperties()
{
    //Order of nodes: truss node first, continuum node second. Both should be at the same point.
    Node *nodeB = this->giveNode(2);

    //Local system in the x-y plane: first axis along the bar, second axis perpendicular to it
    FloatArray normal(2), s(2);

    normal.at(1) = this->directionVector.at(1);
    normal.at(2) = this->directionVector.at(2);
    if ( normal.computeNorm() == 0. ) {
        OOFEM_ERROR("dirvector of element %d has no component in the x-y plane", this->giveNumber() );
    }
    normal.normalize();

    s.at(1) = -normal.at(2);
    s.at(2) = normal.at(1);

    this->localCoordinateSystem.resize(2, 2);
    this->localCoordinateSystem.zero();

    for ( int i = 1; i <= 2; i++ ) {
        this->localCoordinateSystem.at(1, i) = normal.at(i);
        this->localCoordinateSystem.at(2, i) = s.at(i);
    }

    this->globalCentroid.resize(3);
    for ( int i = 1; i <= 3; i++ ) {
        this->globalCentroid.at(i) = nodeB->giveCoordinate(i);
    }

    this->geometryFlag = 1;

    return;
}

void
BondLink2dTruss::saveContext(DataStream &stream, ContextMode mode)
{
    StructuralElement::saveContext(stream, mode);
}


void
BondLink2dTruss::restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralElement::restoreContext(stream, mode);
}


void
BondLink2dTruss::giveInternalForcesVector(FloatArray &answer,
                                          TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix b;
    FloatArray u, stress, strain;

    this->computeVectorOf(VM_Total, tStep, u);
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }
    answer.clear();

    for ( GaussPoint *gp : * this->giveDefaultIntegrationRulePtr() ) {
        this->computeBmatrixAt(gp, b);

        if ( useUpdatedGpRecord == 1 ) {
            auto status = gp->giveMaterialStatus();
            StructuralMaterialStatus *matStat = dynamic_cast< StructuralMaterialStatus * >( status );
            if ( matStat ) {
                stress = matStat->giveStressVector();
            } else {
                StructuralInterfaceMaterialStatus *ms = static_cast< StructuralInterfaceMaterialStatus * >( status );
                stress = ms->giveTraction();
            }
            // the status holds the 3-component traction of the material; keep slip and lateral
            if ( stress.giveSize() > 2 ) {
                stress.resizeWithValues(2);
            }
        } else {
            if ( !this->isActivated(tStep) ) {
                strain.resize( StructuralMaterial::giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
                strain.zero();
            }
            strain.beProductOf(b, u);
            this->computeStressVector(stress, strain, gp, tStep);
        }

        if ( stress.giveSize() == 0 ) {
            break;
        }

        // f = B^T * stress * bond area
        double dV = this->computeVolumeAround(gp) / this->giveLength();
        if ( stress.giveSize() == 6 ) {
            // Reduce a full 3d stress vector to the material mode.
            FloatArray stressTemp;
            StructuralMaterial::giveReducedSymVectorForm( stressTemp, stress, gp->giveMaterialMode() );
            answer.plusProduct(b, stressTemp, dV);
        } else {
            answer.plusProduct(b, stress, dV);
        }
    }

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}


double
BondLink2dTruss::giveLength()
{
    // Bond length
    return this->bondLength;
}

void
BondLink2dTruss::computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    //The material works in 3D (slip, lateral, lateral). Use the slip/lateral block.
    auto d = static_cast< StructuralInterfaceCrossSection * >( this->giveCrossSection() )->giveStiffnessMatrix_Eng_ntt(rMode, gp, tStep);
    answer.resize(2, 2);
    for ( int i = 1; i <= 2; i++ ) {
        for ( int j = 1; j <= 2; j++ ) {
            answer.at(i, j) = d.at(i, j);
        }
    }
}

void
BondLink2dTruss::computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    //Pad the 2D jump (slip, lateral) with a zero out-of-plane component for the 3D material
    FloatArrayF< 3 >jump = { strain.at(1), strain.at(2), 0. };
    auto traction = static_cast< StructuralInterfaceCrossSection * >( this->giveCrossSection() )->giveEngTraction_ntt(jump, gp, tStep);
    answer.resize(2);
    answer.at(1) = traction.at(1);
    answer.at(2) = traction.at(2);
}



#ifdef __OOFEG

void
BondLink2dTruss::drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
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



void BondLink2dTruss::drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    WCRec p[ 2 ];  /* points */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
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

void BondLink2dTruss::drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    double defScale = gc.getDefScale();

    WCRec p[ 2 ];  /* points */

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
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
