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

#include "lattice3dboundary_mt.h"
#include "tm/Materials/transportmaterial.h"
#include "tm/Materials/LatticeMaterials/latticetransmat.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "load.h"
#include "classfactory.h"




#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Element(Lattice3dboundary_mt);

Lattice3dboundary_mt :: Lattice3dboundary_mt(int n, Domain *aDomain, ElementMode em) :
    Lattice3d_mt(n, aDomain, em)
    // Constructor.
{
    numberOfDofMans  = 3;
    length = -1.0;
}

Lattice3dboundary_mt :: ~Lattice3dboundary_mt()
// Destructor
{}

void
Lattice3dboundary_mt ::   giveDofManDofIDMask(int inode, IntArray &answer) const {
    // returns DofId mask array for inode element node.
    // DofId mask array determines the dof ordering requsted from node.
    // DofId mask array contains the DofID constants (defined in cltypes.h)
    // describing physical meaning of particular DOFs.

    //The 3rd node contains gradients in x, y and z direction
    if ( inode == 3 ) {
        answer.resize(3);
        answer.at(1) = D_u;
        answer.at(2) = D_v;
        answer.at(3) = D_w;
    } else {
        answer.resize(1);
        answer.at(1) = P_f;
    }
}



void
Lattice3dboundary_mt :: computeFlow(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1);
    answer.zero();

    IntArray projectionComponent(3);
    projectionComponent.zero();

    FloatArray r;
    IntArray dofid;
    double dV;
    double length = giveLength();
    double k;
    answer.resize(1);

    k = static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicValue(Conductivity, gp, tStep);
    dV = this->computeVolumeAround(gp);

    this->computeVectorOf(VM_Total, tStep, r);

    for ( int i = 0; i < 2; i++ ) {
        if ( i + 1 == 1 ) {
            if ( location.at(1) != 0 ) {
                giveSwitches( projectionComponent, location.at(1) );
            }
        } else {
            if ( location.at(2) != 0 ) {
                giveSwitches( projectionComponent, location.at(2) );
            }
        }
        r.at(i + 1) = r.at(i + 1) + r.at(3) * projectionComponent.at(1)   + r.at(4) * projectionComponent.at(2) + r.at(5) * projectionComponent.at(3);
        projectionComponent.zero();
    }

    double flow;
    double dP;
    dP = r.at(2) - r.at(1);
    flow = k * ( dP ) * dV / pow(length, 2.);
    answer.at(1) = fabs(flow);

    if ( !isActivated(tStep) ) {
        answer.zero();
    }
}



void
Lattice3dboundary_mt :: computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // The 5x5 dof layout is {P_f node1, P_f node2, 3 macro-gradient dofs on
    // the control node}. Mirroring computeConductivityMatrix and computeFlow,
    // the effective endpoint pressure on the periodic-image side is
    //     p̃_i = p_i + ∇p · shift_i,
    // so storage rate ∂p̃/∂t couples through the gradient DOFs. The capacity
    // is therefore C = T^T M T, where M is the standard 2-node mass matrix
    // and T = [I_2 | S] carries the periodic-shift switches in S. Without
    // this projection the gradient DOFs would receive no time term while
    // contributing to the spatial flow — inconsistent with conductivity.
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    const double c = static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicValue(Capacity, gp, tStep);
    const double V = this->computeVolumeAround(gp);

    FloatMatrix M(2, 2);
    if ( this->lumpedCapacity ) {
        // Strut-level lumping: each endpoint carries c·V/(2·dim).
        M.at(1, 1) = 1.;
        M.at(2, 2) = 1.;
        M.times(c * V / ( 2.0 * this->dimension ));
    } else {
        M.at(1, 1) = 2.;
        M.at(1, 2) = 1.;
        M.at(2, 1) = 1.;
        M.at(2, 2) = 2.;
        M.times(c * V / ( 6.0 * this->dimension ));
    }

    IntArray sw1(3), sw2(3);
    sw1.zero();
    sw2.zero();
    if ( location.at(1) != 0 ) {
        giveSwitches(sw1, location.at(1));
    }
    if ( location.at(2) != 0 ) {
        giveSwitches(sw2, location.at(2));
    }

    FloatMatrix T(2, 5);
    T.zero();
    T.at(1, 1) = 1.;
    T.at(2, 2) = 1.;
    T.at(1, 3) = sw1.at(1);
    T.at(1, 4) = sw1.at(2);
    T.at(1, 5) = sw1.at(3);
    T.at(2, 3) = sw2.at(1);
    T.at(2, 4) = sw2.at(2);
    T.at(2, 5) = sw2.at(3);

    FloatMatrix MT;
    MT.beProductOf(M, T);          // 2x5
    answer.beTProductOf(T, MT);    // 5x5  =  T^T (M T)
}


void
Lattice3dboundary_mt :: initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority)
{

    this->Lattice3d_mt :: initializeFrom(ir, priority);

    location.resize(2);
    IR_GIVE_FIELD(ir, location, _IFT_Lattice3dboundary_mt_location); // Macro

}


void
Lattice3dboundary_mt :: computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
{
    double l = giveLength();

    //    Assemble Gradient Matrix used to compute temperature gradient

    answer.resize(1, 3);
    answer.zero();
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 1.;

    answer.times(1. / l);


    return;
}


void
Lattice3dboundary_mt :: computeBCSubVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode, int indx)
{
    FloatArray vec;

    answer.resize(this->giveNumberOfDofManagers() + 2);
    answer.zero();

    // loop over boundary load array
    int nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nLoads; i++ ) {
        int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int id = boundaryLoadArray.at(i * 2);
        Load *load = domain->giveLoad(n);
        bcGeomType ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeBCSubVectorAt(vec, load, id, tStep, mode, indx);
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceBCSubVectorAt(vec, load, id, tStep, mode, indx);
        } else {
            OOFEM_ERROR("unsupported bc type encountered");
        }

        answer.add(vec);
    }
}



void
Lattice3dboundary_mt :: computeInternalForcesVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    FloatArray tmp;
    FloatArray unknowns;
    FloatMatrix s;
    this->computeVectorOf(mode, tStep, unknowns);

    this->computeConductivityMatrix(s, Conductivity, tStep);
    answer.beProductOf(s, unknowns);

    this->computeInternalSourceRhsVectorAt(tmp, tStep, mode);
    answer.subtract(tmp);

    FloatMatrix bc_tangent;
    this->computeBCMtrxAt(bc_tangent, tStep, VM_Total);
    if ( bc_tangent.isNotEmpty() ) {
        tmp.beProductOf(bc_tangent, unknowns);
        answer.add(tmp);
    }

    IntArray projectionComponent(3);
    giveSwitches( projectionComponent, location.at(1) );


    IntArray projectionComponentNodeOne(3);
    projectionComponentNodeOne.zero();
    if ( location.at(1) != 0 ) {
        giveSwitches( projectionComponentNodeOne, location.at(1) );
    }

    IntArray projectionComponentNodeTwo(3);
    projectionComponentNodeTwo.zero();

    if ( location.at(2) != 0 ) {
        giveSwitches( projectionComponentNodeTwo, location.at(2) );
    }

    answer.at(3) = projectionComponentNodeOne.at(1) * answer.at(1) + projectionComponentNodeTwo.at(1) * answer.at(2);
    answer.at(4) = projectionComponentNodeOne.at(2) * answer.at(1) + projectionComponentNodeTwo.at(2) * answer.at(2);
    answer.at(5) = projectionComponentNodeOne.at(3) * answer.at(1) + projectionComponentNodeTwo.at(3) * answer.at(2);
}


void
Lattice3dboundary_mt :: computeHomogenisedInternalForcesVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode, FloatArray &unknowns)
{
    FloatArray tmp;

    FloatMatrix s;


    this->computeConductivityMatrix(s, Conductivity, tStep);
    answer.beProductOf(s, unknowns);

    this->computeInternalSourceRhsVectorAt(tmp, tStep, mode);
    answer.subtract(tmp);

    FloatMatrix bc_tangent;
    this->computeBCMtrxAt(bc_tangent, tStep, VM_Total);
    if ( bc_tangent.isNotEmpty() ) {
        tmp.beProductOf(bc_tangent, unknowns);
        answer.add(tmp);
    }

    IntArray projectionComponent(3);
    giveSwitches( projectionComponent, location.at(1) );


    IntArray projectionComponentNodeOne(3);
    projectionComponentNodeOne.zero();
    if ( location.at(1) != 0 ) {
        giveSwitches( projectionComponentNodeOne, location.at(1) );
    }

    IntArray projectionComponentNodeTwo(3);
    projectionComponentNodeTwo.zero();

    if ( location.at(2) != 0 ) {
        giveSwitches( projectionComponentNodeTwo, location.at(2) );
    }

    answer.at(3) = projectionComponentNodeOne.at(1) * answer.at(1) + projectionComponentNodeTwo.at(1) * answer.at(2);
    answer.at(4) = projectionComponentNodeOne.at(2) * answer.at(1) + projectionComponentNodeTwo.at(2) * answer.at(2);
    answer.at(5) = projectionComponentNodeOne.at(3) * answer.at(1) + projectionComponentNodeTwo.at(3) * answer.at(2);
}

void
Lattice3dboundary_mt :: computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rmode, TimeStep *tStep)
{
    double dV;
    FloatMatrix b, d, db, bi, bj, dbj, dij;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    double length = giveLength();
    double k;

    answer.resize(5, 5);
    answer.zero();

    FloatMatrix answerTemp(2, 2);
    answerTemp.zero();

    IntArray projectionComponentNodeOne(3);
    projectionComponentNodeOne.zero();
    if ( location.at(1) != 0 ) {
        giveSwitches( projectionComponentNodeOne, location.at(1) );
    }

    IntArray projectionComponentNodeTwo(3);
    projectionComponentNodeTwo.zero();
    if ( location.at(2) != 0 ) {
        giveSwitches( projectionComponentNodeTwo, location.at(2) );
    }

    answer.at(1, 1) = 1.;
    answer.at(1, 2) = -1.;
    answer.at(2, 1) = -1;
    answer.at(2, 2) = 1.;

    k = static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicValue(Conductivity, gp, tStep);
    dV      = this->computeVolumeAround(gp);
    double temp = k * dV / pow(length, 2.);

    answer.times(temp);

    answer.at(1, 3) = projectionComponentNodeOne.at(1) * answer.at(1, 1) + projectionComponentNodeTwo.at(1) * answer.at(1, 2);
    answer.at(1, 4) = projectionComponentNodeOne.at(2) * answer.at(1, 1) + projectionComponentNodeTwo.at(2) * answer.at(1, 2);
    answer.at(1, 5) = projectionComponentNodeOne.at(3) * answer.at(1, 1) + projectionComponentNodeTwo.at(3) * answer.at(1, 2);

    answer.at(2, 3) = projectionComponentNodeOne.at(1) * answer.at(2, 1) + projectionComponentNodeTwo.at(1) * answer.at(2, 2);
    answer.at(2, 4) = projectionComponentNodeOne.at(2) * answer.at(2, 1) + projectionComponentNodeTwo.at(2) * answer.at(2, 2);
    answer.at(2, 5) = projectionComponentNodeOne.at(3) * answer.at(2, 1) + projectionComponentNodeTwo.at(3) * answer.at(2, 2);

    answer.at(3, 1) = projectionComponentNodeOne.at(1) * answer.at(1, 1) + projectionComponentNodeTwo.at(1) * answer.at(2, 1);
    answer.at(3, 2) = projectionComponentNodeOne.at(1) * answer.at(1, 2) + projectionComponentNodeTwo.at(1) * answer.at(2, 2);
    answer.at(3, 3) = projectionComponentNodeOne.at(1) * projectionComponentNodeOne.at(1) * answer.at(1, 1) + projectionComponentNodeOne.at(1) * projectionComponentNodeTwo.at(1) * answer.at(1, 2) + projectionComponentNodeOne.at(1) * projectionComponentNodeTwo.at(1) * answer.at(2, 1) + projectionComponentNodeTwo.at(1) * projectionComponentNodeTwo.at(1) * answer.at(2, 2);
    answer.at(3, 4) = projectionComponentNodeOne.at(1) * projectionComponentNodeOne.at(2) * answer.at(1, 1) + projectionComponentNodeOne.at(1) * projectionComponentNodeTwo.at(2) * answer.at(1, 2) + projectionComponentNodeOne.at(2) * projectionComponentNodeTwo.at(1) * answer.at(2, 1) + projectionComponentNodeTwo.at(1) * projectionComponentNodeTwo.at(2) * answer.at(2, 2);
    answer.at(3, 5) = projectionComponentNodeOne.at(1) * projectionComponentNodeOne.at(3) * answer.at(1, 1) + projectionComponentNodeOne.at(1) * projectionComponentNodeTwo.at(3) * answer.at(1, 2) + projectionComponentNodeOne.at(3) * projectionComponentNodeTwo.at(1) * answer.at(2, 1) + projectionComponentNodeTwo.at(1) * projectionComponentNodeTwo.at(3) * answer.at(2, 2);

    answer.at(4, 1) = projectionComponentNodeOne.at(2) * answer.at(1, 1) + projectionComponentNodeTwo.at(2) * answer.at(2, 1);
    answer.at(4, 2) = projectionComponentNodeOne.at(2) * answer.at(1, 2) + projectionComponentNodeTwo.at(2) * answer.at(2, 2);
    answer.at(4, 3) = projectionComponentNodeOne.at(1) * projectionComponentNodeOne.at(2) * answer.at(1, 1) + projectionComponentNodeOne.at(2) * projectionComponentNodeTwo.at(1) * answer.at(1, 2) + projectionComponentNodeOne.at(1) * projectionComponentNodeTwo.at(2) * answer.at(2, 1) + projectionComponentNodeTwo.at(1) * projectionComponentNodeTwo.at(2) * answer.at(2, 2);
    answer.at(4, 4) = projectionComponentNodeOne.at(2) * projectionComponentNodeOne.at(2) * answer.at(1, 1) + projectionComponentNodeOne.at(2) * projectionComponentNodeTwo.at(2) * answer.at(1, 2) + projectionComponentNodeOne.at(2) * projectionComponentNodeTwo.at(2) * answer.at(2, 1) + projectionComponentNodeTwo.at(2) * projectionComponentNodeTwo.at(2) * answer.at(2, 2);
    answer.at(4, 5) =  projectionComponentNodeOne.at(2) * projectionComponentNodeOne.at(3) * answer.at(1, 1) + projectionComponentNodeOne.at(2) * projectionComponentNodeTwo.at(3) * answer.at(1, 2) + projectionComponentNodeOne.at(3) * projectionComponentNodeTwo.at(2) * answer.at(2, 1) + projectionComponentNodeTwo.at(2) * projectionComponentNodeTwo.at(3) * answer.at(2, 2);

    answer.at(5, 1) = projectionComponentNodeOne.at(3) * answer.at(1, 1) + projectionComponentNodeTwo.at(3) * answer.at(2, 1);
    answer.at(5, 2) = projectionComponentNodeOne.at(3) * answer.at(1, 2) +  projectionComponentNodeTwo.at(3) * answer.at(2, 2);
    answer.at(5, 3) = projectionComponentNodeOne.at(1) * projectionComponentNodeOne.at(3) * answer.at(1, 1) + projectionComponentNodeOne.at(3) * projectionComponentNodeTwo.at(1) * answer.at(1, 2) + projectionComponentNodeOne.at(1) * projectionComponentNodeTwo.at(3) * answer.at(2, 1) + projectionComponentNodeTwo.at(1) * projectionComponentNodeTwo.at(3) * answer.at(2, 2);
    answer.at(5, 4) = projectionComponentNodeOne.at(2) * projectionComponentNodeOne.at(3) * answer.at(1, 1) + projectionComponentNodeOne.at(3) * projectionComponentNodeTwo.at(2) * answer.at(1, 2) + projectionComponentNodeOne.at(2) * projectionComponentNodeTwo.at(3) * answer.at(2, 1) + projectionComponentNodeTwo.at(2) * projectionComponentNodeTwo.at(3) * answer.at(2, 2);
    answer.at(5, 5) =  projectionComponentNodeOne.at(3) * projectionComponentNodeOne.at(3) * answer.at(1, 1) + projectionComponentNodeOne.at(3) * projectionComponentNodeTwo.at(3) * answer.at(1, 2) + projectionComponentNodeOne.at(3) * projectionComponentNodeTwo.at(3) * answer.at(2, 1) + projectionComponentNodeTwo.at(3) * projectionComponentNodeTwo.at(3) * answer.at(2, 2);
    return;
}


double
Lattice3dboundary_mt :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    return this->giveArea() * this->giveLength() / 2.;
}

void
Lattice3dboundary_mt :: giveVTKCoordinates(int nodeNumber, FloatArray &coords) {
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
            giveSwitches( projectionComponent, location.at(1) );
        }
    } else if ( nodeNumber == 2 ) {
        node  = this->giveNode(2);
        if ( location.at(2) != 0 ) {
            giveSwitches( projectionComponent, location.at(2) );
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
Lattice3dboundary_mt :: computeGeometryProperties()
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
        giveSwitches( projectionComponentNodeOne, location.at(1) );
    }

    IntArray projectionComponentNodeTwo(3);
    projectionComponentNodeTwo.zero();
    if ( location.at(2) != 0 ) {
        giveSwitches( projectionComponentNodeTwo, location.at(2) );
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

    // Compute midpoint
    this->midPoint.resize(3);
    for ( int i = 0; i < 3; i++ ) {
        this->midPoint.at(i + 1) = 0.5 * ( coordsB.at(i + 1) + coordsA.at(i + 1) );
    }

    this->length  = sqrt( pow(this->normal.at(1), 2.) + pow(this->normal.at(2), 2.) + pow(this->normal.at(3), 2.) );

    //Treat zero length elements
    if ( this->length < this->minLength ) {
        this->length = this->minLength;
        computeSpecialCrossSectionProperties();
        this->geometryFlag = 1;
        return;
    }

    for ( int i = 0; i < 3; i++ ) {
        this->normal.at(i + 1) /= length;
    }

    computeCrossSectionProperties();
    //Set geometry flag to 1 so that this is done only once
    this->geometryFlag = 1;

    return;
}

void
Lattice3dboundary_mt :: giveSwitches(IntArray &answer, int location) {
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


#ifdef __OOFEG

void Lattice3dboundary_mt :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    WCRec p [ 2 ]; /* points */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);


    IntArray projectionComponentNodeOne(3);
    projectionComponentNodeOne.zero();
    if ( location.at(1) != 0 ) {
        giveSwitches( projectionComponentNodeOne, location.at(1) );
    }

    IntArray projectionComponentNodeTwo(3);
    projectionComponentNodeTwo.zero();
    if ( location.at(2) != 0 ) {
        giveSwitches( projectionComponentNodeTwo, location.at(2) );
    }

    FloatArray specimenDimension(3);
    specimenDimension.at(1) =  this->giveNode(3)->giveCoordinate(1);
    specimenDimension.at(2) =  this->giveNode(3)->giveCoordinate(2);
    specimenDimension.at(3) =  this->giveNode(3)->giveCoordinate(3);

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

#endif
} // end namespace oofem
