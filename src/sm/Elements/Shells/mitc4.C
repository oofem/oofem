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

#include "sm/Elements/Shells/mitc4.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "fei2dquadlin.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "sm/CrossSections/variablecrosssection.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "intarray.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "classfactory.h"
#include "connectivitytable.h"

namespace oofem {
REGISTER_Element(MITC4Shell);

FEI2dQuadLin MITC4Shell::interp_lin(1, 2);

MITC4Shell::MITC4Shell(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this),
    SPRNodalRecoveryModelInterface(), SpatialLocalizerInterface(this),
    nPointsXY(4),
    nPointsZ(2)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = nPointsXY * nPointsZ;
}


FEInterpolation *
MITC4Shell::giveInterpolation() const { return & interp_lin; }


FEInterpolation *
MITC4Shell::giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}


Interface *
MITC4Shell::giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >( this );
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >( this );
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >( this );
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >( this );
    }

    return nullptr;
}


void
MITC4Shell::SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(numberOfDofMans);
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}


void
MITC4Shell::SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        if ( this->giveNode(i)->giveNumber() == pap ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("unknown node number %d", pap);
    }
}


int
MITC4Shell::SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
}


SPRPatchType
MITC4Shell::SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiLin;
}


void
MITC4Shell::computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ] = std::make_unique< GaussIntegrationRule >(1, this, 1, 10);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], nPointsXY, nPointsZ, this);
    }
}


std::array< FloatArrayF< 3 >, 4 >
MITC4Shell::giveDirectorVectors()
{
    if ( directorType == 0 ) { // normal to the midplane
        auto e = this->computeLocalBaseVectors();
        return { e [ 2 ], e [ 2 ], e [ 2 ], e [ 2 ] };
    } else if ( directorType == 1 ) {          // nodal average
        int csNum = this->giveCrossSection()->giveNumber();
        std::array< FloatArrayF< 3 >, 4 >directors;

        IntArray neighbours;
        ConnectivityTable *conTable = this->giveDomain()->giveConnectivityTable();
        IntArray node(1);

        for ( int i = 1; i <= 4; i++ ) {
            FloatArrayF< 3 >nodeDir;
            node.at(1) = this->giveNode(i)->giveNumber();
            conTable->giveNodeNeighbourList(neighbours, node);

            for ( int j = 1; j <= neighbours.giveSize(); j++ ) {
                auto neighbour = dynamic_cast< MITC4Shell * >( this->giveDomain()->giveElement(neighbours.at(j) ) );
                if ( neighbour ) {
                    if ( neighbour->giveCrossSection()->giveNumber() == csNum ) {
                        auto e = neighbour->computeLocalBaseVectors();
                        nodeDir += e [ 2 ];
                    }
                }
            }
            directors [ i ] = normalize(nodeDir);
        }
        return directors;
    } else if ( directorType == 2 ) {          // specified at crosssection
        std::array< FloatArrayF< 3 >, 4 >V;
        for ( int i = 0; i < 4; ++i ) {
            const auto &c = this->giveNode(i + 1)->giveCoordinates();
            FloatArrayF< 3 >v = {
                this->giveCrossSection()->give(CS_DirectorVectorX, c, this, false),
                this->giveCrossSection()->give(CS_DirectorVectorY, c, this, false),
                this->giveCrossSection()->give(CS_DirectorVectorZ, c, this, false),
            };
            V [ i ] = normalize(v);
        }
        return V;
    } else {
        throw std::runtime_error("Unsupported directorType");
    }
}


std::array< FloatArrayF< 3 >, 4 >
MITC4Shell::giveLocalDirectorVectors()
{
    auto Vg = this->giveDirectorVectors();
    return {
        dot(GtoLRotationMatrix, Vg [ 0 ]),
        dot(GtoLRotationMatrix, Vg [ 1 ]),
        dot(GtoLRotationMatrix, Vg [ 2 ]),
        dot(GtoLRotationMatrix, Vg [ 3 ]),
    };
}


void
MITC4Shell::computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the [6x24] displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
// Zeroes in rows 4, 5, 6.
{
    auto h = interp_lin.evalN(FloatArrayF< 3 >(iLocCoord) [ { 0, 1 } ]);
    auto a = this->giveThickness();
    auto V = this->giveLocalDirectorVectors();

    FloatArrayF< 3 >e2 = { 0., 1., 0. };

    answer.resize(6, 6 * 4);
    answer.zero();
    for ( int i = 0; i < 4; ++i ) {
        auto Ve = normalize( cross(e2, V [ i ]) );
        auto VVe = cross(V [ i ], Ve);
        answer.at(1, 1 + i * 6) = answer.at(2, 2 + i * 6) = answer.at(3, 3 + i * 6) = h [ i ];
        answer.at(1, 4 + i * 6) = -iLocCoord.at(3) / 2.0 * a [ i ] * h [ i ] * VVe.at(1);
        answer.at(1, 5 + i * 6) =  iLocCoord.at(3) / 2.0 * a [ i ] * h [ i ] * Ve.at(1);
        answer.at(2, 4 + i * 6) = -iLocCoord.at(3) / 2.0 * a [ i ] * h [ i ] * VVe.at(2);
        answer.at(2, 5 + i * 6) =  iLocCoord.at(3) / 2.0 * a [ i ] * h [ i ] * Ve.at(2);
        answer.at(3, 4 + i * 6) = -iLocCoord.at(3) / 2.0 * a [ i ] * h [ i ] * VVe.at(3);
        answer.at(3, 5 + i * 6) =  iLocCoord.at(3) / 2.0 * a [ i ] * h [ i ] * Ve.at(3);
    }
}


void
MITC4Shell::computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    auto scs = dynamic_cast< StructuralCrossSection * >( this->giveCrossSection() );
    //auto cs = dynamic_cast< SimpleCrossSection * >( this->giveCrossSection() );
    answer = scs->give3dDegeneratedShellStiffMtrx(rMode, gp, tStep);
}


std::array< FloatArrayF< 3 >, 4 >
MITC4Shell::giveNodeCoordinates()
{
    std::array< FloatArrayF< 3 >, 4 >c;
    for ( int i = 0; i < 4; ++i ) {
        c [ i ] = this->giveLocalCoordinates( this->giveNode(i + 1)->giveCoordinates() );
    }
    return c;
}

FloatArrayF< 3 >
MITC4Shell::giveLocalCoordinates(const FloatArrayF< 3 > &global)
{
    auto offset = global - FloatArrayF< 3 >( this->giveNode(1)->giveCoordinates() );
    return dot(GtoLRotationMatrix, offset);
}

void
MITC4Shell::initializeFrom(InputRecord &ir)
{
    StructuralElement::initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, nPointsXY, _IFT_Element_nip);
    IR_GIVE_OPTIONAL_FIELD(ir, nPointsZ, _IFT_MITC4Shell_nipZ);
    //@todo: extend for nonlinear geometry
    //IR_GIVE_OPTIONAL_FIELD(ir, nlGeometry, _IFT_NLStructuralElement_nlgeoflag);

    directorType = 0; // default
    IR_GIVE_OPTIONAL_FIELD(ir, directorType, _IFT_MITC4Shell_directorType);
}


void
MITC4Shell::giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

double
MITC4Shell::computeVolumeAround(GaussPoint *gp)
{
    FloatArrayF< 3 >lcoords = {
        gp->giveNaturalCoordinate(1),
        gp->giveNaturalCoordinate(2),
        gp->giveNaturalCoordinate(3),
    };

    auto jacobianMatrix = this->giveJacobian(lcoords);
    return det(jacobianMatrix) * gp->giveWeight();
}


FloatMatrixF< 3, 3 >
MITC4Shell::giveJacobian(const FloatArrayF< 3 > &lcoords)
{
    // derivatives of interpolation functions
    auto dn = interp_lin.evaldNdxi(lcoords [ { 0, 1 } ]);
    auto hk1 = dn.row< 0 >(); // dh(r1,r2)/dr1
    auto hk2 = dn.row< 1 >(); // dh(r1,r2)/dr2

    // interpolation functions - h(r1,r2)
    //auto [h1, h2, h3, h4] = interp_lin.evalN(lcoords[{0, 1}]);
    auto h = interp_lin.evalN(lcoords [ { 0, 1 } ]);
    double h1 = h [ 0 ];
    double h2 = h [ 1 ];
    double h3 = h [ 2 ];
    double h4 = h [ 3 ];

    //auto [a1, a2, a3, a4] = this->giveThickness();
    auto a = this->giveThickness();
    double a1 = a [ 0 ];
    double a2 = a [ 1 ];
    double a3 = a [ 2 ];
    double a4 = a [ 3 ];

    // get node coordinates
    double r3 = lcoords.at(3);
    auto c = this->giveNodeCoordinates();
    //auto [x1, y1, z1] = c[0];
    //auto [x2, y2, z2] = c[1];
    //auto [x3, y3, z3] = c[2];
    //auto [x4, y4, z4] = c[3];
    double x1 = c [ 0 ] [ 0 ], y1 = c [ 0 ] [ 1 ];
    double x2 = c [ 1 ] [ 0 ], y2 = c [ 1 ] [ 1 ];
    double x3 = c [ 2 ] [ 0 ], y3 = c [ 2 ] [ 1 ];
    double x4 = c [ 3 ] [ 0 ], y4 = c [ 3 ] [ 1 ];

    //auto [V1, V2, V3, V4] = this->giveLocalDirectorVectors();
    auto V = this->giveLocalDirectorVectors();
    auto V1 = V [ 0 ], V2 = V [ 1 ], V3 = V [ 2 ], V4 = V [ 3 ];

    // Jacobian Matrix
    FloatMatrixF< 3, 3 >jacobianMatrix;
    jacobianMatrix.at(1, 1) = hk1.at(1) * x1 + hk1.at(2) * x2 + hk1.at(3) * x3 + hk1.at(4) * x4 + r3 / 2. * ( a1 * hk1.at(1) * V1.at(1) + a2 * hk1.at(2) * V2.at(1) + a3 * hk1.at(3) * V3.at(1) + a4 * hk1.at(4) * V4.at(1) );
    jacobianMatrix.at(1, 2) = hk1.at(1) * y1 + hk1.at(2) * y2 + hk1.at(3) * y3 + hk1.at(4) * y4 + r3 / 2. * ( a1 * hk1.at(1) * V1.at(2) + a2 * hk1.at(2) * V2.at(2) + a3 * hk1.at(3) * V3.at(2) + a4 * hk1.at(4) * V4.at(2) );
    jacobianMatrix.at(1, 3) = r3 / 2. * ( a1 * hk1.at(1) * V1.at(3) + a2 * hk1.at(2) * V2.at(3) + a3 * hk1.at(3) * V3.at(3) + a4 * hk1.at(4) * V4.at(3) );
    jacobianMatrix.at(2, 1) = hk2.at(1) * x1 + hk2.at(2) * x2 + hk2.at(3) * x3 + hk2.at(4) * x4 + r3 / 2. * ( a1 * hk2.at(1) * V1.at(1) + a2 * hk2.at(2) * V2.at(1) + a3 * hk2.at(3) * V3.at(1) + a4 * hk2.at(4) * V4.at(1) );
    jacobianMatrix.at(2, 2) = hk2.at(1) * y1 + hk2.at(2) * y2 + hk2.at(3) * y3 + hk2.at(4) * y4 + r3 / 2. * ( a1 * hk2.at(1) * V1.at(2) + a2 * hk2.at(2) * V2.at(2) + a3 * hk2.at(3) * V3.at(2) + a4 * hk2.at(4) * V4.at(2) );
    jacobianMatrix.at(2, 3) = r3 / 2. * ( a1 * hk2.at(1) * V1.at(3) + a2 * hk2.at(2) * V2.at(3) + a3 * hk2.at(3) * V3.at(3) + a4 * hk2.at(4) * V4.at(3) );
    jacobianMatrix.at(3, 1) =  1. / 2. * ( a1 * h1 * V1.at(1) + a2 * h2 * V2.at(1) + a3 * h3 * V3.at(1) + a4 * h4 * V4.at(1) );
    jacobianMatrix.at(3, 2) =  1. / 2. * ( a1 * h1 * V1.at(2) + a2 * h2 * V2.at(2) + a3 * h3 * V3.at(2) + a4 * h4 * V4.at(2) );
    jacobianMatrix.at(3, 3) =  1. / 2. * ( a1 * h1 * V1.at(3) + a2 * h2 * V2.at(3) + a3 * h3 * V3.at(3) + a4 * h4 * V4.at(3) );
    return jacobianMatrix;
}

void
MITC4Shell::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // This element adds an additional stiffness for the so called drilling dofs.
    StructuralElement::computeStiffnessMatrix(answer, rMode, tStep);

    bool drillType = this->giveStructuralCrossSection()->give(CS_DrillingType, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0) );
    if ( drillType == 1 ) {
        double relDrillCoeff = this->giveStructuralCrossSection()->give(CS_RelDrillingStiffness, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0) );
        if ( relDrillCoeff == 0.0 ) {
            relDrillCoeff = 0.001; // default
        }

        int j = 1;
        while ( answer.at(j, j) == 0 ) {
            j++;
        }
        drillCoeff = answer.at(j, j);
        // find the smallest non-zero number on the diagonal
        for ( int i = j; i <= 24; i++ ) {
            if ( drillCoeff > answer.at(i, i) && answer.at(i, i) != 0 ) {
                drillCoeff = answer.at(i, i);
            }
        }

        drillCoeff *= relDrillCoeff;

        IntArray drillDofs = { 6, 12, 18, 24 };
        auto drillStiffness = eye< 4 >() * drillCoeff;
#if 0
        FloatMatrix drillStiffness;
        for ( auto &gp : * integrationRulesArray [ 0 ] ) {
            double dV = this->computeVolumeAround(gp);
            // double drillCoeff = this->giveStructuralCrossSection()->give(CS_DrillingStiffness, gp);
            double coeff = drillCoeff;
            if ( this->giveStructuralCrossSection()->give(CS_DrillingStiffness, gp) > 0 ) {
                drillCoeff *= this->giveStructuralCrossSection()->give(CS_DrillingStiffness, gp);
            }
            // Drilling stiffness is here for improved numerical properties
            this->interp_lin.evalN(n, gp->giveNaturalCoordinates(), FEIVoidCellGeometry() );
            for ( int j = 0; j < 4; j++ ) {
                n [ j ] -= 0.25;
            }
            drillStiffness.plusDyadSymmUpper(n, coeff * dV);
        }
        drillStiffness.symmetrized();
#endif
        answer.assemble(drillStiffness, drillDofs);
    }
}

void
MITC4Shell::giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // This element adds an additional stiffness for the so called drilling dofs.
    StructuralElement::giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);

    bool drillType = this->giveStructuralCrossSection()->give(CS_DrillingType, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0) );

    if ( drillType == 1 ) {
        FloatArray n, tmp;
        FloatArray drillUnknowns, drillMoment;
        IntArray drillDofs = {
            6, 12, 18, 24
        };
        this->computeVectorOf(VM_Total, tStep, tmp);
        drillUnknowns.beSubArrayOf(tmp, drillDofs);

#if 0
        for ( auto &gp : * integrationRulesArray [ 0 ] ) {
            double dV = this->computeVolumeAround(gp);
            // double drillCoeff = this->giveStructuralCrossSection()->give(CS_DrillingStiffness, gp);
            double coeff = drillCoeff;
            if ( this->giveStructuralCrossSection()->give(CS_DrillingStiffness, gp) > 0 ) {
                drillCoeff *= this->giveStructuralCrossSection()->give(CS_DrillingStiffness, gp);
            }

            this->interp_lin.evalN(n, gp->giveNaturalCoordinates(), FEIVoidCellGeometry() );
            for ( int j = 0; j < 4; j++ ) {
                n [ j ] -= 0.25;
            }
            double dtheta = n.dotProduct(drillUnknowns);
            drillMoment.add(coeff * dV * dtheta, n);
        }
#endif

        FloatMatrix drillStiffness;
        drillStiffness.resize(4, 4);
        drillStiffness.zero();
        for ( int i = 1; i <= 4; i++ ) {
            drillStiffness.at(i, i) = drillCoeff;
        }

        drillMoment.beProductOf(drillStiffness, drillUnknowns);

        answer.assemble(drillMoment, drillDofs);
    }
}

void
MITC4Shell::computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x20] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    //auto [c1, c2, c3, c4] = this->giveNodeCoordinates();
    auto c = this->giveNodeCoordinates();
    auto c1 = c [ 0 ], c2 = c [ 1 ], c3 = c [ 2 ], c4 = c [ 3 ];

    double r1 = gp->giveNaturalCoordinate(1);
    double r2 = gp->giveNaturalCoordinate(2);
    double r3 = gp->giveNaturalCoordinate(3);
    FloatArrayF< 3 >lcoords = { r1, r2, r3 };

    //auto [a1, a2, a3, a4] = this->giveThickness();
    auto a = this->giveThickness();
    double a1 = a [ 0 ], a2 = a [ 1 ], a3 = a [ 2 ], a4 = a [ 3 ];

    //auto h = interp_lin.evalN(lcoords);
    //auto dn = interp_lin.evaldNdxi(lcoords);

    //auto [hkx, hky] = this->givedNdx(lcoords);
    auto hk = this->givedNdx(lcoords);
    auto hkx = hk [ 0 ], hky = hk [ 1 ];

    auto J = this->giveJacobian(lcoords);
    auto invJ = inv(J);
    double sb = 2 * invJ.at(1, 1) * invJ.at(3, 3);
    double sa = 2 * invJ.at(1, 2) * invJ.at(3, 3);
    double cb = 2 * invJ.at(2, 1) * invJ.at(3, 3);
    double ca = 2 * invJ.at(2, 2) * invJ.at(3, 3);

    //auto [V1, V2, V3, V4] = this->giveLocalDirectorVectors();
    auto V = this->giveLocalDirectorVectors();
    auto V1 = V [ 0 ], V2 = V [ 1 ], V3 = V [ 2 ], V4 = V [ 3 ];

    FloatArrayF< 3 >e2 = { 0., 1., 0. };

    auto V11 = normalize( cross(e2, V1) );
    auto V12 = normalize( cross(e2, V2) );
    auto V13 = normalize( cross(e2, V3) );
    auto V14 = normalize( cross(e2, V4) );

    auto V21 = cross(V1, V11);
    auto V22 = cross(V2, V12);
    auto V23 = cross(V3, V13);
    auto V24 = cross(V4, V14);

    answer.resize(6, 24);
    answer.zero();

    answer.at(4, 1) = 1. / 32. * ( ( a1 * V1.at(1) + a2 * V2.at(1) ) * ( cb * ( 1. + r2 ) ) + ( a1 * V1.at(1) + a4 * V4.at(1) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 2) = 1. / 32. * ( ( a1 * V1.at(2) + a2 * V2.at(2) ) * ( cb * ( 1. + r2 ) ) + ( a1 * V1.at(2) + a4 * V4.at(2) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 3) = 1. / 32. * ( ( a1 * V1.at(3) + a2 * V2.at(3) ) * ( cb * ( 1. + r2 ) ) + ( a1 * V1.at(3) + a4 * V4.at(3) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 4) = -a1 / 32. * ( ( dot(V21, c1 - c2) * ( cb * ( 1. + r2 ) ) ) + ( dot(V21, c1 - c4) * ( ca * ( 1. + r1 ) ) ) );
    answer.at(4, 5) = a1 / 32. * ( ( dot(V11, c1 - c2) * ( cb * ( 1. + r2 ) ) ) + ( dot(V11, c1 - c4) * ( ca * ( 1. + r1 ) ) ) );

    answer.at(4, 7) = 1. / 32. * ( -( a1 * V1.at(1) + a2 * V2.at(1) ) * ( cb * ( 1. + r2 ) ) + ( a2 * V2.at(1) + a3 * V3.at(1) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 8) = 1. / 32. * ( -( a1 * V1.at(2) + a2 * V2.at(2) ) * ( cb * ( 1. + r2 ) ) + ( a2 * V2.at(2) + a3 * V3.at(2) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 9) = 1. / 32. * ( -( a1 * V1.at(3) + a2 * V2.at(3) ) * ( cb * ( 1. + r2 ) ) + ( a2 * V2.at(3) + a3 * V3.at(3) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 10) = -a1 / 32. * ( ( dot(V21, c1 - c2) * ( cb * ( 1. + r2 ) ) ) + ( dot(V21, c2 - c3) * ( ca * ( 1. - r1 ) ) ) );
    answer.at(4, 11) = a1 / 32. * ( ( dot(V11, c1 - c2) * ( cb * ( 1. + r2 ) ) ) + ( dot(V11, c2 - c3) * ( ca * ( 1. - r1 ) ) ) );

    answer.at(4, 13) = 1. / 32. * ( -( a3 * V3.at(1) + a4 * V4.at(1) ) * ( cb * ( 1. - r2 ) ) - ( a2 * V2.at(1) + a3 * V3.at(1) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 14) = 1. / 32. * ( -( a3 * V3.at(2) + a4 * V4.at(2) ) * ( cb * ( 1. - r2 ) ) - ( a2 * V2.at(2) + a3 * V3.at(2) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 15) = 1. / 32. * ( -( a3 * V3.at(3) + a4 * V4.at(3) ) * ( cb * ( 1. - r2 ) ) - ( a2 * V2.at(3) + a3 * V3.at(3) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 16) = -a1 / 32. * ( ( dot(V21, c4 - c3) * ( cb * ( 1. - r2 ) ) ) + ( dot(V21, c2 - c3) * ( ca * ( 1. - r1 ) ) ) );
    answer.at(4, 17) = a1 / 32. * ( ( dot(V11, c4 - c3) * ( cb * ( 1. - r2 ) ) ) + ( dot(V11, c2 - c3) * ( ca * ( 1. - r1 ) ) ) );

    answer.at(4, 19) = 1. / 32. * ( ( a3 * V3.at(1) + a4 * V4.at(1) ) * ( cb * ( 1. - r2 ) ) - ( a1 * V1.at(1) + a4 * V4.at(1) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 20) = 1. / 32. * ( ( a3 * V3.at(2) + a4 * V4.at(2) ) * ( cb * ( 1. - r2 ) ) - ( a1 * V1.at(2) + a4 * V4.at(2) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 21) = 1. / 32. * ( ( a3 * V3.at(3) + a4 * V4.at(3) ) * ( cb * ( 1. - r2 ) ) - ( a1 * V1.at(3) + a4 * V4.at(3) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 22) = -a1 / 32. * ( ( dot(V21, c4 - c3) * ( cb * ( 1. - r2 ) ) ) + ( dot(V21, c1 - c4) * ( ca * ( 1. + r1 ) ) ) );
    answer.at(4, 23) = a1 / 32. * ( ( dot(V11, c4 - c3) * ( cb * ( 1. - r2 ) ) ) + ( dot(V11, c1 - c4) * ( ca * ( 1. + r1 ) ) ) );

    answer.at(5, 1) = 1. / 32. * ( ( a1 * V1.at(1) + a2 * V2.at(1) ) * ( sb * ( 1. + r2 ) ) + ( a1 * V1.at(1) + a4 * V4.at(1) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 2) = 1. / 32. * ( ( a1 * V1.at(2) + a2 * V2.at(2) ) * ( sb * ( 1. + r2 ) ) + ( a1 * V1.at(2) + a4 * V4.at(2) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 3) = 1. / 32. * ( ( a1 * V1.at(3) + a2 * V2.at(3) ) * ( sb * ( 1. + r2 ) ) + ( a1 * V1.at(3) + a4 * V4.at(3) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 4) = -a1 / 32. * ( ( dot(V21, c1 - c2) * ( sb * ( 1. + r2 ) ) ) + ( dot(V21, c1 - c4) * ( sa * ( 1. + r1 ) ) ) );
    answer.at(5, 5) = a1 / 32. * ( ( dot(V11, c1 - c2) * ( sb * ( 1. + r2 ) ) ) + ( dot(V11, c1 - c4) * ( sa * ( 1. + r1 ) ) ) );

    answer.at(5, 7) = 1. / 32. * ( -( a1 * V1.at(1) + a2 * V2.at(1) ) * ( sb * ( 1. + r2 ) ) + ( a2 * V2.at(1) + a3 * V3.at(1) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 8) = 1. / 32. * ( -( a1 * V1.at(2) + a2 * V2.at(2) ) * ( sb * ( 1. + r2 ) ) + ( a2 * V2.at(2) + a3 * V3.at(2) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 9) = 1. / 32. * ( -( a1 * V1.at(3) + a2 * V2.at(3) ) * ( sb * ( 1. + r2 ) ) + ( a2 * V2.at(3) + a3 * V3.at(3) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 10) = -a1 / 32. * ( ( dot(V21, c1 - c2) * ( sb * ( 1. + r2 ) ) ) + ( dot(V21, c2 - c3) * ( sa * ( 1. - r1 ) ) ) );
    answer.at(5, 11) = a1 / 32. * ( ( dot(V11, c1 - c2) * ( sb * ( 1. + r2 ) ) ) + ( dot(V11, c2 - c3) * ( sa * ( 1. - r1 ) ) ) );

    answer.at(5, 13) = 1. / 32. * ( -( a3 * V3.at(1) + a4 * V4.at(1) ) * ( sb * ( 1. - r2 ) ) - ( a2 * V2.at(1) + a3 * V3.at(1) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 14) = 1. / 32. * ( -( a3 * V3.at(2) + a4 * V4.at(2) ) * ( sb * ( 1. - r2 ) ) - ( a2 * V2.at(2) + a3 * V3.at(2) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 15) = 1. / 32. * ( -( a3 * V3.at(3) + a4 * V4.at(3) ) * ( sb * ( 1. - r2 ) ) - ( a2 * V2.at(3) + a3 * V3.at(3) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 16) = -a1 / 32. * ( ( dot(V21, c4 - c3) * ( sb * ( 1. - r2 ) ) ) + ( dot(V21, c2 - c3) * ( sa * ( 1. - r1 ) ) ) );
    answer.at(5, 17) = a1 / 32. * ( ( dot(V11, c4 - c3) * ( sb * ( 1. - r2 ) ) ) + ( dot(V11, c2 - c3) * ( sa * ( 1. - r1 ) ) ) );

    answer.at(5, 19) = 1. / 32. * ( ( a3 * V3.at(1) + a4 * V4.at(1) ) * ( sb * ( 1. - r2 ) ) - ( a1 * V1.at(1) + a4 * V4.at(1) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 20) = 1. / 32. * ( ( a3 * V3.at(2) + a4 * V4.at(2) ) * ( sb * ( 1. - r2 ) ) - ( a1 * V1.at(2) + a4 * V4.at(2) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 21) = 1. / 32. * ( ( a3 * V3.at(3) + a4 * V4.at(3) ) * ( sb * ( 1. - r2 ) ) - ( a1 * V1.at(3) + a4 * V4.at(3) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 22) = -a1 / 32. * ( ( dot(V21, c4 - c3) * ( sb * ( 1. - r2 ) ) ) + ( dot(V21, c1 - c4) * ( sa * ( 1. + r1 ) ) ) );
    answer.at(5, 23) = a1 / 32. * ( ( dot(V11, c4 - c3) * ( sb * ( 1. - r2 ) ) ) + ( dot(V11, c1 - c4) * ( sa * ( 1. + r1 ) ) ) );

    answer.at(1, 1) = hkx.at(1);
    answer.at(1, 4) = -r3 / 2. * a1 * hkx.at(1) * V21.at(1);
    answer.at(1, 5) = r3 / 2. * a1 * hkx.at(1) * V11.at(1);

    answer.at(1, 7) = hkx.at(2);
    answer.at(1, 10) = -r3 / 2. * a2 * hkx.at(2) * V22.at(1);
    answer.at(1, 11) = r3 / 2. * a2 * hkx.at(2) * V12.at(1);

    answer.at(1, 13) = hkx.at(3);
    answer.at(1, 16) = -r3 / 2. * a3 * hkx.at(3) * V23.at(1);
    answer.at(1, 17) = r3 / 2. * a3 * hkx.at(3) * V13.at(1);

    answer.at(1, 19) = hkx.at(4);
    answer.at(1, 22) = -r3 / 2. * a4 * hkx.at(4) * V24.at(1);
    answer.at(1, 23) = r3 / 2. * a4 * hkx.at(4) * V14.at(1);

    answer.at(2, 2) = hky.at(1);
    answer.at(2, 4) = -r3 / 2. * a1 * hky.at(1) * V21.at(2);
    answer.at(2, 5) = r3 / 2. * a1 * hky.at(1) * V11.at(2);

    answer.at(2, 8) = hky.at(2);
    answer.at(2, 10) = -r3 / 2. * a2 * hky.at(2) * V22.at(2);
    answer.at(2, 11) = r3 / 2. * a2 * hky.at(2) * V12.at(2);

    answer.at(2, 14) = hky.at(3);
    answer.at(2, 16) = -r3 / 2. * a3 * hky.at(3) * V23.at(2);
    answer.at(2, 17) = r3 / 2. * a3 * hky.at(3) * V13.at(2);

    answer.at(2, 20) = hky.at(4);
    answer.at(2, 22) = -r3 / 2. * a4 * hky.at(4) * V24.at(2);
    answer.at(2, 23) = r3 / 2. * a4 * hky.at(4) * V14.at(2);

    answer.at(6, 1) = hky.at(1);
    answer.at(6, 2) = hkx.at(1);
    answer.at(6, 4) = -r3 / 2. * a1 * ( hkx.at(1) * V21.at(2) + hky.at(1) * V21.at(1) );
    answer.at(6, 5) = r3 / 2. * a1 * ( hky.at(1) * V11.at(1) + hky.at(1) * V11.at(2) );

    answer.at(6, 7) = hky.at(2);
    answer.at(6, 8) = hkx.at(2);
    answer.at(6, 10) = -r3 / 2. * a2 * ( hkx.at(2) * V22.at(2) + hky.at(2) * V22.at(1) );
    answer.at(6, 11) = r3 / 2. * a2 * ( hky.at(2) * V12.at(1) + hky.at(2) * V12.at(2) );

    answer.at(6, 13) = hky.at(3);
    answer.at(6, 14) = hkx.at(3);
    answer.at(6, 16) = -r3 / 2. * a3 * ( hkx.at(3) * V23.at(2) + hky.at(3) * V23.at(1) );
    answer.at(6, 17) = r3 / 2. * a3 * ( hky.at(3) * V13.at(1) + hky.at(3) * V13.at(2) );

    answer.at(6, 19) = hky.at(4);
    answer.at(6, 20) = hkx.at(4);
    answer.at(6, 22) = -r3 / 2. * a4 * ( hkx.at(4) * V24.at(2) + hky.at(4) * V24.at(1) );
    answer.at(6, 23) = r3 / 2. * a4 * ( hky.at(4) * V14.at(1) + hky.at(4) * V14.at(2) );
}


std::array< double, 4 >
MITC4Shell::giveThickness()
{
    const auto &c1 = this->giveNode(1)->giveCoordinates();
    const auto &c2 = this->giveNode(2)->giveCoordinates();
    const auto &c3 = this->giveNode(3)->giveCoordinates();
    const auto &c4 = this->giveNode(4)->giveCoordinates();

    return {
        this->giveCrossSection()->give(CS_Thickness, c1, this, false),
        this->giveCrossSection()->give(CS_Thickness, c2, this, false),
        this->giveCrossSection()->give(CS_Thickness, c3, this, false),
        this->giveCrossSection()->give(CS_Thickness, c4, this, false),
    };
}


void
MITC4Shell::postInitialize()
{
    StructuralElement::postInitialize();

    auto e = this->computeLocalBaseVectors();

    for ( int i = 1; i <= 3; i++ ) {
        GtoLRotationMatrix.at(1, i) = e [ 0 ].at(i);
        GtoLRotationMatrix.at(2, i) = e [ 1 ].at(i);
        GtoLRotationMatrix.at(3, i) = e [ 2 ].at(i);
    }
}

std::array< FloatArrayF< 3 >, 3 >
MITC4Shell::computeLocalBaseVectors()
{
    // compute A - (node2+node3)/2
    auto coordA = 0.5 * ( FloatArrayF< 3 >( this->giveNode(2)->giveCoordinates() ) + FloatArrayF< 3 >( this->giveNode(3)->giveCoordinates() ) );
    // compute B - (node1+node4)/2
    auto coordB = 0.5 * ( FloatArrayF< 3 >( this->giveNode(1)->giveCoordinates() ) + FloatArrayF< 3 >( this->giveNode(4)->giveCoordinates() ) );
    // compute e1' = [B-A]
    auto e1 = normalize(coordB - coordA);

    // compute A - (node3+node4)/2
    auto coordC = 0.5 * ( FloatArrayF< 3 >( this->giveNode(4)->giveCoordinates() ) + FloatArrayF< 3 >( this->giveNode(3)->giveCoordinates() ) );
    // compute B - (node2+node1)/2
    auto coordD = 0.5 * ( FloatArrayF< 3 >( this->giveNode(1)->giveCoordinates() ) + FloatArrayF< 3 >( this->giveNode(2)->giveCoordinates() ) );

    // compute help = [D-C]
    auto help = coordD - coordC;
    // compute e3' : vector product of e1' x help
    auto e3 = normalize( cross(e1, help) );
    // now from e3' x e1' compute e2'
    auto e2 = cross(e3, e1);
    return { e1, e2, e3 };
}

std::array< FloatMatrixF< 3, 3 >, 4 >
MITC4Shell::computeLToDirectorRotationMatrix()
// Returns the rotation matrix of the reciever of the size [3,3]
// {alpha_i,beta_i} = Ti * {rotL_xi, rotL_yi, rotL_zi}
//      alpha_i, beta_i - rotations about the components of director vector at node i
//      r1_i, r2_i, r3_i, - rotations about local coordinates e1', e2', e3'
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N4-N3]    Ni - means i - th node
// help   : [N2-N3]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
    //auto [e1, e2, e3, e4] = this->computeLocalBaseVectors();
    auto e = this->computeLocalBaseVectors();
    auto V = this->giveDirectorVectors();

    std::array< FloatMatrixF< 3, 3 >, 4 >answer;
    for ( int i = 0; i < 4; ++i ) {
        auto Ve = normalize( cross(e [ 1 ], V [ i ]) );
        auto VV = cross(V [ i ], Ve);

        answer [ i ].at(1, 1) = dot(Ve, e [ 0 ]);
        answer [ i ].at(1, 2) = dot(Ve, e [ 1 ]);
        answer [ i ].at(1, 3) = dot(Ve, e [ 2 ]);

        answer [ i ].at(2, 1) = dot(VV, e [ 0 ]);
        answer [ i ].at(2, 2) = dot(VV, e [ 1 ]);
        answer [ i ].at(2, 3) = dot(VV, e [ 2 ]);

        answer [ i ].at(3, 1) = dot(V [ i ], e [ 0 ]);
        answer [ i ].at(3, 2) = dot(V [ i ], e [ 1 ]);
        answer [ i ].at(3, 3) = dot(V [ i ], e [ 2 ]);
    }
    return answer;
}


bool
MITC4Shell::computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [24,24]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v,w,alpha,beta} = T * {u,v,w,r1,r2,r3}
{
    auto LtoDir = this->computeLToDirectorRotationMatrix();

    answer.resize(24, 24);
    answer.zero();

    for ( int i = 0; i <= 3; i++ ) {
        answer.setSubMatrix(GtoLRotationMatrix, i * 6 + 1, i * 6 + 1);
        auto help = dot(LtoDir [ i ], GtoLRotationMatrix);
        answer.setSubMatrix(help, i * 6 + 4, i * 6 + 4);
    }

    return 1;
}


void
MITC4Shell::computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveRealStress_3dDegeneratedShell(strain, gp, tStep);
}


FloatMatrix
MITC4Shell::giveCharacteristicTensor(CharTensor type, GaussPoint *gp, TimeStep *tStep)
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveStructuralCrossSection()->giveMaterial(gp) );

    if ( type == GlobalForceTensor ) {
        FloatArray localStress, localStrain;
        this->computeStrainVector(localStrain, gp, tStep);
        this->computeStressVector(localStress, localStrain, gp, tStep);
        auto stress = mat->transformStressVectorTo(GtoLRotationMatrix, localStress, false);
        return from_voigt_stress(stress);
    } else if ( type == GlobalStrainTensor ) {
        FloatArray localStrain;
        this->computeStrainVector(localStrain, gp, tStep);
        auto strain = mat->transformStrainVectorTo(GtoLRotationMatrix, localStrain, false);
        return from_voigt_strain(strain);
    } else {
        throw std::runtime_error("unsupported tensor mode");
    }
}

void
MITC4Shell::printOutputAt(FILE *file, TimeStep *tStep)
{
    FloatArray v;
    fprintf(file, "element %d (%8d):\n", this->giveLabel(), number);

    for ( int i = 0; i < nPointsXY; i++ ) {
        fprintf(file, "  GP %d :", i + 1);

        fprintf(file, "  forces     ");
        for ( auto &val : this->giveMidplaneIPValue(i, IST_ShellForceTensor, tStep) ) {
            fprintf(file, " %.4e", val);
        }

        fprintf(file, "\n          moments    ");
        for ( auto &val : this->giveMidplaneIPValue(i, IST_ShellMomentTensor, tStep) ) {
            fprintf(file, " %.4e", val);
        }

        fprintf(file, "\n          strains    ");
        for ( auto &val : this->giveMidplaneIPValue(i, IST_ShellStrainTensor, tStep) ) {
            fprintf(file, " %.4e", val);
        }

        fprintf(file, "\n          curvatures ");
        for ( auto &val : this->giveMidplaneIPValue(i, IST_CurvatureTensor, tStep) ) {
            fprintf(file, " %.4e", val);
        }

        for ( int j = 0; j < nPointsZ; j++ ) {
            auto gp = integrationRulesArray [ 0 ]->getIntegrationPoint(nPointsZ * i + j);

            fprintf(file, "\n          GP %d.%d :", i + 1, j + 1);

            this->giveIPValue(v, gp, IST_StrainTensor, tStep);
            fprintf(file, "    strains    ");
            for ( auto &val : v ) {
                fprintf(file, " %.4e", val);
            }

            this->giveIPValue(v, gp, IST_StressTensor, tStep);
            fprintf(file, "\n                      stresses   ");
            for ( auto &val : v ) {
                fprintf(file, " %.4e", val);
            }
        }
        fprintf(file, "\n");
    }
}

FloatArray
MITC4Shell::giveMidplaneIPValue(int gpXY, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_ShellMomentTensor || type == IST_ShellForceTensor ) {
        FloatArrayF< 6 >mLocal;

        for ( int i = 0; i < nPointsZ; i++ ) {
            auto gp = integrationRulesArray [ 0 ]->getIntegrationPoint(nPointsZ * gpXY + i);
            double thickness = this->giveCrossSection()->give(CS_Thickness, gp->giveGlobalCoordinates(), this, false);
            double J = thickness / 2.0;
            double z;
            if (  type == IST_ShellMomentTensor ) {
                z = gp->giveNaturalCoordinates().at(3) * ( thickness / 2 );
            } else { /*if (  type == IST_ShellForceTensor )*/
                z = 1;
            }
            double w = gp->giveWeight() * J * z;

            FloatArray localStress, localStrain;
            this->computeStrainVector(localStrain, gp, tStep);
            this->computeStressVector(localStress, localStrain, gp, tStep);
            mLocal += w * FloatArrayF< 6 >(localStress);
        }

        // local to global
        return StructuralMaterial::transformStressVectorTo(GtoLRotationMatrix, mLocal, false);
    } else if ( type == IST_CurvatureTensor ) {
        auto gp = integrationRulesArray [ 0 ]->getIntegrationPoint(nPointsZ * gpXY);
        const auto &coords = gp->giveNaturalCoordinates();

        //auto [hkx, hky] = this->givedNdx(coords);
        auto hk = this->givedNdx(coords);

        FloatArray dofs(24);
        this->computeVectorOf(VM_Total, tStep, dofs);

        FloatArrayF< 4 >rotX, rotY;
        for ( int i = 0; i < 4; i++ ) {
            rotX [ i ] = dofs.at(i * 6 + 4);
            rotY [ i ] = dofs.at(i * 6 + 5);
        }
        FloatArrayF< 6 >cLocal;
        cLocal.at(1) = dot(rotY, hk [ 0 ]);
        cLocal.at(2) = -dot(rotX, hk [ 1 ]);
        cLocal.at(6) = dot(rotY, hk [ 1 ]) - dot(rotX, hk [ 0 ]);

        return StructuralMaterial::transformStrainVectorTo(GtoLRotationMatrix, cLocal, false);
    } else if ( type == IST_ShellStrainTensor ) {
        auto gp = integrationRulesArray [ 0 ]->getIntegrationPoint(nPointsZ * gpXY);
        auto coords = gp->giveNaturalCoordinates();
        coords.at(3) = 0.; //set to midplane
        GaussIntegrationRule iRule(1, this, 1, 10);
        GaussPoint midGP(& iRule, 1, coords, 1, this->giveMaterialMode() );

        FloatArray answer;
        this->giveIPValue(answer, & midGP, IST_StrainTensor, tStep);
        return answer;
    } else {
        throw std::runtime_error("unknown type");
    }
}

std::array< FloatArrayF< 4 >, 2 >
MITC4Shell::givedNdx(const FloatArrayF< 3 > &coords)
{
    auto dn = interp_lin.evaldNdxi(coords [ { 0, 1 } ]);
    auto J = this->giveJacobian(coords);
    auto invJ = inv(J);
    auto invJ2 = invJ({ 0, 1 }, { 0, 1 });
    auto dndx = dot(invJ2, dn);

    auto hkx = dndx.row< 0 >();
    auto hky = dndx.row< 1 >();
    return { hkx, hky };
}

void
MITC4Shell::setupIRForMassMtrxIntegration(IntegrationRule &iRule)
{
    iRule.setUpIntegrationPoints(this->giveIntegrationDomain(), nPointsXY, nPointsZ, this->giveMaterialMode() );
}

int
MITC4Shell::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_StrainTensor ) {
        auto globTensor = this->giveCharacteristicTensor(GlobalStrainTensor, gp, tStep);
        answer = to_voigt_strain(globTensor);
        return 1;
    } else if ( type == IST_StressTensor ) {
        auto globTensor = this->giveCharacteristicTensor(GlobalForceTensor, gp, tStep);
        answer = to_voigt_stress(globTensor);
        return 1;
    } else if ( type == IST_ShellMomentTensor || type == IST_ShellForceTensor || type == IST_CurvatureTensor || type == IST_ShellStrainTensor ) {
        int gpnXY = ( gp->giveNumber() - 1 ) / 2;
        answer = this->giveMidplaneIPValue(gpnXY, type, tStep);

        return 1;
    } else {
        return StructuralElement::giveIPValue(answer, gp, type, tStep);
    }
}


bool
MITC4Shell::computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
//converts global coordinates to local planar area coordinates,
//does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
{
    // rotate the input point Coordinate System into the element CS
    FloatArray llc;
    auto inputCoords_ElCS = this->giveLocalCoordinates(coords);
    std::vector< FloatArray >lc(3);
    for ( int _i = 0; _i < 4; _i++ ) {
        lc [ _i ] = this->giveLocalCoordinates( this->giveNode(_i + 1)->giveCoordinates() );
    }
    bool inplane = interp_lin.global2local(llc, inputCoords_ElCS, FEIVertexListGeometryWrapper(lc) ) > 0;
    answer.resize(2);
    answer.at(1) = inputCoords_ElCS.at(1);
    answer.at(2) = inputCoords_ElCS.at(2);
    GaussPoint _gp(nullptr, 1, answer, 2.0, _2dPlate);
    // now check if the third local coordinate is within the thickness of element
    bool outofplane = ( fabs(inputCoords_ElCS.at(3) ) <= this->giveCrossSection()->give(CS_Thickness, & _gp) / 2. );

    return inplane && outofplane;
}



int
MITC4Shell::computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    FloatMatrix N;
    computeNmatrixAt(lcoords, N);

    answer.resize(3);
    for ( int _i = 1; _i <= 3; _i++ ) {
        answer.at(_i) = N.at(_i, _i) * this->giveNode(1)->giveCoordinate(_i) + N.at(_i, _i + 6) * this->giveNode(2)->giveCoordinate(_i) + N.at(_i, _i + 12) * this->giveNode(3)->giveCoordinate(_i) + N.at(_i, _i + 18) * this->giveNode(4)->giveCoordinate(_i);
    }
    return true;
}


int
MITC4Shell::computeLoadGToLRotationMtrx(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [5,6]
// f(local) = T * f(global)
{
    answer.resize(6, 6);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = answer.at(4, i + 3) = GtoLRotationMatrix.at(1, i);
        answer.at(2, i) = answer.at(5, i + 3) = GtoLRotationMatrix.at(2, i);
        answer.at(3, i) = answer.at(6, i + 3) = GtoLRotationMatrix.at(3, i);
    }

    return 1;
}


void
MITC4Shell::NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                       InternalStateType type, TimeStep *tStep)
{
    FloatMatrixF< 3, 3 >A;
    std::vector< FloatArrayF< 3 > >r;
    FloatArray val;
    int size = 0;
    for ( GaussPoint *gp : * integrationRulesArray [ 0 ] ) {
        giveIPValue(val, gp, type, tStep);
        if ( size == 0 ) {
            size = val.giveSize();
            r.resize(size);
        }

        const auto &coord = gp->giveNaturalCoordinates();
        double u = coord.at(1);
        double v = coord.at(2);

        A.at(1, 1) += 1;
        A.at(1, 2) += u;
        A.at(1, 3) += v;
        A.at(2, 1) += u;
        A.at(2, 2) += u * u;
        A.at(2, 3) += u * v;
        A.at(3, 1) += v;
        A.at(3, 2) += v * u;
        A.at(3, 3) += v * v;

        for ( int j = 1; j <= size; j++ ) {
            double y = val.at(j);
            r [ j ].at(1) += y;
            r [ j ].at(2) += y * u;
            r [ j ].at(3) += y * v;
        }
    }

    auto invA = inv(A);
    std::vector< FloatArrayF< 3 > >b(size);
    for ( int i = 0; i < size; ++i ) {
        b [ i ] = dot(invA, r [ i ]);
    }

    double x1 = 0.0, x2 = 0.0;
    switch ( node ) {
    case 1:
        x1 =  1.0;
        x2 =  1.0;
        break;
    case 2:
        x1 = -1.0;
        x2 =  1.0;
        break;
    case 3:
        x1 = -1.0;
        x2 = -1.0;
        break;
    case 4:
        x1 =  1.0;
        x2 = -1.0;
        break;
    default:
        OOFEM_ERROR("unsupported node");
    }

    answer.resize(size);
    for ( int j = 1; j <= size; j++ ) {
        answer.at(j) = b [ j ].at(1) + x1 * b [ j ].at(2) + x2 * b [ j ].at(3);
    }
}


void
MITC4Shell::giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    if ( iEdge == 1 ) { // edge between nodes 1 2
        answer = {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
        };
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer = {
            7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
        };
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
        answer = {
            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24
        };
    } else if ( iEdge == 4 ) { // edge between nodes 4 1
        answer = {
            19, 20, 21, 22, 23, 24, 1, 2, 3, 4, 5, 6
        };
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}


double
MITC4Shell::computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    auto lcF = this->giveNodeCoordinates();
    std::vector< FloatArray >lc;
    lc [ 0 ] = lcF [ 0 ];
    lc [ 1 ] = lcF [ 1 ];
    lc [ 2 ] = lcF [ 2 ];
    lc [ 3 ] = lcF [ 3 ];
    double detJ = this->interp_lin.edgeGiveTransformationJacobian(iEdge, gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lc) );
    return detJ * gp->giveWeight();
}

#if 0
void
MITC4Shell::computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    auto lc = this->giveNodeCoordinates();

    FloatArray local;
    this->interp_lin.edgeLocal2global(local, iEdge, gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lc) );
    local.resize(3);
    local.at(3) = 0.;

    answer = local; // MITC4 - todo
    // answer.beProductOf(this->lcsMatrix, local);
}
#endif

int
MITC4Shell::computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    auto e = this->computeLocalBaseVectors();

    const auto &edgeNodes = this->interp_lin.computeLocalEdgeMapping(iEdge);

    auto n1 = FloatArrayF< 3 >( this->giveNode(edgeNodes.at(1) )->giveCoordinates() );
    auto n2 = FloatArrayF< 3 >( this->giveNode(edgeNodes.at(2) )->giveCoordinates() );

    auto xl = normalize(n1 - n2);
    auto yl = cross(e [ 2 ], xl);

    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) = answer.at(4, 4) = dot(e [ 0 ], xl);
    answer.at(1, 2) = answer.at(4, 5) = dot(e [ 0 ], yl);
    answer.at(1, 3) = answer.at(4, 6) = dot(e [ 0 ], e [ 2 ]);
    answer.at(2, 1) = answer.at(5, 4) = dot(e [ 1 ], xl);
    answer.at(2, 2) = answer.at(5, 5) = dot(e [ 1 ], yl);
    answer.at(2, 3) = answer.at(5, 6) = dot(e [ 1 ], e [ 2 ]);
    answer.at(3, 1) = answer.at(6, 4) = dot(e [ 2 ], xl);
    answer.at(3, 2) = answer.at(6, 5) = dot(e [ 2 ], yl);
    answer.at(3, 3) = answer.at(6, 6) = dot(e [ 2 ], e [ 2 ]);

    return 1;
}



void
MITC4Shell::computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    const auto &coords2 = sgp->giveNaturalCoordinates();
    FloatArray coords = { coords2 [ 0 ], coords [ 1 ], 0. };
    this->computeNmatrixAt(coords, answer);
}

void
MITC4Shell::giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    if ( iSurf == 1 ) {
        answer.enumerate(24);
    } else {
        OOFEM_ERROR("wrong surface number");
    }
}

double
MITC4Shell::computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    FloatArrayF< 2 >lcoords = {
        gp->giveNaturalCoordinate(1),
        gp->giveNaturalCoordinate(2),
    };

    auto xyz = this->giveNodeCoordinates();
    auto dn = interp_lin.evaldNdxi(lcoords);

    FloatMatrixF< 2, 2 >jacobianMatrix;
    for ( std::size_t i = 0; i < dn.cols(); i++ ) {
        double x = xyz [ i ] [ 0 ];
        double y = xyz [ i ] [ 1 ];

        jacobianMatrix(0, 0) += dn(0, i) * x;
        jacobianMatrix(0, 1) += dn(0, i) * y;
        jacobianMatrix(1, 0) += dn(1, i) * x;
        jacobianMatrix(1, 1) += dn(1, i) * y;
    }

    return det(jacobianMatrix) * gp->giveWeight();
}

void
MITC4Shell::computeEdgeNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords)
{
    FloatArray n_vec;
    this->giveInterpolation()->boundaryEdgeEvalN(n_vec, boundaryID, lcoords, FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n_vec, 6);
}


void
MITC4Shell::computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords)
{
    FloatArray n_vec;
    this->giveInterpolation()->boundarySurfaceEvalN(n_vec, boundaryID, lcoords, FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n_vec, 6);
}
} // end namespace oofem
