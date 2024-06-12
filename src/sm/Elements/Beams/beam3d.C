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

#include "sm/Elements/Beams/beam3d.h"
#include "sm/Materials/structuralms.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "engngm.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fei3dlinelin.h"
#include "classfactory.h"
#include "elementinternaldofman.h"
#include "masterdof.h"
#include "bctracker.h"

#include "bodyload.h"
#include "boundaryload.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Beam3d);

FEI3dLineLin Beam3d :: interp;

Beam3d :: Beam3d(int n, Domain *aDomain) : BeamBaseElement(n, aDomain)
{
    numberOfDofMans = 2;
    referenceNode = 0;

    numberOfGaussPoints = 3;
    length = 0.;
    kappay = kappaz = -1.0;

    ghostNodes [ 0 ] = ghostNodes [ 1 ] = NULL;
    numberOfCondensedDofs = 0;
    subsoilMat = 0;
}

Beam3d :: ~Beam3d()
{
    delete ghostNodes [ 0 ];
    delete ghostNodes [ 1 ];
}

FEInterpolation *Beam3d :: giveInterpolation() const { return & interp; }

void
Beam3d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// eeps = {\eps_x, \gamma_xz, \gamma_xy, \der{phi_x}{x}, \kappa_y, \kappa_z}^T
{
    double l, ksi, kappay, kappaz, c1y, c1z;
    TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();

    l     = this->computeLength();
    ksi   = 0.5 + 0.5 * gp->giveNaturalCoordinate(1);
    kappay = this->giveKappayCoeff(tStep);
    kappaz = this->giveKappazCoeff(tStep);
    c1y = 1. + 2. * kappay;
    c1z = 1. + 2. * kappaz;

    answer.resize(6, 12);
    answer.zero();

    answer.at(1, 1) =  -1. / l;
    answer.at(1, 7) =   1. / l;

    answer.at(2, 3) =   ( -2. * kappay ) / ( l * c1y );
    answer.at(2, 5) =   kappay / ( c1y );
    answer.at(2, 9) =   2. * kappay / ( l * c1y );
    answer.at(2, 11) =  kappay / ( c1y );

    answer.at(3, 2) =   ( -2. * kappaz ) / ( l * c1z );
    answer.at(3, 6) =   -kappaz / ( c1z );
    answer.at(3, 8) =   2. * kappaz / ( l * c1z );
    answer.at(3, 12) =  -kappaz / ( c1z );


    answer.at(4, 4) =  -1. / l;
    answer.at(4, 10) =   1. / l;

    answer.at(5, 3) =   ( 6. - 12. * ksi ) / ( l * l * c1y );
    answer.at(5, 5) =   ( -2. * ( 2. + kappay ) + 6. * ksi ) / ( l * c1y );
    answer.at(5, 9) =   ( -6. + 12. * ksi ) / ( l * l * c1y );
    answer.at(5, 11) =   ( -2. * ( 1. - kappay ) + 6. * ksi ) / ( l * c1y );

    answer.at(6, 2) =   -1.0 * ( 6. - 12. * ksi ) / ( l * l * c1z );
    answer.at(6, 6) =   ( -2. * ( 2. + kappaz ) + 6. * ksi ) / ( l * c1z );
    answer.at(6, 8) =   -1.0 * ( -6. + 12. * ksi ) / ( l * l * c1z );
    answer.at(6, 12) =  ( -2. * ( 1. - kappaz ) + 6. * ksi ) / ( l * c1z );
}


void Beam3d :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        // the gauss point is used only when methods from crosssection and/or material
        // classes are requested
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ] = std :: make_unique< GaussIntegrationRule >(1, this, 1, 2);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}


void
Beam3d :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp. Used for numerical calculation of consistent mass
// matrix. Must contain only interpolation for displacement terms,
// not for any rotations. (Inertia forces do not work on rotations).
// r = {u1,v1,w1,fi_x1,fi_y1,fi_z1,u2,v2,w2,fi_x2,fi_y2,fi_21}^T
{
    double l, ksi, ksi2, ksi3, kappay, kappaz, c1y, c1z;
    TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();

    l     = this->computeLength();
    ksi =   0.5 + 0.5 * iLocCoord.at(1);
    kappay = this->giveKappayCoeff(tStep);
    kappaz = this->giveKappazCoeff(tStep);
    c1y = 1. + 2. * kappay;
    c1z = 1. + 2. * kappaz;
    ksi2 = ksi * ksi;
    ksi3 = ksi2 * ksi;

    answer.resize(6, 12);
    answer.zero();

    answer.at(1, 1) = 1. - ksi;
    answer.at(1, 7) = ksi;
    answer.at(2, 2) = ( c1z - 2. * kappaz * ksi - 3. * ksi2 + 2. * ksi3 ) / c1z;
    answer.at(2, 6) = -l * ( -( 1. + kappaz ) * ksi + ( 2. + kappaz ) * ksi2 - ksi3 ) / c1z;
    answer.at(2, 8) = ( 2. * kappaz * ksi + 3. * ksi2 - 2. * ksi3 ) / c1z;
    answer.at(2, 12) = -l * ( kappaz * ksi + ( 1. - kappaz ) * ksi2 - ksi3 ) / c1z;
    answer.at(3, 3) = ( c1y - 2. * kappay * ksi - 3. * ksi2 + 2. * ksi3 ) / c1y;
    answer.at(3, 5) = l * ( -( 1. + kappay ) * ksi + ( 2. + kappay ) * ksi2 - ksi3 ) / c1y;
    answer.at(3, 9) = ( 2. * kappay * ksi + 3. * ksi2 - 2. * ksi3 ) / c1y;
    answer.at(3, 11) = l * ( kappay * ksi + ( 1. - kappay ) * ksi2 - ksi3 ) / c1y;

    // rotations excluded
    answer.at(4, 4) = 1. - ksi;
    answer.at(4, 10) = ksi;
    answer.at(5, 3) = ( 6. * ksi - 6. * ksi2 ) / ( l * c1y );
    answer.at(5, 5) = ( c1y - 2. * ( 2. + kappay ) * ksi + 3. * ksi2 ) / c1y;
    answer.at(5, 9) = -( 6. * ksi - 6. * ksi2 ) / ( l * c1y );
    answer.at(5, 11) = ( -2. * ( 1. - kappay ) * ksi + 3. * ksi2 ) / c1y;
    answer.at(6, 2) = -( 6. * ksi - 6. * ksi2 ) / ( l * c1z );
    answer.at(6, 6) = ( c1z - 2. * ( 2. + kappaz ) * ksi + 3. * ksi2 ) / c1z;
    answer.at(6, 8) = ( 6. * ksi - 6. * ksi2 ) / ( l * c1z );
    answer.at(6, 12) = ( -2. * ( 1. - kappaz ) * ksi + 3. * ksi2 ) / c1z;
}


void
Beam3d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    double l = this->computeLength();
    FloatMatrix B, DB, d;
    answer.clear();
    for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {
        this->computeBmatrixAt(gp, B);
        this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
        double dV = gp->giveWeight() * 0.5 * l;
        DB.beProductOf(d, B);
        answer.plusProductSymmUpper(B, DB, dV);
    }
    answer.symmetrized();

    if ( subsoilMat ) {
        FloatMatrix k;
        this->computeSubSoilStiffnessMatrix(k, rMode, tStep);
        answer.add(k);
    }
}


void
Beam3d :: computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep, bool global)
{
    answer.clear();

    if ( edge != 1 ) {
        OOFEM_ERROR("Beam3D only has 1 edge (the midline) that supports loads. Attempted to apply load to edge %d", edge);
    }

    if ( type != ExternalForcesVector ) {
        return;
    }

    double l = this->computeLength();
    FloatArray coords, t;
    FloatMatrix N, T;

    for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        this->computeNmatrixAt(lcoords, N);
        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValues(t, tStep, lcoords, { D_u, D_v, D_w, R_u, R_v, R_w }, mode);
        } else {
            this->computeGlobalCoordinates(coords, lcoords);
            load->computeValues(t, tStep, coords, { D_u, D_v, D_w, R_u, R_v, R_w }, mode);
        }

        if ( load->giveCoordSystMode() == Load :: CST_Global ) {
            if ( this->computeLoadGToLRotationMtrx(T) ) {
                t.rotatedWith(T, 'n');
            }
        }

        double dl = gp->giveWeight() * 0.5 * l;
        answer.plusProduct(N, t, dl);
    }

    if ( global ) {
        // Loads from sets expects global c.s.
        this->computeGtoLRotationMatrix(T);
        answer.rotatedWith(T, 't');
    }
}


int
Beam3d :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
{
    FloatMatrix lcs;

    answer.resize(6, 6);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(3 + i, 3 + j) = lcs.at(i, j);
        }
    }
    return 1;
}


bool
Beam3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    int ndofs = computeNumberOfGlobalDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
            answer.at(i + 6, j + 6) = lcs.at(i, j);
            answer.at(i + 9, j + 9) = lcs.at(i, j);
        }
    }

    for ( int i = 13; i <= ndofs; i++ ) {
        answer.at(i, i) = 1.0;
    }

    if ( this->hasDofs2Condense() ) {
        int condensedDofCounter = 0;
        DofIDItem dofids[] = {
            D_u, D_v, D_w, R_u, R_v, R_w
        };
        FloatMatrix l2p(12, ndofs); // local -> primary
        l2p.zero();
        // loop over nodes
        for ( int inode = 0; inode < 2; inode++ ) {
            // loop over DOFs
            for ( int idof = 0; idof < 6; idof++ ) {
                int eq = inode * 6 + idof + 1;
                if ( ghostNodes [ inode ] ) {
                    if ( ghostNodes [ inode ]->hasDofID(dofids [ idof ]) ) {
                        condensedDofCounter++;
                        l2p.at(eq, 12 + condensedDofCounter) = 1.0;
                        continue;
                    }
                }
                l2p.at(eq, eq) = 1.0;
            }
        }

        FloatMatrix g2l(answer);
        answer.beProductOf(l2p, g2l);
    }

    return true;
}



FloatMatrixF<6,6>
Beam3d :: B3SSMI_getUnknownsGtoLRotationMatrix() const
// Returns the rotation matrix for element unknowns
{
    FloatMatrix lcs;

    FloatMatrixF<6,6> answer;

    const_cast<Beam3d*>(this)->giveLocalCoordinateSystem(lcs);
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
        }
    }
    return answer;
}

double
Beam3d :: computeVolumeAround(GaussPoint *gp)
{
    return gp->giveWeight() * 0.5 * this->computeLength();
}


int
Beam3d :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_BeamForceMomentTensor ) {
        answer = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        return 1;
    } else if ( type == IST_BeamStrainCurvatureTensor ) {
        answer = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        return 1;
    } else if ( type == IST_ShellForceTensor || type == IST_ShellStrainTensor ) {
        // Order in generalized strain is:
        // {\eps_x, \gamma_xz, \gamma_xy, \der{phi_x}{x}, \kappa_y, \kappa_z}
        const FloatArray &help = type == IST_ShellForceTensor ?
                                 static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector() :
                                 static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();

        answer.resize(6);
        answer.at(1) = help.at(1); // nx
        answer.at(2) = 0.0;        // ny
        answer.at(3) = 0.0;        // nz
        answer.at(4) = 0.0;        // vyz
        answer.at(5) = help.at(2); // vxz
        answer.at(6) = help.at(3); // vxy
        return 1;
    } else if ( type == IST_ShellMomentTensor || type == IST_CurvatureTensor ) {
        const FloatArray &help = type == IST_ShellMomentTensor ?
                                 static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector() :
                                 static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        answer.resize(6);
        answer.at(1) = help.at(4); // mx
        answer.at(2) = 0.0;        // my
        answer.at(3) = 0.0;        // mz
        answer.at(4) = 0.0;        // myz
        answer.at(5) = help.at(6); // mxz
        answer.at(6) = help.at(5); // mxy
        return 1;
    } else {
        return BeamBaseElement :: giveIPValue(answer, gp, type, tStep);
    }
}


void
Beam3d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}


double
Beam3d :: computeLength()
{
    double dx, dy, dz;
    Node *nodeA, *nodeB;

    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
        dz      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        length  = sqrt(dx * dx + dy * dy + dz * dz);
    }

    return length;
}


void
Beam3d :: computeKappaCoeffs(TimeStep *tStep)
{
    // kappa_y = (6*E*Iy)/(k*G*A*l^2)

    FloatMatrix d;
    double l = this->computeLength();

    this->computeConstitutiveMatrixAt(d, ElasticStiffness, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    //  kappay = 6. * d.at(5, 5) / ( d.at(3, 3) * l * l );
    //  kappaz = 6. * d.at(6, 6) / ( d.at(2, 2) * l * l );
    if ( d.at(3, 3) != 0. ) {
        kappay = 6. * d.at(5, 5) / ( d.at(3, 3) * l * l );
    } else {
        kappay = 0.;
    }
    if ( d.at(2, 2) != 0. ) {
        kappaz = 6. * d.at(6, 6) / ( d.at(2, 2) * l * l );
    } else {
        kappaz = 0.;
    }
}


double
Beam3d :: giveKappayCoeff(TimeStep *tStep)
{
    if ( kappay < 0.0 ) {
        this->computeKappaCoeffs(tStep);
    }

    return kappay;
}


double
Beam3d :: giveKappazCoeff(TimeStep *tStep)
{
    if ( kappaz < 0.0 ) {
        this->computeKappaCoeffs(tStep);
    }

    return kappaz;
}


int
Beam3d :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx, ly, lz, help(3);
    Node *nodeA, *nodeB;
    nodeA = this->giveNode(1);
    nodeB = this->giveNode(2);

    lx.beDifferenceOf( nodeB->giveCoordinates(), nodeA->giveCoordinates() );
    lx.normalize();

    if ( this->referenceNode ) {
        Node *refNode = this->giveDomain()->giveNode(this->referenceNode);
        help.beDifferenceOf( refNode->giveCoordinates(), nodeA->giveCoordinates() );

        lz.beVectorProductOf(lx, help);
        lz.normalize();
    } else if ( this->zaxis.giveSize() > 0 ) {
        lz = this->zaxis;
        lz.normalize();
        lz.add(lz.dotProduct(lx), lx);
        lz.normalize();
    } else if ( this->yaxis.giveSize() > 0 ) {
        ly = this->yaxis;
        ly.normalize();
        ly.add(ly.dotProduct(lx), lx);
        ly.normalize();
        lz.beVectorProductOf(ly, lx);
        lz.normalize();
    } else {
        FloatMatrix rot(3, 3);
        double theta = referenceAngle * M_PI / 180.0;

        rot.at(1, 1) = cos(theta) + pow(lx.at(1), 2) * ( 1 - cos(theta) );
        rot.at(1, 2) = lx.at(1) * lx.at(2) * ( 1 - cos(theta) ) - lx.at(3) * sin(theta);
        rot.at(1, 3) = lx.at(1) * lx.at(3) * ( 1 - cos(theta) ) + lx.at(2) * sin(theta);

        rot.at(2, 1) = lx.at(2) * lx.at(1) * ( 1 - cos(theta) ) + lx.at(3) * sin(theta);
        rot.at(2, 2) = cos(theta) + pow(lx.at(2), 2) * ( 1 - cos(theta) );
        rot.at(2, 3) = lx.at(2) * lx.at(3) * ( 1 - cos(theta) ) - lx.at(1) * sin(theta);

        rot.at(3, 1) = lx.at(3) * lx.at(1) * ( 1 - cos(theta) ) - lx.at(2) * sin(theta);
        rot.at(3, 2) = lx.at(3) * lx.at(2) * ( 1 - cos(theta) ) + lx.at(1) * sin(theta);
        rot.at(3, 3) = cos(theta) + pow(lx.at(3), 2) * ( 1 - cos(theta) );

        help.at(3) = 1.0;         // up-vector
        // here is ly is used as a temp var
        if ( fabs( lx.dotProduct(help) ) > 0.999 ) { // Check if it is vertical
            ly = {
                0., 1., 0.
            };
        } else {
            ly.beVectorProductOf(lx, help);
        }
        lz.beProductOf(rot, ly);
        lz.normalize();
    }

    ly.beVectorProductOf(lz, lx);
    ly.normalize();

    answer.resize(3, 3);
    answer.zero();
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}


void
Beam3d :: initializeFrom(InputRecord &ir)
{
    BeamBaseElement :: initializeFrom(ir);

    referenceNode = 0;
    referenceAngle = 0;
    this->yaxis.clear();
    this->zaxis.clear();
    
    if ( ir.hasField(_IFT_Beam3d_yaxis) ) {
        IR_GIVE_FIELD(ir, this->yaxis, _IFT_Beam3d_yaxis);
    } else if ( ir.hasField(_IFT_Beam3d_zaxis) ) {
        IR_GIVE_FIELD(ir, this->zaxis, _IFT_Beam3d_zaxis);
    } else if ( ir.hasField(_IFT_Beam3d_refnode) ) {
        IR_GIVE_FIELD(ir, referenceNode, _IFT_Beam3d_refnode);
        if ( referenceNode == 0 ) {
            OOFEM_WARNING("wrong reference node specified. Using default orientation.");
        }
    } else if ( ir.hasField(_IFT_Beam3d_refangle) ) {
        IR_GIVE_FIELD(ir, referenceAngle, _IFT_Beam3d_refangle);
    } else {
        throw ValueInputException(ir, _IFT_Beam3d_zaxis, "axis, reference node, or angle not set");
    }

    if ( ir.hasField(_IFT_Beam3d_dofstocondense) ) {
        IntArray val;
        IR_GIVE_FIELD(ir, val, _IFT_Beam3d_dofstocondense);
        if ( val.giveSize() >= 12 ) {
            throw ValueInputException(ir, _IFT_Beam3d_dofstocondense, "wrong input data for condensed dofs");
        }

        //dofsToCondense = new IntArray(val);
        DofIDItem mask[] = {
            D_u, D_v, D_w, R_u, R_v, R_w
        };
        this->numberOfCondensedDofs = val.giveSize();
        for ( int i = 1; i <= val.giveSize(); i++ ) {
            if ( val.at(i) <= 6 ) {
                if ( ghostNodes [ 0 ] == NULL ) {
                    ghostNodes [ 0 ] = new ElementDofManager(1, giveDomain(), this);
                }
                ghostNodes [ 0 ]->appendDof( new MasterDof(ghostNodes [ 0 ], mask [ val.at(i) - 1 ]) );
            } else {
                if ( ghostNodes [ 1 ] == NULL ) {
                    ghostNodes [ 1 ] = new ElementDofManager(2, giveDomain(), this);
                }
                ghostNodes [ 1 ]->appendDof( new MasterDof(ghostNodes [ 1 ], mask [ val.at(i) - 7 ]) );
            }
        }
    } else {
        //dofsToCondense = NULL;
    }

    this->subsoilMat = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->subsoilMat, _IFT_Beam3d_subsoilmat);
}


void
Beam3d :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
#if 0
    FloatMatrix stiffness;
    FloatArray u;

    this->computeStiffnessMatrix(stiffness, SecantStiffness, tStep);
    this->computeVectorOf(VM_Total, tStep, u);
    answer.beProductOf(stiffness, u);
#else
    BeamBaseElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
    if ( subsoilMat ) {
        // add internal forces due to subsoil interaction
        // @todo: linear subsoil assumed here; more general approach should integrate internal forces
        FloatMatrix k;
        FloatArray u, F;
        this->computeSubSoilStiffnessMatrix(k, TangentStiffness, tStep);
        this->computeVectorOf(VM_Total, tStep, u);
        F.beProductOf(k, u);
        answer.add(F);
    }
#endif
}


void
Beam3d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->give3dBeamStiffMtrx(rMode, gp, tStep);
}


void
Beam3d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveGeneralizedStress_Beam3d(strain, gp, tStep);
}


void
Beam3d :: giveEndForcesVector(FloatArray &answer, TimeStep *tStep)
{
    // computes exact global end-forces vector
    FloatArray loadEndForces;

    this->giveInternalForcesVector(answer, tStep);

    // add exact end forces due to nonnodal loading
    this->computeLocalForceLoadVector(loadEndForces, tStep, VM_Total); // will compute only contribution of loads applied directly on receiver (not using sets)
    if ( loadEndForces.giveSize() ) {
        answer.subtract(loadEndForces);
    }
    /*
     * if (subsoilMat) {
     * // @todo: linear subsoil assumed here; more general approach should integrate internal forces
     * FloatMatrix k;
     * FloatArray u, F;
     * this->computeSubSoilStiffnessMatrix(k, TangentStiffness, tStep);
     * this->computeVectorOf(VM_Total, tStep, u);
     * F.beProductOf(k, u);
     * answer.add(F);
     * }
     */
}


void
Beam3d :: printOutputAt(FILE *File, TimeStep *tStep)
{
    FloatArray rl, Fl;

    fprintf( File, "beam element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    // ask for global element displacement vector
    this->computeVectorOf(VM_Total, tStep, rl);
    // ask for global element end forces vector
    this->giveEndForcesVector(Fl, tStep);

    fprintf(File, "  local displacements ");
    for ( auto &val : rl ) {
        fprintf(File, " %.4e", val);
    }

    fprintf(File, "\n  local end forces    ");
    for ( auto &val : Fl ) {
        fprintf(File, " %.4e", val);
    }

    fprintf(File, "\n");

    for ( auto &iRule : integrationRulesArray ) {
        iRule->printOutputAt(File, tStep);
    }
}


void
Beam3d :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
{
    FloatArray lc(1);
    BeamBaseElement :: computeBodyLoadVectorAt(answer, load, tStep, mode);
    answer.times( this->giveCrossSection()->give(CS_Area, lc, this) );
}


void
Beam3d :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity)
{
  if (0) {
    StructuralElement::computeConsistentMassMatrix(answer, tStep, mass, ipDensity);
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    double area = this->giveCrossSection()->give(CS_Area, gp); // constant area assumed
    answer.times(area);
    mass *= area;
    return;
  } else {
    

  FloatMatrix stiff;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    /*
     * SructuralElement::computeMassMatrix(answer, tStep);
     * answer.times(this->giveCrossSection()->give('A'));
     */
    double l = this->computeLength();
    double kappay = this->giveKappayCoeff(tStep);
    double kappaz = this->giveKappazCoeff(tStep);
    double kappay2 = kappay * kappay;
    double kappaz2 = kappaz * kappaz;

    double density = this->giveStructuralCrossSection()->give('d', gp); // constant density assumed
    if ( ipDensity != NULL ) {
        // Override density if desired
        density = * ipDensity;
    }

    double area = this->giveCrossSection()->give(CS_Area, gp); // constant area assumed
    double Ikk = this->giveCrossSection()->give(CS_InertiaMomentY, gp) +
      this->giveCrossSection()->give(CS_InertiaMomentZ, gp); // polar moment of inertia 
    double c2y = ( area * density ) / ( ( 1. + 2. * kappay ) * ( 1. + 2. * kappay ) );
    double c2z = ( area * density ) / ( ( 1. + 2. * kappaz ) * ( 1. + 2. * kappaz ) );
    double c1 = ( area * density );

    answer.resize(12, 12);
    answer.zero();

    answer.at(1, 1) = c1 * l / 3.0;
    answer.at(1, 7) = c1 * l / 6.0;
    answer.at(2, 2) = c2z * l * ( 13. / 35. + 7. * kappaz / 5. + 4. * kappaz2 / 3. );
    answer.at(2, 6) = c2z * l * l * ( 11. / 210. + kappaz * 11. / 60. + kappaz2 / 6. );
    answer.at(2, 8) = c2z * l * ( 9. / 70. + kappaz * 3. / 5. + kappaz2 * 2. / 3. );
    answer.at(2, 12) = -c2z * l * l * ( 13. / 420. + kappaz * 3. / 20. + kappaz2 / 6. );
    answer.at(3, 3) = c2y * l * ( 13. / 35. + 7. * kappay / 5. + 4. * kappay2 / 3. );
    answer.at(3, 5) = -c2y * l * l * ( 11. / 210. + kappay * 11. / 60. + kappay2 / 6. );
    answer.at(3, 9) = c2y * l * ( 9. / 70. + kappay * 3. / 5. + kappay2 * 2. / 3. );
    answer.at(3, 11) = c2y * l * l * ( 13. / 420. + kappay * 3. / 20. + kappay2 / 6. );
    answer.at(4,4) = Ikk*density*l/3.;
    answer.at(4,10) = Ikk*density*l/6.;
    answer.at(5, 5) = c2y * l * l * l * ( 1. / 105. + kappay / 30. + kappay2 / 30. );
    answer.at(5, 9) = -c2y * l * l * ( 13. / 420. + kappay * 3. / 20. + kappay2 / 6. );
    answer.at(5, 11) = -c2y * l * l * l * ( 1. / 140. + kappay / 30. + kappay2 / 30. );
    answer.at(6, 6) = c2z * l * l * l * ( 1. / 105. + kappaz / 30. + kappaz2 / 30. );
    answer.at(6, 8) = c2z * l * l * ( 13. / 420. + kappaz * 3. / 20. + kappaz2 / 6. );
    answer.at(6, 12) = -c2z * l * l * l * ( 1. / 140. + kappaz / 30. + kappaz2 / 30. );


    answer.at(7, 7) = c1 * l / 3.0;
    answer.at(8, 8) = c2z * l * ( 13. / 35. + kappaz * 7. / 5. + kappaz2 * 4. / 3. );
    answer.at(8, 12) = -c2z * l * l * ( 11. / 210. + kappaz * 11. / 60. + kappaz2 / 6. );
    answer.at(9, 9) = c2y * l * ( 13. / 35. + kappay * 7. / 5. + kappay2 * 4. / 3. );
    answer.at(9, 11) = c2y * l * l * ( 11. / 210. + kappay * 11. / 60. + kappay2 / 6. );
    answer.at(10,10) = Ikk*density*l/3.;
    answer.at(11, 11) = c2y * l * l * l * ( 1. / 105. + kappay / 30. + kappay2 / 30. );
    answer.at(12, 12) = c2z * l * l * l * ( 1. / 105. + kappaz / 30. + kappaz2 / 30. );

    answer.symmetrized();

    mass = area * l * density;
  }
}


int
Beam3d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 *this->giveNode(2)->giveCoordinate(1);
    answer.at(2) = n1 * this->giveNode(1)->giveCoordinate(2) + n2 *this->giveNode(2)->giveCoordinate(2);
    answer.at(3) = n1 * this->giveNode(1)->giveCoordinate(3) + n2 *this->giveNode(2)->giveCoordinate(3);

    return 1;
}


void
Beam3d :: computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // computes initial stress matrix of receiver (or geometric stiffness matrix)

    FloatMatrix stiff;
    FloatArray endForces;

    double l = this->computeLength();
    double kappay = this->giveKappayCoeff(tStep);
    double kappaz = this->giveKappazCoeff(tStep);
    double kappay2 = kappay * kappay;
    double kappaz2 = kappaz * kappaz;
    double minVal;
    double denomy = ( 1. + 2. * kappay ) * ( 1. + 2. * kappay ), denomz = ( 1. + 2. * kappaz ) * ( 1. + 2. * kappaz );
    double N;

    answer.resize(12, 12);
    answer.zero();

    answer.at(2, 2) = ( 4. * kappaz2 + 4. * kappaz + 6. / 5. ) / denomz;
    answer.at(2, 6) = ( l / 10. ) / denomz;
    answer.at(2, 8) = ( -4. * kappaz2 - 4. * kappaz - 6. / 5. ) / denomz;
    answer.at(2, 12) = ( l / 10. ) / denomz;

    answer.at(3, 3) = ( 4. * kappay2 + 4. * kappay + 6. / 5. ) / denomy;
    answer.at(3, 5) = ( -l / 10. ) / denomy;
    answer.at(3, 9) = ( -4. * kappay2 - 4. * kappay - 6. / 5. ) / denomy;
    answer.at(3, 11) = ( -l / 10. ) / denomy;

    answer.at(5, 5) = l * l * ( kappay2 / 3. + kappay / 3. + 2. / 15. ) / denomy;
    answer.at(5, 9) = ( l / 10. ) / denomy;
    answer.at(5, 11) = -l * l * ( kappay2 / 3. + kappay / 3. + 1. / 30. ) / denomy;

    answer.at(6, 6) = l * l * ( kappaz2 / 3. + kappaz / 3. + 2. / 15. ) / denomz;
    answer.at(6, 8) = ( -l / 10. ) / denomz;
    answer.at(6, 12) = -l * l * ( kappaz2 / 3. + kappaz / 3. + 1. / 30. ) / denomz;

    answer.at(8, 8) = ( 4. * kappaz2 + 4. * kappaz + 6. / 5. ) / denomz;
    answer.at(8, 12) = ( -l / 10. ) / denomz;

    answer.at(9, 9) = ( 4. * kappay2 + 4. * kappay + 6. / 5. ) / denomy;
    answer.at(9, 11) = ( l / 10. ) / denomy;

    answer.at(11, 11) = l * l * ( kappay2 / 3. + kappay / 3. + 2. / 15. ) / denomy;
    answer.at(12, 12) = l * l * ( kappaz2 / 3. + kappaz / 3. + 2. / 15. ) / denomz;

    minVal = min( fabs( answer.at(2, 2) ), fabs( answer.at(3, 3) ) );
    minVal = min( minVal, fabs( answer.at(5, 5) ) );
    minVal = min( minVal, fabs( answer.at(6, 6) ) );

    answer.at(1, 1) = minVal / 1000.;
    answer.at(1, 7) = -answer.at(1, 1);
    answer.at(7, 7) = answer.at(1, 1);

    answer.at(4, 4)   = minVal / 1000.;
    answer.at(4, 10)  = -answer.at(4, 4);
    answer.at(10, 10) = answer.at(4, 4);



    answer.symmetrized();
    // ask end forces in g.c.s
    this->giveEndForcesVector(endForces, tStep);

    N = ( -endForces.at(1) + endForces.at(7) ) / 2.;
    answer.times(N / l);

    //answer.beLumpedOf (mass);
}


void
Beam3d :: computeLumpedInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // computes initial stress matrix of receiver (or geometric stiffness matrix)

    FloatMatrix stiff;
    FloatArray endForces;

    double l = this->computeLength();
    double N;

    answer.resize(12, 12);
    answer.zero();

    answer.at(2, 2) = 1.;
    answer.at(2, 8) =-1.;
    answer.at(8, 2) =-1.;
    answer.at(8, 8) = 1.;

    answer.at(3, 3) = 1.;
    answer.at(3, 9) =-1.;
    answer.at(9, 3) =-1.;
    answer.at(9, 9) = 1.;


    // ask end forces in g.c.s
      // ask end forces in g.c.s
    this->giveEndForcesVector(endForces, tStep);
    N = ( -endForces.at(1) + endForces.at(7) ) / 2.;
    answer.times(N / l);
    
}


void
Beam3d :: FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, const FloatArray &masterGpStrain,
                                                                  GaussPoint *slaveGp, TimeStep *tStep)
{
    double layerYCoord, layerZCoord;

    layerZCoord = slaveGp->giveNaturalCoordinate(2);
    layerYCoord = slaveGp->giveNaturalCoordinate(1);

    answer.resize(3);  // {Exx,GMzx,GMxy}

    answer.at(1) = masterGpStrain.at(1) + masterGpStrain.at(5) * layerZCoord - masterGpStrain.at(6) * layerYCoord;
    answer.at(2) = masterGpStrain.at(2) + layerYCoord * masterGpStrain.at(4);
    answer.at(3) = masterGpStrain.at(3) - layerZCoord * masterGpStrain.at(4);
}


Interface *
Beam3d :: giveInterface(InterfaceType interface)
{
    if ( interface == FiberedCrossSectionInterfaceType ) {
        return static_cast< FiberedCrossSectionInterface * >(this);
    } else if ( interface == Beam3dSubsoilMaterialInterfaceType ) {
        return static_cast< Beam3dSubsoilMaterialInterface * >(this);
    } else if ( interface == VTKXMLExportModuleElementInterfaceType ) {
        return static_cast< VTKXMLExportModuleElementInterface * >(this);
    }

    return NULL;
}


void
Beam3d :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    BeamBaseElement :: updateLocalNumbering(f);
    if ( this->referenceNode ) {
        this->referenceNode = f(this->referenceNode, ERS_DofManager);
    }
}


#ifdef __OOFEG
void Beam3d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    //  if (!go) { // create new one
    WCRec p [ 2 ];   /* poin */
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


void Beam3d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
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

void
Beam3d :: computeSubSoilNMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // only winkler model supported now (passing only unknown interpolation)
    this->computeNmatrixAt(gp->giveNaturalCoordinates(), answer);
}


void
Beam3d :: computeSubSoilStiffnessMatrix(FloatMatrix &answer,
                                        MatResponseMode rMode, TimeStep *tStep)
{
    double l = this->computeLength();
    FloatMatrix N, DN, d;
    answer.clear();
    for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
        this->computeSubSoilNMatrixAt(gp, N);
        d = static_cast<StructuralMaterial *>(this->domain->giveMaterial(subsoilMat))->give3dBeamSubSoilStiffMtrx(rMode, gp, tStep);
        double dV = gp->giveWeight() * 0.5 * l;
        DN.beProductOf(d, N);
        answer.plusProductSymmUpper(N, DN, dV);
    }
    answer.symmetrized();
}

int
Beam3d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords, const FloatArray &pCoords)
{
    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 *pCoords.at(1);
    answer.at(2) = n1 * this->giveNode(1)->giveCoordinate(2) + n2 *pCoords.at(2);
    answer.at(3) = n1 * this->giveNode(1)->giveCoordinate(3) + n2 *pCoords.at(3);

    return 1;
}

void
Beam3d :: giveInternalForcesVectorAtPoint(FloatArray &answer, TimeStep *tStep, FloatArray &coords)
{
    // computes exact global end-forces vector
    FloatArray loadEndForces, iF;
    IntArray leftIndx = {
        1, 2, 3, 4, 5, 6
    };
    this->giveEndForcesVector(iF, tStep);

    answer.beSubArrayOf(iF, leftIndx);
    Node *nodeA;

    nodeA   = this->giveNode(1);
    double dx      = nodeA->giveCoordinate(1) - coords.at(1);
    double dy      = nodeA->giveCoordinate(2) - coords.at(2);
    double dz      = nodeA->giveCoordinate(3) - coords.at(3);
    double ds  = sqrt(dx * dx + dy * dy + dz * dz);

    answer.at(5) += iF.at(3) * ds;
    answer.at(6) -= iF.at(2) * ds;



    // loop over body load array first
    int nBodyLoads = this->giveBodyLoadArray()->giveSize();
    FloatArray help;

    for ( int i = 1; i <= nBodyLoads; i++ ) {
        int id = bodyLoadArray.at(i);
        Load *load = domain->giveLoad(id);
        bcGeomType ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            this->computeInternalForcesFromBodyLoadVectorAtPoint(help, load, tStep, VM_Total, coords, ds); // this one is local
            answer.add(help);
        } else {
            if ( load->giveBCValType() != TemperatureBVT && load->giveBCValType() != EigenstrainBVT ) {
                // temperature and eigenstrain is handled separately at computeLoadVectorAt subroutine
                OOFEM_ERROR("body load %d is of unsupported type (%d)", id, ltype);
            }
        }
    }

    // loop over boundary load array
    int nBoundaryLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nBoundaryLoads; i++ ) {
        int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int id = boundaryLoadArray.at(i * 2);
        Load *load = domain->giveLoad(n);
        BoundaryLoad *bLoad;
        if ( ( bLoad = dynamic_cast< BoundaryLoad * >(load) ) ) {
            bcGeomType ltype = load->giveBCGeoType();
            if ( ltype == EdgeLoadBGT ) {
                this->computeInternalForcesFromBoundaryEdgeLoadVectorAtPoint(help, bLoad, id,
                                                                             ExternalForcesVector, VM_Total, tStep, coords, ds, false);
                answer.add(help);
            } else {
                OOFEM_ERROR("boundary load %d is of unsupported type (%d)", id, ltype);
            }
        }
    }
    // add exact end forces due to nonnodal loading
    //    this->computeForceLoadVectorAt(loadEndForces, tStep, VM_Total, coords); // will compute only contribution of loads applied directly on receiver (not using sets)
    if ( loadEndForces.giveSize() ) {
        answer.subtract(loadEndForces);
    }

    // add exact end forces due to nonnodal loading applied indirectly (via sets)
    BCTracker *bct = this->domain->giveBCTracker();
    BCTracker :: entryListType bcList = bct->getElementRecords(this->number);

    for ( BCTracker :: entryListType :: iterator it = bcList.begin(); it != bcList.end(); ++it ) {
        GeneralBoundaryCondition *bc = this->domain->giveBc( ( * it ).bcNumber );
        BodyLoad *bodyLoad;
        BoundaryLoad *boundaryLoad;
        if ( bc->isImposed(tStep) ) {
            if ( ( bodyLoad = dynamic_cast< BodyLoad * >(bc) ) ) { // body load
                this->computeInternalForcesFromBodyLoadVectorAtPoint(help, bodyLoad, tStep, VM_Total, coords, ds); // this one is local
                //answer.subtract(help);
            } else if ( ( boundaryLoad = dynamic_cast< BoundaryLoad * >(bc) ) ) {
                // compute Boundary Edge load vector in GLOBAL CS !!!!!!!
                this->computeInternalForcesFromBoundaryEdgeLoadVectorAtPoint(help, boundaryLoad, ( * it ).boundaryId,
                                                                             ExternalForcesVector, VM_Total, tStep, coords, ds, false);
            }
            answer.add(help);
        }
    }




    if ( subsoilMat ) {
        // @todo: linear subsoil assumed here; more general approach should integrate internal forces
        FloatMatrix k;
        FloatArray u, F;
        this->computeSubSoilStiffnessMatrix(k, TangentStiffness, tStep);
        this->computeVectorOf(VM_Total, tStep, u);
        F.beProductOf(k, u);
        answer.add(F);
    }


    answer.times(-1);
}


void
Beam3d :: computeInternalForcesFromBoundaryEdgeLoadVectorAtPoint(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep, FloatArray &pointCoords, double ds, bool global)
{
    answer.clear();

    if ( edge != 1 ) {
        OOFEM_ERROR("Beam3D only has 1 edge (the midline) that supports loads. Attempted to apply load to edge %d", edge);
    }

    if ( type != ExternalForcesVector ) {
        return;
    }

    FloatArray coords, t;
    FloatMatrix T;


    for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        this->computeGlobalCoordinates(coords, lcoords, pointCoords);
        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValues(t, tStep, lcoords, { D_u, D_v, D_w, R_u, R_v, R_w }, mode);
        } else {
            load->computeValues(t, tStep, coords, { D_u, D_v, D_w, R_u, R_v, R_w }, mode);
        }

        if ( load->giveCoordSystMode() == Load :: CST_Global ) {
            if ( this->computeLoadGToLRotationMtrx(T) ) {
                t.rotatedWith(T, 'n');
            }
        }


        double dl = gp->giveWeight() * 0.5 * ds;
        FloatArray f;
        f = t;
        f.at(5) += f.at(3) * ( lcoords.at(1) + 1 ) * ds / 2;
        f.at(6) -= f.at(2) * ( lcoords.at(1) + 1 ) * ds / 2;
        answer.add(dl, f);
    }

    if ( global ) {
        // Loads from sets expects global c.s.
        this->computeGtoLRotationMatrix(T);
        answer.rotatedWith(T, 't');
    }
}

void
Beam3d :: computeInternalForcesFromBodyLoadVectorAtPoint(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode, FloatArray &pointCoords, double ds)
// Computes numerically the load vector of the receiver due to the body
// loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV;
    FloatArray force, ntf;
    FloatMatrix n, T;
    FloatArray lc(1);

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, tStep, mode);
    force.times( this->giveCrossSection()->give(CS_Area, lc, this) );
    // transform from global to element local c.s
    if ( this->computeLoadGToLRotationMtrx(T) ) {
        force.rotatedWith(T, 'n');
    }

    answer.clear();

    if ( force.giveSize() ) {
        for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();
            this->computeNmatrixAt(gp->giveSubPatchCoordinates(), n);
            dV  = gp->giveWeight() * 0.5 * ds;
            dens = this->giveCrossSection()->give('d', gp);
            FloatArray iF;
            iF = force;
            iF.at(5) += force.at(3) *  ( lcoords.at(1) + 1 ) * ds / 2;
            iF.at(6) -= force.at(2) *  ( lcoords.at(1) + 1 ) * ds / 2;
            answer.add(dV * dens, iF);
        }
    } else {
        return;
    }
}


void
Beam3d :: giveCompositeExportData(std :: vector< ExportRegion > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep)
{
    // divide element into several small ones
    vtkPieces.resize(1);
    vtkPieces [ 0 ].setNumberOfCells(Beam3d_nSubBeams);
    int nNodes = 2 * Beam3d_nSubBeams;
    vtkPieces [ 0 ].setNumberOfNodes(nNodes);
    FloatArray nodeXi(nNodes), xi(1);

    Node *nodeA, *nodeB;
    nodeA = this->giveNode(1);
    nodeB   = this->giveNode(2);
    double dx = ( nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1) ) / Beam3d_nSubBeams;
    double dy = ( nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2) ) / Beam3d_nSubBeams;
    double dz = ( nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3) ) / Beam3d_nSubBeams;
    FloatArray coords(3);
    int nodeNumber = 1;
    int val = 1;
    int offset = 0;
    IntArray connectivity(2);
    for ( int i = 0; i < Beam3d_nSubBeams; i++ ) {
        for ( int j = 0; j < 2; j++ ) {
            coords.at(1) = nodeA->giveCoordinate(1) + ( i + j ) * dx;
            coords.at(2) = nodeA->giveCoordinate(2) + ( i + j ) * dy;
            coords.at(3) = nodeA->giveCoordinate(3) + ( i + j ) * dz;
            vtkPieces [ 0 ].setNodeCoords(nodeNumber, coords);
            nodeXi.at(nodeNumber) = -1.0 + ( 2.0 / Beam3d_nSubBeams ) * ( i + j );
            nodeNumber++;
            connectivity.at(j + 1) = val++;
        }
        vtkPieces [ 0 ].setConnectivity(i + 1, connectivity);
        offset += 2;
        vtkPieces [ 0 ].setOffset(i + 1, offset);
        vtkPieces [ 0 ].setCellType(i + 1, 3);
    }

    InternalStateType isttype;
    int n = internalVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfInternalVarsToExport(internalVarsToExport, nNodes);
    for ( int i = 1; i <= n; i++ ) {
        isttype = ( InternalStateType ) internalVarsToExport.at(i);
        for ( int nN = 1; nN <= nNodes; nN++ ) {
            if ( isttype == IST_BeamForceMomentTensor ) {
                FloatArray coords = vtkPieces [ 0 ].giveNodeCoords(nN);
                FloatArray endForces;
                this->giveInternalForcesVectorAtPoint(endForces, tStep, coords);
                vtkPieces [ 0 ].setInternalVarInNode(isttype, nN, endForces);
            } else if ( isttype == IST_X_LCS || isttype == IST_Y_LCS || isttype == IST_Z_LCS ) {
                FloatArray answer;
                this->giveLocalCoordinateSystemVector(isttype, answer);
                vtkPieces[0].setInternalVarInNode(isttype, nN, answer);
            } else {
                fprintf( stderr, "VTKXMLExportModule::exportIntVars: unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
            }
        }
    }

    n = primaryVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfPrimaryVarsToExport(primaryVarsToExport, nNodes);
    for ( int i = 1; i <= n; i++ ) {
        UnknownType utype = ( UnknownType ) primaryVarsToExport.at(i);
        if ( utype == DisplacementVector ) {
            FloatMatrix Tgl, n;
            FloatArray d(3);

            Tgl = this->B3SSMI_getUnknownsGtoLRotationMatrix();
            for ( int nN = 1; nN <= nNodes; nN++ ) {
                FloatArray u, dl, dg;
                this->computeVectorOf(VM_Total, tStep, u);
                xi.at(1) = nodeXi.at(nN);
                this->computeNmatrixAt(xi, n);
                dl.beProductOf(n, u); // local interpolated displacement
                dg.beTProductOf(Tgl, dl); // local displacement tranformed to global c.s.
                d.at(1) = dg.at(1);
                d.at(2) = dg.at(2);
                d.at(3) = dg.at(3);
                vtkPieces [ 0 ].setPrimaryVarInNode(utype, nN, d);
            }
        } else {
            fprintf( stderr, "VTKXMLExportModule::exportPrimaryVars: unsupported variable type %s\n", __UnknownTypeToString(utype) );
        }
    }
}
} // end namespace oofem
