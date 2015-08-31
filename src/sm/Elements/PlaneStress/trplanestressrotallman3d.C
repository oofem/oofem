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

#include "../sm/Elements/PlaneStress/trplanestressrotallman3d.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Materials/structuralms.h"
#include "material.h"
#include "node.h"
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "fei2dtrlin.h"

namespace oofem {
REGISTER_Element(TrPlanestressRotAllman3d);

TrPlanestressRotAllman3d :: TrPlanestressRotAllman3d(int n, Domain *aDomain) : TrPlanestressRotAllman(n, aDomain)
{
    GtoLRotationMatrix = NULL;
}


void
TrPlanestressRotAllman3d :: computeLocalNodalCoordinates(std :: vector< FloatArray > &lxy)
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
    // test if pereviously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }


    lxy.resize(6);
    const FloatArray *nc;
    for ( int i = 0; i < 3; i++ ) {
        nc = this->giveNode(i + 1)->giveCoordinates();
        lxy [ i ].beProductOf(* GtoLRotationMatrix, * nc);
    }
    lxy [ 3 ].resize(3);
    lxy [ 4 ].resize(3);
    lxy [ 5 ].resize(3);
    for ( int i = 1; i <= 3; i++ ) {
        lxy [ 3 ].at(i) = 0.5 * ( lxy [ 0 ].at(i) + lxy [ 1 ].at(i) );
        lxy [ 4 ].at(i) = 0.5 * ( lxy [ 1 ].at(i) + lxy [ 2 ].at(i) );
        lxy [ 5 ].at(i) = 0.5 * ( lxy [ 2 ].at(i) + lxy [ 0 ].at(i) );
    }
}


double
TrPlanestressRotAllman3d :: computeVolumeAround(GaussPoint *gp)
{
    double detJ, weight;

    std :: vector< FloatArray >lc;
    this->computeLocalNodalCoordinates(lc);

    weight = gp->giveWeight();
    detJ = fabs( this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lc) ) );
    return detJ * weight * this->giveStructuralCrossSection()->give(CS_Thickness, gp);
}


void
TrPlanestressRotAllman3d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}


const FloatMatrix *
TrPlanestressRotAllman3d :: computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiver of the size [3,3]
// coords(local) = T * coords(global)
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N2-N1]    Ni - means i - th node
// help   : [N3-N1]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
    if ( GtoLRotationMatrix == NULL ) {
        FloatArray e1(3), e2(3), e3(3), help(3);

        // compute e1' = [N2-N1]  and  help = [N3-N1]
        e1.beDifferenceOf( * this->giveNode(2)->giveCoordinates(),  * this->giveNode(1)->giveCoordinates() );
        help.beDifferenceOf( * this->giveNode(3)->giveCoordinates(),  * this->giveNode(1)->giveCoordinates() );

        // let us normalize e1'
        e1.normalize();

        // compute e3' : vector product of e1' x help
        e3.beVectorProductOf(e1, help);
        // let us normalize
        e3.normalize();

        // now from e3' x e1' compute e2'
        e2.beVectorProductOf(e3, e1);

        //
        GtoLRotationMatrix = new FloatMatrix(3, 3);

        for ( int i = 1; i <= 3; i++ ) {
            GtoLRotationMatrix->at(1, i) = e1.at(i);
            GtoLRotationMatrix->at(2, i) = e2.at(i);
            GtoLRotationMatrix->at(3, i) = e3.at(i);
        }
    }

    return GtoLRotationMatrix;
}


bool
TrPlanestressRotAllman3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [9,18]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v,r3} = T * {u,v,w,r1,r2,r3}
{
    // test if pereviously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize(9, 18);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = answer.at(1 + 3, i  + 6) = answer.at(1 + 6, i  + 12) = GtoLRotationMatrix->at(1, i);
        answer.at(2, i) = answer.at(2 + 3, i  + 6) = answer.at(2 + 6, i  + 12) = GtoLRotationMatrix->at(2, i);
        answer.at(3, i + 3) = answer.at(3 + 3, i + 3 + 6) = answer.at(3 + 6, i + 3 + 12) = GtoLRotationMatrix->at(3, i);
    }

    return 1;
}


void
TrPlanestressRotAllman3d :: giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep)
// returns characteristic tensor of the receiver at given gp and tStep
// strain vector = (Eps_X, Eps_y, Gamma_xy, Kappa_z)
{
    FloatArray charVect;
    StructuralMaterialStatus *ms = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );

    answer.resize(3, 3);
    answer.zero();

    if ( ( type == LocalForceTensor ) || ( type == GlobalForceTensor ) ) {
        //this->computeStressVector(charVect, gp, tStep);
        charVect = ms->giveStressVector();

        answer.at(1, 1) = charVect.at(1);
        answer.at(2, 2) = charVect.at(2);
        answer.at(1, 2) = charVect.at(3);
        answer.at(2, 1) = charVect.at(3);
    } else if ( ( type == LocalMomentumTensor ) || ( type == GlobalMomentumTensor ) ) {} else if ( ( type == LocalStrainTensor ) || ( type == GlobalStrainTensor ) ) {
        //this->computeStrainVector(charVect, gp, tStep);
        charVect = ms->giveStrainVector();

        answer.at(1, 1) = charVect.at(1);
        answer.at(2, 2) = charVect.at(2);
        answer.at(1, 2) = charVect.at(3) / 2.;
        answer.at(2, 1) = charVect.at(3) / 2.;
    } else if ( ( type == LocalCurvatureTensor ) || ( type == GlobalCurvatureTensor ) ) {} else {
        OOFEM_ERROR("unsupported tensor mode");
        exit(1);
    }

    if ( ( type == GlobalForceTensor  ) || ( type == GlobalMomentumTensor  ) ||
         ( type == GlobalStrainTensor ) || ( type == GlobalCurvatureTensor ) ) {
        this->computeGtoLRotationMatrix();
        answer.rotatedWith(* GtoLRotationMatrix);
    }
}


int
TrPlanestressRotAllman3d :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatMatrix globTensor;
    CharTensor cht;

    answer.resize(6);

    if ( type == IST_ShellForceTensor || type == IST_ShellStrainTensor ) {
        double c = 1.0;
        if ( type == IST_ShellForceTensor ) {
            cht = GlobalForceTensor;
        } else {
            cht = GlobalStrainTensor;
            c = 2.0; // tensor components reported
        }

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = c * globTensor.at(2, 3); //yz
        answer.at(5) = c * globTensor.at(1, 3); //xz
        answer.at(6) = c * globTensor.at(2, 3); //yz
        // mutiply stresses by thickness to get forces

        return 1;
    } else if ( ( type == IST_ShellMomentumTensor ) || ( type == IST_ShellCurvatureTensor ) ) {
        answer.clear();
        return 1;
    } else {
        answer.clear();
        return 0;
    }
}

int
TrPlanestressRotAllman3d :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [6,6]
// f(local) = T * f(global)
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize(6, 6);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = answer.at(4, i + 3) = GtoLRotationMatrix->at(1, i);
        answer.at(2, i) = answer.at(5, i + 3) = GtoLRotationMatrix->at(2, i);
        answer.at(3, i) = answer.at(6, i + 3) = GtoLRotationMatrix->at(3, i);
    }

    return 1;
}


void
TrPlanestressRotAllman3d :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV, load;
    GaussPoint *gp = NULL;
    FloatArray force;
    FloatMatrix T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    GaussIntegrationRule irule(0, this);
    irule.SetUpPointsOnTriangle( 1, this->giveMaterialMode() );

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, tStep, mode);

    if ( force.giveSize() ) {
        gp = irule.getIntegrationPoint(0);

        dens = this->giveStructuralCrossSection()->give('d', gp);
        dV   = this->computeVolumeAround(gp);

        answer.resize(18);
        answer.zero();

        load = force.at(1) * dens * dV / 3.0;
        answer.at(1) = load;
        answer.at(7) = load;
        answer.at(13) = load;

        load = force.at(2) * dens * dV / 3.0;
        answer.at(2) = load;
        answer.at(8) = load;
        answer.at(14) = load;

        load = force.at(3) * dens * dV / 3.0;
        answer.at(3) = load;
        answer.at(9) = load;
        answer.at(15) = load;

        // transform result from global cs to local element cs.
        if ( this->computeGtoLRotationMatrix(T) ) {
            answer.rotatedWith(T, 'n');
        }
    } else {
        answer.clear();          // nil resultant
    }
}

void
TrPlanestressRotAllman3d :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    FloatArray v;

    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    for ( int i = 0; i < ( int ) integrationRulesArray.size(); i++ ) {
        for ( GaussPoint *gp : *integrationRulesArray [ i ] ) {
            fprintf( file, "  GP %2d.%-2d :", i + 1, gp->giveNumber() );

            this->giveIPValue(v, gp, IST_ShellStrainTensor, tStep);
            fprintf(file, "  strains    ");
            for ( auto &val : v ) {
                fprintf(file, " %.4e", val);
            }

            this->giveIPValue(v, gp, IST_ShellCurvatureTensor, tStep);
            fprintf(file, "\n              curvatures ");
            for ( auto &val : v ) {
                fprintf(file, " %.4e", val);
            }

            // Forces - Moments
            this->giveIPValue(v, gp, IST_ShellForceTensor, tStep);
            fprintf(file, "\n              stresses   ");
            for ( auto &val : v ) {
                fprintf(file, " %.4e", val);
            }

            this->giveIPValue(v, gp, IST_ShellMomentumTensor, tStep);
            fprintf(file, "\n              moments    ");
            for ( auto &val : v ) {
                fprintf(file, " %.4e", val);
            }

            fprintf(file, "\n");
        }
    }
}
} // end namespace oofem
