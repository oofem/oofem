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

#include "../sm/Elements/Shells/cct3d.h"
#include "../sm/Materials/structuralms.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "load.h"
#include "mathfem.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "classfactory.h"

#include <cstdlib>

namespace oofem {
REGISTER_Element(CCTPlate3d);

CCTPlate3d :: CCTPlate3d(int n, Domain *aDomain) : CCTPlate(n, aDomain)
{
}


void
CCTPlate3d :: giveLocalCoordinates(FloatArray &answer, FloatArray &global)
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
    FloatArray offset;
    // test the parametr
    if ( global.giveSize() != 3 ) {
        OOFEM_ERROR("cannot transform coordinates - size mismatch");
        exit(1);
    }

    // first ensure that receiver's GtoLRotationMatrix[3,3] is defined
    this->computeGtoLRotationMatrix();

    offset = global;
    offset.subtract( * this->giveNode(1)->giveCoordinates() );
    answer.beProductOf(GtoLRotationMatrix, offset);
}


void
CCTPlate3d :: giveNodeCoordinates(double &x1, double &x2, double &x3,
                                  double &y1, double &y2, double &y3,
                                  double &z1, double &z2, double &z3)
{
    FloatArray nc1(3), nc2(3), nc3(3);

    this->giveLocalCoordinates( nc1, * ( this->giveNode(1)->giveCoordinates() ) );
    this->giveLocalCoordinates( nc2, * ( this->giveNode(2)->giveCoordinates() ) );
    this->giveLocalCoordinates( nc3, * ( this->giveNode(3)->giveCoordinates() ) );

    x1 = nc1.at(1);
    x2 = nc2.at(1);
    x3 = nc3.at(1);

    y1 = nc1.at(2);
    y2 = nc2.at(2);
    y3 = nc3.at(2);

    z1 = nc1.at(3);
    z2 = nc2.at(3);
    z3 = nc3.at(3);

}


void
CCTPlate3d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, D_w, R_u, R_v, R_w};
}


bool
CCTPlate3d :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
//converts global coordinates to local planar area coordinates,
//does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
{
    // rotate the input point Coordinate System into the element CS
    FloatArray inputCoords_ElCS;
    std::vector< FloatArray > lc(3);
    FloatArray llc;
    this->giveLocalCoordinates( inputCoords_ElCS, const_cast< FloatArray & >(coords) );
    for ( int _i = 0; _i < 3; _i++ ) {
        this->giveLocalCoordinates( lc [ _i ], * this->giveNode(_i + 1)->giveCoordinates() );
    }
    FEI2dTrLin _interp(1, 2);
    bool inplane = _interp.global2local(llc, inputCoords_ElCS, FEIVertexListGeometryWrapper(lc)) > 0;
    answer.resize(2);
    answer.at(1) = inputCoords_ElCS.at(1);
    answer.at(2) = inputCoords_ElCS.at(2);
    GaussPoint _gp(NULL, 1, answer, 2.0, _2dPlate);
    // now check if the third local coordinate is within the thickness of element
    bool outofplane = ( fabs( inputCoords_ElCS.at(3) ) <= this->giveCrossSection()->give(CS_Thickness, & _gp) / 2. );

    return inplane && outofplane;
}


int
CCTPlate3d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = 1. - l2 - l1;

    answer.resize(3);
    for ( int _i = 1; _i <= 3; _i++ ) {
        answer.at(_i) = l1 * this->giveNode(1)->giveCoordinate(_i) + l2 *this->giveNode(2)->giveCoordinate(_i) + l3 *this->giveNode(3)->giveCoordinate(_i);
    }
    return true;
}


const FloatMatrix *
CCTPlate3d :: computeGtoLRotationMatrix()
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
    if ( !GtoLRotationMatrix.isNotEmpty() ) {
        FloatArray e1, e2, e3, help;

        // compute e1' = [N2-N1]  and  help = [N3-N1]
        e1.beDifferenceOf(*this->giveNode(2)->giveCoordinates(), *this->giveNode(1)->giveCoordinates());
        help.beDifferenceOf(*this->giveNode(3)->giveCoordinates(), *this->giveNode(1)->giveCoordinates());

        // let us normalize e1'
        e1.normalize();

        // compute e3' : vector product of e1' x help
        e3.beVectorProductOf(e1, help);
        // let us normalize
        e3.normalize();

        // now from e3' x e1' compute e2'
        e2.beVectorProductOf(e3, e1);

        //
        GtoLRotationMatrix.resize(3, 3);

        for ( int i = 1; i <= 3; i++ ) {
            GtoLRotationMatrix.at(1, i) = e1.at(i);
            GtoLRotationMatrix.at(2, i) = e2.at(i);
            GtoLRotationMatrix.at(3, i) = e3.at(i);
        }
    }

    return &GtoLRotationMatrix;
}


bool
CCTPlate3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [9,18]
// r(local) = T * r(global)
// for one node (r written transposed): {w,r1,r2} = T * {u,v,w,r1,r2,r3}
{
    this->computeGtoLRotationMatrix();

    answer.resize(9, 18);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = answer.at(1 + 3, i  + 6) = answer.at(1 + 6, i  + 12) = GtoLRotationMatrix.at(3, i);
        answer.at(2, i + 3) = answer.at(2 + 3, i + 3 + 6) = answer.at(2 + 6, i + 3 + 12) = GtoLRotationMatrix.at(1, i);
        answer.at(3, i + 3) = answer.at(3 + 3, i + 3 + 6) = answer.at(3 + 6, i + 3 + 12) = GtoLRotationMatrix.at(2, i);
    }

    return 1;
}

void
CCTPlate3d :: giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep)
// returns characteristic tensor of the receiver at given gp and tStep
// strain vector = (Kappa_x, Kappa_y, Kappa_xy, Gamma_zx, Gamma_zy)
{
    FloatArray charVect;
    StructuralMaterialStatus *ms = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );

    answer.resize(3, 3);
    answer.zero();

    if ( ( type == LocalForceTensor ) || ( type == GlobalForceTensor ) ) {
        //this->computeStressVector(charVect, gp, tStep);
        charVect = ms->giveStressVector();

        answer.at(1, 3) = charVect.at(4);
        answer.at(3, 1) = charVect.at(4);
        answer.at(2, 3) = charVect.at(5);
        answer.at(3, 2) = charVect.at(5);
    } else if ( ( type == LocalMomentumTensor ) || ( type == GlobalMomentumTensor ) ) {
        //this->computeStressVector(charVect, gp, tStep);
        charVect = ms->giveStressVector();

        answer.at(1, 1) = charVect.at(1);
        answer.at(2, 2) = charVect.at(2);
        answer.at(1, 2) = charVect.at(3);
        answer.at(2, 1) = charVect.at(3);
    } else if ( ( type == LocalStrainTensor ) || ( type == GlobalStrainTensor ) ) {
        //this->computeStrainVector(charVect, gp, tStep);
        charVect = ms->giveStrainVector();

        answer.at(1, 3) = charVect.at(4) / 2.;
        answer.at(3, 1) = charVect.at(4) / 2.;
        answer.at(2, 3) = charVect.at(5) / 2.;
        answer.at(3, 2) = charVect.at(5) / 2.;
    } else if ( ( type == LocalCurvatureTensor ) || ( type == GlobalCurvatureTensor ) ) {
        //this->computeStrainVector(charVect, gp, tStep);
        charVect = ms->giveStrainVector();

        answer.at(1, 1) = charVect.at(1);
        answer.at(2, 2) = charVect.at(2);
        answer.at(1, 2) = charVect.at(3) / 2.;
        answer.at(2, 1) = charVect.at(3) / 2.;
    } else {
        OOFEM_ERROR("unsupported tensor mode");
        exit(1);
    }

    if ( ( type == GlobalForceTensor  ) || ( type == GlobalMomentumTensor  ) ||
        ( type == GlobalStrainTensor ) || ( type == GlobalCurvatureTensor ) ) {
        this->computeGtoLRotationMatrix();
        answer.rotatedWith(GtoLRotationMatrix);
    }
}


int
CCTPlate3d :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatMatrix globTensor;
    CharTensor cht;

    answer.resize(6);

    if (  type == IST_ShellCurvatureTensor || type == IST_ShellStrainTensor ) {
        if ( type == IST_ShellCurvatureTensor ) {
            cht = GlobalCurvatureTensor;
        } else {
            cht = GlobalStrainTensor;
        }

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = 2 * globTensor.at(2, 3); //yz
        answer.at(5) = 2 * globTensor.at(1, 3); //xz
        answer.at(6) = 2 * globTensor.at(2, 3); //yz

        return 1;
    } else if ( type == IST_ShellMomentumTensor || type == IST_ShellForceTensor ) {
        if ( type == IST_ShellMomentumTensor ) {
            cht = GlobalMomentumTensor;
        } else {
            cht = GlobalForceTensor;
        }

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = globTensor.at(2, 3); //yz
        answer.at(5) = globTensor.at(1, 3); //xz
        answer.at(6) = globTensor.at(2, 3); //yz

        return 1;
    } else {
        return NLStructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}

int
CCTPlate3d :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [6,6]
// f(local) = T * f(global)
{
    this->computeGtoLRotationMatrix();

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
CCTPlate3d :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV, load;
    FloatArray force;
    FloatMatrix T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    GaussIntegrationRule irule(1, this, 1, 5);
    irule.SetUpPointsOnTriangle(1, _2dPlate);

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, tStep, mode);

    if ( force.giveSize() ) {
        GaussPoint *gp = irule.getIntegrationPoint(0);

        dens = this->giveStructuralCrossSection()->give('d', gp); // constant density assumed
        dV   = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness, gp); // constant thickness assumed

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
CCTPlate3d :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    FloatMatrix ne;
    this->computeNmatrixAt(sgp->giveNaturalCoordinates(), ne);

    answer.resize(6, 18);
    answer.zero();
    int ri[] = {
        2, 3, 4
    };
    int ci[] = {
        2, 3, 4, 8, 9, 10, 14, 15, 16
    };

    for ( int i = 0; i < 3; i++ ) {
        for ( int j = 0; j < 9; j++ ) {
            answer(ri [ i ], ci [ j ]) = ne(i, j);
        }
    }
}

void
CCTPlate3d :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    answer.resize(18);
    answer.zero();
    if ( iSurf == 1 ) {
        answer.at(3) = 1; // node 1
        answer.at(4) = 2;
        answer.at(5) = 3;

        answer.at(9)  = 4; // node 2
        answer.at(10) = 5;
        answer.at(11) = 6;

        answer.at(15) = 7; // node 3
        answer.at(16) = 8;
        answer.at(17) = 9;
    } else {
        OOFEM_ERROR("wrong surface number");
    }
}

IntegrationRule *
CCTPlate3d :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, approxOrder);
    iRule->SetUpPointsOnTriangle(npoints, _Unknown);
    return iRule;
}

double
CCTPlate3d :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    return this->computeVolumeAround(gp);
}


void
CCTPlate3d :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int isurf)
{
    this->computeGlobalCoordinates( answer, gp->giveNaturalCoordinates() );
}


int
CCTPlate3d :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int isurf, GaussPoint *gp)
{
    return 0;
}

void
CCTPlate3d :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    FloatArray v;

    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    for ( int i = 0; i < (int)integrationRulesArray.size(); i++ ) {
        for ( GaussPoint *gp: *integrationRulesArray [ i ] ) {

            fprintf( file, "  GP %2d.%-2d :", i + 1, gp->giveNumber() );

            this->giveIPValue(v, gp, IST_ShellStrainTensor, tStep);
            fprintf(file, "  strains    ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);

            this->giveIPValue(v, gp, IST_ShellCurvatureTensor, tStep);
            fprintf(file, "\n              curvatures ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);

            // Forces - Moments
            this->giveIPValue(v, gp, IST_ShellForceTensor, tStep);
            fprintf(file, "\n              stresses   ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);

            this->giveIPValue(v, gp, IST_ShellMomentumTensor, tStep);
            fprintf(file, "\n              moments    ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);

            fprintf(file, "\n");
        }
    }
}
} // end namespace oofem
