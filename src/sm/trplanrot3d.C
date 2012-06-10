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

#include "trplanrot3d.h"
#include "material.h"
#include "crosssection.h"
#include "node.h"
#include "load.h"
#include "structuralms.h"
#include "mathfem.h"

namespace oofem {
TrPlaneStrRot3d :: TrPlaneStrRot3d(int n, Domain *aDomain) : TrPlaneStrRot(n, aDomain)
{
    GtoLRotationMatrix = NULL;
}


void
TrPlaneStrRot3d :: giveLocalCoordinates(FloatArray &answer, FloatArray &global)
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
    // test the parameter
    if ( global.giveSize() != 3 ) {
        _error("giveLocalCoordinate : cannot transform coordinates - size mismatch");
        exit(1);
    }

    // first ensure that receiver's GtoLRotationMatrix[3,3] is defined
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    answer.beProductOf(* GtoLRotationMatrix, global);
}


void
TrPlaneStrRot3d :: giveNodeCoordinates(FloatArray &x, FloatArray &y)
{
    FloatArray nc1(3), nc2(3), nc3(3);

    this->giveLocalCoordinates( nc1, * ( this->giveNode(1)->giveCoordinates() ) );
    this->giveLocalCoordinates( nc2, * ( this->giveNode(2)->giveCoordinates() ) );
    this->giveLocalCoordinates( nc3, * ( this->giveNode(3)->giveCoordinates() ) );

    x.at(1) = nc1.at(1);
    x.at(2) = nc2.at(1);
    x.at(3) = nc3.at(1);

    y.at(1) = nc1.at(2);
    y.at(2) = nc2.at(2);
    y.at(3) = nc3.at(2);

    //if (z) {
    //  z[0] = nc1->at(3);
    //  z[1] = nc2->at(3);
    //  z[2] = nc3->at(3);
    //}
}


void
TrPlaneStrRot3d :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    answer.setValues(6, D_u, D_v, D_w, R_u, R_v, R_w);
}


const FloatMatrix *
TrPlaneStrRot3d :: computeGtoLRotationMatrix()
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
        int i;
        FloatArray e1(3), e2(3), e3(3), help(3);

        // compute e1' = [N2-N1]  and  help = [N3-N1]
        for ( i = 1; i <= 3; i++ ) {
            e1.at(i) = ( this->giveNode(2)->giveCoordinate(i) - this->giveNode(1)->giveCoordinate(i) );
            help.at(i) = ( this->giveNode(3)->giveCoordinate(i) - this->giveNode(1)->giveCoordinate(i) );
        }

        // let us normalize e1'
        e1.normalize();

        // compute e3' : vector product of e1' x help
        e3.beVectorProductOf(e1,help);
        // let us normalize
        e3.normalize();

        // now from e3' x e1' compute e2'
        e2.beVectorProductOf(e3,e1);

        //
        GtoLRotationMatrix = new FloatMatrix(3, 3);

        for ( i = 1; i <= 3; i++ ) {
            GtoLRotationMatrix->at(1, i) = e1.at(i);
            GtoLRotationMatrix->at(2, i) = e2.at(i);
            GtoLRotationMatrix->at(3, i) = e3.at(i);
        }
    }

    return GtoLRotationMatrix;
}


bool
TrPlaneStrRot3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
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
TrPlaneStrRot3d :: giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep)
// returns characteristic tensor of the receiver at given gp and tStep
// strain vector = (Eps_X, Eps_y, Gamma_xy, Kappa_z)
{
    FloatArray charVect;
    Material *mat = this->giveMaterial();

    answer.resize(3, 3);
    answer.zero();

    if ( ( type == LocalForceTensor ) || ( type == GlobalForceTensor ) ) {
        //this->computeStressVector(charVect, gp, tStep);
        charVect = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )->giveStressVector();

        answer.at(1, 1) = charVect.at(1);
        answer.at(2, 2) = charVect.at(2);
        answer.at(1, 2) = charVect.at(3);
        answer.at(2, 1) = charVect.at(3);
    } else if ( ( type == LocalMomentumTensor ) || ( type == GlobalMomentumTensor ) ) {
        //this->computeStressVector(charVect, gp, tStep);
        charVect = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )->giveStressVector();

        answer.at(3, 3) = charVect.at(4);
    } else if ( ( type == LocalStrainTensor ) || ( type == GlobalStrainTensor ) ) {
        //this->computeStrainVector(charVect, gp, tStep);
        charVect = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )->giveStrainVector();

        answer.at(1, 1) = charVect.at(1);
        answer.at(2, 2) = charVect.at(2);
        answer.at(1, 2) = charVect.at(3) / 2.;
        answer.at(2, 1) = charVect.at(3) / 2.;
    } else if ( ( type == LocalCurvatureTensor ) || ( type == GlobalCurvatureTensor ) ) {
        //this->computeStrainVector(charVect, gp, tStep);
        charVect = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )->giveStrainVector();

        answer.at(3, 3) = charVect.at(4);
    } else {
        _error("GiveCharacteristicTensor: unsupported tensor mode");
        exit(1);
    }

    if ( ( type == GlobalForceTensor  ) || ( type == GlobalMomentumTensor  ) ||
        ( type == GlobalStrainTensor ) || ( type == GlobalCurvatureTensor ) ) {
        this->computeGtoLRotationMatrix();
        answer.rotatedWith(* GtoLRotationMatrix);
    }
}


int
TrPlaneStrRot3d :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *tStep)
{
    FloatMatrix globTensor;
    CharTensor cht;

    answer.resize(12);

    if ( ( type == IST_ShellForceMomentumTensor ) || ( type == IST_ShellStrainCurvatureTensor ) ) {
        if ( type == IST_ShellForceMomentumTensor ) {
            cht = GlobalForceTensor;
        } else {
            cht = GlobalStrainTensor;
        }

        this->giveCharacteristicTensor(globTensor, cht, aGaussPoint, tStep);

        answer.at(1) = globTensor.at(1, 1); //sxForce
        answer.at(2) = globTensor.at(2, 2); //syForce
        answer.at(3) = globTensor.at(3, 3); //szForce
        answer.at(4) = globTensor.at(2, 3); //syzForce
        answer.at(5) = globTensor.at(1, 3); //qxzForce
        answer.at(6) = globTensor.at(1, 2); //qxyForce

        if ( type == IST_ShellForceMomentumTensor ) {
            cht = GlobalMomentumTensor;
        } else {
            cht = GlobalCurvatureTensor;
        }

        this->giveCharacteristicTensor(globTensor, cht, aGaussPoint, tStep);

        answer.at(7)  = globTensor.at(1, 1); //mxForce
        answer.at(8)  = globTensor.at(2, 2); //myForce
        answer.at(9)  = globTensor.at(3, 3); //mzForce
        answer.at(10) = globTensor.at(2, 3); //myzForce
        answer.at(11) = globTensor.at(1, 3); //mxzForce
        answer.at(12) = globTensor.at(1, 2); //mxyForce

        return 1;
    } else {
        answer.resize(0);
        return 0;
    }
}


int
TrPlaneStrRot3d :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
    if ( ( type == IST_ShellForceMomentumTensor ) || ( type == IST_ShellStrainCurvatureTensor ) ) {
        answer.resize(12);
        for ( int i = 1; i <= 12; i++ ) {
            answer.at(i) = i;
        }

        return 1;
    } else {
        return TrPlaneStrRot :: giveIntVarCompFullIndx(answer, type);
    }
}


void
TrPlaneStrRot3d :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body loads, at stepN.
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
        _error("computeBodyLoadVectorAt: unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, stepN, mode);

    if ( force.giveSize() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        dens = this->giveMaterial()->give('d', gp);
        dV   = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness);

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
        answer.resize(0);          // nil resultant
    }
}


void
TrPlaneStrRot3d :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    int i, j;
    GaussPoint *gp;
    FloatArray v;

#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        for ( j = 0; j < integrationRulesArray [ i ]->getNumberOfIntegrationPoints(); j++ ) {
            gp = integrationRulesArray [ i ]->getIntegrationPoint(j);

            // gp   -> printOutputAt(file,stepN) ;

            fprintf( file, "  GP %2d.%-2d :", i + 1, gp->giveNumber() );

            this->giveIPValue(v, gp, IST_ShellStrainCurvatureTensor, tStep);
            fprintf(file, "  strains ");
            fprintf( file,
                    " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
                    v.at(1), v.at(2), v.at(3),  2. * v.at(4), 2. * v.at(5), 2. * v.at(6),
                    v.at(7), v.at(8), v.at(9),  2. * v.at(10), 2. * v.at(11), 2. * v.at(12) );

            this->giveIPValue(v, gp, IST_ShellForceMomentumTensor, tStep);
            fprintf(file, "\n              stresses");
            fprintf( file,
                    " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
                    v.at(1), v.at(2), v.at(3),  v.at(4), v.at(5), v.at(6),
                    v.at(7), v.at(8), v.at(9),  v.at(10), v.at(11), v.at(12) );

            fprintf(file, "\n");
        }
    }
}


GaussPoint *
TrPlaneStrRot3d :: giveMiddleGaussPoint()
// return gausspoint with local coordinates (0.333,0.333,0.333)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(3);

    if ( ( fabs(gp->giveCoordinate(1) - 0.333333333333) > 0.00000000001 ) ||
        ( fabs(gp->giveCoordinate(2) - 0.333333333333) > 0.00000000001 ) ) {
        _error("the gausspoint has not local coordinate (0.333,0.333,0.333)");
    }

    return gp;
}
} // end namespace oofem
