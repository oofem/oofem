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

#include "../sm/Elements/PlaneStress/linquad3d_planestress.h"
#include "classfactory.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "gausspoint.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "connectivitytable.h"
 #include "Materials/rcm2.h"
#endif

namespace oofem {
REGISTER_Element(LinQuad3DPlaneStress);

LinQuad3DPlaneStress :: LinQuad3DPlaneStress(int n, Domain *aDomain) :
    PlaneStress2d(n, aDomain)
    // Constructor.
{
    this->GtoLRotationMatrix = NULL;
}

LinQuad3DPlaneStress :: ~LinQuad3DPlaneStress()
// Destructor
{
    delete this->GtoLRotationMatrix;
}

Interface *
LinQuad3DPlaneStress :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    }

    return NULL;
}


FEICellGeometry* 
LinQuad3DPlaneStress::giveCellGeometryWrapper()
{
    if (cellGeometryWrapper) {
        return cellGeometryWrapper;
    } else {
        this->computeLocalNodalCoordinates(lc);
        return (cellGeometryWrapper = new FEIVertexListGeometryWrapper(lc));
    }
}


void
LinQuad3DPlaneStress :: computeLocalNodalCoordinates(std::vector< FloatArray > &lxy)
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }


    lxy.resize(4);
    const FloatArray *nc;
    for ( int i = 0; i < 4; i++ ) {
        nc = this->giveNode(i + 1)->giveCoordinates();
        lxy[i].beProductOf(* GtoLRotationMatrix, *nc);
    }
}

void
LinQuad3DPlaneStress :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, D_w};
}




const FloatMatrix *
LinQuad3DPlaneStress :: computeGtoLRotationMatrix()
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
        FloatArray e1, e2, e3, help;

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
LinQuad3DPlaneStress :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [8,12]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v} = T * {u,v,w}
{
    // test if pereviously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize(8, 12);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
      answer.at(1, i) = answer.at(3, i  + 3) = answer.at(5, i  + 6) = answer.at(7, i+9) = GtoLRotationMatrix->at(1, i);
      answer.at(2, i) = answer.at(4, i  + 3) = answer.at(6, i  + 6) = answer.at(8, i+9) = GtoLRotationMatrix->at(2, i);
    }

    return 1;
}
int
LinQuad3DPlaneStress :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
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
LinQuad3DPlaneStress :: giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep)
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
    } else if ( ( type == LocalMomentumTensor ) || ( type == GlobalMomentumTensor ) ) {
    } else if ( ( type == LocalStrainTensor ) || ( type == GlobalStrainTensor ) ) {
        //this->computeStrainVector(charVect, gp, tStep);
        charVect = ms->giveStrainVector();

        answer.at(1, 1) = charVect.at(1);
        answer.at(2, 2) = charVect.at(2);
        answer.at(1, 2) = charVect.at(3) / 2.;
        answer.at(2, 1) = charVect.at(3) / 2.;
    } else if ( ( type == LocalCurvatureTensor ) || ( type == GlobalCurvatureTensor ) ) {
    } else {
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
LinQuad3DPlaneStress :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
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
            c = 2.0;
        }

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = c * globTensor.at(2, 3); //yz
        answer.at(5) = c * globTensor.at(1, 3); //xz
        answer.at(6) = c * globTensor.at(2, 3); //yz

        return 1;
    } else if ( type == IST_ShellMomentumTensor || type == IST_ShellCurvatureTensor ) {
        answer.clear();
        return 1;
    } else {
        answer.clear();
        return 0;
    }
}


void
LinQuad3DPlaneStress :: printOutputAt(FILE *file, TimeStep *tStep)
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

            /*
            this->giveIPValue(v, gp, IST_ShellCurvatureTensor, tStep);
            fprintf(file, "\n              curvatures ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);
            */
            // Forces - Moments
            this->giveIPValue(v, gp, IST_ShellForceTensor, tStep);
            fprintf(file, "\n              forces     ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);
            /*
            this->giveIPValue(v, gp, IST_ShellMomentumTensor, tStep);
            fprintf(file, "\n              moments    ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);
            */
            fprintf(file, "\n");
        }
    }
}


} // end namespace oofem
