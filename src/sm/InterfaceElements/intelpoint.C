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

#include "intelpoint.h"
#include "structuralinterfacecrosssection.h"

#include "domain.h"
#include "node.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"



namespace oofem {
REGISTER_Element(IntElPoint);

IntElPoint :: IntElPoint(int n, Domain *aDomain) :
    StructuralInterfaceElement(n, aDomain)
{
    numberOfDofMans = 2;
    referenceNode = 0;
    normal.resize(3);
    normal.zero();
}


void
IntElPoint :: setCoordMode()
{
    switch ( domain->giveNumberOfSpatialDimensions() ) {
    case 1:
        this->mode = ie1d_1d;
        break;
    case 2:
        this->mode = ie1d_2d;
        break;
    case 3:
        this->mode = ie1d_3d;
        break;
    default:
        OOFEM_ERROR("Unsupported domain type")
    }
}


MaterialMode
IntElPoint :: giveMaterialMode()
{
    setCoordMode();
    switch ( mode ) {
    case ie1d_1d: return _1dInterface;

    case ie1d_2d: return _2dInterface;

    case ie1d_3d: return _3dInterface;

    default: OOFEM_ERROR("Unsupported coord mode");
    }
    return _1dInterface; // to make the compiler happy
}






void
IntElPoint :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    setCoordMode();
    

    FloatArray N;
    switch ( mode ) {
    case ie1d_1d:
        answer.resize(1, 2);
        answer.at(1, 1) = -1.0;
        answer.at(1, 2) = +1.0;
        break;
    case ie1d_2d:
        answer.resize(2, 4);
        answer.zero();
        answer.at(1, 1) = -1.0;
        answer.at(1, 3) = +1.0;
        answer.at(2, 2) = -1.0;
        answer.at(2, 4) = +1.0;
        break;
    case ie1d_3d:
        answer.resize(3, 6);
        answer.zero();
        answer.at(1, 1) = -1.0;
        answer.at(1, 4) = +1.0;
        answer.at(2, 2) = -1.0;
        answer.at(2, 5) = +1.0;
        answer.at(3, 3) = -1.0;
        answer.at(3, 6) = +1.0;
        break;
    default:
        OOFEM_ERROR("Unsupported mode");
    }

}


void
IntElPoint :: computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Computes transformation matrix to local coordinate system.
    setCoordMode();
    switch ( mode ) {
    case ie1d_1d:
        answer.resize(1, 1);
        answer.at(1, 1) = 1.;
        return;

    case ie1d_2d:
        answer.resize(2, 2);
        answer.at(1, 1) =  normal.at(1);
        answer.at(1, 2) =  normal.at(2);
        answer.at(2, 1) = -normal.at(2);
        answer.at(2, 2) =  normal.at(1);
        return;

    case ie1d_3d:
    {
        //FloatMatrix test;
        //test.beLocalCoordSys(normal.normalize());

        FloatArray ly(3), lz(3);
        normal.normalize();
        ly.zero();
        if ( fabs( normal.at(1) ) > fabs( normal.at(2) ) ) {
            ly.at(2) = 1.;
        } else {
            ly.at(1) = 1.;
        }

        lz.beVectorProductOf(normal, ly);
        lz.normalize();
        ly.beVectorProductOf(lz, normal);
        ly.normalize();

        answer.resize(3, 3);
        int i;
        for ( i = 1; i <= 3; i++ ) {
            answer.at(1, i) = normal.at(i);
            answer.at(2, i) = ly.at(i);
            answer.at(3, i) = lz.at(i);
        }

        return;
    }

    default:
        OOFEM_ERROR("Unsupported mode");
    }
}


void
IntElPoint :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        int numberOfIntegrationRules = 1;
        integrationRulesArray.resize(numberOfIntegrationRules);
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints( _Line, 1, this->giveMaterialMode() );
    }
}


int
IntElPoint :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    answer.resize(3);
    answer.at(1) = 0.5* ( this->giveNode(1)->giveCoordinate(1) + this->giveNode(2)->giveCoordinate(1) );
    answer.at(2) = 0.5* ( this->giveNode(1)->giveCoordinate(2) + this->giveNode(2)->giveCoordinate(2) );
    answer.at(3) = 0.5* ( this->giveNode(1)->giveCoordinate(3) + this->giveNode(2)->giveCoordinate(3) );

    return 1;
}



double
IntElPoint :: computeAreaAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
        //double area = this->giveCrossSection()->give(CS_Area);
    return 1.0;  ///@todo read input from cs
}


IRResultType
IntElPoint :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralInterfaceElement :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, referenceNode, _IFT_IntElPoint_refnode);
    IR_GIVE_OPTIONAL_FIELD(ir, normal, _IFT_IntElPoint_normal);
    if ( referenceNode == 0 && normal.at(1) == 0 && normal.at(2) == 0 && normal.at(1) == 0 && normal.at(3) == 0 ) {
        OOFEM_ERROR("Wrong reference node or normal specified");
    }

    this->computeLocalSlipDir(normal); ///@todo Move into postInitialize ?
    return IRRT_OK;
}


int
IntElPoint :: computeNumberOfDofs()
{
    setCoordMode();
    switch ( mode ) {
    case ie1d_1d:
        return 2;

    case ie1d_2d:
        return 4;

    case ie1d_3d:
        return 6;

    default:
        OOFEM_ERROR("Unsupported mode");
    }

    return 0; // to suppress compiler warning
}


void
IntElPoint :: giveDofManDofIDMask(int inode, IntArray &answer) const
{

    switch ( domain->giveNumberOfSpatialDimensions() ) {
    case 1:
        answer = IntArray({ D_u });
        break;
    case 2:
        answer = { D_u, D_v };
        break;
    case 3:
        answer = { D_u, D_v, D_w };
        break;
    default:
        OOFEM_ERROR("Unsupported mode");
    }

}


void
IntElPoint :: computeLocalSlipDir(FloatArray &normal)
{
    normal.resizeWithValues(3);
    if ( this->referenceNode ) {
        // normal
        normal.at(1) = domain->giveNode(this->referenceNode)->giveCoordinate(1) - this->giveNode(1)->giveCoordinate(1);
        normal.at(2) = domain->giveNode(this->referenceNode)->giveCoordinate(2) - this->giveNode(1)->giveCoordinate(2);
        normal.at(3) = domain->giveNode(this->referenceNode)->giveCoordinate(3) - this->giveNode(1)->giveCoordinate(3);
    } else {
        if ( normal.at(1) == 0 && normal.at(2) == 0 && normal.at(3) == 0 ) {
            OOFEM_ERROR("Normal is not defined (referenceNode=0,normal=(0,0,0))");
        }
    }

    normal.normalize();
}

} // end namespace oofem
