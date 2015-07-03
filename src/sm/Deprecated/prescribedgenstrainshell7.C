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

#include "prescribedgenstrainshell7.h"
#include "dofiditem.h"
#include "dofmanager.h"
#include "dof.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "function.h"
#include "engngm.h"
#include "set.h"
#include "node.h"
#include "element.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "feinterpol.h"
#include "sparsemtrx.h"
#include "sparselinsystemnm.h"
#include "Elements/Shells/shell7base.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGenStrainShell7);

double PrescribedGenStrainShell7 :: give(Dof *dof, ValueModeType mode, double time)
{
    DofIDItem id = dof->giveDofID();
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();

    if ( coords->giveSize() != this->centerCoord.giveSize() ) {
        OOFEM_ERROR("PrescribedGenStrainShell7 :: give - Size of coordinate system different from center coordinate in b.c.");
    }

    double factor = 0;
    if ( mode == VM_Total ) {
        factor = this->giveTimeFunction()->evaluateAtTime(time);
    } else if ( mode == VM_Velocity ) {
        factor = this->giveTimeFunction()->evaluateVelocityAtTime(time);
    } else if ( mode == VM_Acceleration ) {
        factor = this->giveTimeFunction()->evaluateAccelerationAtTime(time);
    } else {
        OOFEM_ERROR("Should not be called for value mode type then total, velocity, or acceleration.");
    }

    // Reminder: u_i = F_ij . (x_j - xb_j) = H_ij . dx_j
    FloatArray dx;
    dx.beDifferenceOf(* coords, this->centerCoord);

    // Assuming the coordinate system to be local, dx(3) = z
    this->setDeformationGradient( dx.at(3) );

    FloatArray u, u2;
    u.beProductOf(gradient, dx);

    // Add second order contribution, note only higher order in the thickness direction
    this->evaluateHigherOrderContribution(u2, dx.at(3), dx);
    u.add(u2);

    u.times( factor );

    switch ( id ) {
    case D_u:
    case V_u:
        return u.at(1);

    case D_v:
    case V_v:
        return u.at(2);

    case D_w:
    case V_w:
        return u.at(3);

    default:
        return 0.0;
    }
}



void
PrescribedGenStrainShell7 :: evalCovarBaseVectorsAt(FloatMatrix &gcov, FloatArray &genEps, double zeta)
{
    // Evaluates the covariant base vectors in the current configuration
    FloatArray g1; FloatArray g2; FloatArray g3;

    FloatArray dxdxi1, dxdxi2, m, dmdxi1, dmdxi2;
    double dgamdxi1, dgamdxi2, gam;
    Shell7Base :: giveGeneralizedStrainComponents(genEps, dxdxi1, dxdxi2, dmdxi1, dmdxi2, m, dgamdxi1, dgamdxi2, gam);
    double fac1 = ( zeta + 0.5 * gam * zeta * zeta );
    double fac2 = ( 0.5 * zeta * zeta );
    double fac3 = ( 1.0 + zeta * gam );

    g1 = dxdxi1 + fac1*dmdxi1 + fac2*dgamdxi1*m;
    g2 = dxdxi2 + fac1*dmdxi2 + fac2*dgamdxi2*m;
    g3 = fac3*m;
    gcov.resize(3,3);
    gcov.setColumn(g1,1); gcov.setColumn(g2,2); gcov.setColumn(g3,3);
}


void
PrescribedGenStrainShell7 :: evalInitialCovarBaseVectorsAt(FloatMatrix &Gcov, FloatArray &genEps,  double zeta)
{
    // Evaluates the initial base vectors given the array of generalized strain
    FloatArray G1(3), G2(3), G3(3); 
    
    G1.at(1) = genEps.at(1) + zeta * genEps.at(7);
    G1.at(2) = genEps.at(2) + zeta * genEps.at(8);
    G1.at(3) = genEps.at(3) + zeta * genEps.at(9);

    G2.at(1) = genEps.at(4) + zeta * genEps.at(10);
    G2.at(2) = genEps.at(5) + zeta * genEps.at(11);
    G2.at(3) = genEps.at(6) + zeta * genEps.at(12);

    G3.at(1) = genEps.at(13);
    G3.at(2) = genEps.at(14);
    G3.at(3) = genEps.at(15);


    Gcov.resize(3,3);
    Gcov.setColumn(G1,1); Gcov.setColumn(G2,2); Gcov.setColumn(G3,3);
}


void
PrescribedGenStrainShell7 :: setDeformationGradient(double zeta)
{
    // Computes the deformation gradient in matrix form as open product(g_i, G^i) = gcov*Gcon^T
    FloatMatrix gcov, Gcon, Gcov;

    this->evalCovarBaseVectorsAt(gcov, this->genEps, zeta);
    this->evalInitialCovarBaseVectorsAt(Gcov, this->initialGenEps, zeta);
    Shell7Base :: giveDualBase(Gcov, Gcon);

    this->gradient.beProductTOf(gcov, Gcon);
    this->gradient.at(1,1) -= 1.0;
    this->gradient.at(2,2) -= 1.0;
    this->gradient.at(3,3) -= 1.0;
}

void
PrescribedGenStrainShell7 :: evaluateHigherOrderContribution(FloatArray &answer, double zeta, FloatArray &dx)
{
    // Computes the higher order contribution from the second gradient F2_ijk = (g3,3)_i * (G^3)_j * (G^3)_k 
    // Simplified version with only contribtion in the xi-direction
    FloatMatrix gcov, Gcon, Gcov;

    this->evalCovarBaseVectorsAt(gcov, this->genEps, zeta);
    this->evalInitialCovarBaseVectorsAt(Gcov, this->initialGenEps, zeta);
    Shell7Base :: giveDualBase(Gcov, Gcon);

    FloatArray G3(3), g3prime(3), m(3);
    G3.at(1) = Gcon.at(1,3);
    G3.at(2) = Gcon.at(2,3);
    G3.at(3) = Gcon.at(3,3);

    double factor = G3.dotProduct(dx);
    double gamma = this->genEps.at(18);
    m.at(1) = this->genEps.at(13);
    m.at(2) = this->genEps.at(14);
    m.at(3) = this->genEps.at(15);
    g3prime = gamma*m;
    
    answer = 0.5*factor*factor * g3prime;


}

#if 0
void PrescribedGenStrainShell7 :: updateCoefficientMatrix(FloatMatrix &C)
// This is written in a very general way, supporting both fm and sm problems.
// v_prescribed = C.d = (x-xbar).d;
// C = [x 0 y]
//     [0 y x]
//     [ ... ] in 2D, voigt form [d_11, d_22, d_12]
// C = [x 0 0 y z 0]
//     [0 y 0 x 0 z]
//     [0 0 z 0 x y]
//     [ ......... ] in 3D, voigt form [d_11, d_22, d_33, d_23, d_13, d_12]
{
    Domain *domain = this->giveDomain();
    int nNodes = domain->giveNumberOfDofManagers();

    int nsd = domain->giveNumberOfSpatialDimensions();
    int npeq = domain->giveEngngModel()->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
    C.resize(npeq, nsd * ( nsd + 1 ) / 2);
    C.zero();

    FloatArray &cCoords = this->giveCenterCoordinate();
    double xbar = cCoords.at(1), ybar = cCoords.at(2), zbar = 0.0;
    if ( nsd == 3 ) {
        zbar = cCoords.at(3);
    }

    for ( int i = 1; i <= nNodes; i++ ) {
        Node *n = domain->giveNode(i);
        FloatArray *coords = n->giveCoordinates();
        Dof *d1 = n->giveDofWithID( this->dofs(0) );
        Dof *d2 = n->giveDofWithID( this->dofs(1) );
        int k1 = d1->__givePrescribedEquationNumber();
        int k2 = d2->__givePrescribedEquationNumber();
        if ( nsd == 2 ) {
            if ( k1 ) {
                C.at(k1, 1) = coords->at(1) - xbar;
                C.at(k1, 3) = coords->at(2) - ybar;
            }

            if ( k2 ) {
                C.at(k2, 2) = coords->at(2) - ybar;
                C.at(k2, 3) = coords->at(1) - xbar;
            }
        } else { // nsd == 3
            OOFEM_ERROR("PrescribedGenStrainShell7 :: updateCoefficientMatrix - 3D Not tested yet!");
            Dof *d3 = n->giveDofWithID( this->dofs(2) );
            int k3 = d3->__givePrescribedEquationNumber();

            if ( k1 ) {
                C.at(k1, 1) = coords->at(1) - xbar;
                C.at(k1, 4) = coords->at(2) - ybar;
                C.at(k1, 5) = coords->at(3) - zbar;
            }
            if ( k2 ) {
                C.at(k2, 2) = coords->at(2) - ybar;
                C.at(k2, 4) = coords->at(1) - xbar;
                C.at(k2, 6) = coords->at(3) - zbar;
            }
            if ( k3 ) {
                C.at(k3, 3) = coords->at(3) - zbar;
                C.at(k3, 5) = coords->at(1) - xbar;
                C.at(k3, 6) = coords->at(2) - ybar;
            }
        }
    }
}
#endif

double PrescribedGenStrainShell7 :: domainSize()
{
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
    }
    return domain_size / nsd;
}




IRResultType PrescribedGenStrainShell7 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->initialGenEps, _IFT_PrescribedGenStrainShell7_initialgeneralizedstrain);
    IR_GIVE_FIELD(ir, this->genEps, _IFT_PrescribedGenStrainShell7_generalizedstrain);

    this->centerCoord.resize( this->gradient.giveNumberOfColumns() );
    this->centerCoord.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, this->centerCoord, _IFT_PrescribedGenStrainShell7_centercoords)

    return GeneralBoundaryCondition :: initializeFrom(ir);
}


void PrescribedGenStrainShell7 :: giveInputRecord(DynamicInputRecord &input)
{
    BoundaryCondition :: giveInputRecord(input);
    input.setField(this->initialGenEps, _IFT_PrescribedGenStrainShell7_initialgeneralizedstrain);
    input.setField(this->genEps, _IFT_PrescribedGenStrainShell7_generalizedstrain);
    input.setField(this->centerCoord, _IFT_PrescribedGenStrainShell7_centercoords);
}
} // end namespace oofem
