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

#include "../sm/Elements/Perleta/pplate.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei2dquadlin.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(PPlate);

FEI2dQuadLin PPlate :: interp_lin(1, 2);

PPlate :: PPlate(int n, Domain *aDomain) : StructuralElement(n, aDomain),


    lnodes(4)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = 4;
    E=0;
    t=0;
    ni=0;
    a=0;
    b=0;
}


FEInterpolation *
PPlate :: giveInterpolation() const { return & interp_lin; }


FEInterpolation *
PPlate :: giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}


void
PPlate :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 5) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
PPlate :: crossSectioncharacteristics()
        //Computes E,ni,t
{

            if( E == 0. ) {
                GaussPoint *gp;
                gp= integrationRulesArray[0]->getIntegrationPoint(0);
                StructuralCrossSection * csec=this->giveStructuralCrossSection();
                StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >(csec->giveMaterial(gp));
                t = csec->give(CS_Thickness, gp);
                E=mat->give('E', gp);
                ni=mat->give('n', gp);

            }
           ;
}

void
PPlate :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [5x12] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    answer.resize(5, 12);
    answer.zero();



}

void
PPlate::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    FloatMatrix Ka,Kb,Kc,Kd,K;
    double D,d;

    this->crossSectioncharacteristics();
    this->computeLength();

    K.resize(12,3);
    K.zero();

    this->computeKa(Ka);
    this->computeKb(Kb);
    this->computeKc(Kc);
    this->computeKd(Kd);

    K.add(Ka);
    K.add(Kb);
    K.add(Kc);
    K.add(Kd);
    d=(E*t*t*t)/(12*(1-ni*ni));
    D=d/(a*b);

    answer.resize(12,12);
    answer.zero();

    answer.at(1,1)= D*K.at(1,1);

    answer.at(2,1)= D*K.at(2,1);
    answer.at(2,2)= D*K.at(2,2);

    answer.at(3,1)= D*K.at(3,1);
    answer.at(3,2)= D*K.at(3,2);
    answer.at(3,3)= D*K.at(3,3);

    answer.at(4,1)= D*K.at(4,1);
    answer.at(4,2)= D*K.at(4,2);
    answer.at(4,3)= D*K.at(4,3);

    answer.at(5,1)= D*K.at(5,1);
    answer.at(5,2)= D*K.at(5,2);
    answer.at(5,3)= D*K.at(5,3);

    answer.at(6,1)= D*K.at(6,1);
    answer.at(6,2)= D*K.at(6,2);
    answer.at(6,3)= D*K.at(6,3);

    answer.at(7,1)= D*K.at(7,1);
    answer.at(7,2)= D*K.at(7,2);
    answer.at(7,3)= D*K.at(7,3);

    answer.at(8,1)= D*K.at(8,1);
    answer.at(8,2)= D*K.at(8,2);
    answer.at(8,3)= D*K.at(8,3);

    answer.at(9,1)= D*K.at(9,1);
    answer.at(9,2)= D*K.at(9,2);
    answer.at(9,3)= D*K.at(9,3);

    answer.at(10,1)= D*K.at(10,1);
    answer.at(10,2)= D*K.at(10,2);
    answer.at(10,3)= D*K.at(10,3);

    answer.at(11,1)= D*K.at(11,1);
    answer.at(11,2)= D*K.at(11,2);
    answer.at(11,3)= D*K.at(11,3);

    answer.at(12,1)= D*K.at(12,1);
    answer.at(12,2)= D*K.at(12,2);
    answer.at(12,3)= D*K.at(12,3);

    answer.at(4,4)= answer.at(1,1);

    answer.at(5,4)= -answer.at(2,1);
    answer.at(5,5)= answer.at(2,2);

    answer.at(6,4)= answer.at(3,1);
    answer.at(6,5)= -answer.at(3,2);
    answer.at(6,6)= answer.at(3,3);

    answer.at(7,4)= answer.at(10,1);
    answer.at(7,5)= -answer.at(10,2);
    answer.at(7,6)= answer.at(10,3);
    answer.at(7,7)= answer.at(1,1);

    answer.at(8,4)= -answer.at(11,1);
    answer.at(8,5)= answer.at(11,2);
    answer.at(8,6)= -answer.at(11,3);
    answer.at(8,7)= answer.at(2,1);
    answer.at(8,8)= answer.at(2,2);

    answer.at(9,4)= answer.at(12,1);
    answer.at(9,5)= -answer.at(12,2);
    answer.at(9,6)= answer.at(12,3);
    answer.at(9,7)= -answer.at(3,1);
    answer.at(9,8)= -answer.at(3,2);
    answer.at(9,9)= answer.at(3,3);

    answer.at(10,4)= answer.at(7,1);
    answer.at(10,5)= -answer.at(7,2);
    answer.at(10,6)= answer.at(7,3);
    answer.at(10,7)= answer.at(4,1);
    answer.at(10,8)= answer.at(4,2);
    answer.at(10,9)= -answer.at(4,3);
    answer.at(10,10)= answer.at(1,1);

    answer.at(11,4)= -answer.at(8,1);
    answer.at(11,5)= answer.at(8,2);
    answer.at(11,6)= -answer.at(8,3);
    answer.at(11,7)= answer.at(5,1);
    answer.at(11,8)= answer.at (5,2);
    answer.at(11,9)= -answer.at(5,3);
    answer.at(11,10)= -answer.at(2,1);
    answer.at(11,11)= answer.at(2,2);

    answer.at(12,4)= answer.at(9,1);
    answer.at(12,5)= -answer.at(9,2);
    answer.at(12,6)= answer.at(9,3);
    answer.at(12,7)= -answer.at(6,1);
    answer.at(12,8)= -answer.at(6,2);
    answer.at(12,9)= answer.at(6,3);
    answer.at(12,10)= -answer.at(3,1);
    answer.at(12,11)= answer.at(3,2);
    answer.at(12,12)= answer.at(3,3);


    for ( int j = 2; j <= 12; j++ ) {
        for ( int i = 1; i < j; i++ ) {
            answer.at(i, j) = answer.at(j, i);
        }
    }


}



void
PPlate::computeKa(FloatMatrix &answer)
{//computes Ka matrix


    double k;//,z;

    //this->crossSectioncharacteristics();
    //this->computeLength();


    k= ((b/a)*(b/a)/70); //(1/3)*(b/a)*(b/a);

    answer.resize(12,3);
    answer.zero();

    answer.at(1,1)= 312*k;
    answer.at(1,2)= 156*k;
    answer.at(1,3)= 44*k;

    answer.at(2,1)= 156*k;
    answer.at(2,2)= 104*k;
    answer.at(2,3)= 22*k;

    answer.at(3,1)= 44*k;
    answer.at(3,2)= 22*k;
    answer.at(3,3)= 8*k;

    answer.at(4,1)= -312*k;
    answer.at(4,2)= -156*k;
    answer.at(4,3)= -44*k;

    answer.at(5,1)= 156*k;
    answer.at(5,2)= 52*k;
    answer.at(5,3)= 22*k;

    answer.at(6,1)= -44*k;
    answer.at(6,2)= -22*k;
    answer.at(6,3)= -8*k;

    answer.at(7,1)= 108*k;
    answer.at(7,2)= 54*k;
    answer.at(7,3)= 26*k;

    answer.at(8,1)= 54*k;
    answer.at(8,2)= 36*k;
    answer.at(8,3)= 13*k;

    answer.at(9,1)= -26*k;
    answer.at(9,2)= -13*k;
    answer.at(9,3)= -6*k;

    answer.at(10,1)= -108*k;
    answer.at(10,2)= -54*k;
    answer.at(10,3)= -26*k;

    answer.at(11,1)= 54*k;
    answer.at(11,2)= 18*k;
    answer.at(11,3)= 13*k;

    answer.at(12,1)= 26*k;
    answer.at(12,2)= 13*k;
    answer.at(12,3)= 6*k;



}

void
PPlate::computeKb(FloatMatrix &answer)
{//computes Kb matrix


    double k;//,z;
    //this->crossSectioncharacteristics();
    //this->computeLength();

    //z=(1/70);

    k=((a/b)*(a/b)/70); //(1/3)*(a/b)*(a/b);

    answer.resize(12,3);
    answer.zero();


    answer.at(1,1)= 312*k;
    answer.at(1,2)= 44*k;
    answer.at(1,3)= 156*k;

    answer.at(2,1)= 44*k;
    answer.at(2,2)= 8*k;
    answer.at(2,3)= 22*k;

    answer.at(3,1)= 156*k;
    answer.at(3,2)= 22*k;
    answer.at(3,3)= 104*k;

    answer.at(4,1)= 108*k;
    answer.at(4,2)= 26*k;
    answer.at(4,3)= 54*k;

    answer.at(5,1)= -26*k;
    answer.at(5,2)= -6*k;
    answer.at(5,3)= -13*k;

    answer.at(6,1)= 54*k;
    answer.at(6,2)= 13*k;
    answer.at(6,3)= 36*k;

    answer.at(7,1)= -312*k;
    answer.at(7,2)= -44*k;
    answer.at(7,3)= -156*k;

    answer.at(8,1)= -44*k;
    answer.at(8,2)= -8*k;
    answer.at(8,3)= -22*k;

    answer.at(9,1)= 156*k;
    answer.at(9,2)= 22*k;
    answer.at(9,3)= 52*k;

    answer.at(10,1)= -108*k;
    answer.at(10,2)= -26*k;
    answer.at(10,3)= -54*k;

    answer.at(11,1)= 26*k;
    answer.at(11,2)= 6*k;
    answer.at(11,3)= 13*k;

    answer.at(12,1)= 54*k;
    answer.at(12,2)= 13*k;
    answer.at(12,3)= 18*k;
}

void
PPlate::computeKc(FloatMatrix &answer)
{//computes Kc matrix


    double k;
    this->crossSectioncharacteristics();
    this->computeLength();

    k=(ni/50);

    answer.resize(12,3);
    answer.zero();

    answer.at(1,1)= 144*k;
    answer.at(1,2)= 72*k;
    answer.at(1,3)= 72*k;

    answer.at(2,1)= 72*k;
    answer.at(2,2)= 16*k;
    answer.at(2,3)= 61*k;

    answer.at(3,1)= 72*k;
    answer.at(3,2)= 61*k;
    answer.at(3,3)= 16*k;

    answer.at(4,1)= -144*k;
    answer.at(4,2)= -12*k;
    answer.at(4,3)= -72*k;

    answer.at(5,1)= 12*k;
    answer.at(5,2)= -4*k;
    answer.at(5,3)= 6*k;

    answer.at(6,1)= -72*k;
    answer.at(6,2)= -6*k;
    answer.at(6,3)= -16*k;

    answer.at(7,1)= -144*k;
    answer.at(7,2)= -72*k;
    answer.at(7,3)= -12*k;

    answer.at(8,1)= -72*k;
    answer.at(8,2)= -16*k;
    answer.at(8,3)= -6*k;

    answer.at(9,1)= 12*k;
    answer.at(9,2)= 6*k;
    answer.at(9,3)= -4*k;

    answer.at(10,1)= 144*k;
    answer.at(10,2)= 12*k;
    answer.at(10,3)= 12*k;

    answer.at(11,1)= -12*k;
    answer.at(11,2)= 4*k;
    answer.at(11,3)= -1*k;

    answer.at(12,1)= -12*k;
    answer.at(12,2)= -1*k;
    answer.at(12,3)= 4*k;



}

void
PPlate::computeKd(FloatMatrix &answer)
{//computes Kd matrix


    double k;
    //this->crossSectioncharacteristics();
    //this->computeLength();

    k=(1-ni)/50;

    answer.resize(12,3);
    answer.zero();

    answer.at(1,1)= 144*k;
    answer.at(1,2)= 12*k;
    answer.at(1,3)= 12*k;

    answer.at(2,1)= 12*k;
    answer.at(2,2)= 16*k;
    answer.at(2,3)= 1*k;

    answer.at(3,1)= 12*k;
    answer.at(3,2)= 1*k;
    answer.at(3,3)= 16*k;

    answer.at(4,1)= -144*k;
    answer.at(4,2)= -12*k;
    answer.at(4,3)= -12*k;

    answer.at(5,1)= 12*k;
    answer.at(5,2)= -4*k;
    answer.at(5,3)= 1*k;

    answer.at(6,1)= -12*k;
    answer.at(6,2)= -1*k;
    answer.at(6,3)= -16*k;

    answer.at(7,1)= -144*k;
    answer.at(7,2)= -12*k;
    answer.at(7,3)= -12*k;

    answer.at(8,1)= -12*k;
    answer.at(8,2)= -16*k;
    answer.at(8,3)= -1*k;

    answer.at(9,1)= 12*k;
    answer.at(9,2)= 1*k;
    answer.at(9,3)= -4*k;

    answer.at(10,1)= 144*k;
    answer.at(10,2)= 12*k;
    answer.at(10,3)= 12*k;

    answer.at(11,1)= -12*k;
    answer.at(11,2)= 4*k;
    answer.at(11,3)= -1*k;

    answer.at(12,1)= -12*k;
    answer.at(12,2)= -1*k;
    answer.at(12,3)= 4*k;


}

double
PPlate :: computeLength()
// Returns the edge lengths of the receiver.
{
        this->computeLCS();
        double  x12 , y12, x14, y14, x23, y23, x34, y34 ;

        if(a+b==0){
            x12=lnodes[1][0]-lnodes[0][0];
            y12=lnodes[1][1]-lnodes[0][1];
            x14=lnodes[3][0]-lnodes[0][0];
            y14=lnodes[3][1]-lnodes[0][1];
            x23=lnodes[2][0]-lnodes[1][0];
            y23=lnodes[2][1]-lnodes[1][1];
            x34=lnodes[2][0]-lnodes[3][0];
            y34=lnodes[2][1]-lnodes[3][1];

            a = (sqrt(x12*x12+y12*y12)+sqrt(x34*x34+y34*y34))/2;
            b = (sqrt(x14*x14+y14*y14)+sqrt(x23*x23+y23*y23))/2;
    }

  return a;
}

void
PPlate :: computeLCS()
{
    lcsMatrix.resize(3, 3); // Note! G -> L transformation matrix
    FloatArray e1, e2, e3, help;

    // compute e1' = [N2-N1]  and  help = [N4-N1]

	// Pointer operators removed because the project could not build. Test later to check if the method still works properly.
    e1.beDifferenceOf( this->giveNode(2)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );
    help.beDifferenceOf( this->giveNode(4)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );

    e1.normalize();
    e3.beVectorProductOf(e1, help);
    e3.normalize();
    e2.beVectorProductOf(e3, e1);
    for ( int i = 1; i <= 3; i++ ) {
        this->lcsMatrix.at(1, i) = e1.at(i);
        this->lcsMatrix.at(2, i) = e2.at(i);
        this->lcsMatrix.at(3, i) = e3.at(i);
    }
    for ( int i = 1; i <= 4; i++ ) {
		// Pointer operator removed because the project could not build. Test later to check if the method still works properly.
        this->lnodes [ i - 1 ].beProductOf( this->lcsMatrix, this->giveNode(i)->giveCoordinates() );
    }

}

bool
PPlate :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    this->computeLCS();
    answer.resize(12, 12);
    answer.zero();

    for ( int i = 0; i <= 3; i++ ) { // Loops over nodes
        // In each node, transform global c.s. {D_w,R_u,R_v} into local c.s.
        answer.setSubMatrix(this->lcsMatrix, 1 + i * 3, 1 + i * 3);     // Displacements
    }

    return true;
}

int
PPlate :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [5,6]
// f(local) = T * f(global)
{
    this->computeLCS();
   // this->computeGtoLRotationMatrix();

    answer.resize(3, 3);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lcsMatrix.at(1, i);
        answer.at(2, i) = lcsMatrix.at(2, i);
        answer.at(3, i) = lcsMatrix.at(3, i);
    }

    return 1;
}

void
PPlate :: giveNodeCoordinates(double &x1, double &x2, double &x3, double &x4,
                                 double &y1, double &y2, double &y3, double &y4,
                                 double &z1, double &z2, double &z3, double &z4)
{
	// Pointer operators removed because the project could not build. Test later to check if the method still works properly.
    FloatArray nc1, nc2, nc3, nc4;
    nc1 = this->giveNode(1)->giveCoordinates();
    nc2 = this->giveNode(2)->giveCoordinates();
    nc3 = this->giveNode(3)->giveCoordinates();
    nc4 = this->giveNode(4)->giveCoordinates();

    x1 = nc1.at(1);
    x2 = nc2.at(1);
    x3 = nc3.at(1);
    x4 = nc4.at(1);

    y1 = nc1.at(2);
    y2 = nc2.at(2);
    y3 = nc3.at(2);
    y4 = nc4.at(2);

    z1 = nc1.at(3);
    z2 = nc2.at(3);
    z3 = nc3.at(3);
    z4 = nc4.at(3);

}


void
PPlate :: initializeFrom(InputRecord &ir)
{
    StructuralElement :: initializeFrom(ir);
}


void
PPlate :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_w, R_u, R_v};
}


void
PPlate :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
// returns normal vector to midPlane in GaussPoinr gp of receiver
{
    FloatArray u, v;
	// Pointer operators removed because the project could not build. Test later to check if the method still works properly.
    u.beDifferenceOf( this->giveNode(2)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );
    v.beDifferenceOf( this->giveNode(3)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );

    answer.beVectorProductOf(u, v);
    answer.normalize();
}




#define POINT_TOL 1.e-3

bool
PPlate :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
//converts global coordinates to local planar area coordinates,
//does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
{
    // get node coordinates
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    this->giveNodeCoordinates(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);

    // Fetch local coordinates.
    bool ok = this->interp_lin.global2local( answer, coords, FEIElementGeometryWrapper(this) ) > 0;

    //check that the point is in the element and set flag
    for ( int i = 1; i <= 4; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            return false;
        }

        if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            return false;
        }
    }

    //get midplane location at this point
    double midplZ;
    midplZ = z1 * answer.at(1) + z2 * answer.at(2) + z3 * answer.at(3) + z4 * answer.at(4);

    //check that the z is within the element
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    double elthick = cs->give(CS_Thickness, answer, this);

    if ( elthick / 2.0 + midplZ - fabs( coords.at(3) ) < -POINT_TOL ) {
        answer.zero();
        return false;
    }


    return ok;
}


int
PPlate :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatArray help;
    answer.resize(6);
    if ( type == IST_ShellForceTensor || type == IST_ShellStrainTensor ) {
        if ( type == IST_ShellForceTensor ) {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        } else {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        }
        answer.at(1) = 0.0; // nx
        answer.at(2) = 0.0; // ny
        answer.at(3) = 0.0; // nz
        answer.at(4) = help.at(5); // vyz
        answer.at(5) = help.at(4); // vxz
        answer.at(6) = 0.0; // vxy
        return 1;
    } else if ( type == IST_ShellMomentTensor || type == IST_CurvatureTensor ) {
        if ( type == IST_ShellMomentTensor ) {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        } else {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        }
        answer.at(1) = help.at(1); // mx
        answer.at(2) = help.at(2); // my
        answer.at(3) = 0.0;        // mz
        answer.at(4) = 0.0;        // myz
        answer.at(5) = 0.0;        // mxz
        answer.at(6) = help.at(3); // mxy
        return 1;
    } else {
        return StructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}


} // end namespace oofem


