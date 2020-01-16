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

#include "../sm/Elements/Perleta/basiclsr3d.h"
#include "../sm/Materials/structuralms.h"
#include "Materials/structuralmaterial.h"
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
#include "fei2dquadlin.h"
#include "classfactory.h"
#include "elementinternaldofman.h"
#include "masterdof.h"
#include "../sm/CrossSections/structuralcrosssection.h"



#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(BasicLSR3d);

FEI2dQuadLin BasicLSR3d :: interp(1, 2);

BasicLSR3d :: BasicLSR3d(int n, Domain *aDomain) : StructuralElement(n, aDomain),

    lnodes(4)
{
    numberOfDofMans = 4;

    numberOfGaussPoints = 4;
    E=0;
    t=0;
    ni=0;
    a=0;
    b=0;
    referenceNode=0;
}

BasicLSR3d :: ~BasicLSR3d()
{ }

FEInterpolation *BasicLSR3d :: giveInterpolation() const { return & interp; }

void
BasicLSR3d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns the [3x8] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{

    FloatArray Gcoords;
    double x,y;

    Gcoords   = gp->giveGlobalCoordinates();
    x   = Gcoords(0);
    y   = Gcoords(1);

    answer.resize(3, 8);
    answer.zero();

    answer.at(1, 1) =  ( y / (a*b))-( 1/a ) ;
    answer.at(1, 2) =  0;
    answer.at(1, 3) =  ( 1/a )- (  y / (a*b) );
    answer.at(1, 4) =  0;
    answer.at(1, 5) =  y / (a*b);
    answer.at(1, 6) =  0;
    answer.at(1, 7) =  -( y / (a*b));
    answer.at(1, 8) =   0;

    answer.at(2, 1) =  0;
    answer.at(2, 2) =  ( x / (a*b))-( 1/b ) ;
    answer.at(2, 3) =  0;
    answer.at(2, 4) =  -( x / (a*b));
    answer.at(2, 5) =  0;
    answer.at(2, 6) =   ( x / (a*b));
    answer.at(2, 7) =  0;
    answer.at(2, 8) =   ( 1/b )- (  x / (a*b) );

    answer.at(3, 1) =  ( x / (a*b))-( 1/b ) ;
    answer.at(3, 2) =  ( y / (a*b))-( 1/a );
    answer.at(3, 3) =  -( x / (a*b)) ;
    answer.at(3, 4) =   ( 1/a )- (  y / (a*b) );
    answer.at(3, 5) =   ( x / (a*b));
    answer.at(3, 6) =   y / (a*b);
    answer.at(3, 7) =   ( 1/b )- (  x / (a*b) );
    answer.at(3, 8) =  -( y / (a*b)) ;
}


void BasicLSR3d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        // the gauss point is used only when methods from crosssection and/or material
        // classes are requested
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
        this->crossSectioncharacteristics();
        this->computeLength();
    }
}


void
BasicLSR3d :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp. Used for numerical calculation of consistent mass
// matrix. Must contain only interpolation for displacement terms,
// not for any rotations. (Inertia forces do not work on rotations).
// r = {u1,v1,u2,v2,u3,v3,u4,v4}^T
{
    FloatArray Gcoords;
    GaussPoint *gp = integrationRulesArray[0]->getIntegrationPoint(0);
    double x,y;

    Gcoords   = gp->giveGlobalCoordinates();

    x   = lnodes[0][0];
    y   = lnodes[0][1];

    answer.resize(2, 8);
    answer.zero();

    answer.at(1, 1) = ( (x-a)*(y-b))/(a*b);
    answer.at(1, 2) = 0;
    answer.at(1, 3) = -(x*(y-b))/(a*b);
    answer.at(1, 4) = 0;
    answer.at(1, 5) =(x*y)/(a*b);
    answer.at(1, 6) = 0;
    answer.at(1, 7) = -((x-a)*y)/(a*b);
    answer.at(1, 8) = 0;

    answer.at(2, 1) = 0;
    answer.at(2, 2) =( (x-a)*(y-b))/(a*b);
    answer.at(2, 3) = 0;
    answer.at(2, 4) = -(x*(y-b))/(a*b);
    answer.at(2, 5) = 0;
    answer.at(2, 6) =(x*y)/(a*b);
    answer.at(2, 7) = 0;
    answer.at(2, 8) = -((x-a)*y)/(a*b);
}

void
BasicLSR3d :: crossSectioncharacteristics()
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
BasicLSR3d :: computeLocalStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    // compute clamped stifness
    this->computeClampedStiffnessMatrix(answer, rMode, tStep);
}


void
BasicLSR3d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    // compute clamped stiffness
    this->computeLocalStiffnessMatrix(answer, rMode, tStep);
}


void
BasicLSR3d :: computeClampedStiffnessMatrix(FloatMatrix &answer,
                                        MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the local
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
     double a2,b2;


     a2=a*a;
     b2=b*b;


     answer.resize(8, 8);
     answer.zero();

     answer.at(1,1)= (a2*ni-2*b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));
     answer.at(1,2)= -(t*E)/(8*(ni-1));
     answer.at(1,3)= (a2*ni+4*b2-a2)*t*E/(12*a*b*(ni-1)*(ni+1));
     answer.at(1,4)= -(3*ni-1)*t*E/(8*(ni-1)*(ni+1));
     answer.at(1,5)= -(a2*ni-2*b2-a2)*t*E/(12*a*b*(ni-1)*(ni+1));
     answer.at(1,6)= (t*E)/(8*(ni-1));
     answer.at(1,7)= -(a2*ni+b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));
     answer.at(1,8)= (3*ni-1)*t*E/(8*(ni-1)*(ni+1));

     answer.at(2,2)= (b2*ni-b2-2*a2)*t*E/(6*a*b*(ni-1)*(ni+1));
     answer.at(2,3)= (3*ni-1)*t*E/(8*(ni-1)*(ni+1));
     answer.at(2,4)= -(b2*ni-b2+a2)*t*E/(6*a*b*(ni-1)*(ni+1));
     answer.at(2,5)= (t*E)/(8*(ni-1));
     answer.at(2,6)= -(b2*ni-b2-2*a2)*t*E/(12*a*b*(ni-1)*(ni+1));
     answer.at(2,7)= -(3*ni-1)*t*E/(8*(ni-1)*(ni+1));
     answer.at(2,8)= (b2*ni-b2+4*a2)*t*E/(12*a*b*(ni-1)*(ni+1));

     answer.at(3,3)= (a2*ni-2*b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));
     answer.at(3,4)= (t*E)/(8*(ni-1));
     answer.at(3,5)= -(a2*ni+b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));
     answer.at(3,6)= -(3*ni-1)*t*E/(8*(ni-1)*(ni+1));
     answer.at(3,7)= -(a2*ni-2*b2-a2)*t*E/(12*a*b*(ni-1)*(ni+1));
     answer.at(3,8)= -(t*E)/(8*(ni-1));

     answer.at(4,4)= (b2*ni-b2-2*a2)*t*E/(6*a*b*(ni-1)*(ni+1));
     answer.at(4,5)= (3*ni-1)*t*E/(8*(ni-1)*(ni+1));
     answer.at(4,6)= (b2*ni-b2+4*a2)*t*E/(12*a*b*(ni-1)*(ni+1));
     answer.at(4,7)= -(t*E)/(8*(ni-1));
     answer.at(4,8)= -(b2*ni-b2-2*a2)*t*E/(12*a*b*(ni-1)*(ni+1));

     answer.at(5,5)= (a2*ni-2*b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));
     answer.at(5,6)= -(t*E)/(8*(ni-1));
     answer.at(5,7)= (a2*ni+4*b2-a2)*t*E/(12*a*b*(ni-1)*(ni+1));
     answer.at(5,8)= -(3*ni-1)*t*E/(8*(ni-1)*(ni+1));

     answer.at(6,6)= (b2*ni-b2-2*a2)*t*E/(6*a*b*(ni-1)*(ni+1));
     answer.at(6,7)= (3*ni-1)*t*E/(8*(ni-1)*(ni+1));
     answer.at(6,8)= -(b2*ni-b2+a2)*t*E/(6*a*b*(ni-1)*(ni+1));

     answer.at(7,7)= (a2*ni-2*b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));
     answer.at(7,8)= (t*E)/(8*(ni-1));

     answer.at(8,8)= (b2*ni-b2-2*a2)*t*E/(6*a*b*(ni-1)*(ni+1));




    answer.symmetrized();
}




double
BasicLSR3d :: computeLength()
// Returns the edge lengths of the receiver.
{
    this->computeLCS();
    double  x12 , y12, x14, y14 ;

    if(a+b==0){
        x12=lnodes[1][0]-lnodes[0][0];
        y12=lnodes[1][1]-lnodes[0][1];
        x14=lnodes[3][0]-lnodes[0][0];
        y14=lnodes[3][1]-lnodes[0][1];
        a = sqrt(x12*x12+y12*y12);
        b = sqrt(x14*x14+y14*y14);

    }

  return a;
}


void
BasicLSR3d :: initializeFrom(InputRecord &ir)
{
    numberOfGaussPoints = 4;
    StructuralElement :: initializeFrom(ir);

    if ( numberOfGaussPoints != 1 && numberOfGaussPoints != 4 && numberOfGaussPoints != 9 && numberOfGaussPoints != 16 && numberOfGaussPoints != 25 ) {
        numberOfGaussPoints = 4;
        OOFEM_WARNING("Number of Gauss points enforced to 1");
    }
}



void
BasicLSR3d :: computeLCS()
{
    lcsMatrix.resize(2, 3); // Note! G -> L transformation matrix
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
        //this->lcsMatrix.at(3, i) = e3.at(i);
    }
    for ( int i = 1; i <= 4; i++ ) {
		// Pointer operator removed because the project could not build. Test later to check if the method still works properly.
        this->lnodes [ i - 1 ].beProductOf( this->lcsMatrix, this->giveNode(i)->giveCoordinates() );
    }

}


bool
BasicLSR3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    this->computeLCS();
    answer.resize(8, 12);
    answer.zero();

    for ( int i = 0; i <= 3; i++ ) { // Loops over nodes
        // In each node, transform global c.s. {D_u, D_v,D_w} into local c.s.
        answer.setSubMatrix(this->lcsMatrix, 1 + i * 2, 1 + i * 3);     // Displacements
    }

    return true;
}

int
BasicLSR3d :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [5,6]
// f(local) = T * f(global)
{
    this->computeLCS();

    answer.resize(2, 3);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lcsMatrix.at(1, i);
        answer.at(2, i) = lcsMatrix.at(2, i);
        //answer.at(3, i) = answer.at(6, i + 3) = GtoLRotationMatrix.at(3, i);
    }

    return 1;
}

void
BasicLSR3d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v,D_w
    };
}

void
BasicLSR3d :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    StructuralElement :: updateLocalNumbering(f);
    if ( this->referenceNode ) {
        this->referenceNode = f(this->referenceNode, ERS_DofManager);
    }
}


} // end namespace oofem
