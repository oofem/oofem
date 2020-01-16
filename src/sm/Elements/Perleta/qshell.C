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

#include "../sm/Elements/Perleta/qshell.h"
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
REGISTER_Element(qshell);

FEI2dQuadLin qshell :: interp(1, 2);

qshell :: qshell(int n, Domain *aDomain) : StructuralElement(n, aDomain),

    lnodes(4)
{
    numberOfDofMans = 4;

    numberOfGaussPoints = 4;
    E=0;
    t=0;
    ni=0;
    a=0;
    b=0;
    K=0;
    referenceNode=0;

}

qshell :: ~qshell()
{ }

FEInterpolation *qshell :: giveInterpolation() const { return & interp; }

void
qshell :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns the [8x24] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{

    answer.resize(8,24);
    answer.zero();

}


void qshell :: computeGaussPoints()
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
qshell :: crossSectioncharacteristics()
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
qshell :: computeLocalStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    // compute clamped stifness
    this->computeClampedStiffnessMatrix(answer, rMode, tStep);
}


void
qshell :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    // compute clamped stiffness
    this->computeLocalStiffnessMatrix(answer, rMode, tStep);
}


void
qshell :: computeClampedStiffnessMatrix(FloatMatrix &answer,
                                        MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the local
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
     double a2,a3,a4,b2,b3,b4,b5,t3,ni2;


       this->crossSectioncharacteristics();
       this->computeLength();

     a2=a*a;
     a3=a*a*a;
     a4=a*a*a*a;
     b2=b*b;
     b3=b*b*b;
     b4=b*b*b*b;
     b5=b*b*b*b*b;
     ni2=ni*ni;
     t3=t*t*t;

     answer.resize(24, 24);
     answer.zero();
// membrane part of stiffness matrix:

     answer.at(1,1)= (a2*ni-2*b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));     // (1,1)
     answer.at(1,2)= -(t*E)/(8*(ni-1));                             // (1,2)
     answer.at(1,7)= (a2*ni+4*b2-a2)*t*E/(12*a*b*(ni-1)*(ni+1));    // (1,3)
     answer.at(1,8)= -(3*ni-1)*t*E/(8*(ni-1)*(ni+1));               // (1,4)
     answer.at(1,13)= -(a2*ni-2*b2-a2)*t*E/(12*a*b*(ni-1)*(ni+1));   // (1,5)
     answer.at(1,14)= (t*E)/(8*(ni-1));                              // (1,6)
     answer.at(1,19)= -(a2*ni+b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));      // (1,7)
     answer.at(1,20)= (3*ni-1)*t*E/(8*(ni-1)*(ni+1));                // (1,8)

     answer.at(2,2)= (b2*ni-b2-2*a2)*t*E/(6*a*b*(ni-1)*(ni+1));     // (2,2)
     answer.at(2,7)= (3*ni-1)*t*E/(8*(ni-1)*(ni+1));                // (2,3)
     answer.at(2,8)= -(b2*ni-b2+a2)*t*E/(6*a*b*(ni-1)*(ni+1));      // (2,4)
     answer.at(2,13)= (t*E)/(8*(ni-1));                              // (2,5)
     answer.at(2,14)= -(b2*ni-b2-2*a2)*t*E/(12*a*b*(ni-1)*(ni+1));   // (2,6)
     answer.at(2,19)= -(3*ni-1)*t*E/(8*(ni-1)*(ni+1));               // (2,7)
     answer.at(2,20)= (b2*ni-b2+4*a2)*t*E/(12*a*b*(ni-1)*(ni+1));    // (2,8)


     answer.at(7,7)= (a2*ni-2*b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));     // (3,3)
     answer.at(7,8)= (t*E)/(8*(ni-1));                              // (3,4)
     answer.at(7,13)= -(a2*ni+b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));      // (3,5)
     answer.at(7,14)= -(3*ni-1)*t*E/(8*(ni-1)*(ni+1));               // (3,6)
     answer.at(7,19)= -(a2*ni-2*b2-a2)*t*E/(12*a*b*(ni-1)*(ni+1));   // (3,7)
     answer.at(7,20)= -(t*E)/(8*(ni-1));                             // (3,8)

     answer.at(8,8)= (b2*ni-b2-2*a2)*t*E/(6*a*b*(ni-1)*(ni+1));     // (4,4)
     answer.at(8,13)= (3*ni-1)*t*E/(8*(ni-1)*(ni+1));                // (4,5)
     answer.at(8,14)= (b2*ni-b2+4*a2)*t*E/(12*a*b*(ni-1)*(ni+1));    // (4,6)
     answer.at(8,19)= -(t*E)/(8*(ni-1));                             // (4,7)
     answer.at(8,20)= -(b2*ni-b2-2*a2)*t*E/(12*a*b*(ni-1)*(ni+1));   // (4,8)



     answer.at(13,13)= (a2*ni-2*b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));     // (5,5)
     answer.at(13,14)= -(t*E)/(8*(ni-1));                             // (5,6)
     answer.at(13,19)= (a2*ni+4*b2-a2)*t*E/(12*a*b*(ni-1)*(ni+1));    // (5,7)
     answer.at(13,20)= -(3*ni-1)*t*E/(8*(ni-1)*(ni+1));               // (5,8)

     answer.at(14,14)= (b2*ni-b2-2*a2)*t*E/(6*a*b*(ni-1)*(ni+1));     // (6,6)
     answer.at(14,19)= (3*ni-1)*t*E/(8*(ni-1)*(ni+1));                // (6,7)
     answer.at(14,20)= -(b2*ni-b2+a2)*t*E/(6*a*b*(ni-1)*(ni+1));      // (6,8)



     answer.at(19,19)= (a2*ni-2*b2-a2)*t*E/(6*a*b*(ni-1)*(ni+1));     // (7,7)
     answer.at(19,20)= (t*E)/(8*(ni-1));                              // (7,8)

     answer.at(20,20)= (b2*ni-b2-2*a2)*t*E/(6*a*b*(ni-1)*(ni+1));     // (8,8)

//plate part of stiffness matrix:


     answer.at(3,3)=((2*a2*b2*ni-10*b4-7*a2*b2-10*a4)*t3*E)/(30*a3*b3*ni2-30*a3*b3);    //(1,1)
     answer.at(3,4)=-((4*b5*ni+b5+10*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));             //(1,2)
     answer.at(3,5)=-((4*a2*b3*ni+10*b5+a2*b3)*t3*E)/(60*a2*b4*ni2-60*a2*b4);           //(1,3)
     answer.at(3,9)=-((2*a2*b2*ni-10*b4-7*a2*b2+5*a4)*t3*E)/(30*a3*b3*ni2-30*a3*b3);    //(1,4)
     answer.at(3,10)=((4*b5*ni+b5-5*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));               //(1,5)
     answer.at(3,11)=-((-a2*b*ni+10*b3+a2*b)*t3*E)/(60*a2*b2*ni2-60*a2*b2);              //(1,6)
     answer.at(3,15)=((2*a2*b2*ni+5*b4-7*a2*b2+5*a4)*t3*E)/(30*a3*b3*ni2-30*a3*b3);      //(1,7)
     answer.at(3,16)=-((b5*ni-b5+5*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));                //(1,8)
     answer.at(3,17)=((-a2*b*ni-5*b3+a2*b)*t3*E)/(60*a2*b2*ni2-60*a2*b2);                //(1,9)
     answer.at(3,21)=-((2*a2*b2*ni+5*b4-7*a2*b2-10*a4)*t3*E)/(30*a3*b3*ni2-30*a3*b3);   //(1,10)
     answer.at(3,22)=((b5*ni-b5-10*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));               //(1,11)
     answer.at(3,23)=((4*a2*b3*ni-5*b5+a2*b3)*t3*E)/(60*a2*b4*ni2-60*a2*b4);            //(1,12)


     answer.at(4,4)=((2*b5*ni-2*b5-10*a2*b3)*t3*E)/(5*(18*a*b4*ni2-18*a*b4));           //(2,2)
     answer.at(4,5)=-(b3*ni*t3*E)/(2*(6*b3*ni2-6*b3));                                  //(2,3)
     answer.at(4,9)=((4*b5*ni+b5-5*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));               //(2,4)
     answer.at(4,10)=-((2*b5*ni-2*b5+5*a2*b3)*t3*E)/(5*(18*a*b4*ni2-18*a*b4));           //(2,5)
     answer.at(4,11)=0;                                                                  //(2,6)
     answer.at(4,15)=((b5*ni-b5+5*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));                 //(2,7)
     answer.at(4,16)=((b5*ni-b5-5*a2*b3)*t3*E)/(10*(18*a*b4*ni2-18*a*b4));               //(2,8)
     answer.at(4,17)=0;                                                                  //(2,9)
     answer.at(4,21)=-((b5*ni-b5-10*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));              //(2,10)
     answer.at(4,22)=-((b5*ni-b5+10*a2*b3)*t3*E)/(10*(18*a*b4*ni2-18*a*b4));            //(2,11)
     answer.at(4,23)=0;                                                                 //(2,12)


     answer.at(5,5)=-((-a2*b*ni+5*b3+a2*b)*t3*E)/(45*a*b2*ni2-45*a*b2);                 //(3,3)
     answer.at(5,9)=((-a2*b*ni+10*b3+a2*b)*t3*E)/(60*a2*b2*ni2-60*a2*b2);               //(3,4)
     answer.at(5,10)=0;                                                                  //(3,5)
     answer.at(5,11)=-((a2*b*ni+10*b3-a2*b)*t3*E)/(180*a*b2*ni2-180*a*b2);               //(3,6)
     answer.at(5,15)=-((-a2*b*ni-5*b3+a2*b)*t3*E)/(60*a2*b2*ni2-60*a2*b2);               //(3,7)
     answer.at(5,16)=0;                                                                  //(3,8)
     answer.at(5,17)=((a2*b*ni-5*b3-a2*b)*t3*E)/(180*a*b2*ni2-180*a*b2);                 //(3,9)
     answer.at(5,21)=((4*a2*b3*ni-5*b5+a2*b3)*t3*E)/(60*a2*b4*ni2-60*a2*b4);            //(3,10)
     answer.at(5,22)=0;                                                                 //(3,11)
     answer.at(5,23)=-((2*a2*b*ni+5*b3-2*a2*b)*t3*E)/(2*(45*a*b2*ni2-45*a*b2));         //(3,12)


     answer.at(9,9)=((2*a2*b2*ni-10*b4-7*a2*b2-10*a4)*t3*E)/(30*a3*b3*ni2-30*a3*b3);    //(4,4)
     answer.at(9,10)=-((4*b5*ni+b5+10*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));             //(4,5)
     answer.at(9,11)=((4*a2*b3*ni+10*b5+a2*b3)*t3*E)/(60*a2*b4*ni2-60*a2*b4);            //(4,6)
     answer.at(9,15)=-((2*a2*b2*ni+5*b4-7*a2*b2-10*a4)*t3*E)/(30*a3*b3*ni2-30*a3*b3);    //(4,7)
     answer.at(9,16)=((b5*ni-b5-10*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));                //(4,8)
     answer.at(9,17)=-((4*a2*b3*ni-5*b5+a2*b3)*t3*E)/(60*a2*b4*ni2-60*a2*b4);            //(4,9)
     answer.at(9,21)=((2*a2*b2*ni+5*b4-7*a2*b2+5*a4)*t3*E)/(30*a3*b3*ni2-30*a3*b3);     //(4,10)
     answer.at(9,22)=-((b5*ni-b5+5*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));               //(4,11)
     answer.at(9,23)=-((-a2*b*ni-5*b3+a2*b)*t3*E)/(60*a2*b2*ni2-60*a2*b2);              //(4,12)


     answer.at(10,10)=((2*b5*ni-2*b5-10*a2*b3)*t3*E)/(5*(18*a*b4*ni2-18*a*b4));           //(5,5)
     answer.at(10,11)=(b3*ni*t3*E)/(2*(6*b3*ni2-6*b3));                                   //(5,6)
     answer.at(10,15)=-((b5*ni-b5-10*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));               //(5,7)
     answer.at(10,16)=-((b5*ni-b5+10*a2*b3)*t3*E)/(10*(18*a*b4*ni2-18*a*b4));             //(5,8)
     answer.at(10,17)=0;                                                                  //(5,9)
     answer.at(10,21)=((b5*ni-b5+5*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));                //(5,10)
     answer.at(10,22)=((b5*ni-b5-5*a2*b3)*t3*E)/(10*(18*a*b4*ni2-18*a*b4));              //(5,11)
     answer.at(10,23)=0;                                                                 //(5,12)


     answer.at(11,11)=-((-a2*b*ni+5*b3+a2*b)*t3*E)/(45*a*b2*ni2-45*a*b2);                 //(6,6)
     answer.at(11,15)=-((4*a2*b3*ni-5*b5+a2*b3)*t3*E)/(60*a2*b4*ni2-60*a2*b4);            //(6,7)
     answer.at(11,16)=0;                                                                  //(6,8)
     answer.at(11,17)=-((2*a2*b*ni+5*b3-2*a2*b)*t3*E)/(2*(45*a*b2*ni2-45*a*b2));          //(6,9)
     answer.at(11,21)=((-a2*b*ni-5*b3+a2*b)*t3*E)/(60*a2*b2*ni2-60*a2*b2);               //(6,10)
     answer.at(11,22)=0;                                                                 //(6,11)
     answer.at(11,23)=((a2*b*ni-5*b3-a2*b)*t3*E)/(180*a*b2*ni2-180*a*b2);                //(6,12)


     answer.at(15,15)=((2*a2*b2*ni-10*b4-7*a2*b2-10*a4)*t3*E)/(30*a3*b3*ni2-30*a3*b3);    //(7,7)
     answer.at(15,16)=((4*b5*ni+b5+10*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));              //(7,8)
     answer.at(15,17)=((4*a2*b3*ni+10*b5+a2*b3)*t3*E)/(60*a2*b4*ni2-60*a2*b4);            //(7,9)
     answer.at(15,21)=-((2*a2*b2*ni-10*b4-7*a2*b2+5*a4)*t3*E)/(30*a3*b3*ni2-30*a3*b3);   //(7,10)
     answer.at(15,22)=-((4*b5*ni+b5-5*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));             //(7,11)
     answer.at(15,23)=((-a2*b*ni+10*b3+a2*b)*t3*E)/(60*a2*b2*ni2-60*a2*b2);              //(7,12)


     answer.at(16,16)=((2*b5*ni-2*b5-10*a2*b3)*t3*E)/(5*(18*a*b4*ni2-18*a*b4));           //(8,8)
     answer.at(16,17)=-(b3*ni*t3*E)/(2*(6*b3*ni2-6*b3));                                  //(8,9)
     answer.at(16,21)=-((4*b5*ni+b5-5*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));             //(8,10)
     answer.at(16,22)=-((2*b5*ni-2*b5+5*a2*b3)*t3*E)/(5*(18*a*b4*ni2-18*a*b4));          //(8,11)
     answer.at(16,23)=0;                                                                 //(8,12)


     answer.at(17,17)=-((-a2*b*ni+5*b3+a2*b)*t3*E)/(45*a*b2*ni2-45*a*b2);                 //(9,9)
     answer.at(17,21)=-((-a2*b*ni+10*b3+a2*b)*t3*E)/(60*a2*b2*ni2-60*a2*b2);             //(9,10)
     answer.at(17,22)=0;                                                                 //(9,11)
     answer.at(17,23)=-((a2*b*ni+10*b3-a2*b)*t3*E)/(180*a*b2*ni2-180*a*b2);              //(9,12)


     answer.at(21,21)=((2*a2*b2*ni-10*b4-7*a2*b2-10*a4)*t3*E)/(30*a3*b3*ni2-30*a3*b3);  //(10,10)
     answer.at(21,22)=((4*b5*ni+b5+10*a2*b3)*t3*E)/(10*(6*a*b5*ni2-6*a*b5));            //(10,11)
     answer.at(21,23)=-((4*a2*b3*ni+10*b5+a2*b3)*t3*E)/(60*a2*b4*ni2-60*a2*b4);         //(10,12)


     answer.at(22,22)=((2*b5*ni-2*b5-10*a2*b3)*t3*E)/(5*(18*a*b4*ni2-18*a*b4));         //(11,11)
     answer.at(22,23)=(b3*ni*t3*E)/(2*(6*b3*ni2-6*b3));                                 //(11,12)


     answer.at(23,23)=-((-a2*b*ni+5*b3+a2*b)*t3*E)/(45*a*b2*ni2-45*a*b2);               //(12,12)

// extra stiffness, R_w
          double max = answer.at(1,1);

          for (int i=1;i<24;i++){
              if (answer.at(i,i)>max){
                  max=answer.at(i,i);
              }
              K= max/10000;
          }

     answer.at(6,6)= K;
     answer.at(12,12)= K;
     answer.at(18,18)= K;
     answer.at(24,24)= K;

     answer.symmetrized();

}




double
qshell :: computeLength()
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
qshell :: initializeFrom(InputRecord &ir)
{
    numberOfGaussPoints = 4;
    StructuralElement :: initializeFrom(ir);

    if ( numberOfGaussPoints != 1 && numberOfGaussPoints != 4 && numberOfGaussPoints != 9 && numberOfGaussPoints != 16 && numberOfGaussPoints != 25 ) {
        numberOfGaussPoints = 4;
        OOFEM_WARNING("Number of Gauss points enforced to 1");
    }
}



void
qshell :: computeLCS()
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
qshell :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    this->computeLCS();
    answer.resize(24, 24);
    answer.zero();

    for ( int i = 0; i <=3; i++ ) { // Loops over nodes
        // In each node, transform global c.s. {D_u, D_v, D_w, R_u, R_v, R_w} into local c.s.
        answer.setSubMatrix(this->lcsMatrix, 1 + i * 6, 1 + i * 6);     // Displacements
        answer.setSubMatrix(this->lcsMatrix, 1 + i * 6 + 3, 1 + i * 6 + 3); // Rotations
    }

    return true;
}

int
qshell :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [5,6]
// f(local) = T * f(global)
{
    this->computeLCS();
   // this->computeGtoLRotationMatrix();

    answer.resize(6,6);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = answer.at(4, i + 3) = lcsMatrix.at(1, i);
        answer.at(2, i) = answer.at(5, i + 3) = lcsMatrix.at(2, i);
        answer.at(3, i) = answer.at(6, i + 3) = lcsMatrix.at(3, i);
    }

    return 1;
}

void
qshell :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v,D_w,R_u,R_v,R_w
    };
}

void
qshell :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    StructuralElement :: updateLocalNumbering(f);
    if ( this->referenceNode ) {
        this->referenceNode = f(this->referenceNode, ERS_DofManager);
    }
}


} // end namespace oofem

