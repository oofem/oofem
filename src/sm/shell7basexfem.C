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

#include "shell7basexfem.h"
#include "shell7base.h"
#include "enrichmentitem.h"
#include "xfemmanager.h"
#include "constantpressureload.h"

namespace oofem {
    IntArray Shell7BaseXFEM :: dofId_Midplane(3);
    IntArray Shell7BaseXFEM :: dofId_Director(3);
    IntArray Shell7BaseXFEM :: dofId_InhomStrain(1); 
    bool Shell7BaseXFEM :: __initializedFieldDofId = Shell7BaseXFEM :: initDofId();
Shell7BaseXFEM :: Shell7BaseXFEM(int n, Domain *aDomain) : Shell7Base(n, aDomain), XfemElementInterface(this) 
{
    //this->xMan =  this->giveDomain()->giveXfemManager(1);
}

IRResultType Shell7BaseXFEM :: initializeFrom(InputRecord *ir)
{
    Shell7Base :: initializeFrom(ir);
    return IRRT_OK; 
};

Interface
*Shell7BaseXFEM :: giveInterface(InterfaceType it)
{
    if ( it != XfemElementInterfaceType ) {
        return Shell7Base :: giveInterface(it);
    } else if ( it == XfemElementInterfaceType ) {
        return ( XfemElementInterface * ) this;
    } else {
        return Shell7Base :: giveInterface(it);
    }
}


double 
Shell7BaseXFEM :: giveGlobalZcoord(GaussPoint *gp) 
{

    this->setupDelaminationXiCoordList();
    this->setupGPDelaminationGroupList();

    double xiRef = gp->giveCoordinate(3);
    int dGroup   = this->giveDelaminationGroupAt( xiRef );
    double xiMid = this->giveDelaminationGroupMidXi(dGroup);
    
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * > (this->element->giveCrossSection());
    
    return (xiRef - xiMid)*layeredCS->computeIntegralThick(); // new xi-coord measured from dGroup c.s. 
    
}




void
Shell7BaseXFEM :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    // Continuous part
    Shell7Base ::giveDofManDofIDMask(inode, ut, answer);

    // Discontinuous part
    DofManager *dMan = this->giveDofManager(inode);
    XfemManager *xMan =  this->giveDomain()->giveXfemManager(1);
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
        for ( int j = 1; j <= ei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( ei->isDofManEnrichedByEnrichmentDomain(dMan,j) ) {
                IntArray eiDofIdArray;
                ei->giveEIDofIdArray(eiDofIdArray, j);
                answer.followedBy(eiDofIdArray);
            }
        }
    }
}


void
Shell7BaseXFEM :: evalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, FloatArray &genEpsC)
{
    // Continuous part
    FloatArray g1c, g2c, g3c;
    Shell7Base :: evalCovarBaseVectorsAt(gp, g1, g2, g3, genEpsC);

    // Discontinuous part - ///@todo naive imlementation should be changed
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    xMan =  this->giveDomain()->giveXfemManager(1);
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) ); // should check success
        EnrichmentFunction *ef = dei->giveEnrichmentFunction(1);

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                double xi0 = dei->enrichmentDomainXiCoords.at(j);
                double H = dei->heaviside(gp->giveCoordinate(3), xi0);        
                
                FloatArray g1d_temp, g2d_temp, g3d_temp, dGenEps;
                computeDiscGeneralizedStrainComponents(dGenEps, gp, dei, j, tStep);
                Shell7Base :: evalCovarBaseVectorsAt(gp, g1d_temp, g2d_temp, g3d_temp, dGenEps);

                g1.add(H,g1d_temp); g2.add(H,g2d_temp); g3.add(H,g3d_temp);
            }
        }
    }
}


void
Shell7BaseXFEM :: computeDiscGeneralizedStrainComponents(FloatArray &answer, GaussPoint *gp, EnrichmentItem *ei, int enrichmentDomainNumber, TimeStep *tStep)
{
    FloatArray dSolVec;
    IntArray eiDofIdArray;
    ei->giveEIDofIdArray(eiDofIdArray, enrichmentDomainNumber);
    this->discGiveUpdatedSolutionVector(dSolVec, eiDofIdArray, tStep);
    FloatMatrix B11, B22, B32, B43, B53;
    Shell7Base :: computeBmatricesAt(gp, B11, B22, B32, B43, B53);
    Shell7Base :: computeGeneralizedStrainVector(answer, dSolVec, B11, B22, B32, B43, B53);
}



void
Shell7BaseXFEM :: discGiveUpdatedSolutionVector(FloatArray &answer, IntArray &eiDofIdArray, TimeStep *tStep)
{
    // Returns the solution vector of discontinuous dofs dx_d & dm_d
    temp_computeVectorOf(eiDofIdArray, VM_Total, tStep, answer);
    //eiDofIdArray.printYourself();
    //answer.printYourself();
}



int
Shell7BaseXFEM :: giveNumberOfDofs()
{
    // Continuous part
    int nDofs = Shell7Base ::giveNumberOfDofs();

    // Discontinuous part
    XfemManager *xMan =  this->giveDomain()->giveXfemManager(1);

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        for ( int j = 1; j <= xMan->giveNumberOfEnrichmentItems(); j++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(j);
            for ( int k = 1; k <= ei->giveNumberOfEnrichmentDomains(); k++ ) {
                if ( ei->isDofManEnrichedByEnrichmentDomain(dMan,k) ) {
                    IntArray eiDofIdArray;
                    ei->giveEIDofIdArray(eiDofIdArray, k);
                    nDofs += eiDofIdArray.giveSize();
                }
            }
        }
    }
    return nDofs;
}




void
Shell7BaseXFEM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
//
// Computes internal forces as a summation of: sectional forces + convective mass force
{

    // compute number of xDofs
    xMan =  this->giveDomain()->giveXfemManager(1);
//    EnrichmentItem *ei = xMan->giveEnrichmentItem(1);
//    int nXdofs = ei->giveNumberOfEnrichmentDomains() * ei->giveEnrichesDofsWithIdArray()->giveSize() * this->giveNumberOfDofManagers(); 
    

    answer.resize( this->giveNumberOfDofs() );
    answer.zero();
    
    FloatArray solVec;

    // Continuous part
    this->giveUpdatedSolutionVector(solVec, tStep);
    solVec.printYourself();
    this->computeSectionalForces(answer, tStep, solVec, useUpdatedGpRecord);

    // Disccontinuous part
    
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) ); // should check success

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                IntArray eiDofIdArray;
                dei->giveEIDofIdArray(eiDofIdArray, j);
                this->discGiveUpdatedSolutionVector(solVec, eiDofIdArray, tStep);
                double xi0 = dei->enrichmentDomainXiCoords.at(j);
                this->discComputeSectionalForces(answer, tStep, solVec, useUpdatedGpRecord, xi0, dei, j);                      

            }
        }
    }
    answer.printYourself();
}


void 
Shell7BaseXFEM :: computeOrderingArray( IntArray &orderingArray, IntArray &activeDofsArray,  EnrichmentItem *ei, int enrichmentDomainNumber, SolutionField field)
{
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value
    IntArray eiDofIdArray;
    ei->giveEIDofIdArray(eiDofIdArray, enrichmentDomainNumber);

    IntArray ordering_cont = this->giveOrdering(field);
    IntArray fieldDofId    = this->giveFieldDofId(field);

    IntArray ordering_temp, temp;
    ordering_temp.resize(ordering_cont.giveSize());
    temp.resize(ordering_cont.giveSize());

    int k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        int pos = (i-1)*fieldDofId.giveSize();
        int pos2 = (i-1)*eiDofIdArray.giveSize();
        if ( ei->isDofManEnrichedByEnrichmentDomain(dMan,enrichmentDomainNumber) ){

            for (int j = 1; j <= fieldDofId.giveSize(); j++ ) {
                if ( dMan->hasDofID( (DofIDItem) fieldDofId.at(j) ) ) {
                    k++;
                    ordering_temp.at(k) = ordering_cont.at(k);
                    temp.at(k) = j + pos;
                }
            }
        }
    }
      
    IntArray ordering; orderingArray.resize(k), activeDofsArray.resize(k);
    ///@todo will not work if there are several ei
    int shift = Shell7Base :: giveNumberOfDofs(); 
    for ( int i = 1; i <= k; i++ ) {
        orderingArray.at(i) = ordering_temp.at(i) + shift;
        activeDofsArray.at(i) = temp.at(i);
    }

}



void
Shell7BaseXFEM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    int ndofs = this->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    
    FloatMatrix temp;
    FloatArray solVec;
    // Continuous part
    this->giveUpdatedSolutionVector(solVec, tStep);
    this->computeBulkTangentMatrix(answer, solVec, rMode, tStep);

    // Disccontinuous part
    
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) ); // should check success

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                IntArray eiDofIdArray;
                dei->giveEIDofIdArray(eiDofIdArray, j);
                this->discGiveUpdatedSolutionVector(solVec, eiDofIdArray, tStep);
                double xi0 = dei->enrichmentDomainXiCoords.at(j);

                discComputeBulkTangentMatrix(answer, solVec, rMode, tStep, xi0, dei, j);

            }
        }
    }
}


void
Shell7BaseXFEM :: discComputeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, MatResponseMode rMode, TimeStep *tStep,
                    double xi0, EnrichmentItem *ei, int enrichmentDomainNumber)
{
    FloatMatrix A [ 3 ] [ 3 ];
    FloatMatrix F1 [ 3 ];
    FloatArray t1(3), t2(3), t3(3);
    FloatArray f3 [ 2 ], f4, f5;

    FloatMatrix L11(6, 6), L12(6, 6), L13(6, 3), L14(6, 2), L15(6, 1),
    L22(6, 6), L23(6, 3), L24(6, 2), L25(6, 1),
    L33(3, 3), L34(3, 2), L35(3, 1),
    L44(2, 2), L45(2, 1),
    L55(1, 1);
    FloatMatrix B11, B22, B32, B43, B53;
    FloatArray S1g(3), S2g(3), S3g(3), m(3), dm1(3), dm2(3), temp1, temp2;
    FloatMatrix K11(18, 18), K12(18, 18), K13(18, 6), K22(18, 18), K23(18, 6), K33(6, 6);
    K11.zero(), K12.zero(), K13.zero(), K22.zero(), K23.zero(), K33.zero();

    FloatMatrix K11temp1, K12temp1, K13temp1, K22temp1, K22temp2, K23temp1, K23temp2, K33temp1, K33temp2;

    //FloatArray solVec;
    //this->giveUpdatedSolutionVector(solVec, tStep);

    //int ndofs = this->giveNumberOfDofs();
    //answer.resize(ndofs, ndofs);
    //answer.zero();
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    int numberOfLayers = layeredCS->giveNumberOfLayers();     // conversion from double to int - fix!

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = layerIntegrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);

            double gam, dg1, dg2;

            this->computeBmatricesAt(gp, B11, B22, B32, B43, B53);
            FloatArray genEps;
            this->computeGeneralizedStrainVector(genEps, solVec, B11, B22, B32, B43, B53);
            this->giveGeneralizedStrainComponents(genEps, temp1, temp1, dm1, dm2, m, dg1, dg2, gam);


            // Material stiffness
            Shell7Base :: computeLinearizedStiffness(gp, mat, tStep, S1g, S2g, S3g, A, genEps);

            // Tangent stiffness
 #if 1
            // thickness coefficients
            double zeta = giveGlobalZcoord(gp);
            double a = zeta + 0.5 * gam * zeta * zeta;
            double b = 0.5 * zeta * zeta;
            double c = 1. + gam * zeta;

            // f1(alpha) = b*A(alpha,beta)*dg(beta) + c*A(alpha,3);
            F1 [ 0 ].resize(0, 0);
            F1 [ 1 ].resize(0, 0);
            F1 [ 2 ].resize(0, 0);
            F1 [ 0 ].add(dg1 * b, A [ 0 ] [ 0 ]);
            F1 [ 0 ].add(dg2 * b, A [ 0 ] [ 1 ]);
            F1 [ 0 ].add(c, A [ 0 ] [ 2 ]);
            F1 [ 1 ].add(dg1 * b, A [ 1 ] [ 0 ]);
            F1 [ 1 ].add(dg2 * b, A [ 1 ] [ 1 ]);
            F1 [ 1 ].add(c, A [ 1 ] [ 2 ]);
            F1 [ 2 ].add(dg1 * b, A [ 2 ] [ 0 ]);
            F1 [ 2 ].add(dg2 * b, A [ 2 ] [ 1 ]);
            F1 [ 2 ].add(c, A [ 2 ] [ 2 ]);

            //f3(alpha) = b*A(alpha,beta)*dm(beta) + zeta*A(alpha,3)*m;
            f3 [ 0 ].resize(0);
            f3 [ 1 ].resize(0);
            t1.beProductOf(A [ 0 ] [ 0 ], dm1);
            t2.beProductOf(A [ 0 ] [ 1 ], dm2);
            t3.beProductOf(A [ 0 ] [ 2 ], m);
            f3 [ 0 ].add(b, t1);
            f3 [ 0 ].add(b, t2);
            f3 [ 0 ].add(zeta, t3);

            //f32 = b*( A.at(2,1)*dm1 + A.at(2,2)*dm2 )  + zeta*A.at(2,3)*m;
            t1.beProductOf(A [ 1 ] [ 0 ], dm1);
            t2.beProductOf(A [ 1 ] [ 1 ], dm2);
            t3.beProductOf(A [ 1 ] [ 2 ], m);
            f3 [ 1 ].add(b, t1);
            f3 [ 1 ].add(b, t2);
            f3 [ 1 ].add(zeta, t3);


            //f4 = b*dm(alpha)*A(alpha,3) + zeta*A(3,3)*m;
            f4.resize(0);
            t1.beTProductOf(A [ 0 ] [ 2 ], dm1);
            t2.beTProductOf(A [ 1 ] [ 2 ], dm2);
            t3.beTProductOf(A [ 2 ] [ 2 ], m);
            f4.add(b, t1);
            f4.add(b, t2);
            f4.add(zeta, t3);

            // f5 = b*F1(alpha)*dm(alpha) + zeta*F2*m + zeta*N3
            f5.resize(0);
            t1.beTProductOf(F1 [ 0 ], dm1);
            t2.beTProductOf(F1 [ 1 ], dm2);
            t3.beTProductOf(F1 [ 2 ], m);
            f5.add(b, t1);
            f5.add(b, t2);
            f5.add(zeta, t3);
            f5.add(zeta, S3g);


            /*
             *
             *     A     a*A            F1           b*A*m             f3
             *    a*A    a^2*A          a*F1         a*b*A*m         a*f3+b*N
             *    F1     a*F1     b*F1*dgam+c*F2     b*(F1*m+N)        f5
             *    b*A*m  a*b*A*m    b*(F1*m+N)       b^2*m*A*m     b*m*f3
             *    f3     a*f3+b*N       f5           b*m*f3    b*dm*f3+zeta*m*f4
             */

            // L(1,1) = A(alpha,beta)
            L11.setSubMatrix(A [ 0 ] [ 0 ], 1, 1);
            L11.setSubMatrix(A [ 0 ] [ 1 ], 1, 4);
            L11.setSubMatrix(A [ 1 ] [ 0 ], 4, 1);
            L11.setSubMatrix(A [ 1 ] [ 1 ], 4, 4);

            // L(1,2) = a*A(alpha,beta)
            L12 = L11;
            L12.times(a);

            // L(1,3) = F1(alpha)
            L13.setSubMatrix(F1 [ 0 ], 1, 1);
            L13.setSubMatrix(F1 [ 1 ], 4, 1);

            // L(1,4) = b*A*m
            L14.zero();
            t1.beProductOf(A [ 0 ] [ 0 ], m);
            t2.beProductOf(A [ 0 ] [ 1 ], m);
            L14.addSubVectorCol(t1, 1, 1);
            L14.addSubVectorCol(t2, 1, 2);
            t1.beProductOf(A [ 1 ] [ 0 ], m);
            t2.beProductOf(A [ 1 ] [ 1 ], m);
            L14.addSubVectorCol(t1, 4, 1);
            L14.addSubVectorCol(t2, 4, 2);
            L14.times(b);

            // L(1,5) = f3(alpha)
            L15.zero();
            L15.addSubVectorCol(f3 [ 0 ], 1, 1);
            L15.addSubVectorCol(f3 [ 1 ], 4, 1);


            // L(2,2) = a^2*A(alpha,beta)= a*L12
            L22 = L12;
            L22.times(a);

            // L(2,3) = a*F1(alpha) = a*L13
            L23 = L13;
            L23.times(a);

            // L(2,4) = a*F1(alpha)=a*L14
            L24 = L14;
            L24.times(a);

            // L(2,5) = a*f3(alpha) + b*N(alpha)
            L25.zero();
            L25.addSubVectorCol(S1g, 1, 1);
            L25.addSubVectorCol(S2g, 4, 1);
            L25.times(b);
            L25.add(a, L15);


            // L(3,3) = b*F1(beta)*dgam(beta) + c*F2
            L33.resize(0, 0);
            L33.add(dg1 * b, F1 [ 0 ]);
            L33.add(dg2 * b, F1 [ 1 ]);
            L33.add(c, F1 [ 2 ]);


            // L(3,4) = b*( F1(beta)*m + N(beta) )
            t1.beTProductOf(F1 [ 0 ], m);
            t1.add(S1g);
            t1.times(b);
            t2.beTProductOf(F1 [ 1 ], m);
            t2.add(S2g);
            t2.times(b);
            L34.setColumn(t1, 1);
            L34.setColumn(t2, 2);

            // L(3,5) = f5
            L35.setColumn(f5, 1);

            // L(4,4) = b^2*m*A*m (2x2)
            t1.beProductOf(A [ 0 ] [ 0 ], m);
            L44.at(1, 1) = t1.dotProduct(m);
            t1.beProductOf(A [ 0 ] [ 1 ], m);
            L44.at(1, 2) = t1.dotProduct(m);
            t1.beProductOf(A [ 1 ] [ 0 ], m);
            L44.at(2, 1) = t1.dotProduct(m);
            t1.beProductOf(A [ 1 ] [ 1 ], m);
            L44.at(2, 2) = t1.dotProduct(m);
            L44.times(b * b);

            // L(4,5) = b*m*f3(alpha)
            L45.at(1, 1) = b * m.dotProduct(f3 [ 0 ]);
            L45.at(2, 1) = b * m.dotProduct(f3 [ 1 ]);

            // L(5,5) = b*m*f3(alpha)
            L55.at(1, 1) = b * dm1.dotProduct(f3 [ 0 ]) + b * dm2.dotProduct(f3 [ 1 ]) + zeta * m.dotProduct(f4);
 #endif
            // K11 = BT11*L11*B11
            // K11 = BT11*K11temp1
            K11temp1.beProductOf(L11, B11);

            // K12 = BT11*(L12*B22+L13*B32)
            // K12 = BT11*K12temp1
            K12temp1.beProductOf(L12, B22);
            K12temp1.addProductOf(L13, B32);

            // K13 = BT11*(L14*B43 + L15*B53)
            // K13 = BT11*K13temp1
            K13temp1.beProductOf(L14, B43);
            K13temp1.addProductOf(L15, B53);

            // K22 = BT22*(L22*B22 + L23*B32) + BT32*(L32*B22 + L33*B32)
            // K22 = BT22*K22temp1            + BT32*K22temp2
            K22temp1.beProductOf(L22, B22);
            K22temp1.addProductOf(L23, B32);
            K22temp2.beTProductOf(L23, B22);
            K22temp2.addProductOf(L33, B32);

            // K23 = BT22*(L24*B43 + L25*B53) + BT32*(L34*B43 + L35*B53)
            // K23 = BT22*K23temp1            + BT32*K23temp2
            K23temp1.beProductOf(L24, B43);
            K23temp1.addProductOf(L25, B53);
            K23temp2.beProductOf(L34, B43);
            K23temp2.addProductOf(L35, B53);
            // K23 = (BT22*L24 + BT32*L34)*B43 + (BT22*L25 + BT32*L35)*B53 // TODO: This order would be faster
            //K23temp1.TbeProductOf(B22,L24);
            //K23temp1.addTProductOf(B32,L34);
            //K23temp2.beTProductOf(L22,L25);
            //K23temp2.addTProductOf(B32,L35);


            // K33 = BT43*(L44*B43 + L45*B53) + BT53*(L54*B43 + L55*B53)
            // K33 = BT43*K33temp1            + BT53*K33temp2
            K33temp1.beProductOf(L44, B43);
            K33temp1.addProductOf(L45, B53);
            K33temp2.beTProductOf(L45, B43);
            K33temp2.add(L55.at(1, 1), B53);            // L55 is just a scalar (and K33temp2 is just an array)

            double dV = this->computeVolumeAroundLayer(gp, layer);

            K11.plusProductSymmUpper(B11, K11temp1, dV);
            K12.plusProductUnsym(B11, K12temp1, dV);
            K13.plusProductUnsym(B11, K13temp1, dV);

            K22.plusProductSymmUpper(B22, K22temp1, dV);
            K22.plusProductSymmUpper(B32, K22temp2, dV);

            K23.plusProductUnsym(B22, K23temp1, dV);
            K23.plusProductUnsym(B32, K23temp2, dV);
            // Potential minor optimization, actually use FloatArray for anything that only has a single column
            //tmpA.beProductOf(K23temp1, B43); K23.add(dV, tmpA);
            //K23.plusDyadUnsymm(K33temp2, B53, dV);

            K33.plusProductSymmUpper(B43, K33temp1, dV);
            K33.plusProductSymmUpper(B53, K33temp2, dV);
            //K33.plusDyadSymmUpper(B53, K33temp2, dV); // Potential minor optimization (probably not worth it)
        }
    }

    K11.symmetrized();
    K22.symmetrized();
    K33.symmetrized();

    IntArray ordering_phibar, ordering_m, ordering_gam;
    IntArray activeDofs_phibar, activeDofs_m, activeDofs_gam;

    computeOrderingArray(ordering_phibar, activeDofs_phibar, ei, enrichmentDomainNumber, Midplane);
    computeOrderingArray(ordering_m, activeDofs_m, ei, enrichmentDomainNumber, Director);
    computeOrderingArray(ordering_gam, activeDofs_gam, ei, enrichmentDomainNumber, InhomStrain);


    FloatMatrix K11Ass, K12Ass, K13Ass, K22Ass, K23Ass, K33Ass, mat1, mat2(3,3);
    IntArray index;
    mat2.at(1,1) = 1.; mat2.at(1,2) = 2.; mat2.at(1,3) = 3.;
    mat2.at(2,1) = 4.; mat2.at(2,2) = 5.; mat2.at(2,3) = 6.;
    mat2.at(3,1) = 7.; mat2.at(3,2) = 8.; mat2.at(3,3) = 9.;
    index.setValues(3, 1, 3, 5);
    mat1.beSubMatrixOfSizeOf(mat2, index, 6);
    mat2.printYourself();
    index.printYourself();
    mat1.printYourself();
    
    mat1.beSubMatrixOf(mat2, index);
    mat1.printYourself();


    answer.assemble(K11, ordering_phibar, ordering_phibar);
    answer.assemble(K12, ordering_phibar, ordering_m);
    answer.assemble(K13, ordering_phibar, ordering_gam);
    answer.assemble(K22, ordering_m,      ordering_m);
    answer.assemble(K23, ordering_m,      ordering_gam);
    answer.assemble(K33, ordering_gam,    ordering_gam);

    FloatMatrix K21, K31, K32;
    K21.beTranspositionOf(K12);
    K31.beTranspositionOf(K13);
    K32.beTranspositionOf(K23);
    answer.assemble(K21, ordering_m,      ordering_phibar);
    answer.assemble(K31, ordering_gam,    ordering_phibar);
    answer.assemble(K32, ordering_gam,    ordering_m);
}



void
Shell7BaseXFEM :: discComputeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord, 
                    double xi0, EnrichmentItem *ei, int enrichmentDomainNumber)
//
{
    FloatMatrix B;
    FloatArray BtF, f, genEps;
    //answer.resize( this->giveNumberOfDofs() );
    //answer.zero();

    
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    int numberOfLayers = layeredCS->giveNumberOfLayers();     // conversion of types
    FloatArray f1(18), f2(18), f3(6);
    f1.zero();
    f2.zero();
    f3.zero();
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = layerIntegrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->getNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);

            if ( gp->giveCoordinate(3) > xi0 ) { // Should be enriched ///@todo not general!

                FloatMatrix B11, B22, B32, B43, B53;
                this->computeBmatricesAt(gp, B11, B22, B32, B43, B53);
                this->computeGeneralizedStrainVector(genEps, solVec, B11, B22, B32, B43, B53);

                double zeta = giveGlobalZcoord(gp);
                FloatArray N, M, T, Ms;
                double Ts = 0.;
                this->computeSectionalForcesAt(N, M, T, Ms, Ts, gp, mat, tStep, genEps, zeta);

                // Computation of sectional forces: f = B^t*[N M T Ms Ts]^t
                FloatArray f1temp(18), f2temp(18), f3temp(6), temp;
                // f1 = BT11*N
                f1temp.beTProductOf(B11, N);

                // f2 = BT22*M + BT32*T
                f2temp.beTProductOf(B22, M);
                temp.beTProductOf(B32, T);
                f2temp.add(temp);

                // f3 = BT43*Ms + BT53*Ts
                f3temp.beTProductOf(B43, Ms);
                for ( int i = 1; i <= 6; i++ ) {
                    f3temp.at(i) += B53.at(1, i) * Ts;
                }

                double dV = this->computeVolumeAroundLayer(gp, layer);
                f1.add(dV, f1temp);
                f2.add(dV, f2temp);
                f3.add(dV, f3temp);
            
            }

        }
    }

    // Should assemble to xfem dofs
    IntArray ordering_phibar, ordering_m, ordering_gam;
    IntArray activeDofs_phibar, activeDofs_m, activeDofs_gam;

    computeOrderingArray(ordering_phibar, activeDofs_phibar, ei, enrichmentDomainNumber, Midplane);
    computeOrderingArray(ordering_m, activeDofs_m, ei, enrichmentDomainNumber, Director);
    computeOrderingArray(ordering_gam, activeDofs_gam, ei, enrichmentDomainNumber, InhomStrain);

/*
    ordering_phibar.printYourself();
    ordering_m.printYourself();
    ordering_gam.printYourself();

    activeDofs_phibar.printYourself();
    activeDofs_m.printYourself();
    activeDofs_gam.printYourself();
    */


    FloatArray f1Ass, f2Ass, f3Ass;
    f1Ass.beSubArrayOf(f1,activeDofs_phibar);
    f2Ass.beSubArrayOf(f2,activeDofs_m);
    f3Ass.beSubArrayOf(f3,activeDofs_gam);
    answer.assemble(f1Ass, ordering_phibar);
    answer.assemble(f2Ass, ordering_m);
    answer.assemble(f3Ass, ordering_gam);

}






void
Shell7BaseXFEM :: computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep) {
    // Num refers in this case to  numerical integration in both in-plane and through the thickness.
    // For analytically integrated throught he thickness, see computeMassMatrix


    FloatMatrix N, Nt, Ntm, NtmN, mass, temp;
    FloatArray solVec, unknowns;
    this->giveUpdatedSolutionVector(solVec, tStep);
    int ndofs = this->giveNumberOfDofs();
    temp.resize(ndofs, ndofs);
    temp.zero();


    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    int numberOfLayers = layeredCS->giveNumberOfLayers();     // conversion of data

    FloatMatrix M11(18, 18), M12(18, 18), M13(18, 6), M22(18, 18), M23(18, 6), M33(6, 6);
    M11.zero();
    M12.zero();
    M13.zero();
    M22.zero();
    M23.zero();
    M33.zero();

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = layerIntegrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->getNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);

            FloatMatrix N11, N22, N33;
            this->computeNmatricesAt(gp, N11, N22, N33);
            FloatArray xbar, m;
            double gam = 0.;
            this->computeSolutionFields(xbar, m, gam, solVec, N11, N22, N33);
            //this->computeNmatrixAt(gp, N);
            //unknowns.beProductOf(N,a); // [xbar, m, gam]^T
            //m.setValues(3, unknowns.at(4), unknowns.at(5), unknowns.at(6) );
            //double gam = unknowns.at(7);


            /* Consistent mass matrix M = int{N^t*mass*N}
             *
             *         3    3    1
             *         3 [a*I  b*I   c*m      [A  B  C
             * mass =   3       d*I   e*m    =     D  E
             *         1  sym       f*m.m]     sym   F]
             */


            double zeta = giveGlobalZcoord(gp);
            double fac1 = 4;
            double fac2 = 2.0 * zeta * ( 2.0 + gam * zeta );
            double fac3 = 2.0 * zeta * zeta;
            double fac4 = zeta * zeta * ( 2.0 + gam * zeta ) * ( 2.0 + gam * zeta );
            double fac5 = zeta * zeta * zeta * ( 2.0 + gam * zeta );
            double fac6 = zeta * zeta * zeta * zeta;
            FloatMatrix mass11(3, 3), mass12(3, 3), mass13(3, 1), mass21(3, 3), mass22(3, 3), mass23(3, 1), mass31(1, 3), mass32(1, 3), mass33(1, 1);
            mass11.zero();
            mass12.zero();
            mass13.zero();
            mass21.zero();
            mass22.zero();
            mass23.zero();
            mass31.zero();
            mass32.zero();
            mass33.zero();
            mass.resize(7, 7);
            mass11.at(1, 1) = mass11.at(2, 2) = mass11.at(3, 3) = fac1;          // A
            mass12.at(1, 1) = mass12.at(2, 2) = mass12.at(3, 3) = fac2;          // B
            mass13.at(1, 1) = fac3 * m.at(1);
            mass13.at(2, 1) = fac3 * m.at(2);
            mass13.at(3, 1) = fac3 * m.at(3);            // C
            mass22.at(1, 1) = mass22.at(2, 2) = mass22.at(3, 3) = fac4;          // D
            mass23.at(1, 1) = fac5 * m.at(1);
            mass23.at(2, 1) = fac5 * m.at(2);
            mass23.at(3, 1) = fac5 * m.at(3);            // E
            mass33.at(1, 1) = fac6 * m.dotProduct(m);            // F
            mass21.beTranspositionOf(mass12);
            mass31.beTranspositionOf(mass13);
            mass32.beTranspositionOf(mass23);
            //mass.symmetrized();

            double dV = this->computeVolumeAroundLayer(gp, layer);
            double rho = mat->give('d', gp);

            FloatMatrix M11temp, M12temp, M13temp, M22temp, M23temp, M33temp;
            this->computeTripleProduct(M11temp, N11, mass11, N11);
            this->computeTripleProduct(M12temp, N11, mass12, N22);
            this->computeTripleProduct(M13temp, N11, mass13, N33);
            this->computeTripleProduct(M22temp, N22, mass22, N22);
            this->computeTripleProduct(M23temp, N22, mass23, N33);
            this->computeTripleProduct(M33temp, N33, mass33, N33);
            M11.add(0.25 * rho * dV, M11temp);
            M12.add(0.25 * rho * dV, M12temp);
            M13.add(0.25 * rho * dV, M13temp);
            M22.add(0.25 * rho * dV, M22temp);
            M23.add(0.25 * rho * dV, M23temp);
            M33.add(0.25 * rho * dV, M33temp);
        }

        //M33.printYourself();
        answer.resize(ndofs, ndofs);
        answer.zero();

        IntArray ordering_phibar = giveOrdering(Midplane);
        IntArray ordering_m = giveOrdering(Director);
        IntArray ordering_gam = giveOrdering(InhomStrain);
        answer.assemble(M11, ordering_phibar, ordering_phibar);
        answer.assemble(M12, ordering_phibar, ordering_m);
        answer.assemble(M13, ordering_phibar, ordering_gam);
        answer.assemble(M22, ordering_m,      ordering_m);
        answer.assemble(M23, ordering_m,      ordering_gam);
        answer.assemble(M33, ordering_gam,    ordering_gam);

        FloatMatrix M21, M31, M32;
        M21.beTranspositionOf(M12);
        M31.beTranspositionOf(M13);
        M32.beTranspositionOf(M23);
        answer.assemble(M21, ordering_m,      ordering_phibar);
        answer.assemble(M31, ordering_gam,    ordering_phibar);
        answer.assemble(M32, ordering_gam,    ordering_m);
        answer.symmetrized();

    }
}






# if 0
void
Shell7BaseXFEM :: discEvalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1d, FloatArray &g2d, FloatArray &g3d, FloatArray &dGenEps)
{

    FloatArray lcoords = * gp->giveCoordinates();
    double zeta = giveGlobalZcoord(gp);

    FloatArray dxdxi1, dxdxi2, dmdxi1, dmdxi2, m;
    this->discGiveGeneralizedStrainComponents(dGenEps, dxdxi1, dxdxi2, dmdxi1, dmdxi2, m);

    g1d = dxdxi1;
    g1d.add(zeta, dmdxi1);
    g2d = dxdxi2;
    g2d.add(zeta, dmdxi2);
    g3d = m;
    
}



void
Shell7BaseXFEM :: discComputeGeneralizedStrainVector(FloatArray &answer, const FloatArray &solVec, const FloatMatrix &B11,
                                                     const FloatMatrix &B22, const FloatMatrix &B32) {
    answer.resize(15);
    answer.zero();
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    int ndofs_gam = this->giveNumberOfFieldDofs(InhomStrain);
    for ( int i = 1; i <= ndofs_xm; i++ ) {
        for ( int j = 1; j <= 6; j++ ) {
            answer.at(j)  += B11.at(j, i) * solVec.at(i);            // dx/dxi
            answer.at(6 + j)  += B22.at(j, i) * solVec.at(i + ndofs_xm);      // dm/dxi
        }
    }

    for ( int i = 1; i <= ndofs_xm; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(12 + j) += B32.at(j, i) * solVec.at(i + ndofs_xm);      // m
        }
    }

}


void
Shell7BaseXFEM :: discGiveGeneralizedStrainComponents(FloatArray &genEps, FloatArray &dphidxi1, FloatArray &dphidxi2, FloatArray &dmdxi1, 
         FloatArray &dmdxi2, FloatArray &m) {
    // generealized strain vector for discontinuous part  [dxdxi, dmdxi, m]^T
    dphidxi1.setValues( 3, genEps.at(1), genEps.at(2), genEps.at(3) );
    dphidxi2.setValues( 3, genEps.at(4), genEps.at(5), genEps.at(6) );
    dmdxi1.setValues( 3, genEps.at(7), genEps.at(8), genEps.at(9) );
    dmdxi2.setValues( 3, genEps.at(10), genEps.at(11), genEps.at(12) );
    m.setValues( 3, genEps.at(13), genEps.at(14), genEps.at(15) );
}
#endif









IntArray
Shell7BaseXFEM :: giveFieldDofId(SolutionField fieldType) const {
    if ( fieldType == Midplane ) {
        return this->dofId_Midplane;
    } else if ( fieldType == Director  )   {
        return this->dofId_Director;
    } else if ( fieldType == InhomStrain  )   {
        return this->dofId_InhomStrain;
    } else {
        _error("giveOrdering: unknown fieldType");
    }
}





// Delamination specific

#if 1

void
Shell7BaseXFEM :: setupDelaminationXiCoordList()
{
    if ( this->delaminationXiCoordList.size()==0 ) {
    // Stores a paired list with the EnrichmentDomain# and the corresponding xi-coord of the delamination.
    xMan =  this->giveDomain()->giveXfemManager(1); // When does one need several xfemman?
    int numEI = xMan->giveNumberOfEnrichmentItems();
    for ( int i = 1; i <= numEI; i++ ) {
        Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) );
        if ( dei ) {
            int numED = dei->giveNumberOfEnrichmentDomains(); // numEnrDomains max possible number
            int pos = 1;
            for ( int j = 1; j <= numED; j++ ) {
                if( dei->isElementEnrichedByEnrichmentDomain(this->element, j) ) { 
                    std::pair<int, double> pid;
                    pid.first  = pos;
                    pid.second = dei->enrichmentDomainXiCoords.at(j); 
                    this->delaminationXiCoordList.push_back(pid); 
                    pos++;
                }
            }
        } 
    }

    // Sort xi-coords in acending order.
    this->delaminationXiCoordList.sort(sortFunc);
    
    }
}

void 
Shell7BaseXFEM :: setupGPDelaminationGroupList() 
{   // Creates a list wich stores the gp# and the dGroup# for quick access later.
  if ( this->gpDelaminationGroupList.size()==0 ) {
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    int numberOfLayers = layeredCS->giveNumberOfLayers();  

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = layerIntegrationRulesArray [ layer - 1 ];

        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);
            std::pair<int, int> pid;
            pid.first  = gp->giveNumber();
            pid.second = giveDelaminationGroupAt( gp->giveCoordinate(3));
            this->gpDelaminationGroupList.push_back(pid); 
        }
    }
  }
}


int
Shell7BaseXFEM :: giveDelaminationGroupAt(double xi) 
{   
    // Starts ordering from 0
    std::list< std::pair<int, double> >::const_iterator iter;
    iter = this->delaminationXiCoordList.begin();

    int nDelam = this->delaminationXiCoordList.size();
    for ( int j = 1; j <= nDelam; j++ ) {
        double xiDelam = (*iter).second;
        if ( xi < xiDelam ) { //belong to the delamination group just below delamination #j. How to deal with poins that lie onthe boundary?
            return j-1;            
        }
        iter++;
    }
    return nDelam;
            
}


void 
Shell7BaseXFEM :: giveDelaminationGroupXiLimits(int &dGroup, double &xiTop, double &xiBottom)
{
    
    int nDelam = this->delaminationXiCoordList.size();
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * > (this->element->giveCrossSection());
    
    std::list< std::pair<int, double> >::const_iterator iter;
    iter = this->delaminationXiCoordList.begin(); 

    if ( dGroup == 0 ) {
        xiBottom = - layeredCS->giveMidSurfaceXiCoordFromBottom();
        xiTop = (*iter).second;
    } else if (dGroup == nDelam) {
        std::advance(iter, dGroup-1);
        xiBottom = (*iter).second;
        xiTop    = -layeredCS->giveMidSurfaceXiCoordFromBottom() + 2.0;
    } else {
        std::advance(iter, dGroup-1);
        xiBottom = (*iter).second;
        iter++;
        xiTop    = (*iter).second;
    }

#if DEBUG
    if ( xiBottom > xiTop ) {
        OOFEM_ERROR2("giveDelaminationGroupZLimits: Bottom xi-coord is larger than top xi-coord in dGroup. (%i)", dGroup);
    }
#endif
}


double 
Shell7BaseXFEM :: giveDelaminationGroupMidXi(int dGroup)
{
    double xiTop=0., xiBottom=0.;
    this->giveDelaminationGroupXiLimits(dGroup, xiTop, xiBottom);
    return 0.5 * ( xiTop + xiBottom );
}

#endif














void 
Shell7BaseXFEM :: setupDelaminationXiCoordsAtGP() 
{
    std::pair<int, double> pid;
    std::list<std::pair<int, double> > *delaminationXiCoordList;

    xMan = this->giveDomain()->giveXfemManager(1);
    int numEI = xMan->giveNumberOfEnrichmentItems();
    for ( int i = 1; i <= numEI; i++ ) {
        Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) );
        if ( dei ) {
            if ( dei->isElementEnriched(this) ) {


                int nDelam = dei->giveNumberOfEnrichmentDomains(); // numEnrDomains max possible number
                int pos = 1;
                for ( int j = 1; j <= nDelam; j++ ) {
                    //if( this->isElementEnriched(element) ) {
                    //pid.first  = pos;
                    //pid.second = this->delaminationZCoords.at(i); 
                    //xiCoordList->push_back(pid); 
                    //pos++;
                }
            }
        } 
    }
}








} // end namespace oofem

