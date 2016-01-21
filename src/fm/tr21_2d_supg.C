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

#include "tr21_2d_supg.h"
#include "fei2dtrquad.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fluiddynamicmaterial.h"
#include "fluidcrosssection.h"
#include "timestep.h"
#include "contextioerr.h"
#include "crosssection.h"
#include "classfactory.h"
#include "engngm.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Element(TR21_2D_SUPG);

FEI2dTrQuad TR21_2D_SUPG :: velocityInterpolation(1, 2);
FEI2dTrLin TR21_2D_SUPG :: pressureInterpolation(1, 2);


TR21_2D_SUPG :: TR21_2D_SUPG(int n, Domain *aDomain) :
    SUPGElement2(n, aDomain), ZZNodalRecoveryModelInterface(this)
    // Constructor.
{
    numberOfDofMans  = 6;
}

TR21_2D_SUPG :: ~TR21_2D_SUPG()
// Destructor
{ }

FEInterpolation *
TR21_2D_SUPG :: giveInterpolation() const
{
    return & this->velocityInterpolation;
}

FEInterpolation *
TR21_2D_SUPG :: giveInterpolation(DofIDItem id) const
{
    if ( id == P_f ) {
        return & this->pressureInterpolation;
    } else {
        return & this->velocityInterpolation;
    }
}

int
TR21_2D_SUPG :: computeNumberOfDofs()
{
    return 15;
}

void
TR21_2D_SUPG :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode < 4 ) {
        answer = {V_u, V_v, P_f};
    } else {
        answer = {V_u, V_v};
    }
}

void
TR21_2D_SUPG :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(3);

        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], 3, this);

        //seven point Gauss integration
        integrationRulesArray [ 1 ].reset( new GaussIntegrationRule(2, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 1 ], 7, this);

        integrationRulesArray [ 2 ].reset( new GaussIntegrationRule(3, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 2 ], 13, this);


        //integrationRulesArray [ 3 ] = new GaussIntegrationRule(4, this, 1, 3);
        //this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[3], 27, this );
    }
}


void
TR21_2D_SUPG :: computeNuMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray n;

    this->velocityInterpolation.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    
    answer.beNMatrixOf(n, 2);
}

void
TR21_2D_SUPG :: computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix n, dn;
    FloatArray u, un;
    this->velocityInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    this->computeNuMatrix(n, gp);
    this->computeVectorOfVelocities(VM_Total, tStep, un);

    u.beProductOf(n, un);

    answer.resize(2, 12);
    answer.zero();
    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1, 2 * i - 1) = dn.at(i, 1) * u.at(1) + dn.at(i, 2) * u.at(2);
        answer.at(2, 2 * i)   = dn.at(i, 1) * u.at(1) + dn.at(i, 2) * u.at(2);
    }
}

void
TR21_2D_SUPG :: computeBMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dn;
    this->velocityInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(3, 12);
    answer.zero();

    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1, 2 * i - 1) = dn.at(i, 1);
        answer.at(2, 2 * i)   = dn.at(i, 2);
        answer.at(3, 2 * i - 1) = dn.at(i, 2);
        answer.at(3, 2 * i)   = dn.at(i, 1);
    }
}

void
TR21_2D_SUPG :: computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dn;
    velocityInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(1, 12);
    answer.zero();

    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1, 2 * i - 1) = dn.at(i, 1);
        answer.at(1, 2 * i) = dn.at(i, 2);
    }
}

void
TR21_2D_SUPG :: computeNpMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray n;
    this->pressureInterpolation.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(1, 3);
    answer.zero();

    answer.at(1, 1)  = n.at(1);
    answer.at(1, 2)  = n.at(2);
    answer.at(1, 3)  = n.at(3);
}


void
TR21_2D_SUPG :: computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    FloatMatrix dn, um(2, 6);
    this->computeVectorOfVelocities(VM_Total, tStep, u);

    velocityInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int i = 1; i <= 6; i++ ) {
        um.at(1, i) = u.at(2 * i - 1);
        um.at(2, i) = u.at(2 * i);
    }

    answer.beProductOf(um, dn);
}

void
TR21_2D_SUPG :: computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dn;
    pressureInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.beTranspositionOf(dn);
}



void
TR21_2D_SUPG :: computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix D, d2n;

    answer.resize(2, 12);
    answer.zero();


    this->velocityInterpolation.evald2Ndx2( d2n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial()->giveDeviatoricStiffnessMatrix(D, TangentStiffness, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1, 2 * i - 1) = D.at(1, 1) * d2n.at(i, 1) + D.at(1, 3) * d2n.at(i, 3) + D.at(3, 1) * d2n.at(i, 3) + D.at(3, 3) * d2n.at(i, 2);
        answer.at(1, 2 * i)     = D.at(1, 2) * d2n.at(i, 3) + D.at(1, 3) * d2n.at(i, 1) + D.at(3, 2) * d2n.at(i, 2) + D.at(3, 3) * d2n.at(i, 3);
        answer.at(2, 2 * i - 1) = D.at(2, 1) * d2n.at(i, 3) + D.at(2, 3) * d2n.at(i, 2) + D.at(3, 1) * d2n.at(i, 1) + D.at(3, 3) * d2n.at(i, 3);
        answer.at(2, 2 * i)     = D.at(2, 2) * d2n.at(i, 2) + D.at(2, 3) * d2n.at(i, 3) + D.at(3, 2) * d2n.at(i, 3) + D.at(3, 3) * d2n.at(i, 1);
    }
}


void
TR21_2D_SUPG :: updateStabilizationCoeffs(TimeStep *tStep)
{
    double mu_min, norm_N, norm_N_d, norm_M_d, norm_LSIC;
    FloatMatrix N, N_d, M_d, LSIC;
    FloatArray u;
    mu_min = 1;
    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {
        double mu = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial()->giveEffectiveViscosity(gp, tStep);
        if ( mu_min > mu ) {
            mu_min = mu;
        }
    }

    //this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), u);

    this->computeAdvectionTerm(N, tStep);
    this->computeAdvectionDeltaTerm(N_d, tStep);
    this->computeMassDeltaTerm(M_d, tStep);
    this->computeLSICTerm(LSIC, tStep);

    norm_N = N.computeFrobeniusNorm();
    norm_N_d = N_d.computeFrobeniusNorm();
    norm_M_d = M_d.computeFrobeniusNorm();
    norm_LSIC = LSIC.computeFrobeniusNorm();

    if ( ( norm_N == 0 ) || ( norm_N_d == 0 ) || ( norm_M_d == 0 ) ) {
        t_supg = 0;
    } else {
        //t_supg =  1. / sqrt( 1. / ( t_s1 * t_s1 ) + 1. / ( t_s2 * t_s2 ) + 1. / ( t_s3 * t_s3 ) );
        t_supg = 0;
    }

    if ( norm_LSIC == 0 ) {
        t_lsic = 0;
    } else {
        //t_lsic = norm_N / norm_LSIC;

        t_lsic = 0;
    }

    t_pspg = 0;
}



void
TR21_2D_SUPG :: computeAdvectionTerm(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix n, b;

    answer.clear();

    /* consistent part + supg stabilization term */
    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix(b, gp, tStep);
        double dV  = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
        answer.plusProductUnsym(n, b, rho * dV);
    }
}


void
TR21_2D_SUPG :: computeAdvectionDeltaTerm(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix n, b;

    answer.clear();

    /* consistent part + supg stabilization term */
    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix(b, gp, tStep);
        double dV  = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);

        answer.plusProductUnsym(b, b, rho * dV);
    }
}



void
TR21_2D_SUPG :: computeMassDeltaTerm(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix n, b;

    answer.clear();

    /* mtrx for computing t_supg, norm of this mtrx is computed */
    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix(b, gp, tStep);
        double dV  = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);

        answer.plusProductUnsym(b, n, rho * dV);
    }
}

void
TR21_2D_SUPG :: computeLSICTerm(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix b;

    answer.clear();

    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {
        double dV  = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
        this->computeDivUMatrix(b, gp);

        answer.plusProductSymmUpper(b, b, dV * rho);
    }

    answer.symmetrized();
}


int
TR21_2D_SUPG :: giveNumberOfSpatialDimensions()
{
    return 2;
}


double
TR21_2D_SUPG :: computeCriticalTimeStep(TimeStep *tStep)
{
    return 1.e6;
}


double
TR21_2D_SUPG :: LS_PCS_computeF(LevelSetPCS *ls, TimeStep *tStep)
{
    double answer = 0.0, vol = 0.0;
    FloatMatrix n, dn;
    FloatArray fi(6), u, un, gfi;

    this->computeVectorOfVelocities(VM_Total, tStep, un);

    for ( int i = 1; i <= 6; i++ ) {
        fi.at(i) = ls->giveLevelSetDofManValue( dofManArray.at(i) );
    }

    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {
        double dV = this->computeVolumeAround(gp);
        velocityInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        this->computeNuMatrix(n, gp);
        u.beProductOf(n, un);
        gfi.beTProductOf(dn, fi);

        vol += dV;
        answer += dV * u.dotProduct(gfi) / gfi.computeNorm();
    }

    return answer / vol;
}

void
TR21_2D_SUPG :: LS_PCS_computedN(FloatMatrix &answer)
{
    FloatMatrix dn;

    answer.clear();

    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {

        velocityInterpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        answer.add(dn); ///@todo This code makes no sense to me. Without the weight and jacobian this is dependant on the number of gauss points.
    }
}

void
TR21_2D_SUPG :: LS_PCS_computeVolume(double &answer, const FloatArray **coordinates)
{
    //double answer = 0.0;
    answer = 0.0;

    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {
        //answer += this->computeVolumeAround(gp);

        double determinant, weight, volume;

        determinant = fabs( this->velocityInterpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );

        weight = gp->giveWeight();
        volume = determinant * weight;

        answer += volume;
    }
}


double
TR21_2D_SUPG :: LS_PCS_computeVolume()
{
    double answer = 0.0;

    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {
        answer += this->computeVolumeAround(gp);
    }

    return answer;
}

double
TR21_2D_SUPG :: LS_PCS_computeS(LevelSetPCS *ls, TimeStep *tStep)
{
    FloatArray fi(6), un, n;

    double vol = 0.0, eps = 0.0, _fi, S = 0.0;

    for ( int i = 1; i <= 6; i++ ) {
        fi.at(i) = ls->giveLevelSetDofManValue( dofManArray.at(i) );
    }

    for ( GaussPoint *gp: *integrationRulesArray [ 1 ] ) {
        double dV = this->computeVolumeAround(gp);
        velocityInterpolation.evalN( n,  gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        vol += dV;
        _fi = n.dotProduct(fi);
        S +=  _fi / ( _fi * _fi + eps * eps ) * dV;
    }

    return S / vol;
}



void
TR21_2D_SUPG :: LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi)
{
    int neg = 0, pos = 0, zero = 0, si = 0, sqi = 0;


    answer.resize(2);

    for ( double ifi: fi ) {  //comparing values of fi in vertices
        if ( ifi > 0. ) {
            pos++;
        } else if ( ifi < 0.0 ) {
            neg++;
        } else {
            zero++;
        }
    }

    //control !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int first_control = 0;
    first_control = this->giveNumber();
    OOFEM_LOG_INFO("TR21_2D_SUPG :: LS_PCS_computeVOFFractions - First control, sign(fi) in vertices, element no. %d", first_control);



    if ( neg == 0 ) {  // all level set values positive
        answer.at(1) = 1.0;
        answer.at(2) = 0.0;
    } else if ( pos == 0 ) {  // all level set values negative
        answer.at(1) = 0.0;
        answer.at(2) = 1.0;
    } else if ( zero == 3 ) {
        // ???????
        answer.at(1) = 1.0;
        answer.at(2) = 0.0;
    } else {
        // zero level set inside
        // three main cases of crossing the triangle are possible

        int inter_case, negq = 0, posq = 0, zeroq = 0; //inter_case is variable that differs type of which fi crosses the triangle

        for ( double ifi: fi ) { //loop over edge nodes
            if ( ifi > 0. ) {
                posq++;
            } else if ( ifi < 0.0 ) {
                negq++;
            } else {
                zeroq++;
            }
        }

        for ( int i = 1; i <= 3; i++ ) { //loop over vertex nodes, searching for which vertex has different value of fi
            if ( ( pos > neg ) && ( fi.at(i) < 0.0 ) ) {
                si = i;
                break;
            }

            if ( ( pos < neg ) && ( fi.at(i) >= 0.0 ) ) {
                si = i;
                break;
            }
        }



        for ( int i = 4; i <= 6; i++ ) { //loop over edge nodes, searching for which edge node has different value of fi
            if ( ( posq > negq ) && ( fi.at(i) < 0.0 ) ) {
                sqi = i;
                break;
            }

            if ( ( posq < negq ) && ( fi.at(i) >= 0.0 ) ) {
                sqi = i;
                break;
            }
        }

        if ( negq == 0 ) {
            inter_case = 11; //only one node has different level set value from other five
        } else if ( posq == 0 ) {
            inter_case = 12; //only one node has different level set value from other five
        } else {       // level set devides element 2:4 or 3:3 nodes
            if ( fi.at(si) * fi.at(sqi) > 0 ) {
                inter_case = 2;
            } else {
                inter_case = 3;
            }
        }


        int edge1 = 0, edge2 = 0;
        FloatArray helplcoords(3);


        if ( inter_case > 3 ) { //only one node (vertex in this case) of triangle has diffenert level set value, than others
            FloatArray inter1(2), inter2(2);

            //kontrola!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            int second_control1;
            second_control1 = this->giveNumber();
            OOFEM_LOG_INFO("TR21_2D_SUPG :: LS_PCS_computeVOFFractions - case 1 - first type of element deviation by LS, element no. %d", second_control1);




            if ( si == 1 ) { // computing intersection points in order to vertex with different sign of level set funct
                this->computeIntersection(1, inter1, fi);


                this->computeIntersection(3, inter2, fi);

                edge1 = 1;
                edge2 = 3;
            } else if ( si == 2 ) {
                this->computeIntersection(2, inter1, fi);


                this->computeIntersection(1, inter2, fi);

                edge1 = 2;
                edge2 = 1;
            } else if ( si == 3 ) {
                this->computeIntersection(3, inter1, fi);


                this->computeIntersection(2, inter2, fi);

                edge1 = 3;
                edge2 = 2;
            }




            OOFEM_LOG_INFO("case 1 - after intersections of LS and edges, element no. %d", second_control1);



            //computing point on zero level set curve: [xM, yM]
            // this point lies on line from the "si" vertex, the condition is fi([xM, yM]) = 0, M = [xM, yM]

            FloatArray _l(4), M(2), X_si(2), _Mid(2), line(6), _X1(2), Mid1(2), Mid2(2), Coeff(3), loc_Mid(3), loc_X1(3), N_Mid, N_X1;
            double x1, xsi, y1, ysi, t, fi_X1, fi_Mid, r1, r11, r12;

            _l.at(1) = inter1.at(1);
            _l.at(2) = inter1.at(2);
            _l.at(3) = inter2.at(1);
            _l.at(4) = inter2.at(2);

            this->computeCenterOf(_Mid, _l, 1);

            xsi = this->giveNode(si)->giveCoordinate(1);
            ysi = this->giveNode(si)->giveCoordinate(2);

            X_si.at(1) = xsi;
            X_si.at(2) = ysi;

            x1 = xsi + 2 * ( _Mid.at(1) - xsi );
            y1 = ysi + 2 * ( _Mid.at(2) - ysi );

            _X1.at(1) = x1;

            _X1.at(2) = y1;

            this->velocityInterpolation.global2local( loc_Mid, _Mid, FEIElementGeometryWrapper(this) );
            this->velocityInterpolation.global2local( loc_X1, _X1, FEIElementGeometryWrapper(this) );

            this->velocityInterpolation.evalN( N_Mid, loc_Mid, FEIElementGeometryWrapper(this) );
            this->velocityInterpolation.evalN( N_X1, loc_X1, FEIElementGeometryWrapper(this) );

            fi_Mid = N_Mid.dotProduct(fi);
            fi_X1 = N_X1.dotProduct(fi);

            line.at(1) = 0;
            line.at(2) = fi.at(si);
            line.at(3) = sqrt( ( _Mid.at(1) - xsi ) * ( _Mid.at(1) - xsi ) + ( _Mid.at(2) - ysi ) * ( _Mid.at(2) - ysi ) );
            line.at(4) = fi_Mid;
            line.at(5) = sqrt( ( _Mid.at(1) - _X1.at(1) ) * ( _Mid.at(1) - _X1.at(1) ) + ( _Mid.at(2) - _X1.at(2) ) * ( _Mid.at(2) - _X1.at(2) ) );
            line.at(6) = fi_X1;

            this->computeQuadraticFunct(Coeff, line);

            this->computeQuadraticRoots(Coeff, r11, r12);

            if ( r11 > line.at(6) || r11 < line.at(1) ) {
                r1 = r12;
            } else {
                r1 = r11;
            }

            t = r1 / sqrt( ( _Mid.at(1) - xsi ) * ( _Mid.at(1) - xsi ) + ( _Mid.at(2) - ysi ) * ( _Mid.at(2) - ysi ) );

            M.at(1) = xsi + t * ( _Mid.at(1) - xsi );
            M.at(2) = ysi + t * ( _Mid.at(2) - ysi );


            OOFEM_LOG_INFO("case 1 - after computing third point on zero level set curve inside element, element no. %d", second_control1);



            //computing middle points to get six-point triangle
            FloatArray borderpoints(4);

            borderpoints.at(1) = xsi;
            borderpoints.at(2) = ysi;
            borderpoints.at(3) = inter1.at(1);
            borderpoints.at(4) = inter1.at(2);

            this->computeMiddlePointOnParabolicArc(Mid1, edge1, borderpoints);

            borderpoints.at(1) = xsi;
            borderpoints.at(2) = ysi;
            borderpoints.at(3) = inter2.at(1);
            borderpoints.at(4) = inter2.at(2);


            this->computeMiddlePointOnParabolicArc(Mid2, edge2, borderpoints);


            const FloatArray *c1 [ 6 ];




            double vol_1, vol;

            c1 [ 0 ] = & X_si;
            c1 [ 1 ] = & inter1;
            c1 [ 2 ] = & inter2;
            c1 [ 3 ] = & Mid1;
            c1 [ 4 ] = & M;
            c1 [ 5 ] = & Mid2;


            //volume of small triangle
            this->LS_PCS_computeVolume(vol_1, c1);
            //volume of whole triangle
            vol = LS_PCS_computeVolume();




            if ( inter_case == 11 ) {
                answer.at(2) = vol_1 / vol;
                answer.at(1) = 1.0 - answer.at(2);
            } else {
                answer.at(1) = vol_1 / vol;
                answer.at(2) = 1.0 - answer.at(1);
            }
        } //end case inter_case > 3
        else if ( inter_case == 2 ) {
            //kontrola!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            int second_control2 = 0;
            second_control2 = this->giveNumber();
            OOFEM_LOG_INFO("case 2 - second type of element deviation by LS, element no. %d", second_control2);


            FloatArray inter1(2), inter2(2), crosssect(4);
            crosssect.zero();

            if ( si == 1 ) { // computing intersection points in order to vertex with different sign of level set funct
                this->computeIntersection(1, inter1, fi);

                this->computeIntersection(3, inter2, fi);

                edge1 = 1;
                edge2 = 3;
            } else if ( si == 2 ) {
                this->computeIntersection(2, inter1, fi);

                this->computeIntersection(1, inter2, fi);

                edge1 = 2;
                edge2 = 1;
            } else if ( si == 3 ) {
                this->computeIntersection(3, inter1, fi);

                this->computeIntersection(2, inter2, fi);

                edge1 = 3;
                edge2 = 2;
            }

            //computing point on zero level set curve: [xM, yM]
            // this point lies on line from the "si" vertex, the condition is fi([xM, yM]) = 0
            FloatArray X_qsi(2), X_si(2), _Mid(2), line(6), _X1(2), Mid1(2), Coeff(3), loc_Mid, loc_X1, N_Mid, N_X1, M(2);
            double x1, xsi, y1, ysi, t, fi_X1, fi_Mid, r1, r11, r12;
            this->computeCenterOf(_Mid, crosssect, 1);

            xsi = this->giveNode(si)->giveCoordinate(1);
            ysi = this->giveNode(si)->giveCoordinate(2);

            X_si.at(1) = xsi;
            X_si.at(2) = ysi;

            X_qsi.at(1) = this->giveNode(sqi)->giveCoordinate(1);
            X_qsi.at(2) = this->giveNode(sqi)->giveCoordinate(2);


            x1 = xsi + 1.3 * ( _Mid.at(1) - xsi );
            y1 = ysi + 1.3 * ( _Mid.at(2) - ysi );

            _X1.at(1) = x1;
            _X1.at(2) = y1;

            this->velocityInterpolation.global2local( loc_Mid, _Mid, FEIElementGeometryWrapper(this) );
            this->velocityInterpolation.global2local( loc_X1, _X1, FEIElementGeometryWrapper(this) );

            this->velocityInterpolation.evalN( N_Mid, loc_Mid, FEIElementGeometryWrapper(this) );
            this->velocityInterpolation.evalN( N_X1, loc_X1, FEIElementGeometryWrapper(this) );

            fi_Mid = N_Mid.dotProduct(fi);
            fi_X1 = N_X1.dotProduct(fi);

            line.at(1) = 0;
            line.at(2) = fi.at(si);
            line.at(3) = sqrt( ( _Mid.at(1) - xsi ) * ( _Mid.at(1) - xsi ) + ( _Mid.at(2) - ysi ) * ( _Mid.at(2) - ysi ) );
            line.at(4) = fi_Mid;
            line.at(5) = sqrt( ( _Mid.at(1) - _X1.at(1) ) * ( _Mid.at(1) - _X1.at(1) ) + ( _Mid.at(2) - _X1.at(2) ) * ( _Mid.at(2) - _X1.at(2) ) );
            line.at(6) = fi_X1;

            this->computeQuadraticFunct(Coeff, line);

            this->computeQuadraticRoots(Coeff, r11, r12);

            if ( r11 > line.at(6) || r11 < line.at(1) ) {
                r1 = r12;
            } else {
                r1 = r11;
            }

            t = r1 / sqrt( ( _Mid.at(1) - xsi ) * ( _Mid.at(1) - xsi ) + ( _Mid.at(2) - ysi ) * ( _Mid.at(2) - ysi ) );

            M.at(1) = xsi + t * ( _Mid.at(1) - xsi );
            M.at(2) = ysi + t * ( _Mid.at(2) - ysi );

            //computing middle points to get six-point triangle
            FloatArray borderpoints(4);
            int sub_case = 0;
            if ( ( si == 1 ) && ( sqi == 4 ) ) {
                borderpoints.at(1) = xsi;
                borderpoints.at(2) = ysi;
                borderpoints.at(3) = inter2.at(1);
                borderpoints.at(4) = inter2.at(2);
                sub_case = 2;

                this->computeMiddlePointOnParabolicArc(Mid1, edge2, borderpoints);
            }

            if ( ( si == 1 ) && ( sqi == 6 ) ) {
                borderpoints.at(1) = xsi;
                borderpoints.at(2) = ysi;
                borderpoints.at(3) = inter1.at(1);
                borderpoints.at(4) = inter1.at(2);
                sub_case = 1;
                this->computeMiddlePointOnParabolicArc(Mid1, edge1, borderpoints);
            }

            if ( ( si == 2 ) && ( sqi == 4 ) ) {
                borderpoints.at(1) = xsi;
                borderpoints.at(2) = ysi;
                borderpoints.at(3) = inter1.at(1);
                borderpoints.at(4) = inter1.at(2);
                sub_case = 1;
                this->computeMiddlePointOnParabolicArc(Mid1, edge1, borderpoints);
            }

            if ( ( si == 2 ) && ( sqi == 5 ) ) {
                borderpoints.at(1) = xsi;
                borderpoints.at(2) = ysi;
                borderpoints.at(3) = inter2.at(1);
                borderpoints.at(4) = inter2.at(2);
                sub_case = 2;
                this->computeMiddlePointOnParabolicArc(Mid1, edge2, borderpoints);
            }

            if ( ( si == 3 ) && ( sqi == 6 ) ) {
                borderpoints.at(1) = xsi;
                borderpoints.at(2) = ysi;
                borderpoints.at(3) = inter2.at(1);
                borderpoints.at(4) = inter2.at(2);
                sub_case = 2;
                this->computeMiddlePointOnParabolicArc(Mid1, edge2, borderpoints);
            }

            if ( ( si == 3 ) && ( sqi == 5 ) ) {
                borderpoints.at(1) = xsi;
                borderpoints.at(2) = ysi;
                borderpoints.at(3) = inter1.at(1);
                borderpoints.at(4) = inter1.at(2);
                sub_case = 1;
                this->computeMiddlePointOnParabolicArc(Mid1, edge1, borderpoints);
            }

            const FloatArray *c1 [ 6 ];



            double vol_1, vol;



            if ( sub_case == 1 ) {
                c1 [ 0 ] = & X_si;
                c1 [ 1 ] = & inter1;
                c1 [ 2 ] = & inter2;
                c1 [ 3 ] = & Mid1;
                c1 [ 4 ] = & M;
                c1 [ 5 ] = & X_qsi;
            } else {
                c1 [ 0 ] = & X_si;
                c1 [ 1 ] = & inter1;
                c1 [ 2 ] = & inter2;
                c1 [ 3 ] = & X_qsi;
                c1 [ 4 ] = & M;
                c1 [ 5 ] = & Mid1;



                //volume of small triangle
                this->LS_PCS_computeVolume(vol_1, c1);
                //volume of whole triangle
                vol = LS_PCS_computeVolume();



                if ( fi(si) < 0 ) {
                    answer.at(2) = vol_1 / vol;
                    answer.at(1) = 1.0 - answer.at(2);
                } else {
                    answer.at(1) = vol_1 / vol;
                    answer.at(2) = 1.0 - answer.at(1);
                }
            } //end case inter_case == 2
        } else {  //inter_case == 3
            //kontrola!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            int second_control3 = 0;
            second_control3 = this->giveNumber();
            OOFEM_LOG_INFO("case 3 - third type of element deviation by LS, element no. %d", second_control3);


            FloatArray inter1(2), inter2(2), crosssect(4);
            crosssect.zero();

            if ( si == 1 ) { // computing intersection points in order to vertex with different sign of level set funct
                this->computeIntersection(1, inter1, fi);

                this->computeIntersection(3, inter2, fi);

                //edge1 = 1;
                //edge2 = 3;
            } else if ( si == 2 ) {
                this->computeIntersection(2, inter1, fi);


                this->computeIntersection(1, inter2, fi);

                //edge1 = 2;
                //edge2 = 1;
            } else if ( si == 3 ) {
                this->computeIntersection(3, inter1, fi);


                this->computeIntersection(2, inter2, fi);

                //edge1 = 3;
                //edge2 = 2;
            }

            //computing point on zero level set curve: [xM, yM]
            // this point lies on line from the "si" vertex, the condition is fi([xM, yM]) = 0
            FloatArray X_si(2), M(2), _Mid(2), line(6), _X1(2), Mid1(2), Coeff(3), loc_Mid, loc_X1, N_Mid, N_X1;
            double x1, xsi, y1, ysi, t, fi_X1, fi_Mid, r1, r11, r12;
            this->computeCenterOf(_Mid, crosssect, 1);

            xsi = this->giveNode(si)->giveCoordinate(1);
            ysi = this->giveNode(si)->giveCoordinate(2);

            X_si.at(1) = xsi;
            X_si.at(2) = ysi;

            x1 = xsi + 0.5 * ( _Mid.at(1) - xsi );
            y1 = ysi + 0.5 * ( _Mid.at(2) - ysi );

            _X1.at(1) = x1;
            _X1.at(2) = y1;

            this->velocityInterpolation.global2local( loc_Mid, _Mid, FEIElementGeometryWrapper(this) );
            this->velocityInterpolation.global2local( loc_X1, _X1, FEIElementGeometryWrapper(this) );

            this->velocityInterpolation.evalN( N_Mid, loc_Mid, FEIElementGeometryWrapper(this) );
            this->velocityInterpolation.evalN( N_X1, loc_X1, FEIElementGeometryWrapper(this) );

            fi_Mid = N_Mid.dotProduct(fi);
            fi_X1 = N_X1.dotProduct(fi);

            line.at(1) = 0;
            line.at(2) = fi.at(si);
            line.at(3) = sqrt( ( _Mid.at(1) - xsi ) * ( _Mid.at(1) - xsi ) + ( _Mid.at(2) - ysi ) * ( _Mid.at(2) - ysi ) );
            line.at(4) = fi_Mid;
            line.at(5) = sqrt( ( _Mid.at(1) - _X1.at(1) ) * ( _Mid.at(1) - _X1.at(1) ) + ( _Mid.at(2) - _X1.at(2) ) * ( _Mid.at(2) - _X1.at(2) ) );
            line.at(6) = fi_X1;

            this->computeQuadraticFunct(Coeff, line);

            this->computeQuadraticRoots(Coeff, r11, r12);

            if ( r11 > line.at(6) || r11 < line.at(1) ) {
                r1 = r12;
            } else {
                r1 = r11;
            }

            t = r1 / sqrt( ( _Mid.at(1) - xsi ) * ( _Mid.at(1) - xsi ) + ( _Mid.at(2) - ysi ) * ( _Mid.at(2) - ysi ) );

            M.at(1) = xsi + t * ( _Mid.at(1) - xsi );
            M.at(2) = ysi + t * ( _Mid.at(2) - ysi );


            FloatArray X_e1(2), X_e2(2); //points on edges close to the vertex si, "quadratic" nodes
            double vol_1, vol;

            X_e1.at(1) = this->giveNode(si + 3)->giveCoordinate(1);
            X_e1.at(2) = this->giveNode(si + 3)->giveCoordinate(2);

            X_e2.at(1) = this->giveNode( ( ( si + 4 ) % 3 ) + 4 )->giveCoordinate(1);
            X_e2.at(2) = this->giveNode( ( ( si + 4 ) % 3 ) + 4 )->giveCoordinate(2);




            const FloatArray *c1 [ 6 ];


            c1 [ 0 ] = & X_si;
            c1 [ 1 ] = & inter1;
            c1 [ 2 ] = & inter2;
            c1 [ 3 ] = & X_e1;
            c1 [ 4 ] = & M;
            c1 [ 5 ] = & X_e2;


            //volume of small triangle
            this->LS_PCS_computeVolume(vol_1, c1);
            //volume of whole triangle
            vol = LS_PCS_computeVolume();




            if ( fi(si) < 0 ) {
                answer.at(2) = vol_1 / vol;
                answer.at(1) = 1.0 - answer.at(2);
            } else {
                answer.at(1) = vol_1 / vol;
                answer.at(2) = 1.0 - answer.at(1);
            }
        }
    }
}

void
TR21_2D_SUPG :: computeIntersection(int iedge, FloatArray &intcoords, FloatArray &fi)
{
    FloatArray Coeff(3), helplcoords(3);
    double fi1, fi2, fi3, r1, r11, r12;
    IntArray edge(3);
    intcoords.resize(2);
    intcoords.zero();

    this->velocityInterpolation.computeLocalEdgeMapping(edge, iedge);
    fi1 = fi.at( edge.at(1) );
    fi2 = fi.at( edge.at(2) );
    fi3 = fi.at( edge.at(3) );

    Coeff.at(1) = fi1 + fi2 - 2 * fi3;
    Coeff.at(2) = fi2 - fi1;
    Coeff.at(3) = 2 * fi3;

    this->computeQuadraticRoots(Coeff, r11, r12);

    if ( r11 > 1.0 || r11 < -1.0 ) {
        r1 = r12;
    } else {
        r1 = r11;
    }

    helplcoords.zero();
    helplcoords.at(1) = r1;

    this->velocityInterpolation.edgeLocal2global( intcoords, 3, helplcoords, FEIElementGeometryWrapper(this) );

    //this->velocityInterpolation.evaldNdx(dn, this->giveDomain(), dofManArray, gp->giveNaturalCoordinates(),tStep->giveTime());
}


#if 0
void
TR21_2D_SUPG :: computeIntersection(int iedge, FloatArray &intcoords, FloatArray &fi)
{
    FloatArray Coeff(3), helplcoords(3);
    double fi1, fi2, fi3, r1, r11, r12;

    intcoords.resize(2);
    intcoords.zero();

    this->velocityInterpolation.computeLocalEdgeMapping(edge, iedge);
    fi1 = fi.at( edge.at(1) );
    fi2 = fi.at( edge.at(2) );
    fi3 = fi.at( edge.at(3) );

    Coeff.at(1) = fi1 + fi2 - 2 * fi3;
    Coeff.at(2) = fi2 - fi1;
    Coeff.at(3) = 2 * fi3;

    this->computeQuadraticRoots(Coeff, r11, r12);

    if ( r11 > 1.0 || r11 < -1.0 ) {
        r1 = r12;
    } else {
        r1 = r11;
    }

    helplcoords.zero();
    helplcoords.at(1) = r1;

    this->velocityInterpolation.edgeLocal2global( intcoords, iedge, this->giveDomain(), dofManArray, helplcoords, tStep->giveTime() );

    //this->velocityInterpolation.evaldNdx(dn, this->giveDomain(), dofManArray, gp->giveNaturalCoordinates(),tStep->giveTime());
}
#endif

void
TR21_2D_SUPG :: computeMiddlePointOnParabolicArc(FloatArray &answer, int iedge, FloatArray borderpoints)
{
    double a, b, c;
    FloatArray Coeff, C(2);




    this->computeQuadraticFunct(Coeff, iedge);


    this->computeCenterOf(C, borderpoints, 1);

    a = Coeff.at(1);
    b = Coeff.at(2);
    c = Coeff.at(3);

    answer.at(1) = C.at(1);
    answer.at(2) = a * C.at(1) * C.at(1) + b *C.at(1) + c;
}




void
TR21_2D_SUPG :: computeCenterOf(FloatArray &C, FloatArray c, int dim)
{
    switch ( dim ) {
    case 1:

        C.at(1) = ( c.at(1) + c.at(3) ) / 2;
        C.at(2) = ( c.at(2) + c.at(4) ) / 2;

        break;

    case 2:

        C.at(1) = ( c.at(1) + c.at(3) + c.at(5) ) / 3;
        C.at(2) = ( c.at(2) + c.at(4) + c.at(6) ) / 3;

        break;
    }
}


void
TR21_2D_SUPG :: computeQuadraticRoots(FloatArray Coeff, double &r1, double &r2)
{
    double a, b, c;

    a = Coeff.at(1);
    b = Coeff.at(2);
    c = Coeff.at(3);

    r1 = ( -b + sqrt(b * b - 4 * a * c) ) / ( 2 * a );
    r2 = ( -b - sqrt(b * b - 4 * a * c) ) / ( 2 * a );
}



void
TR21_2D_SUPG :: computeCoordsOfEdge(FloatArray &answer, int iedge)

{
    IntArray edge;

    velocityInterpolation.computeLocalEdgeMapping(edge, iedge);

    answer.at(1) = this->giveNode( edge.at(1) )->giveCoordinate(1);
    answer.at(2) = this->giveNode( edge.at(1) )->giveCoordinate(2);

    answer.at(3) = this->giveNode( edge.at(2) )->giveCoordinate(1);
    answer.at(4) = this->giveNode( edge.at(2) )->giveCoordinate(2);

    answer.at(5) = this->giveNode( edge.at(3) )->giveCoordinate(1);
    answer.at(6) = this->giveNode( edge.at(3) )->giveCoordinate(2);
}



//this function computes coefficients of quadratic function a,b,c in y = a*x^2 + b*x + c
//iedge is number of i-th edge of quadratic triangle
void
TR21_2D_SUPG :: computeQuadraticFunct(FloatArray &answer, int iedge)
{
    double x1, y1, x2, y2, x3, y3, detA, detA1;
    FloatMatrix A(3, 3), A1(3, 3);
    FloatArray edge, crds(6);


    answer.resize(3);



    this->computeCoordsOfEdge(crds, iedge);


    x1 = crds.at(1);
    y1 = crds.at(2);
    x2 = crds.at(5);
    y2 = crds.at(6);
    x3 = crds.at(3);
    y3 = crds.at(4);

    A.at(1, 1) = x1 * x1;
    A.at(2, 1) = x2 * x2;
    A.at(3, 1) = x3 * x3;
    A.at(1, 2) = x1;
    A.at(2, 2) = x2;
    A.at(3, 2) = x3;
    A.at(1, 3) = 1.0;
    A.at(2, 3) = 1.0;
    A.at(3, 3) = 1.0;

    FloatArray b(3);

    detA = A.giveDeterminant();

    b.at(1) = y1;
    b.at(2) = y2;
    b.at(3) = y3;

    for ( int k = 1; k <= 3; k++ ) {
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                A1.at(i, j) = A.at(i, j);
            }

            A1.at(i, k) = b.at(i);
        }

        detA1 = A1.giveDeterminant();
        answer.at(k) = detA1 / detA;
    }
}

void
TR21_2D_SUPG :: computeQuadraticFunct(FloatArray &answer, FloatArray line)
{
    double x1, y1, x2, y2, x3, y3, detA, detA1;
    FloatMatrix A(3, 3), A1(3, 3);
    FloatArray edge;


    answer.resize(3);


    x1 = line.at(1);
    y1 = line.at(2);
    x2 = line.at(3);
    y2 = line.at(4);
    x3 = line.at(5);
    y3 = line.at(6);

    A.at(1, 1) = x1 * x1;
    A.at(2, 1) = x2 * x2;
    A.at(3, 1) = x3 * x3;
    A.at(1, 2) = x1;
    A.at(2, 2) = x2;
    A.at(3, 2) = x3;
    A.at(1, 3) = 1.0;
    A.at(2, 3) = 1.0;
    A.at(3, 3) = 1.0;

    FloatArray b(3);

    detA = A.giveDeterminant();

    b.at(1) = y1;
    b.at(2) = y2;
    b.at(3) = y3;

    for ( int k = 1; k <= 3; k++ ) {
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                A1.at(i, j) = A.at(i, j);
            }

            A1.at(i, k) = b.at(i);
        }

        detA1 = A1.giveDeterminant();
        answer.at(k) = detA1 / detA;
    }
}


void
TR21_2D_SUPG :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                           InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveIPValue(answer, gp, type, tStep);
}


void
TR21_2D_SUPG :: initGeometry()
{ }


int
TR21_2D_SUPG :: checkConsistency()
{
    return SUPGElement :: checkConsistency();
}


void
TR21_2D_SUPG :: updateYourself(TimeStep *tStep)
{
    SUPGElement :: updateYourself(tStep);
}

int
TR21_2D_SUPG :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    return SUPGElement2 :: giveIPValue(answer, gp, type, tStep);
}

contextIOResultType TR21_2D_SUPG :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = SUPGElement :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}



contextIOResultType TR21_2D_SUPG :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = SUPGElement :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


double
TR21_2D_SUPG :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;

    determinant = fabs( this->velocityInterpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );


    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}


Interface *
TR21_2D_SUPG :: giveInterface(InterfaceType interface)
{
    if ( interface == LevelSetPCSElementInterfaceType ) {
        return static_cast< LevelSetPCSElementInterface * >(this);
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    return NULL;
}


void
TR21_2D_SUPG :: giveLocalVelocityDofMap(IntArray &map)
{
    map = {1, 2, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15};
}

void
TR21_2D_SUPG :: giveLocalPressureDofMap(IntArray &map)
{
    map = {3, 6, 9};
}





#ifdef __OOFEG
int
TR21_2D_SUPG :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int node, TimeStep *tStep)
{
    return SUPGElement :: giveInternalStateAtNode(answer, type, mode, node, tStep);
}



void
TR21_2D_SUPG :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void TR21_2D_SUPG :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    FloatArray v1, v2, v3;
    double s [ 3 ];

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    // if ((gc.giveIntVarMode() == ISM_local) && (gc.giveIntVarType() ==  IST_VOFFraction)) {
    if ( ( gc.giveIntVarType() ==  IST_VOFFraction ) && ( gc.giveIntVarMode() == ISM_local ) ) {
        Polygon matvolpoly;
        //this->formMaterialVolumePoly(matvolpoly, NULL, temp_normal, temp_p, false);
        EASValsSetColor( gc.getStandardSparseProfileColor() );
        //GraphicObj *go = matvolpoly.draw(gc,true,OOFEG_VARPLOT_PATTERN_LAYER);
        matvolpoly.draw(gc, true, OOFEG_VARPLOT_PATTERN_LAYER);
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, gc.giveIntVarType(), gc.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, gc.giveIntVarType(), gc.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, gc.giveIntVarType(), gc.giveIntVarMode(), 3, tStep);
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        result += giveIPValue(v1, gp, gc.giveIntVarType(), tStep);
        v2 = v1;
        v3 = v1;
        result *= 3;
    }

    if ( result != 3 ) {
        return;
    }

    indx = gc.giveIntVarIndx();

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = 0.;
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        gc.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    } else if ( ( gc.getScalarAlgo() == SA_ZPROFILE ) || ( gc.getScalarAlgo() == SA_COLORZPROFILE ) ) {
        double landScale = gc.getLandScale();

        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = s [ i ] * landScale;
        }

        if ( gc.getScalarAlgo() == SA_ZPROFILE ) {
            EASValsSetColor( gc.getDeformedElementColor() );
            EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangle3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
        } else {
            gc.updateFringeTableMinMax(s, 3);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
            EGWithMaskChangeAttributes(FILL_MASK | LAYER_MASK, tr);
        }

        EMAddGraphicsToModel(ESIModel(), tr);
    }
}



#endif
} // end namespace oofem
