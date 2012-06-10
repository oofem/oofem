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

#include "homogenize.h"
#include "mathfem.h"

namespace oofem {
Homogenize :: Homogenize()
{ }

//parallel configuration of Voigt
void Homogenize :: voigt(FloatMatrix &PhaseMatrix) {
    double k, mu;
    double volTot = 0.;
    int r;
    int NumPhases = PhaseMatrix.giveNumberOfRows();

    checkVolFraction(PhaseMatrix);

    k_hmg = 0.;
    mu_hmg = 0.;
    for ( r = 0; r < NumPhases; r++ ) {
        ENuToKMu(PhaseMatrix(r, 1), PhaseMatrix(r, 2), k, mu);
        k_hmg += k * PhaseMatrix(r, 0);
        mu_hmg += mu * PhaseMatrix(r, 0);
        volTot += PhaseMatrix(r, 0);
    }

    kMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);
}


//serial configuration of Reuss
void Homogenize :: reuss(FloatMatrix &PhaseMatrix) {
    double k, mu;
    int r;
    int NumPhases = PhaseMatrix.giveNumberOfRows();

    checkVolFraction(PhaseMatrix);

    k_hmg = 0.;
    mu_hmg = 0.;
    for ( r = 0; r < NumPhases; r++ ) {
        ENuToKMu(PhaseMatrix(r, 1), PhaseMatrix(r, 2), k, mu);
        k_hmg += PhaseMatrix(r, 0) / k;
        mu_hmg += PhaseMatrix(r, 0) / mu;
    }

    k_hmg = 1. / k_hmg;
    mu_hmg = 1. / mu_hmg;
    kMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);
}

/*
 * Hashin-Shtrikman-Walpole lower and upper bounds, arbitrary number of phases
 * If phases satisfy ascending order of bulk or shear moduli simultaneously, i.e. k_0<k_1<...k_N together
 * with mu_0<mu_1<...mu_N, Hashin-Shtrikman bounds are obtained
 * If the order of bulk and shear moduli is not satisfied, wider Hashin-Shtrikman-Walpole bounds are computed
 * Implementation according to J.G. Berryman: Mixture Theories for Rock Properties, Physics and Phase Relations,
 * A Handbook of Physical Constants, Am. Geophys, 1995 or download directly from the web
 */
void Homogenize :: hashinShtrikmanWalpole(FloatMatrix &PhaseMatrix) {
    int phase, i, r;
    int phaseMuMin, phaseMuMax;
    int counter = 0;
    int numRows = PhaseMatrix.giveNumberOfRows();
    double k, k_min, mu, mu_phase, muMin, muMax;
    double dummy;
    FloatMatrix PhaseMatrix1;
    PhaseMatrix1 = PhaseMatrix;
    FloatMatrix SortedPhaseMatrix(numRows, 3);


    checkVolFraction(PhaseMatrix1);

    //Sort phases in ascending order of bulk moduli
    for ( i = 0; i < numRows; i++ ) {
        k_min = 1.e+40;
        for ( r = 0; r < numRows; r++ ) {
            if ( PhaseMatrix1(r, 0) > 0. + 1.e-40 ) {
                ENuToKMu(PhaseMatrix1(r, 1), PhaseMatrix1(r, 2), k, mu);
                if ( k < k_min ) {
                    phase = r;
                    k_min = k;
                    mu_phase = mu;
                }
            }
        }

        SortedPhaseMatrix(counter, 0) = PhaseMatrix1(phase, 0);
        SortedPhaseMatrix(counter, 1) = k_min;
        SortedPhaseMatrix(counter, 2) = mu_phase;

        PhaseMatrix1(phase, 0) = 0.;
        counter++;
    }

    //SortedPhaseMatrix contains (vol. fraction, bulk modulus, shear modulus), sorted in ascending order of bulk modulus
    //Find phase numbers for maximum and minimum shear moduli
    muMin = 1.e+40;
    muMax = 0.;
    for ( i = 0; i < numRows; i++ ) {
        if ( SortedPhaseMatrix(i, 2) > muMax ) {
            phaseMuMax = i;
            muMax = SortedPhaseMatrix(i, 2);
        }

        if ( SortedPhaseMatrix(i, 2) < muMin ) {
            phaseMuMin = i;
            muMin = SortedPhaseMatrix(i, 2);
        }
    }

    // find bounds for bulk moduli
    k_hmg = 0; //lower bound
    k_hmg_2 = 0; //upper bound
    for ( i = 0; i < numRows; i++ ) {
        k_hmg += SortedPhaseMatrix(i, 0) / ( SortedPhaseMatrix(i, 1) + 4. / 3. * muMin );
        k_hmg_2 += SortedPhaseMatrix(i, 0) / ( SortedPhaseMatrix(i, 1) + 4. / 3. * muMax );
    }

    k_hmg = 1. / k_hmg;
    k_hmg -= 4. / 3. * muMin;
    k_hmg_2 = 1. / k_hmg_2;
    k_hmg_2 -= 4. / 3. * muMax;

    //find bounds for shear moduli
    mu_hmg = gamma( SortedPhaseMatrix, zeta(SortedPhaseMatrix(0, 1), muMin) ); //lower bound
    mu_hmg_2 = gamma( SortedPhaseMatrix, zeta(SortedPhaseMatrix(numRows - 1, 1), muMax) ); //upper bound

    //Young's moduli are determined from k_min, mu_min and k_max an mu_max
    kMuToENu(k_hmg, mu_hmg, E_hmg, dummy);
    kMuToENu(k_hmg_2, mu_hmg_2, E_hmg_2, dummy);
    //Poisson's ratii are determined from k_min, mu_max and k_max, mu_min
    kMuToENu(k_hmg, mu_hmg_2, dummy, nu_hmg);
    kMuToENu(k_hmg_2, mu_hmg, dummy, nu_hmg_2);
}

/*
 * Mori-Tanaka homogenization method
 * isotropic strain concentration tensor = 3*k_S[r]*J + 2*mu_S[r]*I
 * refRow is the reference phase of this scheme (0 is the first phase)
 * for derivation see Bernard et al.: A multiscale micromechanics-hydration model... CCR,pp. 1293-1309, 2003
 */
void Homogenize :: moriTanaka(FloatMatrix &PhaseMatrix, int refRow) {
    double f_r, E_r, nu_r, k_r, mu_r;
    double f_m, E_m, nu_m, k_m, mu_m;
    double fr_tot = 0.;
    double nom_k_MT = 0., denom_k_MT = 0., nom_mu_MT = 0., denom_mu_MT = 0.;
    double k_S_denom = 0., mu_S_denom = 0.;
    double alpha_m, beta_m;
    int numRows = PhaseMatrix.giveNumberOfRows();
    int r;
    double *k_S, *mu_S;
    k_S = new double [ numRows ];
    mu_S = new double [ numRows ];

    checkVolFraction(PhaseMatrix);

    //determine matrix parameters
    f_m = PhaseMatrix(refRow, 0);
    E_m = PhaseMatrix(refRow, 1);
    nu_m = PhaseMatrix(refRow, 2);

    ENuToKMu(E_m, nu_m, k_m, mu_m);

    //auxiliary parameters
    alpha_m = 3. * k_m / ( 3. * k_m + 4 * mu_m );
    beta_m  = 6. * ( k_m + 2. * mu_m ) / ( 5. * ( 3. * k_m + 4. * mu_m ) );

    //cycle through all phases, including the matrix phase
    for ( r = 0; r < numRows; r++ ) {
        f_r = PhaseMatrix(r, 0);
        E_r = PhaseMatrix(r, 1);
        nu_r = PhaseMatrix(r, 2);

        fr_tot += f_r;
        ENuToKMu(E_r, nu_r, k_r, mu_r);

        nom_k_MT += f_r * k_r / ( 1. + alpha_m * ( k_r / k_m - 1. ) );
        denom_k_MT += f_r / ( 1. + alpha_m * ( k_r / k_m - 1. ) );

        nom_mu_MT += f_r * mu_r / ( 1. + beta_m * ( mu_r / mu_m - 1. ) );
        denom_mu_MT += f_r / ( 1. + beta_m * ( mu_r / mu_m - 1. ) );

        k_S [ r ] = 1. / ( 1. + alpha_m * ( k_r / k_m - 1. ) );
        mu_S [ r ] = 1. / ( 1. + beta_m * ( mu_r / mu_m - 1. ) );
    }


    k_hmg = nom_k_MT / denom_k_MT;
    mu_hmg = nom_mu_MT / denom_mu_MT;

    kMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);

    //calculate denominator of strain concentration tensor (in terms of k_S and mu_S) and final values
    for ( r = 0; r < numRows; r++ ) {
        f_r = PhaseMatrix(r, 0);
        E_r = PhaseMatrix(r, 1);
        nu_r = PhaseMatrix(r, 2);
        ENuToKMu(E_r, nu_r, k_r, mu_r);
        k_S_denom += f_r * 1. / ( 1. + alpha_m * ( k_r / k_m - 1. ) );
        mu_S_denom += f_r * 1. / ( 1. + beta_m * ( mu_r / mu_m - 1. ) );
    }

    for ( r = 0; r < numRows; r++ ) {
        E_r = PhaseMatrix(r, 1);
        nu_r = PhaseMatrix(r, 2);
        ENuToKMu(E_r, nu_r, k_r, mu_r);
        k_S [ r ] /= k_S_denom;
        mu_S [ r ] /= mu_S_denom;
    }

    delete [] k_S;
    delete [] mu_S;
}

// for derivation see Bernard et al.:A multiscale micromechanics-hydration model... CCR,pp. 1293-1309, 2003
void Homogenize :: selfConsistent(FloatMatrix &PhaseMatrix) {
    double f_r, E_r, nu_r, k_r, mu_r;
    double k_SCS, mu_SCS, nom_k_SCS, denom_k_SCS, nom_mu_SCS, denom_mu_SCS;
    double fr_tot;
    double alpha_m, beta_m;
    double k_S_denom = 0., mu_S_denom = 0.;
    int numRows = PhaseMatrix.giveNumberOfRows();
    double *k_S, *mu_S;
    int i, r;

    checkVolFraction(PhaseMatrix);

    k_S = new double [ numRows ];
    mu_S = new double [ numRows ];

    /*first estimation*/
    k_SCS = 10.;
    mu_SCS = .3;

    /*iteration for nonlinear equations*/
    for ( i = 1; i < 100; i++ ) {
        fr_tot = 0.;
        nom_k_SCS = 0;
        denom_k_SCS = 0;
        nom_mu_SCS = 0;
        denom_mu_SCS = 0;
        for ( r = 0; r < numRows; r++ ) {
            f_r = PhaseMatrix(r, 0);
            E_r = PhaseMatrix(r, 1);
            nu_r = PhaseMatrix(r, 2);

            fr_tot += f_r;

            ENuToKMu(E_r, nu_r, k_r, mu_r);

            alpha_m = 3. * k_SCS / ( 3. * k_SCS + 4 * mu_SCS );
            beta_m = 6. * ( k_SCS + 2. * mu_SCS ) / ( 5. * ( 3. * k_SCS + 4. * mu_SCS ) );

            nom_k_SCS += f_r * k_r / ( 1 + alpha_m * ( k_r / k_SCS - 1 ) );
            denom_k_SCS += f_r / ( 1 + alpha_m * ( k_r / k_SCS - 1 ) );

            nom_mu_SCS += f_r * mu_r / ( 1 + beta_m * ( mu_r / mu_SCS - 1 ) );
            denom_mu_SCS += f_r / ( 1 + beta_m * ( mu_r / mu_SCS - 1 ) );
        }

        k_SCS = nom_k_SCS / denom_k_SCS;
        mu_SCS = nom_mu_SCS / denom_mu_SCS;
    }

    k_hmg = k_SCS;
    mu_hmg = mu_SCS;

    kMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);

    //calculate nominator and denominator of strain concentration tensor (in terms of k_S and mu_S) and final values
    for ( r = 0; r < numRows; r++ ) {
        f_r = PhaseMatrix(r, 0);
        E_r = PhaseMatrix(r, 1);
        nu_r = PhaseMatrix(r, 2);
        ENuToKMu(E_r, nu_r, k_r, mu_r);
        k_S [ r ] = 1. / ( 1. + alpha_m * ( k_r / k_hmg - 1. ) );
        mu_S [ r ] = 1. / ( 1. + beta_m * ( mu_r / mu_hmg - 1. ) );
        k_S_denom += f_r * 1. / ( 1. + alpha_m * ( k_r / k_hmg - 1. ) );
        mu_S_denom += f_r * 1. / ( 1. + beta_m * ( mu_r / mu_hmg - 1. ) );
    }

    for ( r = 0; r < numRows; r++ ) {
        E_r = PhaseMatrix(r, 1);
        nu_r = PhaseMatrix(r, 2);
        ENuToKMu(E_r, nu_r, k_r, mu_r);
        k_S [ r ] /= k_S_denom;
        mu_S [ r ] /= mu_S_denom;
    }

    delete [] k_S;
    delete [] mu_S;
}

/*
 * Homogenization scheme of Herve and Zaoui for n-spherical isotropic domains
 * see Herve E., Zaoui A. n-layered Inclusion-based Micromechanical Modelling,
 * Int. J. Engng Sci. 31, pp 1-10, 1993.
 *
 * The sum of n phases must yield a unity, (n+1) may have arbitrary radius
 *
 * The scheme is derived from stress and displacement equivalence at boundaries among phases
 * homogenized bulk and shear response is calculated
 * NumPhases is the number of phases including homogeneous medium (n+1)
 * In homogenized values, (n+1) phase may have arbitrary values - if it has the right
 * homogenized values of the assembly, the equations (49) and (51) are satisfied exactly without solution
 *
 * In case of n=2, the Hashin-Shtrikman bounds are calculated if the well-ordered moduli are used
 * Homogenized bulk modulus is 100 % OK (corresponds to HS bounds for 2 phases)
 *
 * Shear modulus is with small difference when compared to HS bound for 2 phases
 * typical usage for mortars and concrete is
 * phase 0 - aggregate
 * phase 1 - ITZ
 * phase 2 (n) - cement paste
 * phase 3 (n+1) - equivalent concrete medium
 *
 * The bulk modulus of this scheme coincides with HS bound. This is not true for shear modulus, see Christensen: Two Theoretical Elasticity Micromechanics Models, Journal of Elasticity 50: 15–25, 1998.
 */

void Homogenize :: herveZaoui(FloatMatrix &PhaseMatrix) {
    int i, j, phi;
    int NumPhases = PhaseMatrix.giveNumberOfRows() + 1; //need to extend for an equivalent homogeneous medium
    FloatMatrix PhaseMatrix1;
    checkVolFraction(PhaseMatrix);
    PhaseMatrix1 = PhaseMatrix;
    PhaseMatrix1.resizeWithData( NumPhases, PhaseMatrix.giveNumberOfColumns() );
    //add arbitrary values for a reference medium
    PhaseMatrix1(NumPhases - 1, 0) = 0.1;
    PhaseMatrix1(NumPhases - 1, 1) = 10.;
    PhaseMatrix1(NumPhases - 1, 2) = 0.3;

    // a)HYDROSTATIC pressure - conditions of displacement and radial stress compatibility between adjacent phases
    FloatArray r(NumPhases - 1), k(NumPhases), mu(NumPhases);
    FloatArray Q11(NumPhases - 1), Q21(NumPhases - 1), F(NumPhases), G(NumPhases);
    FloatMatrix J(2, 2), Jinv(2, 2), N(2, 2), Nhelp(2, 2), Q(2, 2);
    //set lambda_0 for strain control - arbitrary value for homogenization
    double lambda_0 = 0.735;
    double temp1;

    if ( NumPhases < 3 ) {
        OOFEM_ERROR1("Number of considered phases must be at least 3 (including hom. medium)\n");
    }

    F(NumPhases - 1) = lambda_0 / 3.;

    //radii of spheres, assume 1 as the maximal radius
    temp1 = 0.;
    for ( i = 0; i < NumPhases - 1; i++ ) {
        temp1 += PhaseMatrix1(i, 0);
        r(i) = 1. * pow(temp1, 1. / 3.);
        //printf( "r(%d)=%f\n", i, r(i) );
    }

    for ( i = 0; i < NumPhases; i++ ) {
        ENuToKMu( PhaseMatrix1(i, 1), PhaseMatrix1(i, 2), k(i), mu(i) );
        //printf("k(%d)=%f mu(%d)=%f\n", i, k(i), i, mu(i));
    }

    //create unit matrix for the first multiplication
    Nhelp(0, 0) = 1.;
    Nhelp(1, 1) = 1.;

    //calculate Q
    for ( phi = 0; phi <= NumPhases - 2; phi++ ) {
        //calculate inverse of J_{phi+1}
        fillJ(J, r(phi), mu, k, phi + 1);
        Jinv.beInverseOf(J);
        //calculate J_{phi}
        fillJ(J, r(phi), mu, k, phi);
        //calulate N^{phi}=J_{phi+1}^{-1}*J_{phi}
        N.beProductOf(Jinv, J);
        //calculate a part of Q^{phi}
        Q.beProductOf(N, Nhelp);
        //store part of matrix Q in vector Q11 and Q21
        Q11(phi) = Q(0, 0);
        Q21(phi) = Q(1, 0);
        //copy Q to Nhelp and use in the next cycle
        Nhelp = Q;
    }

    //calculate V (contains vector components F,G), F(0) and G(0) make no sense
    for ( phi = 1; phi <= NumPhases - 1; phi++ ) {
        F(phi) = F(NumPhases - 1) * Q11(phi - 1) / Q11(NumPhases - 2);
        G(phi) = F(NumPhases - 1) * Q21(phi - 1) / Q11(NumPhases - 2);
        //printf("F(%d)=%f G(%d)=%f\n", phi, F(phi), phi, G(phi));
    }

    //check displacement and stress continuity, phi is the more distant phase
    /*
     * for ( int phi = 1; phi <= NumPhases - 1; phi++ ) {
     *  printf("r=%f displ(in)=%f displ(out)=%f\n", r(phi-1), F(phi-1)*r(phi-1)+G(phi-1)/(r(phi-1)*r(phi-1)), F(phi)*r(phi-1)+G(phi)/(r(phi-1)*r(phi-1)));
     *  printf("r=%f stress(in)=%f stress(out)=%f\n", r(phi-1), 3.*k(phi-1)*F(phi-1)-4.*mu(phi-1)*G(phi-1)/(r(phi-1)*r(phi-1)*r(phi-1)), 3.*k(phi)*F(phi)-4.*mu(phi)*G(phi)/(r(phi-1)*r(phi-1)*r(phi-1)));
     * }
     */

    /*replace n+1 layer with equivalent homogeneous media calculate homogenized k
     * checked with Mori-Tanaka method for two phases, bulk modulus is exact
     */
    k_hmg = ( 3. * k(NumPhases - 2) * pow(r(NumPhases - 2), 3.) * Q11(NumPhases - 3) - 4. * mu(NumPhases - 2) * Q21(NumPhases - 3) ) / ( 3. * ( pow(r(NumPhases - 2), 3.) * Q11(NumPhases - 3) + Q21(NumPhases - 3) ) );

    //printf("Homogenized k=%f\n", k_hmg);

    //b)SIMPLE SHEAR - continuity of radial and phi displacements and two shear stresses
    FloatArray B(NumPhases), C(NumPhases), D(NumPhases), A(NumPhases);
    FloatArray P_v(4), Phelp(4);
    FloatMatrix L(4, 4), Linv(4, 4), M(4, 4), Mhelp(4, 4), Z(4, 4), M_test(4, 4);
    FloatMatrix *P_arr; //array to store individual P matrices
    //set gamma for displacement approach - arbitrary value for homogenization
    //variables for homogenized values
    double a, b, c, nu_n, sqr = 0.;
    double gamma = 10.; //no effect for homogenziation

    A.zero();
    A(NumPhases - 1) = gamma;

    //create unit matrix for the first multiplication
    Mhelp.zero();
    Mhelp.beUnitMatrix();

    //allocate memory for P matrices, store there P11, P12, P21, P22
    P_arr = new FloatMatrix [ NumPhases - 1 ];

    //calculate P
    for ( phi = 0; phi < NumPhases - 1; phi++ ) {
        //calculate inverse of L_{phi+1}
        fillL(L, r(phi), mu, k, phi + 1);
        //printf("%f", L(0,1));
        Linv.beInverseOf(L);
        //calculate L_{phi}
        fillL(L, r(phi), mu, k, phi);
        //calulate M^{phi}=L_{phi+1}^{-1}*L_{phi}
        M.beProductOf(Linv, L);

        /*test routine directly for M
         * this part calculate inverse of M directly
         * M=M_test so the L, Linv and M assembly are OK
         *
         *  double ak, bk, ck, dk, ek, fk, alphak, nuk, nukp1, koef;
         *
         *  nuk=PhaseMatrix1(phi,2);
         *  nukp1=PhaseMatrix1(phi+1,2);
         *
         *  ak=mu(phi)/mu(phi+1)*(7.+5.*nuk)*(7.-10.*nukp1)-(7.-10*nuk)*(7.+5.*nukp1);
         *  bk=4.*(7.-10.*nuk)+mu(phi)/mu(phi+1)*(7.+5.*nuk);
         *  ck=(7.-5.*nukp1)+2.*mu(phi)/mu(phi+1)*(4.-5.*nukp1);
         *  dk=(7.+5.*nukp1)+4.*mu(phi)/mu(phi+1)*(7.-10.*nukp1);
         *  ek=2*(4.-5.*nuk)+mu(phi)/mu(phi+1)*(7.-5.*nuk);
         *  fk=(4.-5.*nuk)*(7.-5.*nukp1)-mu(phi)/mu(phi+1)*(4.-5.*nukp1)*(7.-5.*nuk);
         *  alphak=mu(phi)/mu(phi+1)-1.;
         *  koef=1./(5.*(1-nukp1));
         *
         *  M_test(0,0)=koef*ck/3.;
         *  M_test(0,1)=koef*pow(r(phi),2)*(3*bk-7*ck)/(5.*(1.-2*nuk));
         *  M_test(0,2)=koef*(-12.)*alphak/pow(r(phi),5);
         *  M_test(0,3)=koef*4.*(fk-27.*alphak)/(15.*(1.-2*nuk)*pow(r(phi),3));
         *
         *  M_test(1,0)=koef*0.;
         *  M_test(1,1)=koef*(1.-2*nukp1)*bk/(7.*(1.-2.*nuk));
         *  M_test(1,2)=koef*(-20.)*(1.-2.*nukp1)*alphak/(7*pow(r(phi),7));
         *  M_test(1,3)=koef*(-12.)*alphak*(1.-2.*nukp1)/(7.*(1.-2.*nuk)*(pow(r(phi),5)));
         *
         *  M_test(2,0)=koef*pow(r(phi),5)*alphak/2.;
         *  M_test(2,1)=koef*(-pow(r(phi),7))*(2*ak+147*alphak)/(70.*(1.-2.*nuk));
         *  M_test(2,2)=koef*dk/7.;
         *  M_test(2,3)=koef*pow(r(phi),2)*(105*(1.-nukp1)+12.*alphak*(7.-10.*nukp1)-7.*ek)/(35.*(1.-2*nuk));
         *
         *  M_test(3,0)=koef*(-5.)/6.*(1.-2*nukp1)*alphak*pow(r(phi),3);
         *  M_test(3,1)=koef*7.*(1.-2*nukp1)*alphak*pow(r(phi),5)/(2.*(1.-2.*nuk));
         *  M_test(3,2)=koef*0.;
         *  M_test(3,3)=koef*ek*(1.-2.*nukp1)/(3.*(1.-2*nuk));
         *
         *  for(int ii=0;ii<4;ii++)
         *    for(int jj=0;jj<4;jj++)
         * M(ii,jj)=M_test(ii,jj);
         */
        P_arr [ phi ].resize(4, 4);
        P_arr [ phi ].beProductOf(M, Mhelp);
        Mhelp = P_arr [ phi ];
    }

    //printf("PPP %f\n", P_arr[1](3,0)*P_arr[1](1,1) - P_arr[1](3,1)*P_arr[1](1,0));

    //calculate vector W (contains scalar components A,B,C,D for each step)

    P_v(0) = P_arr [ NumPhases - 2 ](1, 1); //matrix stored column wise (P11, P21, P31....), P_22
    P_v(1) = -P_arr [ NumPhases - 2 ](1, 0); // P_21
    P_v(2) = P_v(3) = 0.;

    for ( phi = 1; phi <= NumPhases - 1; phi++ ) {
        Phelp.beProductOf(P_arr [ phi - 1 ], P_v);
        Phelp.times( A(NumPhases - 1) / sqrt(2.) / ( P_arr [ NumPhases - 2 ](0, 0) * P_arr [ NumPhases - 2 ](1, 1) - P_arr [ NumPhases - 2 ](0, 1) * P_arr [ NumPhases - 2 ](1, 0) ) );
        A(phi) = Phelp(0);
        B(phi) = Phelp(1);
        C(phi) = Phelp(2);
        D(phi) = Phelp(3);
        //printf("A(%d]=%f B(%d)=%f C(%d)=%f D(%d)=%f\n", phi, A(phi), phi, B(phi), phi, C(phi), phi, D(phi));
    }

    //check displacement and stress continuity, phi is the less distant phase
    //this cycle for for checking purposes only
    for ( phi = 0; phi <= NumPhases - 2; phi++ ) {
        //phase closer to the center
        fillL(L, r(phi), mu, k, phi);
        Phelp(0) = A(phi);
        Phelp(1) = B(phi);
        Phelp(2) = C(phi);
        Phelp(3) = D(phi);
        P_v.beProductOf(L, Phelp);
        //printf("(in) n=%d r=%f ur=%f uphi=%f srp=%f srf=%f\n",phi, r(phi), P_v(0),P_v(1), P_v(2), P_v(3));

        fillL(L, r(phi), mu, k, phi + 1);
        Phelp(0) = A(phi + 1);
        Phelp(1) = B(phi + 1);
        Phelp(2) = C(phi + 1);
        Phelp(3) = D(phi + 1);
        P_v.beProductOf(L, Phelp);
        //printf("(out) n=%d r=%f ur=%f uphi=%f srp=%f srf=%f\n", phi,r(phi), P_v(0),P_v(1), P_v(2), P_v(3));
    }

    // replace n+1 layer with equivalent homogeneous media, calculate homogenized mu
    for ( i = 0; i <= 3; i++ ) { //alpha, rows in Z
        for ( j = 0; j <= 3; j++ ) { //beta, columns in Z
            Z(i, j) = P_arr [ NumPhases - 3 ](i, 0) * P_arr [ NumPhases - 3 ](j, 1) - P_arr [ NumPhases - 3 ](j, 0) * P_arr [ NumPhases - 3 ](i, 1);
            //      printf("Z=%f ", Z(i,j));
        }
    }

    nu_n = ( 3. * k(NumPhases - 2) - 2. * mu(NumPhases - 2) ) / ( 6. * k(NumPhases - 2) + 2. * mu(NumPhases - 2) );

    a = 4. * pow(r(NumPhases - 2), 10.) * ( 1. - 2. * nu_n ) * ( 7. - 10. * nu_n ) * Z(0, 1) +
        20. * pow(r(NumPhases - 2), 7.) * ( 7. - 12. * nu_n + 8. * nu_n * nu_n ) * Z(3, 1) +
        12. * pow(r(NumPhases - 2), 5.) * ( 1. - 2. * nu_n ) * ( Z(0, 3) - 7. * Z(1, 2) ) +
        20. * pow(r(NumPhases - 2), 3.) * ( 1. - 2. * nu_n ) * ( 1. - 2. * nu_n ) * Z(0, 2) +
        16. * ( 4. - 5. * nu_n ) * ( 1. - 2. * nu_n ) * Z(3, 2);

    b = 3. * pow(r(NumPhases - 2), 10.) * ( 1. - 2. * nu_n ) * ( 15. * nu_n - 7. ) * Z(0, 1) +
        60. * pow(r(NumPhases - 2), 7.) * ( nu_n - 3. ) * nu_n * Z(3, 1) -
        24. * pow(r(NumPhases - 2), 5.) * ( 1. - 2. * nu_n ) * ( Z(0, 3) - 7. * Z(1, 2) ) -
        40. * pow(r(NumPhases - 2), 3.) * ( 1. - 2. * nu_n ) * ( 1. - 2. * nu_n ) * Z(0, 2) -
        8. * ( 1. - 5. * nu_n ) * ( 1. - 2. * nu_n ) * Z(3, 2);

    c = -pow(r(NumPhases - 2), 10.) * ( 1. - 2. * nu_n ) * ( 7. + 5 * nu_n ) * Z(0, 1) +
        10. * pow(r(NumPhases - 2), 7.) * ( 7. - nu_n * nu_n ) * Z(3, 1) +
        12. * pow(r(NumPhases - 2), 5.) * ( 1. - 2. * nu_n ) * ( Z(0, 3) - 7. * Z(1, 2) ) +
        20. * pow(r(NumPhases - 2), 3.) * ( 1. - 2. * nu_n ) * ( 1. - 2. * nu_n ) * Z(0, 2) -
        8. * ( 7. - 5. * nu_n ) * ( 1. - 2. * nu_n ) * Z(3, 2);

    //solve quadratic equation, report higher number
    sqr = b * b - 4. * a * c;
    if ( sqr < 0 ) {
        OOFEM_ERROR1("Shear modulus does not yield a real number\n");
    }

    if ( ( -b - pow(sqr, 0.5) ) >= 0 ) {
        OOFEM_WARNING1("Two solutions for the shear modulus were found, continuing with the higher value\n");
    }

    //usually this higher value is reported
    mu_hmg = mu(NumPhases - 2) * ( -b + pow(sqr, 0.5) ) / ( 2. * a );
    //uncomment if lower value of mu should be used, only when positive
    //mu_hmg=mu(NumPhases-2)*(-b+pow(sqr, 0.5))/(2.*a);

    //printf("Homogenized mu=%f\n", mu_hmg);

    kMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);

    delete [] P_arr;
}

/*Hirsh model as a combination of Voigt and Reuss (parameter chi)*/
void Homogenize :: hirsch(FloatMatrix &PhaseMatrix, double chi) {
    double k_Voigt, mu_Voigt, k_Reuss, mu_Reuss;
    voigt(PhaseMatrix);
    k_Voigt = k_hmg;
    mu_Voigt = mu_hmg;
    reuss(PhaseMatrix);
    k_Reuss = k_hmg;
    mu_Reuss = mu_hmg;
    k_hmg = ( 1. - chi ) * k_Reuss + chi * k_Voigt;
    mu_hmg = ( 1. - chi ) * mu_Reuss + chi * mu_Voigt;
    kMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);
}


void Homogenize :: hansen(FloatMatrix &PhaseMatrix) {
    double f_i, E_i, nu_i, f_m, E_m, nu_m;

    if ( PhaseMatrix.giveNumberOfRows() != 2 ) {
        OOFEM_ERROR1("Only two phases are allowed\n");
    }

    checkVolFraction(PhaseMatrix);

    f_i = PhaseMatrix(0, 0);
    E_i = PhaseMatrix(0, 1);
    nu_i = PhaseMatrix(0, 2);
    f_m = PhaseMatrix(1, 0);
    E_m = PhaseMatrix(1, 1);
    nu_m = PhaseMatrix(1, 2);

    E_hmg = E_i * ( ( f_i * E_i + ( 1 + f_m ) * E_m ) / ( ( 1. + f_m ) * E_i + f_i * E_m ) );
    nu_hmg = 0.2;
    ENuToKMu(E_hmg, nu_hmg, k_hmg, mu_hmg);
}


void Homogenize :: counto(FloatMatrix &PhaseMatrix) {
    double f_i, E_i, nu_i, f_m, E_m, nu_m;
    if ( PhaseMatrix.giveNumberOfRows() != 2 ) {
        OOFEM_ERROR1("Only two phases are allowed\n");
    }

    checkVolFraction(PhaseMatrix);

    f_i = PhaseMatrix(0, 0);
    E_i = PhaseMatrix(0, 1);
    nu_i = PhaseMatrix(0, 2);
    f_m = PhaseMatrix(1, 0);
    E_m = PhaseMatrix(1, 1);
    nu_m = PhaseMatrix(1, 2);

    E_hmg = 1. / ( ( 1. - sqrt(f_m) ) / E_i + 1. / ( E_i * ( ( 1. - sqrt(f_m) ) / sqrt(f_m) ) + E_m ) );
    nu_hmg = nu_m;
    ENuToKMu(E_hmg, nu_hmg, k_hmg, mu_hmg);
}

/* Kuster-Toksoz model from Monteiro: Elastic Moduli Lightweight Aggregate Model,
 * Cement and Concrete Research, 1995, No.2, pp. 276-280.
 */
void Homogenize :: kusterToksoz(FloatMatrix &PhaseMatrix) {
    double k_KT, mu_KT, nom, denom;
    double f_i, E_i, nu_i, k_i, mu_i, f_m, E_m, nu_m, mu_m, k_m;
    if ( PhaseMatrix.giveNumberOfRows() != 2 ) {
        OOFEM_ERROR1("Only two phases are allowed\n");
    }

    checkVolFraction(PhaseMatrix);

    f_i = PhaseMatrix(0, 0);
    E_i = PhaseMatrix(0, 1);
    nu_i = PhaseMatrix(0, 2);
    ENuToKMu(E_i, nu_i, k_i, mu_i);
    f_m = PhaseMatrix(1, 0);
    E_m = PhaseMatrix(1, 1);
    nu_m = PhaseMatrix(1, 2);
    ENuToKMu(E_m, nu_m, k_m, mu_m);

    nom = 1 + ( 4 * mu_m * ( k_i - k_m ) / ( ( 3 * k_i * +4 * mu_m ) * k_m ) ) * f_i;
    denom = 1 - ( 3 * ( k_i - k_m ) / ( 3 * k_i + 4 * mu_m ) ) * f_i;
    k_KT = k_m * nom / denom;

    nom = ( 6 * k_m + 12 * mu_m ) * mu_i + ( 9 * k_m + 8 * mu_m ) * ( f_m * mu_m + f_i * mu_i );
    denom = ( 9 * k_m + 8 * mu_m ) * mu_m + ( 6 * k_m + 12 * mu_m ) * ( f_m * mu_i + f_i * mu_m );
    mu_KT = mu_m * nom / denom;

    k_hmg = k_KT;
    mu_hmg = mu_KT;

    kMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);
}

void Homogenize :: printYourself(int out)
{
    if ( out == 0 ) {
        printf("E %.4e\t nu %.4e\t k %.4e\t mu %.4e\n", this->E_hmg, this->nu_hmg, this->k_hmg, this->mu_hmg);
    } else {
        printf("E %.4e\t nu %.4e\t k %.4e\t mu %.4e\n", this->E_hmg_2, this->nu_hmg_2, this->k_hmg_2, this->mu_hmg_2);
    }
}


void Homogenize :: ENuToKMu(const double E, const double nu, double &k, double &mu)
{
    mu = E / ( 2. * ( 1. + nu ) );
    k = E / ( 3. * ( 1. - 2. * nu ) );
}

void Homogenize :: kMuToENu(const double k, const double mu, double &E, double &nu)
{
    E  = 9. * k * mu / ( 3. * k + mu );
    nu = ( 3. * k - 2. * mu ) / ( 6. * k + 2 * mu );
}

//Auxiliary function for Hashin-Shtrikman-Walpole shear bounds
double Homogenize :: zeta(double k, double mu) {
    return mu / 6. * ( 9. * k + 8. * mu ) / ( k + 2. * mu );
}

//Auxiliary function for Hashin-Shtrikman-Walpole shear bounds
double Homogenize :: gamma(FloatMatrix &SortedPhaseMatrix, double zeta) {
    double gamma = 0;
    int r;
    int NumPhases = SortedPhaseMatrix.giveNumberOfRows();

    for ( r = 0; r < NumPhases; r++ ) {
        gamma += SortedPhaseMatrix(r, 0) / ( SortedPhaseMatrix(r, 2) + zeta );
    }

    gamma = 1. / gamma;
    gamma -= zeta;
    return gamma;
}

void Homogenize :: checkVolFraction(FloatMatrix &PhaseMatrix) {
    double volTot = 0.;
    int NumPhases = PhaseMatrix.giveNumberOfRows();
    int r;
    for ( r = 0; r < NumPhases; r++ ) {
        volTot += PhaseMatrix(r, 0);
    }

    if ( volTot < 0.99 || volTot > 1.01 ) {
        OOFEM_ERROR3("Volumetric fraction of phases 0-%d is %f, which not the unity\n", NumPhases, volTot);
    }
}

//help function for filling the J matrix
void Homogenize :: fillJ(FloatMatrix &J, double r, const FloatArray &mu, const FloatArray &k, int phase) {
    J(0, 0) = r;
    J(0, 1) = 1. / ( pow(r, 2.) );
    J(1, 0) = 3 * k(phase);
    J(1, 1) = -4. * mu(phase) / ( pow(r, 3.) );
}

//help function for filling the L matrix
void Homogenize :: fillL(FloatMatrix &L, double r, const FloatArray &mu, const FloatArray &k, int phase) {
    double mu_l = mu(phase);
    double k_l = k(phase);
    double nu_l = ( 3. * k_l - 2. * mu_l ) / ( 6. * k_l + 2. * mu_l ); //OK

    L(0, 0) = r;
    L(0, 1) = -6. *nu_l *pow(r, 3.) / ( 1. - 2. * nu_l );
    L(0, 2) = 3. / pow(r, 4.);
    L(0, 3) = ( 5. - 4. * nu_l ) / ( ( 1. - 2. * nu_l ) * pow(r, 2.) );

    L(1, 0) = r;
    L(1, 1) = -( 7. - 4. * nu_l ) * pow(r, 3.) / ( 1. - 2. * nu_l );
    L(1, 2) = -2. / pow(r, 4.);
    L(1, 3) = 2. / pow(r, 2.);

    L(2, 0) = mu_l;
    L(2, 1) = 3. *nu_l *mu_l *pow(r, 2.) / ( 1. - 2. * nu_l );
    L(2, 2) = -12. * mu_l / pow(r, 5.);
    L(2, 3) = 2. * ( nu_l - 5. ) * mu_l / ( ( 1. - 2. * nu_l ) * pow(r, 3.) );

    L(3, 0) = mu_l;
    L(3, 1) = -( 7. + 2 * nu_l ) * mu_l * pow(r, 2.) / ( 1. - 2. * nu_l );
    L(3, 2) = 8. * mu_l / pow(r, 5.);
    L(3, 3) = 2. * ( 1 + nu_l ) * mu_l / ( ( 1. - 2. * nu_l ) * pow(r, 3.) );
}
} // end namespace oofem
