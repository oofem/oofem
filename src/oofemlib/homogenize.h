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

#ifndef homogenize_h
#define homogenize_h

#include "oofemcfg.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "error.h"

namespace oofem {
/**
 * Class for elastic homogenization. Following homogenization estimates and bounds are implemented : Voigt, Reuss, Hirsch,
 * Counto, Mori-Tanaka, Self-consitent, Hashin-Shtrikman-Walpole, Kuster-Toksoz, and Herve-Zaoui's n-layered spherical assemblage.
 * Only isotropic materials and spherical inclusions are considered for the homogenization. Many schemes work with unlimited number of elastic phases.
 */

class OOFEM_EXPORT Homogenize
{
public:

    /// Constructor
    Homogenize();

    /// Destructor
    ~Homogenize(void) { }

    /**  Parallel scheme of Voigt for any number of isotropic phases
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase
     */
    void voigt(FloatMatrix &PhaseMatrix);

    /**  Serial scheme of Reuss for any number of isotropic phases
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase
     */
    void reuss(FloatMatrix &PhaseMatrix);

    /** Hashin-Shtrikman-Walpole lower and upper bounds for arbitrary number of isotropic phases
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase
     * implemented according to J.G. Berryman: Mixture Theories for Rock Properties, Physics and Phase Relations, A Handbook of Physical Constants, Am. Geophys, 1995
     */
    void hashinShtrikmanWalpole(FloatMatrix &PhaseMatrix);

    /** Mori-Tanaka homogenization method for spherical isotropic inclusions
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase
     * @param refRow row of the reference matrix phase (0 is the first)
     */
    void moriTanaka(FloatMatrix &PhaseMatrix, int refRow);

    /** Self-consistent homogenization method of Hill and Budiansky for spherical isotropic inclusions
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase
     */
    void selfConsistent(FloatMatrix &PhaseMatrix);

    /** Herve and Zaoui's homogenization scheme for n-spherical isotropic domains and arbitrary number of isotropic phases
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase
     */
    void herveZaoui(FloatMatrix &PhaseMatrix);

    /**  Hirsch's scheme combining Voigt and Reuss bounds
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase.
     * @param chi a weight parameter between Reuss (=0) and Voigt (=1) bounds
     */
    void hirsch(FloatMatrix &PhaseMatrix, double chi);

    /** Hansen's model for a spherical inclusion in a spherical matrix, for Poisson's ratio 0.2 and two phases
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase. The first row is an inclusion, the second row corresponds to the matrix.
     */
    void hansen(FloatMatrix &PhaseMatrix);

    /** Counto's model for a prismatic inclusion in a prismatic matrix, two phases
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase. The first row is an inclusion, the second row corresponds to the matrix.
     */
    void counto(FloatMatrix &PhaseMatrix);

    /** Kuster-Toksoz model for two phases
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase. The first row is an inclusion, the second row corresponds to the matrix.
     */
    void kusterToksoz(FloatMatrix &PhaseMatrix);

    /** Print homogenized output
     * @param out 0 (default) prints homogenized values or a lower bound, 1 prints an upper bound
     */
    void printYourself(int out = 0);

    /// Convert Young's modulus and Poisson's ratio to bulk and shear moduli, isotropic material
    void ENuToKMu(const double E, const double nu, double &k, double &mu);

    /// Convert bulk and shear moduli to Young's modulus and Poisson's ratio, isotropic material
    void kMuToENu(const double k, const double mu, double &E, double &nu);

    /// Effective Young's modulus or the lower bound
    double E_hmg;

    /// Upper bound of Young's modulus if applicable
    double E_hmg_2;

    /// Effective Poisson's ratio
    double nu_hmg;

    /// Effective Poisson's ratio or the lower bound
    double nu_hmg_2;

    /// Effective bulk modulus or the lower bound
    double k_hmg;

    /// Upper bound of bulk modulus if applicable
    double k_hmg_2;

    /// Effective shear modulus or the lower bound
    double mu_hmg;

    /// Upper shear modulus if applicable
    double mu_hmg_2;

private:
    /// Auxiliary function
    double lambda(double *PhaseMatrixKMu, int NumRows, double mu);

    /// Auxiliary function for Hashin-Shtrikman bounds
    double zeta(double k, double mu);

    /// Auxiliary function
    double gamma(FloatMatrix &SortedPhaseMatrix, double zeta);

    /** Check that the total volumetric fraction is the unity
     * @param PhaseMatrix matrix containing in each row the volume fraction, the Young modulus and the Poisson ratio for each phase
     */
    void checkVolFraction(FloatMatrix &PhaseMatrix);

    /// Auxiliary function for Herve-Zaoui scheme
    void fillJ(FloatMatrix &J, double r, const FloatArray &mu, const FloatArray &k, int phase);

    /// Auxiliary function for Herve-Zaoui scheme
    void fillL(FloatMatrix &L, double r, const FloatArray &mu, const FloatArray &k, int phase);
};
} // end namespace oofem
#endif // homogenize_h
