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

#ifndef microplanematerial_h
#define microplanematerial_h

#include "sm/Materials/structuralmaterial.h"
#include "matconst.h"
#include "floatarrayf.h"

///@name Input fields for MicroplaneMaterial
//@{
#define _IFT_MicroplaneMaterial_nmp "nmp"
#define _IFT_MicroplaneMaterial_e "e"
#define _IFT_MicroplaneMaterial_n "n"
//@}

namespace oofem {

#define MAX_NUMBER_OF_MICROPLANES 61

/// Defines the stress or strain state in a micro plane
struct MicroplaneState
{
    double n = 0.;
    double v = 0.;
    double m = 0.;
    double l = 0.;
};


/**
 * Abstract base class for all microplane models.
 *
 * This class provides methods for setting up the integration weights and normals
 * for particular microplanes. Also projection tensors are computed in advance
 * (see method initializeData) and stored.
 * Methods for computing macro strain projections onto particular microplanes and
 * for homogenization of stress vector are provided.
 */
class MicroplaneMaterial : public StructuralMaterial
{
protected:
    /// Number of microplanes.
    int numberOfMicroplanes = 0;

    /// Integration weights of microplanes.
    FloatArray microplaneWeights;
    /// Normals of microplanes.
    std::vector<FloatArrayF<3>> microplaneNormals;

    /// Kronecker's delta.
    FloatArrayF<6> Kronecker;

    /**
     * Normal projection tensors for all microplanes.
     * Due to symmetry, compressed form is stored.
     */
    std::vector<FloatArrayF<6>> N;
    /**
     * Shear projection tensors (m direction) for all microplanes.
     * Due to symmetry, compressed form is stored.
     */
    std::vector<FloatArrayF<6>> M;
    /**
     * Shear projection tensors (l direction) for all microplanes.
     * Due to symmetry, compressed form is stored.
     */
    std::vector<FloatArrayF<6>> L;

    /// Young's modulus
    double E = 0.;

    /// Poisson's ratio
    double nu = 0.;

public:

    /**
     * Constructor. Creates Microplane Material belonging to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    MicroplaneMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }

    /**
     * Computes the length of normal strain vector on given microplane.
     */
    double computeNormalStrainComponent(int mnumber, const FloatArray &macroStrain) const;
    /**
     * Computes the normal volumetric component of macro strain on given microplane.
     */
    double computeNormalVolumetricStrainComponent(const FloatArray &macroStrain) const;
    /**
     * Computes the normal deviatoric component of macro strain on given microplane.
     */
    double computeNormalDeviatoricStrainComponent(int mnumber, const FloatArray &macroStrain) const;
    /**
     * Computes the shear component (in m direction) of macro strain on given microplane.
     */
    double computeShearMStrainComponent(int mnumber, const FloatArray &macroStrain) const;
    /**
     * Computes the shear component (in l direction) of macro strain on given microplane.
     */
    double computeShearLStrainComponent(int mnumber, const FloatArray &macroStrain) const;
    /**
     * Computes the vector of all micro stress components (Ev, En, Em, El) of macro strain
     * vector on given microplane.
     */
    MicroplaneState computeStrainVectorComponents(int mnumber, const FloatArray &macroStrain) const;


    /**
     * Returns microplane integration weight.
     * @param mplane Microplane.
     * @return Integration weight of given microplane.
     */
    double giveMicroplaneIntegrationWeight(int mnumber) const;

    /**
     * Initializes internal data (integration weights,
     * microplane normals and computes projection tensors).
     * @param numberOfMicroplanes Number of required microplanes.
     */
    virtual void initializeData(int numberOfMicroplanes);

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
};
} // end namespace oofem
#endif // microplanematerial_h
