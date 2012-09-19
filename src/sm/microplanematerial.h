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

#ifndef microplanematerial_h
#define microplanematerial_h

#include "structuralmaterial.h"
#include "gausspnt.h"
#include "matconst.h"

namespace oofem {
class Microplane;

#define MAX_NUMBER_OF_MICROPLANES 61

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
    int numberOfMicroplanes;

    /// Integration weights of microplanes.
    double microplaneWeights [ MAX_NUMBER_OF_MICROPLANES ];
    /// Normals of microplanes.
    double microplaneNormals [ MAX_NUMBER_OF_MICROPLANES ] [ 3 ];

    /// Kronecker's delta.
    double Kronecker [ 6 ];

    /**
     * Normal projection tensors for all microplanes.
     * Due to symmetry, compressed form is stored.
     */
    double N [ MAX_NUMBER_OF_MICROPLANES ] [ 6 ];
    /**
     * Shear projection tensors (m direction) for all microplanes.
     * Due to symmetry, compressed form is stored.
     */
    double M [ MAX_NUMBER_OF_MICROPLANES ] [ 6 ];
    /**
     * Shear projection tensors (l direction) for all microplanes.
     * Due to symmetry, compressed form is stored.
     */
    double L [ MAX_NUMBER_OF_MICROPLANES ] [ 6 ];

public:

    /**
     * Constructor. Creates Microplane Material belonging to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    MicroplaneMaterial(int n, Domain *d) : StructuralMaterial(n, d) { numberOfMicroplanes = 0; }
    /// Destructor.
    virtual ~MicroplaneMaterial() { }

    /**
     * Computes real stress vector on given microplane
     * (the meaning of  values depends on particular implementation,
     * e.g, can contain volumetric, deviatoric normal stresses and shear stresses on microplane)
     * for given increment of microplane strains.
     * @param answer Computed result.
     * @param mplane Pointer to microplane object, for which response is computed.
     * @param strain Strain vector.
     * @param tStep Time step.
     */
    virtual void giveRealMicroplaneStressVector(FloatArray &answer, Microplane *mplane,
                                                const FloatArray &strain, TimeStep *tStep) = 0;

    /**
     * Computes the length of normal strain vector on given microplane.
     */
    double computeNormalStrainComponent(Microplane *mplane, const FloatArray &macroStrain);
    /**
     * Computes the normal volumetric component of macro strain on given microplane.
     */
    double computeNormalVolumetricStrainComponent(Microplane *mplane, const FloatArray &macroStrain);
    /**
     * Computes the normal deviatoric component of macro strain on given microplane.
     */
    double computeNormalDeviatoricStrainComponent(Microplane *mplane, const FloatArray &macroStrain);
    /**
     * Computes the shear component (in m direction) of macro strain on given microplane.
     */
    double computeShearMStrainComponent(Microplane *mplane, const FloatArray &macroStrain);
    /**
     * Computes the shear component (in l direction) of macro strain on given microplane.
     */
    double computeShearLStrainComponent(Microplane *mplane, const FloatArray &macroStrain);
    /**
     * Computes the vector of all micro stress components (Ev, En, Em, El) of macro strain
     * vector on given microplane.
     */
    void computeStrainVectorComponents(FloatArray &answer, Microplane *mplane,
                                       const FloatArray &macroStrain);


    /**
     * Computes normal of given microplane.
     * @param answer Normal of given microplane.
     * @param mplane Microplane, which normal will be computed.
     */
    virtual void giveMicroplaneNormal(FloatArray &answer, Microplane *mplane);
    /**
     * Returns microplane integration weight.
     * @param mplane Microplane.
     * @return Integration weight of given microplane.
     */
    virtual double giveMicroplaneIntegrationWeight(Microplane *mplane);

    /// Returns i-th microplane belonging to master-macro-integration point. )-based indexing.
    Microplane *giveMicroplane(int i, GaussPoint *masterGp);

    /**
     * Initializes internal data (integration weights,
     * microplane normals and computes projection tensors).
     * @param numberOfMicroplanes Number of required microplanes.
     */
    virtual void initializeData(int numberOfMicroplanes);

    /**
     * Returns corresponding material mode for microplane according to macro integration mode
     */
    virtual MaterialMode giveCorrespondingSlaveMaterialMode(MaterialMode masterMode);

    virtual contextIOResultType saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);
    virtual contextIOResultType restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "MicroplaneMaterial"; }
    virtual classType giveClassID() const { return MicroplaneMaterialClass; }

    virtual IntegrationPointStatus *giveMicroplaneStatus(GaussPoint *gp);

protected:
    virtual MaterialStatus *CreateMicroplaneStatus(GaussPoint *gp) = 0;
    virtual void initTempStatus(GaussPoint *gp);
};
} // end namespace oofem
#endif // microplanematerial_h
