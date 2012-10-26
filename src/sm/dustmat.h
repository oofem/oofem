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
#ifndef dustmat_h
#define dustmat_h

#include "flotarry.h"
#include "flotmtrx.h"

#include "structuralms.h"
#include "strainvector.h"
#include "stressvector.h"
#include "structuralmaterial.h"
#include "isolinearelasticmaterial.h"

namespace oofem {
/**
 * This class implements TODO
 * @author Jan Stransky
 */
class DustMaterialStatus : public StructuralMaterialStatus
{
protected:
	/// Volumetric plastic strain
	double volumetricPlasticStrain;
	double tempVolumetricPlasticStrain;

	/// Deviatoric of plastic strain
    StrainVector plasticStrainDeviator;
    StrainVector tempPlasticStrainDeviator;
	 double q;
	 double tempQ;



public:
    /// Constructor
    DustMaterialStatus(int n, Domain *d, GaussPoint *gp);

    /// Destructor
    virtual ~DustMaterialStatus();

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "DustMaterialStatus"; }
    virtual classType giveClassID() const { return DustMaterialStatusClass; }
    /**
     * Get the full plastic strain vector from the material status.
     * @param answer Plastic strain vector.
     */
    void  giveFullPlasticStrainVector(StrainVector &answer) const
    {
        StrainVector plasticStrainVector(_3dMat);
        plasticStrainDeviator.computeDeviatoricVolumetricSum(plasticStrainVector,
                                                             volumetricPlasticStrain);
        plasticStrainVector.convertToFullForm(answer);
    }
    /**
     * Get the plastic strain deviator from the material status.
     * @param answer Plastic strain deviator.
     */
    void givePlasticStrainDeviator(StrainVector &answer) const { answer = plasticStrainDeviator; }
    /**
     * Get the volumetric plastic strain from the material status.
     * @return Volumetric plastic strain.
     */
    double giveVolumetricPlasticStrain() const { return volumetricPlasticStrain; }
	 double giveQ() const { return q; }

    /**
     * Get the temp value of the full plastic strain vector from the material status.
     * @param answer Temp value of plastic strain vector.
     */
    void giveTempPlasticStrainVector(StrainVector &answer) const
    {
        tempPlasticStrainDeviator.computeDeviatoricVolumetricSum(answer,
                                                                 tempVolumetricPlasticStrain);
    }
    /**
     * Get the temp value of the plastic strain deviator from the material status.
     * @param answer Temp value of plastic strain deviator.
     */
    void giveTempPlasticStrainDeviator(StrainVector &answer) const { answer = tempPlasticStrainDeviator; }
    /**
     * Get the temp value of the volumetric strain deviator from the material status.
     * @return Temp value of volumetric plastic strain
     */
    double giveTempVolumetricPlasticStrain() const { return tempVolumetricPlasticStrain; }
    /**
     * Get the temp value of the hardening variable from the material status.
     * @return Temp value of hardening variable kappa.
     */
	  double giveTempQ() const { return tempQ; }

    /**
     * Assign the temp value of deviatoric plastic strain.
     * @param v New temp value of deviatoric plastic strain.
     */
    void letTempPlasticStrainDeviatorBe(const StrainVector &v) { tempPlasticStrainDeviator = v; }
    /**
     * Assign the temp value of volumetric plastic strain.
     * @param v New temp value of volumetric plastic strain.
     */
    void letTempVolumetricPlasticStrainBe(double v) { tempVolumetricPlasticStrain = v; }
	 void letTempQBe(double q) { tempQ = q; }
};

/**
 * This class implements TODO
 * @author Jan Stransky
 */
class DustMaterial : public StructuralMaterial
{
protected:
    IsotropicLinearElasticMaterial *LEMaterial;
	 double bulkModulus;
	 double shearModulus;

	 double alpha;
	 double beta;
	 double lambda;
	 double theta;
	 double ft;
	 double rEllipse;
	 int hardeningType;
	 double mHard;
	   
public:
    /// Constructor
    DustMaterial(int n, Domain *d);
    /// Destructor
    virtual ~DustMaterial();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasMaterialModeCapability(MaterialMode mMode);

    virtual const char *giveClassName() const { return "DustMaterial"; }
    virtual classType giveClassID() const { return DustMaterialClass; }

    virtual void giveRealStressVector(FloatArray &answer,
                              MatResponseForm form,
                              GaussPoint *gp,
                              const FloatArray &strainVector,
                              TimeStep *atTime);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                       MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime);

    virtual int giveIPValueSize(InternalStateType type,
                        GaussPoint *gp);

    virtual int giveIntVarCompFullIndx(IntArray &answer,
                               InternalStateType type,
                               MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *atTime)
    {
        LEMaterial->giveThermalDilatationVector(answer, gp, atTime);
    }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

	 double functionFe(double i1);
	 double functionFc(double sn, double i1, double q);
	 double functionX(double q);
	 double yieldFunction1(double sn, double i1);
	 double yieldFunction2(double sn, double i1, double q);
	 double yieldFunction3(double i1);
};
} // end namespace oofem
#endif // dustmat_h
