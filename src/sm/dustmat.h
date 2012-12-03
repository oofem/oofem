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
public:
    enum stateFlagValues { DM_Elastic, DM_Unloading, DM_Yielding1, DM_Yielding2, DM_Yielding3 };

protected:
    /// Deviatoric of plastic strain
    StrainVector plasticStrain;
    StrainVector tempPlasticStrain;

    double q;
    double tempQ;
    double mVol;

    double bulkModulus;
    double shearModulus;
    double youngsModulus;

    /// Indicates the state (i.e. elastic, yielding, vertex, unloading) of the Gauss point
    int stateFlag;
    int tempStateFlag;

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

    void givePlasticStrain(StrainVector &answer) const { answer = plasticStrain; }
	 double giveVolumetricPlasticStrain() const { return plasticStrain.computeVolumetricPart(); }
    double giveQ() const { return q; }
    int giveStateFlag() const { return stateFlag; }

    void giveTempPlasticStrain(StrainVector &answer) const { answer = tempPlasticStrain; }
    double giveTempQ() const { return tempQ; }
    int giveTempStateFlag() const { return tempStateFlag; }

    void letTempPlasticStrainBe(const StrainVector &v) { tempPlasticStrain = v; }
    void letTempQBe(double v) { tempQ = v; }
    void letTempStateFlagBe(int v) { tempStateFlag = v; }

    void letPlasticStrainBe(const StrainVector &v) { plasticStrain = v; }
    void letQBe(double v) { q = v; }

    void setBulkModulus(double m) { bulkModulus = m; }
    void setShearModulus(double m) { shearModulus = m; }
    void setYoungsModulus(double m) { youngsModulus = m; }
    double giveBulkModulus() { return bulkModulus; }
    double giveShearModulus() { return shearModulus; }
    double giveYoungsModulus() { return youngsModulus; }

    void setMVol(double m) { mVol = m; }
};

/**
 * This class implements TODO
 * @author Jan Stransky
 */
class DustMaterial : public StructuralMaterial
{
protected:
    IsotropicLinearElasticMaterial *LEMaterial;

    double alpha;
    double beta;
    double lambda;
    double theta;
    double ft;
    double rEll;
    int hardeningType;
    double mStiff;
    double x0;
    double q0;
    double wHard;
    double dHard;
    double newtonTol;
    int newtonIter;

    double fFe(double i1);
    double fFeDI1(double i1);
    double fFeDI1DI1(double i1);
    double fFc(double rho, double i1, double q);
    double fX(double q);
    double fXdQ(double q);
    double yieldf1(double rho, double i1);
    double yieldf2(double rho, double i1, double q);
    double yieldf3(double i1);
    void solveQ0(double &q0);
    void computeAndSetBulkAndShearModuli(double &bulkModulus, double &shearModulus, GaussPoint *gp);
    void performStressReturn(GaussPoint *gp, StrainVector strain);
    void fM1(StrainVector &answer, const StressVector &stressDeviator, double rho, double i1, double q);
    void fM2(StrainVector &answer, const StressVector &stressDeviator, double rho, double i1, double q);
    void fM3(StrainVector &answer, const StressVector &stressDeviator, double rho, double i1, double q);
    double fH(double q, double tempQ);
    double fHdq(double tempQ);
    double fI1(double q, double tempQ, double i1, double bulkModulus);
    double fI1DQ(double tempQ, double bulkModulus);
    void performF1return(double i1, double rho, GaussPoint *gp);
    void performF2return(double i1, double rho, GaussPoint *gp);
    void computeQFromPlastVolEps(double &answer, double q, double volumetricPlasticStrain);
    double dGamma2(double tempQ, double q, double i1, double bulkModulus);
    double dGamma2DQ(double tempQ, double q, double i1, double bulkModulus);
    double fTempR2(double tempQ, double q, double i1, double rho, double bulkModulus, double shearModulus);
    double fTempR2Dq(double tempQ, double q, double i1, double rho, double bulkModulus, double shearModulus);

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

    virtual int setIPValue(const FloatArray value, GaussPoint *gp, InternalStateType type);

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

    double giveQ0() { return q0; }
};
} // end namespace oofem
#endif // dustmat_h
