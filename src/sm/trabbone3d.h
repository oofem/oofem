/*
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

//   **************************************
//   *** CLASS TRABECULAR BONE MATERIAL ***
//   **************************************

#ifndef trabbone3d_h
#define trabbone3d_h

#include "structuralmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "matconst.h"
#include "matstatus.h"
#include "strainvector.h"


#include "linearelasticmaterial.h"
#include "dictionr.h"

#include "structuralms.h"
#include "cltypes.h"

namespace oofem {
/*
 * This class implements associated Material Status to TrabBone3D.
 * It is atribute of matStatusDictionary at every GaussPoint, for which this material
 * is active.
 */


/////////////////////////////////////////////////////////////////
//////////////////TRABECULAR BONE STATUS/////////////////////////
/////////////////////////////////////////////////////////////////

class TrabBone3DStatus : public StructuralMaterialStatus
{
protected:

    /////////////////////////////////////////////////////////////////
    // STATE VARIABLE DECLARATION
    //

    double kappa, tempKappa, dam, tempDam, tempPSED, tempTSED, tsed, beta;
    FloatArray tempPlasDef, plasDef, effectiveStress, tempEffectiveStress, plasFlowDirec, tempStrain;;
    FloatMatrix smtrx, tangentMatrix, SSaTensor;
    /// number of substeps in the last iteration 
    int nss;
    //densificator criterion
    double densG;
   

public:

    /////////////////////////////////////////////////////////////////
    // CONSTRUCTOR
    //

    TrabBone3DStatus(int n, Domain *d, GaussPoint *g);


    /////////////////////////////////////////////////////////////////
    // DESTRUCTOR
    //

    ~TrabBone3DStatus();


    /////////////////////////////////////////////////////////////////
    // OUTPUT PRINT
    // Prints the receiver state to stream

    void   printOutputAt(FILE *file, TimeStep *tStep);


    /////////////////////////////////////////////////////////////////
    // STATE VARIABLE
    // declare state variable access and modification methods

    double giveKappa();
    double giveTempKappa();
    double giveDam();
    double giveTempDam();
    double giveTempPSED();
    double giveTSED();
    double giveTempTSED();
    double giveBeta();
    int giveNsubsteps() {return nss;}
    double giveDensG(){return densG;}

    const FloatArray *givePlasDef();
    const FloatArray *giveTempPlasDef();
    const FloatArray *giveTempEffectiveStress();
    const FloatArray *givePlasFlowDirec();
    const FloatMatrix *giveTangentMatrix();
    const FloatMatrix *giveSmtrx();
    const FloatMatrix *giveSSaTensor();
  

    void setTempKappa(double al) { tempKappa = al; }
    void setKappa(double values){kappa = values;}
    void setTempDam(double da) { tempDam = da; }
    void setTempPSED(double pse) { tempPSED = pse; }
    void setTempTSED(double tse) { tempTSED = tse; }
    void setBeta(double be) { beta = be; }
    void setTempEffectiveStress(FloatArray &sc) { tempEffectiveStress = sc; }
    void setTempPlasDef(FloatArray &epsip) { tempPlasDef = epsip; }
    void setPlasFlowDirec(FloatArray &pfd) { plasFlowDirec = pfd; }
    void setSmtrx(FloatMatrix &smt) { smtrx = smt; }
    void setTangentMatrix(FloatMatrix &tmm) { tangentMatrix = tmm; }
    void setSSaTensor(FloatMatrix &ssa) { SSaTensor = ssa; }
    void setNsubsteps(int n)  { nss = n; }

    void setDensG(double g) { densG = g; }
    

    /////////////////////////////////////////////////////////////////
    // DEFINITION
    //

    const char *giveClassName() const { return "TrabBone3DStatus"; }
    classType   giveClassID() const { return TrabBone3DStatusClass; }


    /////////////////////////////////////////////////////////////////
    // INITIALISATION OF TEMPORARY VARIABLES
    // Initializes the temporary internal variables, describing the current state according to
    // previously reached equilibrium internal variables.

    virtual void initTempStatus();


    /////////////////////////////////////////////////////////////////
    // UPDATE VARIABLES
    // Update equilibrium history variables according to temp-variables.
    // Invoked, after new equilibrium state has been reached.

    virtual void updateYourself(TimeStep *);


    /////////////////////////////////////////////////////////////////
    // SAVE CONTEXT - INTERUPT/RESTART
    // saves current context(state) into stream
    // Only non-temp internal history variables are stored.
    // @param stream stream where to write data
    // @param obj pointer to integration point, which invokes this method
    // @return contextIOResultType.

    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);


    /////////////////////////////////////////////////////////////////
    // RESTORE CONTEXT - INTERUPT/RESTART
    // Restores context of receiver from given stream.
    // @param stream stream where to read data
    // @param obj pointer to integration point, which invokes this method
    // @return contextIOResultType.

    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};


/////////////////////////////////////////////////////////////////
////////////////TRABECULAR BONE 3D///////////////////////////////
/////////////////////////////////////////////////////////////////


class TrabBone3D : public StructuralMaterial
{
protected:

    /////////////////////////////////////////////////////////////////
    // STATE VARIABLE DECLARATION
    // declare material properties here

  double m1, m2, rho, eps0, nu0, mu0, expk, expl, sig0Pos, sig0Neg, chi0Pos,chi0, chi0Neg, tau0, expq, expp;
  double plasHardFactor, expPlasHard, expDam, critDam,pR;
    int printflag, abaqus, max_num_iter, max_num_substeps;
    double rel_yield_tol, strain_tol;
    // local coordinate system
    double x1,x2,x3,y1,y2,y3,z1,z2,z3;
    //densificator properties
    double  gammaL0, gammaP0, tDens, densCrit, rL,rP, gammaL, gammaP;
    //viscosity parameter
    double viscosity;
    //Hadi post-yield function
    double yR,kappaMax,kappaMin,kappaSlope,N,gMin, formulation;
    double hardFactor;

public:

    /////////////////////////////////////////////////////////////////
    // CONSTRUCTOR
    //

    TrabBone3D(int n, Domain *d);


    /////////////////////////////////////////////////////////////////
    // DESTRUCTOR
    //

    // NO DESTRUCTOR

    /////////////////////////////////////////////////////////////////
    // INITIALIZATION OF FUNCTION/SUBROUTINE
    //
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }
    double evaluateCurrentYieldStress(const double kappa);
    double evaluateCurrentPlasticModulus(const double kappa);
    double evaluateCurrentViscousStress(const double deltaKappa, TimeStep* atTime);
    double evaluateCurrentViscousModulus(const double deltaKappa, TimeStep* atTime);

    bool projectOnYieldSurface(double &tempKappa, FloatArray &tempEffectiveStress, FloatArray &tempPlasDef, const FloatArray &trialEffectiveStress, const FloatMatrix &elasticity, const FloatMatrix &compliance, TrabBone3DStatus *status,TimeStep *atTime, GaussPoint* gp, int lineSearchFlag);

  

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain,TimeStep* atTime, MaterialMode mode);

    void  constructPlasFlowDirec(FloatArray &answer,double &norm, FloatMatrix &fabric, FloatArray &F, FloatArray &S);
    void  constructDerivativeOfPlasFlowDirec(FloatMatrix &answer, FloatMatrix &fabric, FloatArray &F, FloatArray &S);
    double evaluatePlasCriterion(FloatMatrix &fabric, FloatArray &F, FloatArray &stress,double kappa, double deltaKappa, TimeStep* atTime);

    double computeDamageParam(double kappa);
    double computeDamageParamPrime(double kappa);


    double computeDamage(GaussPoint *gp, TimeStep *atTime);

    virtual void computeCumPlastStrain(double& kappa, GaussPoint *gp, TimeStep *atTime);

    void computePlasStrainEnerDensity(GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &totalStress);

    void computeDensificationStress(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime);

    // Construct anisotropic compliance tensor.
    void          constructAnisoComplTensor(FloatMatrix &answer);
    // Construct anisotropic stiffness tensor.
    void          constructAnisoStiffnessTensor(FloatMatrix &answer);

    // Construct anisotropic fabric tensor.
    void          constructAnisoFabricTensor(FloatMatrix &answer);
    void          constructAnisoFtensor(FloatArray &answer);

    void          constructStiffnessTransformationMatrix(FloatMatrix &answer);
 void          constructFabricTransformationMatrix(FloatMatrix &answer);
    // Construct Tensor to adjust Norm.
    void         constructNormAdjustTensor(FloatMatrix &answer);


    //Default implementation computes 3d stifness matrix using give3dMaterialStiffnessMatrix and
    //reduces it to 1d stiffness using reduce method described above.
    //Howewer, this reduction is quite time consuming and if it is possible,
    //it is recomended to overload this method and provide direct method for computing
    //particular stifness matrix.
    //@param answer stifness matrix
    //@param form material response form
    //@param mode material response mode
    //@param gp integration point, which load history is used
    //@param atTime time step (most models are able to respond only when atTime is current time step)

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm, MatResponseMode, GaussPoint * gp,
                                               TimeStep * atTime);

    //Computes the real stress vector for given total strain and integration point.
    //The total strain is defined as strain computed directly from displacement field at given time.
    //The stress independent parts (temperature, eigen strains) are subtracted in constitutive
    //driver.
    //The service should use previously reached equilibrium history variables. Also
    //it should update temporary history variables in status according to newly reached state.
    //The temporary history variables are moved into equilibrium ones after global structure
    //equlibrium has been reached by iteration process.
    //@param answer contains result
    //@param form material response form
    //@param gp integration point
    //@param reducedStrain strain vector in reduced form
    //@param tStep current time step (most models are able to respond only when atTime is current time step)

    virtual void giveRealStressVector(FloatArray & answer, MatResponseForm, GaussPoint *,
                                      const FloatArray &, TimeStep *);

    //Requests material mode capability.
    //@param mode material mode requested
    //@return nonzero if available

    virtual int hasMaterialModeCapability(MaterialMode);

    // Returns class name of the receiver.

    const char *giveClassName() const { return "TrabBone3D"; }

    // Returns classType id of receiver.

    classType giveClassID() const { return TrabBone3DClass; }

    //Initializes receiver acording to object description stored in input record.
    //The density of material is read into property dictionary (keyword 'd')

    IRResultType initializeFrom(InputRecord *ir);

    MaterialStatus *CreateStatus(GaussPoint *gp) const;

    /**
     * Returns the integration point corresponding value in Reduced form.
     * @param answer contain corresponding ip value, zero sized if not available
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);

    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);

    /**
     * Returns the type of internal variable (scalar, vector, tensor,...).
     * @param type determines the type of internal variable
     * @returns type of internal variable
     */
    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    /**
     * Returns the corresponding integration point  value size in Reduced form.
     * @param type determines the type of internal variable
     * @returns var size, zero if var not supported
     */
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
};
} //end namespace oofem
#endif
