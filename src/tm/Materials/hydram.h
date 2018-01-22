/* TODO:
 *    4/2/2004
 *     need to return a matrix of derivatives
 *      dHeat/dT, dHeat/dw
 *      dWater/dT, dWater/dw
 *     (or something like that) to enable quadratic convergence of the global iteration with respect to the heat and water internal source.
 *    now I have only double dHeatdT(ksi, T, dt) - internal, no water dependence. To master material, the hydration model must provide the derivatives based on the internal state of hydration.
 *
 *    12/1:
 *    - track versions - save intermediate versions of modified files in separate folders!
 *    11/1:
 *    - where to choose, which HydrationModel to create? HydrationModelInterface->constructor can't know the type of hydram
 *      but HydrationModelInterface->initializeFrom(ir) can select and create the appropriate hydration model
 *      maybe it should be added to usrdefsub.C
 *    x try to eliminate HydrationModelStatusInterface, to enable use of different status classes without redefining (sm / tm)
 *      will not do, master status must have reference to HydrationModelStatus
 *      tm status will be defined by simply inheriting from TransportMaterialStatus && HydrationModelStatusInterface
 *    - where to put gp hydration degree access routines? Master material->hmInterface->giveTempKsi or hydrationModel->giveTempHydrationDegree(gp)
 *    10/1:
 *    - cleanly: base class HydrationModelStatus + HydrationModel, defining interface data and services, derive HellmichHydrationModel etc.
 *    - water content influence on hydration -
 *      moisture consumption (stoichiometry of c-s-h?) -> computeInternalSourceVector(+heat,-water content)
 *    - affinity function from CEMHYD3D results: ksi(t)-> dksidt(ksi)->A(ksi) - \
 *    - real-time from CEMHYD3D (better not all IP, only selected range)      -/ needs CEMHYD3D input preprocessor for generating microstructure (VCCTL?) or input microstructure data
 *    - or add alcali reaction
 *    - hydration interpolated from CEMHYD3D with moisture and temperature
 *      needs several hydration model / master material managed statuses and services for interpolating hydration degree for given coordinates
 */

/*
 * USAGE:
 *
 * use multiple inheritance from HydrationModelInterface and HydrationModelStatusInterface.
 *
 * add HydrationModelInterface to material:
 * constructor
 * initializeFrom
 * saveContext / restoreContext
 * TransportMaterial::giveCharacteristicValue
 * use HMI::updateInternalState(tempStateVector, GaussPoint, TimeStep) to update the hydration status
 *
 * add HydrationModelStatusInterface to status:
 * constructor
 * updateYourself
 * printOutputAt
 * add giveInterface
 *
 * then the following services are available:
 * giveHydrationDegree(GaussPoint, TimeStep, ValueModeType)
 * computeInternalSourceVector(sourceVector, GaussPoint, TimeStep, ValueModeType)
 * computeIntSourceLHSCoeff(MAtResponseMode, GaussPoint, TimeStep)
 *
 * missing - consistent tangent, see TODO
 *
 */

// class HydrationModel - returns hydration degree ... dksi(ksi, dt, T) ... localResidual, Affinity
// ( - for non-isothermal analysis               hydrationHeat, capacitycoeff(ksi, dt, T) ... dksidT, dksidh, dAdksi

// Concrete hydration model - Hellmich
// Material extension, derived from Material to use MateriaStatus services

// material->computeInternalSourceVector(val, gp, tStep, mode)
//                                       answer,     chartypemode - total, velocity, acceleration, incremental
// = le*dksi

// capacity ... material->giveCharacteristicValue(rmode, gp, tStep)
//              rmode: material responsemode = capacity
// = rc - le*dksidT


#ifndef hydram_h
#define hydram_h

#include "floatarray.h"
#include "floatmatrix.h"
#include "timestep.h"
#include "material.h"
#include "interface.h"
#include "valuemodetype.h"

#include <memory>

///@name Input fields for HydrationModel
//@{
#define _IFT_HydrationModel_Name "hydrationmodel" ///@todo Is this right? double check. And isn't REGISTER_Material() missing?
#define _IFT_HydrationModel_hydration "hydration"
#define _IFT_HydrationModel_c60mix "c60mix"
#define _IFT_HydrationModel_timeScale "timescale"
#define _IFT_HydrationModel_hheat "hheat"
#define _IFT_HydrationModel_cv "cv"
#define _IFT_HydrationModel_water "water"
//@}

///@name Input fields for HydrationModelInterface
//@{
#define _IFT_HydrationModelInterface_hydration "hydration"
#define _IFT_HydrationModelInterface_castAt "castat"
//@}

namespace oofem {
#define ROOT_PRECISION_DKSI 1e-14
#define BINARY_TREE_STEPS 2

// default time steps for evaluating the hydration degree increment
#define HYDRATION_MAXSTEP0 3600 //600
#define HYDRATION_MAXSTEP1 86400 //3600


// =========== Hydration Model Status class ============
/**
 * This class implements associated Status to HydrationModel.
 * It is attribute of owner material status for each GaussPoint, for which that material
 * is active.
 */
class HydrationModelStatus : public MaterialStatus
{
protected:
    // hydration degree at beginning of current time step
    double hydrationDegree;
    // hydration degree at end of current time step (or during equilibrium iteration)
    double tempHydrationDegree;

public:
    HydrationModelStatus(int n, Domain * d, GaussPoint * g);
    virtual ~HydrationModelStatus() { }

    /// Returns the temp hydration degree.
    double giveTempHydrationDegree() { return tempHydrationDegree; }
    /// Returns the non-temp hydration degree. Used for step restart and postprocessing.
    double giveHydrationDegree() { return hydrationDegree; }
    void setHydrationDegree(double v) { hydrationDegree = v; }
    void setTempHydrationDegree(double v) { tempHydrationDegree = v; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    // --- identification and auxiliary functions ---
    virtual const char *giveClassName() const { return "HydrationModelStatus"; }
};


// =========== Hydration Model class ============

enum FindRootMethod { frRegula = 1, frBinTree = 2, frMixed = 3 };
enum MixtureType { mtLafarge = 1, mtHuber = 2, mtC60, mtC100 };

// derived from material to use status services
class HydrationModel : public Material
{
protected:
    /// Used concrete mixture
    MixtureType mixture;
    /// Time step lenghts at zero and complete hydration
    double hydrationStartMaxStep, hydrationEndMaxStep;
    //!!! initial hydration degree - set in initialize From, but not used
    double initialHydrationDegree;
    /// time scale - used for time input in other units than seconds
    double timeScale;

    // === Material parameters ===
    double aa, ///< Normalized chemical affinity regression function coefficients.
           ba,
           ca,
           da,

           e0, ///< ksi_0.
           ear, ///< Activation term [K].
           le, ///< Latent heat [kJ/m3].
           cv, ///< Input cement content kg/m3 for evaluation of total water consumption.
           we; ///< Total water consumption for hydration [kg/m3].

    // === Hydration degree increment evaluation ===
    // auxiliary values to enable external root finding method without passing them as parameters in each call
    //!!! possible problem for parallel computation, performance???
    double auxksi, auxdt, auxT, auxh;
    double localResidual(double dks); // G(dksi) 4.12

    // Auxiliary functions - root finding
    ///
    double regulafindroot();
    double bintreefindroot();
    double mixedfindroot();

    // === Material functions ===
    /// Returns the normalized chemical affinity A~(ksi) [1/s].
    double affinity(double ksi);
    /// Returns the derivation of chemical affinity dA~/dksi(ksi).
    double dAdksi(double ksi);
    double dksidT(double ksi, double T, double h, double dt);
    double dksidh(double ksi, double T, double h, double dt);

    /// Computes and returns the derivatives of the material-generated Internal Source with respect to the tm state vector.
    double computeIntSource(const FloatArray &vec, GaussPoint *gp, TimeStep *tStep, MatResponseMode rmode);
    /**
     * Computes and returns hydration degree increment for given ksi, T [K], dt [s].
     * Called by updateInternalState(val, gp, tStep)
     * @note Formerly dksi, maybe should be changed to use gp & timestep.
     */
    double computeHydrationDegreeIncrement(double ksi, double T, double h, double dt);

public:
    FindRootMethod useFindRoot;
    /// Constructor
    HydrationModel();
    /// Constructor setting the mixture type and root-finding method
    HydrationModel(MixtureType mix, FindRootMethod usefr);
    /// Destructor
    virtual ~HydrationModel() { }


    // === Hydration model services ===

    /// Sets the mixture type and appropriate material parameters
    void setMixture(MixtureType mix);

    /**
     * Initializes the hydration model according to object description stored in input record.
     * Called from master material initializeFrom, hydrationModelInterface initializeFrom selects the appropriate hydration model type.
     *
     * Not a standard material - initializes from master material record, doesn't call parent initializeFrom
     * Use hm_ prefix in parameter names to avoid confusion with master material parameters
     */
    virtual IRResultType initializeFrom(InputRecord *ir);

    /**
     * Returns the hydration degree at end of TimeStep tStep in given integraion point.
     * The value is obtained from hydration status, in integration point, on material level, or ... interpolated from several statuses
     * @param gp integration point
     * @param tStep solution step
     * @param mode value mode VM_Incremental or VM_Total
     * @return hydration degree or increment in given gp
     */
    double giveHydrationDegree(GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    // how to determine whether hydration degree is uptodate?
    // - in tm, this can be ensured by calling computeHydrationDegreeIncrement when updateInternalState is used
    // - in sm, the tm state changes only at the start of each step (now updateAuxState - 'auxiliary', because the data is copied from tm analysis)
    // where to get temperature / moisture? from gp->status->giveTempStateVector(val) => val[T,w]
    //   => implement updateInternalState here, that is called from master material updateInternalState/updateAuxState
    // use giveTempStateVector in HellmichMaterial
    /**
     * Updates internal state of material according to new state vector - computes the hydration degree for time tStep.
     * @param vec new state vector
     * @param gp integration point
     * @param tStep solution step
     */
    virtual void updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *tStep);

    /**
     * Returns material status of receiver in given integration point.
     * In hydration model, the status is obtained from master gp status (status->giveHydrationModelStatus)
     * If status does not exist, a new one is created.
     * @param gp Returns reference to hydration model status belonging to integration
     * point gp.
     * @return hydration model status associated with given integration point.
     */
    virtual MaterialStatus *giveStatus(GaussPoint *gp) const;

    /// Returns generated heat for given gp [kJ/m3], eventually water consumption
    void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode);
    // } end new 5.1.2004
    /// Returns coefficients for LHS contribution from internal sources (dHeat/dT, dWaterSource/dw) for given temp state vector.
    virtual double giveCharacteristicValue(const FloatArray &vec, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep);
    // --- identification and auxiliary functions ---
    virtual const char *giveInputRecordName() const { return _IFT_HydrationModel_Name; }
    virtual const char *giveClassName() const { return "HydrationModel"; }

protected:
    /// Creates and returns new HydrationModelStatus instance
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
};

// =========== Interfaces for materials using the hydration model ============
class HydrationModelStatusInterface : public Interface
{
protected:
    /// Reference to associated hydration model status
    std :: unique_ptr< HydrationModelStatus > hydrationModelStatus;
public:
    /// Constructor. Nulls the hydrationModelStatus pointer.
    HydrationModelStatusInterface() {}
    /// Destructor. Deletes the associated hydration model status.
    virtual ~HydrationModelStatusInterface() {}

    /// Returns the associated hydration model status.
    HydrationModelStatus *giveHydrationModelStatus() { return hydrationModelStatus.get(); }
    /// Sets the associated hydration model status. Analogue to gp->setMaterialStatus.
    void setHydrationModelStatus(HydrationModelStatus *s) { hydrationModelStatus.reset(s); }

    /// Updates the equilibrium variables to temporary values.
    void updateYourself(TimeStep *tStep);
    /// Outputs the status variables
    void printOutputAt(FILE *file, TimeStep *tStep);
};

class HydrationModelInterface : public Interface
{
    // interface for material access to HydrationModel
protected:
    /// Reference to the associated hydrationModel instance
    std :: unique_ptr< HydrationModel > hydrationModel;
    /// Material cast time - start of hydration
    double castAt;
    /// Constant hydration degree for analysis without hydration model
    double constantHydrationDegree;

public:
    /// Returns the associated hydration model.
    HydrationModel *giveHydrationModel() { return hydrationModel.get(); }
    /// Destructor. Deletes the associated hydration model.
    virtual ~HydrationModelInterface() {}
    /**
     * Creates and initializes the hydration model according to object description stored in input record.
     * The parent class instanciateFrom method is not called here.
     */
    IRResultType initializeFrom(InputRecord *ir);

    contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL)
    {
        if ( hydrationModel ) {
            hydrationModel->saveContext(stream, mode, obj);
        }

        return CIO_OK;
    }
    contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL)
    {
        if ( hydrationModel ) {
            hydrationModel->restoreContext(stream, mode, obj);
        }

        return CIO_OK;
    }

    /**
     * Calls hydrationModel->updateInternalState, if the material is already cast.
     * In case the cast time lies within the span of current timestep, the timestep increment is set to (time-castAt).
     * @param vec New state vector.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *tStep);
    /**
     * Returns the hydration degree at end of TimeStep tStep in given integration point.
     * The value is obtained from gp hydration status via the hydration model or the constantHydrationDegree value is returned.
     * @param gp Integration point.
     * @param tStep Solution step.
     * @param mode Value mode VM_Incremental or VM_Total.
     * @return Hydration degree or increment in given gp.
     */
    double giveHydrationDegree(GaussPoint *gp, TimeStep *tStep, ValueModeType mode);
};
} // end namespace oofem
#endif // hydram_h
