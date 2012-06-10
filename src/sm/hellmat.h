// (ChemoViscoPlastic material model (short-term creep, long-term creep, microcracking - multisurface plasticity (Drucker-Prager with hardening + tension cut-off)

// Chemical hardening capability included (for constant material temperature)
//  Can be used in coupled analysis, reading only temperature from tm domain,
//    computing hydration degree increment here
//  Or copying both T and hydration degree from thermochemical analysis

#ifndef hellmat_h
#define hellmat_h

// oofem includes
#include "flotarry.h"
#include "flotmtrx.h"
#include "structuralms.h"
#include "structuralmaterial.h"
#include "isolinearelasticmaterial.h"

#ifdef __TM_MODULE
// Hellmat includes
 #include "../tm/hemotkmat.h"     // for inverse sorption isotherm function
 #include "hydram.h"                     // hydration model
#endif

namespace oofem {
#ifdef __TM_MODULE

class GaussPoint;

// Output of all key functions
// #define VERBOSE_HELLMAT
// routine-specific verbose output
//#define VERBOSEFINDROOT // plastic increment findroot statistic
//#define VERBOSE_TANGENT // consistent tangent statistic

// substitution for 0 hydration degree for determining non-singular stiffness matrix
 #define HYDRATION_MINDEGREE 1e-3 // 1e-6 deteriorates the iterative solver efficiency for consequent cast

// === additional output file ===
 #define OUTPUTFILE_SS "out/plastic.out"

// settings for local plastic Newton iteration
// approxNewton solution parameter - dx for derivative evaluation
// also used for dRcdchi1
 #define NEWTON_DERIVATIVEDX 1e-12
// relative error bound for plastic local iteration (d f1)
 #define ROOT_PRECISION_DLAMBDA 1e-10
// absolute error bound for plastic local iteration (D f1) - has same scale as stress ... 1e-8/1e6 = 1e-14
 #define FINDROOT_SMALLNUM 1e-8

// active plastic surface
enum ActiveSurface { asNone = 0, asDP = 1, asTC = 2, asCorner = 3 };
// material options mask {
typedef unsigned int MaterialOptions;
// options
 #define moIsothermal 1         // DEFAULT constant temperature in space
 #define moHydration 2          // DEFAULT compute hydration
 #define moShrinkage 4          // DEFAULT autogenous shrinkage
 #define moHardening 8          // DEFAULT plastic hardening of Drucker-Prager surface
 #define moSqrtHardeningLaw 16  // DEFAULT for smooth stress-strain curve, otherwise quadratic hardening is used
 #define moCreep 32             // DEFAULT creep
 #define moApproxNewton 64      //         approximates the derivative in local plastic iteration with finite dx
 #define moPlasticity 128       // DEFAULT use plasticity
 #define moPlotStressStrain 256 //         output stress-strain data
 #define moPlotStressStrainIter 512 //        for each iteration
 #define moComputedl 1024       //         direct computation of plastic increment (quadratic equation); if not set, newton iteration is used
 #define moDeviatoricCreepE 2048 //         only deviatoric part of stress is taken into account for flow creep strain increment evaluation
 #define moDeviatoricCreepF 4096 //         -"- for viscous creep strain increment evaluation
 #define moLinearEModulus   8192 //        force linear relation E(ksi)
 #define moHumidityStrain  16384 //        reversible humidity dilatation
 #define moDryingShrinkage 32768 //        include drying shrinkage and stress-induced temperature dilatation in the flow creep strain / that's a shit!
// option sets
 #define moHasShrinkage 6       //         moHydration + moShrinkage
 #define moDefault 191          //         1+2+4+8+16+32+128 - hydrating plasticity with sqrt hardening, shrinkage, creep.
 #define moMaterialLevel 0      //         nonIsothermal, no creep, no plasticity, ??? hydration - depends
// }

// ==== Aging isotropic linear elastic material class allowing to set E ====
class AgingIsoLEMaterial : public IsotropicLinearElasticMaterial
{
public:
    /**
     * Constructor. Creates a new AgingIsoLEMaterial class instance
     * with given number belonging to domain d.
     * @param n material model number in domain
     * @param d domain which receiver belongs to
     */
    AgingIsoLEMaterial(int n, Domain *d, double E, double nu);
    void setE(double newE);
};

// ================== STATUS =======================

// auxiliary data classes for status modules - noniso, creep, plast; noniso includes creep, thats OK
/// Class for saving plasticity status data
class PlastData
{
public:
    // equilibrated
    FloatArray plasticStrainVector;
    double hardeningVar;
    // temporary
    FloatArray tempPlasticStrainVector;
    double tempHardeningVar;

    // auxiliary - can be determined from stress & hardeningVar - only temp,
    // must be saved in context or maybe evaluated in restoreContext to enable full postprocessing
    ActiveSurface activeSurface;
    double flowIncrement; // dl1 - stored for tangent stiffness matrix evaluation
    FloatArray trialStressVector; // stored to enable stress return output

    PlastData();
    void initTempStatus(GaussPoint *gp);
    void updateYourself();

    contextIOResultType saveContext(DataStream *stream, ContextMode mode);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode);
};

/// Class for storing creep status data
class CreepData
{
public:
    // viscousStrain - short-term creep, flowStrain - long-term creep
    FloatArray viscousStrainVector, flowStrainVector, tempViscousStrainVector, tempFlowStrainVector;
    CreepData() : viscousStrainVector(), flowStrainVector(), tempViscousStrainVector(), tempFlowStrainVector() { }
    void initTempStatus(GaussPoint *gp);
    void updateYourself();
    contextIOResultType saveContext(DataStream *stream, ContextMode mode);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode);
};
/// Class for storing non-isothermal or material-level status data
class NonisoData
{
public:
    // temperature and relative humidity auxiliary here, mapped at step start from tm domain
    double temperature;  // at step start
    double humidity;
    double tempTemperature; // at step end
    double tempHumidity;
    // gp initial temperature for determining thermal strains - must be included in Material's EigenStrains,
    // because StructuralElement's computeTemperatureStrainVector accounts only
    // for prescribed thermal field loads.
    /* May be copied from thermal analysis' initial conditions, currently is set at first solution step,
     * i.e. to the equilibrated temperature at cast. Maybe should be initialized at material setting point,
     * because no stresses are induced by expansion of the concrete mixture
     */
    // x StructuralMaterial::referenceTemperature initialized to 0.0, protected, no set method!
    double initialTemperature;
    // Temperature-dependent creep status
    // allocated even in case that creep is not used
    double gamma0;      // base value of microprestress
    double viscousSlip;
    double tempViscousSlip;
    double viscosity;
    /**
     * Time at which stored auxiliary non-isothermal status values were computed. This concerns only nonisothermal status values -
     * - temperature, viscousslip, viscosity; hydration
     * Used to avoid multiple computation of hydration degree, flow creep viscosity and other
     * time-dependent values during equilibrium iteration.
     * Need to be updated from previous values when iteration is restarted with shorter step length!
     * 25/01/2004 changed to flag - see .C for explanation
     */
    int auxStatusUpdateFlag;
    NonisoData();
    // status update - temp->equilib, sets update flag to 0
    void updateYourself();
    contextIOResultType saveContext(DataStream *stream, ContextMode mode);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode);
};

class HellmichMaterialStatus : public StructuralMaterialStatus, public HydrationModelStatusInterface
{
protected:
    // 10.1.2004 hydration degree not stored here, but in HydrationModelStatus, accessed via HydrationModelStatusInterface
    virtual Interface *giveInterface(InterfaceType);
    /**
     * Non-isothermal analysis
     * non-iso values are stored in material integration point in isothermal case
     * includes also the time stamp when the variables are evaluated to avoid computing them several times
     */
    NonisoData *nonisoData;
    /// Plasticity
    PlastData *plastData;
    /// Creep strain vectors
    CreepData *creepData;
public:
    HellmichMaterialStatus(int n, Domain *d, GaussPoint *g);
    virtual ~HellmichMaterialStatus();
    // === Status variables access routines ===
    // all access routines are defined, appropriate data class is checked and used
 #ifdef VERBOSE_HELLMAT // error-checking status access
    // Non-isothermal
    double giveTemperature() {
        if ( !nonisoData ) { _error("giveTemperature: missing nonisoData."); }

        return nonisoData->temperature;
    }
    double giveTempTemperature() {
        if ( !nonisoData ) { _error("giveTempTemperature: missing nonisoData."); }

        return nonisoData->tempTemperature;
    }
    double giveHumidity() {
        if ( !nonisoData ) { _error("giveHumidity: missing nonisoData."); }

        return nonisoData->humidity;
    }
    double giveTempHumidity() {
        if ( !nonisoData ) { _error("giveTempHumidity: missing nonisoData."); }

        return nonisoData->tempHumidity;
    }
    double giveInitialTemperature() {
        if ( !nonisoData ) { _error("giveInitialTemperature: missing nonisoData."); }

        return nonisoData->initialTemperature;
    }
    void setInitialTemperature(double v);
    void setTempTemperature(double v) {
        if ( !nonisoData ) { _error("setTempTemperature: missing nonisoData."); }

        nonisoData->tempTemperature = v;
    }
    void setTempHumidity(double v) {
        if ( !nonisoData ) { _error("setTempHumidity: missing nonisoData."); }

        nonisoData->tempHumidity = v;
    }
    // Non-isothermal creep
    double giveGamma0() {
        if ( !nonisoData ) { _error("giveGamma0: missing nonisoData."); }

        return nonisoData->gamma0;
    }
    double giveViscousSlip() {
        if ( !nonisoData ) { _error("giveViscousSlip: missing nonisoData."); }

        return nonisoData->viscousSlip;
    }
    double giveTempViscousSlip() {
        if ( !nonisoData ) { _error("giveTempViscousSlip: missing nonisoData."); }

        return nonisoData->tempViscousSlip;
    }
    double giveViscosity() {
        if ( !nonisoData ) { _error("giveViscosity: missing nonisoData."); }

        return nonisoData->viscosity;
    }
    void setGamma0(double v) {
        if ( !nonisoData ) { _error("setGamma0: missing nonisoData."); }

        nonisoData->gamma0 = v;
    }
    void setViscousSlip(double v) {
        if ( !nonisoData ) { _error("setViscousSlip: missing nonisoData."); }

        nonisoData->viscousSlip = v;
    }
    void setTempViscousSlip(double v) {
        if ( !nonisoData ) { _error("setTempViscousSlip: missing nonisoData."); }

        nonisoData->tempViscousSlip = v;
    }
    void setViscosity(double v) {
        if ( !nonisoData ) { _error("setViscosity: missing nonisoData."); }

        nonisoData->viscosity = v;
    }

    // Plasticity
    void givePlasticStrainVector(FloatArray &answer) const {
        if ( !plastData ) { _error("givePlasticStrainVector: missing plastData."); }

        answer = plastData->plasticStrainVector;
    }
    void giveTempPlasticStrainVector(FloatArray &answer) const {
        if ( !plastData ) { _error("giveTempPlasticStrainVector: missing plastData."); }

        answer = plastData->tempPlasticStrainVector;
    }
    double giveHardeningVar() { if ( plastData ) { return plastData->hardeningVar; } else { return 0.; } }
    double giveTempHardeningVar() { if ( plastData ) { return plastData->tempHardeningVar; } else { return 0.; } }
    void giveTrialStressVector(FloatArray &answer) const {
        if ( !plastData ) { _error("giveTrialStressVector: missing plastData."); }

        answer = plastData->trialStressVector;
    }
    ActiveSurface giveActiveSurface() { if ( !plastData ) { return asNone; } else { return plastData->activeSurface; } }
    double giveFlowIncrement() {
        if ( !plastData ) { _error("giveFlowIncrement: missing plastData."); }

        return plastData->flowIncrement;
    }

    void setPlasticStrainVector(const FloatArray &v) {
        if ( !plastData ) { _error("setPlasticStrainVector: missing plastData."); }

        plastData->plasticStrainVector = v;
    }
    void setTempPlasticStrainVector(const FloatArray &v) {
        if ( !plastData ) { _error("setTempPlasticStrainVector: missing plastData."); }

        plastData->tempPlasticStrainVector = v;
    }
    void setHardeningVar(double v) {
        if ( !plastData ) { _error("setHardeningVar: missing plastData."); }

        plastData->hardeningVar = v;
    }
    void setTempHardeningVar(double v) {
        if ( !plastData ) { _error("setTempHardeningVar: missing plastData."); }

        plastData->tempHardeningVar = v;
    }
    void setTrialStressVector(const FloatArray &v) {
        if ( !plastData ) { _error("setTrialStressVector: missing plastData."); }

        plastData->trialStressVector = v;
    }
    void setActiveSurface(ActiveSurface as) {
        if ( !plastData ) { _error("setActiveSurface: missing plastData."); }

        plastData->activeSurface = as;
    }
    void setFlowIncrement(double v) {
        if ( !plastData ) { _error("setFlowIncrement: missing plastData."); }

        plastData->flowIncrement = v;
    }

    // Creep strain vectors
    void giveViscousStrainVector(FloatArray &answer) const {
        if ( !creepData ) { _error("giveViscousStrainVector: missing creepData."); }

        answer = creepData->viscousStrainVector;
    }
    void giveTempViscousStrainVector(FloatArray &answer) const {
        if ( !creepData ) { _error("giveTempViscousStrainVector: missing creepData."); }

        answer = creepData->tempViscousStrainVector;
    }
    void giveFlowStrainVector(FloatArray &answer) const {
        if ( !creepData ) { _error("giveFlowStrainVector: missing creepData."); }

        answer = creepData->flowStrainVector;
    }
    void giveTempFlowStrainVector(FloatArray &answer) const {
        if ( !creepData ) { _error("giveTempFlowStrainVector: missing creepData."); }

        answer = creepData->tempFlowStrainVector;
    }

    void setViscousStrainVector(const FloatArray &v) {
        if ( !creepData ) { _error("setViscousStrainVector: missing creepData."); }

        creepData->viscousStrainVector = v;
    }
    void setTempViscousStrainVector(const FloatArray &v) {
        if ( !creepData ) { _error("setTempViscousStrainVector: missing creepData."); }

        creepData->tempViscousStrainVector = v;
    }
    void setFlowStrainVector(const FloatArray &v) {
        if ( !creepData ) { _error("setFlowStrainVector: missing creepData."); }

        creepData->flowStrainVector = v;
    }
    void setTempFlowStrainVector(const FloatArray &v) {
        if ( !creepData ) { _error("setTempFlowStrainVector: missing creepData."); }

        creepData->tempFlowStrainVector = v;
    }
 #else // standard status access
       // Non-isothermal
    double giveTemperature() { if ( nonisoData ) { return nonisoData->temperature; } else { return 0; } }
    double giveTempTemperature() { if ( nonisoData ) { return nonisoData->tempTemperature; } else { return 0; } }
    double giveHumidity() { if ( nonisoData ) { return nonisoData->humidity; } else { return 0; } }
    double giveTempHumidity() { if ( nonisoData ) { return nonisoData->tempHumidity; } else { return 0; } }
    double giveInitialTemperature() { if ( nonisoData ) { return nonisoData->initialTemperature; } else { return 0; } }
    void setInitialTemperature(double v);
    void setTempTemperature(double v) { if ( nonisoData ) { nonisoData->tempTemperature = v; } }
    void setTempHumidity(double v) { if ( nonisoData ) { nonisoData->tempHumidity = v; } }

    // Non-isothermal creep
    double giveGamma0() { if ( nonisoData ) { return nonisoData->gamma0; } else { return 0; } }
    double giveViscousSlip() { if ( nonisoData ) { return nonisoData->viscousSlip; } else { return 0; } }
    double giveTempViscousSlip() { if ( nonisoData ) { return nonisoData->tempViscousSlip; } else { return 0; } }
    double giveViscosity() { if ( nonisoData ) { return nonisoData->viscosity; } else { return 0; } }
    void setGamma0(double v) { if ( nonisoData ) { nonisoData->gamma0 = v; } }
    void setViscousSlip(double v) { if ( nonisoData ) { nonisoData->viscousSlip = v; } }
    void setTempViscousSlip(double v) { if ( nonisoData ) { nonisoData->tempViscousSlip = v; } }
    void setViscosity(double v) { if ( nonisoData ) { nonisoData->viscosity = v; } }

    // Plasticity
    void givePlasticStrainVector(FloatArray &answer) const { if ( plastData ) { answer = plastData->plasticStrainVector; } }
    void giveTempPlasticStrainVector(FloatArray &answer) const { if ( plastData ) { answer = plastData->tempPlasticStrainVector; } }
    double giveHardeningVar() { if ( plastData ) { return plastData->hardeningVar; } else { return 0.; } }
    double giveTempHardeningVar() { if ( plastData ) { return plastData->tempHardeningVar; } else { return 0.; } }
    void giveTrialStressVector(FloatArray &answer) const { if ( plastData ) { answer = plastData->trialStressVector; } }
    ActiveSurface giveActiveSurface() { if ( !plastData ) { return asNone; } else { return plastData->activeSurface; } }
    double giveFlowIncrement() { if ( plastData ) { return plastData->flowIncrement; } else { return 0.; } }

    void setPlasticStrainVector(const FloatArray &v) { if ( plastData ) { plastData->plasticStrainVector = v; } }
    void setTempPlasticStrainVector(const FloatArray &v) { if ( plastData ) { plastData->tempPlasticStrainVector = v; } }
    void setHardeningVar(double v) { if ( plastData ) { plastData->hardeningVar = v; } }
    void setTempHardeningVar(double v) { if ( plastData ) { plastData->tempHardeningVar = v; } }
    void setTrialStressVector(const FloatArray &v) { if ( plastData ) { plastData->trialStressVector = v; } }
    void setActiveSurface(ActiveSurface as) { if ( plastData ) { plastData->activeSurface = as; } }
    void setFlowIncrement(double v) { if ( plastData ) { plastData->flowIncrement = v; } }

    // Creep strain vectors
    void giveViscousStrainVector(FloatArray &answer) const { if ( creepData ) { answer = creepData->viscousStrainVector; } }
    void giveTempViscousStrainVector(FloatArray &answer) const { if ( creepData ) { answer = creepData->tempViscousStrainVector; } }
    void giveFlowStrainVector(FloatArray &answer) const { if ( creepData ) { answer = creepData->flowStrainVector; } }
    void giveTempFlowStrainVector(FloatArray &answer) const { if ( creepData ) { answer = creepData->tempFlowStrainVector; } }

    void setViscousStrainVector(const FloatArray &v) { if ( creepData ) { creepData->viscousStrainVector = v; } }
    void setTempViscousStrainVector(const FloatArray &v) { if ( creepData ) { creepData->tempViscousStrainVector = v; } }
    void setFlowStrainVector(const FloatArray &v) { if ( creepData ) { creepData->flowStrainVector = v; } }
    void setTempFlowStrainVector(const FloatArray &v) { if ( creepData ) { creepData->tempFlowStrainVector = v; } }
 #endif // status access

    /// Returns the material options mask. Retrieves options from material or uses moMaterialLevel for gp number 0
    MaterialOptions giveMaterialOptions();

    /**
     * Returns and sets the auxiliary (non-isothermal) status values update flag.
     * Used to ensure that initAuxStatus is run only once per step.
     * Causes error if nonisoData is not initialized.
     */
    int giveUpdateFlag() { if ( nonisoData ) { return nonisoData->auxStatusUpdateFlag; } else { return 1; } } // no initialization needed in case noniso data is not present
    void setUpdateFlag(int v) { if ( nonisoData ) { nonisoData->auxStatusUpdateFlag = v; } } // no error in case of isothermal analysis - don't need init

    // === Status services ===
    // Call data init/update according to material options
    virtual void initTempStatus(); // equilibrated -> temp - start of new iteration
    virtual void updateYourself(TimeStep *atTime);

    virtual void printOutputAt(FILE *file, TimeStep *atTime);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // identification
    virtual const char *giveClassName() const { return "HellmichMaterialStatus"; }
    virtual classType giveClassID() const { return ( classType ) HellmichMaterialStatusClass; }
};


// ================== MATERIAL =====================
class HellmichMaterial : public StructuralMaterial, public HydrationModelInterface
{
    // === Material parameters ===
    //
    // Lafarge mixture (pg. 48)
    // c = 380kg/m3 cement content
    // a/c = 4.79 (1820kg/m3 calcareous aggregate 0-8)
    // w/c = 0.6
    // f_c,28 = 39.6 MPa
    // (loading at 28 days)

protected:
    MixtureType mixture;
    // material options mask
    MaterialOptions options;

    // index of output stress/strain component and gp
    int pssIndex, pssElement, pssGaussPoint;
    /// time scale - used for time input in hours instead of seconds
    //!!! should be implemented on domain level?
    double timeScale;

    // === input material parameters ===
    double ny, // Hellmich 3.2.4.2 (pg 38)
           kappa, //   f_b-biaxial c/ f_c-uniaxial
           omega, //   f_y-yield c / f_c
           delta, //   f_t-tensile / f_c
           epscu, //   strain at f_c,8
           fc, //   f_c,8
           ae, //   E,8
           kshr, // moisture content shrinkage coefficient Kappa: eshr = kshr * (1-h^3)
           ashr, // shrinkage A.15
           bshr,
           ksi0; //   percolation threshold

    // creep 3.5
    double modulusH, // microprestress force softening modulus
           ur, // U/R (activation energy / universal gas const)
           refT, // reference temperature for long-term creep
           jv, // viscous creep compliance J,v,8
           tw, //   short-term creep characteristic time t,w,8
               //   (p=2, q=1; c = c' ... no effect)
           c;  // drying shrinkage coefficient - need to take from dT scale to gamma scale

    // auxiliary constant material parameters - computed at init
    double alpha, // Drucker-Prager surface slope
           kDP0, // initial D-P surface parameter
           auxk, //
           auxkdp, // kDP = auxkdp*Rc
           auxd, // auxiliary coeff for df1/dlambda
           chi1u, // chi1 at fc
           tc0; // ft,8

    // === Prestress ===
    // !!! should not be here, HellmichMaterial used as linear elastic steel prestressed cable, no relaxation...
    // prestress stress in x direction (should add DOF components array) [Pa]
    double prestress;
    // time of start and end of prestress application (may be identical)
    // stored in seconds (intrinsic time)
    double prestressFrom, prestressTo;

    // === Staggered analysis ===
    // Bit mask specifying flattening of gp coordinates for obtaining temperature from temperature field
    // 1: x=0, 2: y=0, 4: z=0
    int flatTemperature;
    /// Base value of the temperatures obtained from temperatureField (273.15 for C->K)
    double temperatureFieldBase;
    /// Initial temperature of the material at cast time (for analysis without thermal field)
    double initialTemperature;

    /// Associated time function for temperature history input
    int tTimeFunction;
    /// Associated time function for humidity history input
    int hTimeFunction;


    /**
     * Isothermal analysis status. Contains reference to material. Initialized in material->initializeFrom...createMaterialGp
     * for isothermal analysis only. It's number is set to 0 to let status services know it's the material-level gp.
     */
    GaussPoint *materialGp;
    /**
     * Flag specifying whether the material-level status needs update.
     * Similar to noniso status auxStateUpdateFlag - set 1->0 in updateYourself, 0->1 in initAuxState.
     * Not absolutely necessary, materialGp->status->nonisoData->auxStateFlag would suffice,
     * but this avoids checking the material-level status flag for each gp.
     */
    int materialGpUpdateFlag;
    /// Stamp specifying the solution state for which the material-level state is initialized. Necessary to avoid multiple computation of InitGpForNewStep after step restart.
    StateCounterType materialGpInitAt;
    // saveContext
    StateCounterType materialGpSaveAt;
    // restoreContext
    StateCounterType materialGpRestoreAt;
    // printOutputAt
    StateCounterType materialGpOutputAt;

    /// Base linear elastic material (aging - with E set method)
    AgingIsoLEMaterial *linMat;
    /// Number of material for heat and moisture transport to use sorption isotherm functions
    int hemoMaterial;

    //!!! === Auxiliary variables ===
    //??? nonsense for parallel computation (
    // temporary coefficients for projection algorithm - simplifying f1 and dfdx functions
    double tempa, tempb;
    // auxiliary values updated for each gp, so that they needn't be passed to f1 and dfdx
    // in the findroot newton iteration
    double auxchi1, auxksi, auxKv, auxGv; // kv - for volumetric stress-strain, Gv - for deviatoric stress-strain

    // === Auxiliary functions ===
    double f1(double dlambda1);
    double dfdx(double dlambda1);
    double computedl(double trial, double EKv);

    double approxnewtonfindroot();
    double newtonfindroot();

    /// aging elasticity E(ksi)
    double agingE(double ksi);
    double agingK(double ksi) { return agingE(ksi) / ( 3. * ( 1 - 2 * ny ) ); } //3.3
    double agingG(double ksi) { return agingE(ksi) / ( 2. * ( 1 + ny ) ); }

    /**
     * Tensorial elastic stress-strain operators
     * @param coeff determines whether creep viscosity coefficients should be taken into account (true)
     */
    void elasticStiffness(FloatArray &stress, FloatArray &strain, GaussPoint *gp, TimeStep *atTime, MatResponseForm form, int coeff);
    void elasticCompliance(FloatArray &strain, FloatArray &stress, GaussPoint *gp, TimeStep *atTime, MatResponseForm form, int coeff);

    /// autogenous shrinkage strains
    double autoShrinkageCoeff(double ksi);
    /// drying shrinkage strains
    double dryingShrinkageCoeff(double h);
    /**
     * !!! Prestress value
     * For use of HellMat as prestressing cable material
     */
    double prestressValue(double time);

    /// short-term creep characteristic time
    double twTime(double ksi);
    /// microprestress force Gamma0
    double computeGamma0(double t, double T);
    /// compute viscous slip (gamma) increment
    double computeViscousSlipIncrement(GaussPoint *gp, TimeStep *atTime);
    /// flow creep viscosity 1/ny,f(gamma, temperature)
    double computeViscosity(GaussPoint *gp, TimeStep *atTime);

    /// strength evolution due to hydration
    double fcStrength(double ksi);
    /// uniaxial compressive threshold - plastic hardening fc...fy
    double rcThreshold(double chi1, double ksi);
    /// derivative d Rc/d chi1
    double dRcdchi1(double chi1, double ksi);

    //maybe should use StructuralMaterial/StressVector services
    double invariantI1(FloatArray &stress);
    void   deviator(FloatArray &dev, FloatArray &stress);
    double deviatorNorm(FloatArray &stress);

    /// returns Delta lambda1,2_n+1 for chi1_n, ksi_n+1
    // uses auxiliary variables stored in material!!!
    void projection(ActiveSurface &active, double &dlambda1, double &dlambda2, FloatArray &stress);
    void stressReturn(FloatArray &stress, FloatArray &trialStress, GaussPoint *gp, TimeStep *atTime);

    virtual void give1dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);

    /**
     * Returns time function for temperature given as direct input.
     */
    LoadTimeFunction *giveTTimeFunction();
    /**
     * Returns time function for relative humidity as direct input
     */
    LoadTimeFunction *givehTimeFunction();


public:
    // === initialization ===
    /// Constructor
    HellmichMaterial(int n, Domain *d);

    /**
     * Initializes auxiliary time-independent material constants according to material parameters.
     */
    void initializeParameters();

    /**
     * Creates isothermal analysis (material-level) gp that contains the auxiliary status values constant in material.
     * Uses virtual Element instance to have reference to material.
     * Sets gp number to 0 to let status services know it's the material-level gp.
     * Material mode set to _1dMat - no stress/strain vectors are used, better _No (_Unknown causes error)
     */
    void createMaterialGp();
    /**
     * Returns isothermal analysis (material-level) gp that contains the auxiliary status values constant in material.
     * Creates the material gp if necessary.
     */
    GaussPoint *giveMaterialGp();

    virtual IRResultType initializeFrom(InputRecord *ir);

    /// destructor
    virtual ~HellmichMaterial();

    /// Returns the material options
    MaterialOptions giveOptions() { return options; }
    int givePssIndex() { return pssIndex; }
    /// Returns the heat&moisture transport material
    HeMoTKMaterial *giveHeMoMaterial();

    // saves current context(state) into stream
    contextIOResultType saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);
    contextIOResultType restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);

    // stress return visualization output
    void   plotReturn(FILE *outputStream, GaussPoint *gp, TimeStep *atTime);
    // stress-strain graph output
    void   plotStressStrain(FILE *outputStream, GaussPoint *gp, TimeStep *atTime, int idx, int id, double err);
    // I1-|s| stress path output
    void   plotStressPath(FILE *outputStream, GaussPoint *gp, TimeStep *atTime, int id, bool trial);

    /**
     * Returns the hydration degree at end of TimeStep atTime in given integraion point.
     * Overriden from HydrationModelInterface to enable Isothermal analysis switch.
     * The value is obtained from gp/material hydration status via the hydration model.
     * @param gp integration point
     * @param atTime solution step
     * @param mode value mode VM_Incremental or VM_Total
     * @return hydration degree or increment in given gp
     */
    double giveHydrationDegree(GaussPoint *gp, TimeStep *atTime, ValueModeType mode);

    // returns current temperature in given gp
    double giveTemperature(GaussPoint *gp);
    // returns the total/incremental temperature change (from IC/previousStep) in given gp
    double giveTemperatureChange(GaussPoint *gp, ValueModeType mode);
    // returns the total/incremental/previous=VM_Velocity humidity in given gp
    double giveHumidity(GaussPoint *gp, ValueModeType mode);
    // returns the current viscous slip value in given gp
    double giveViscousSlip(GaussPoint *gp);
    // returns the base microprestress value for given gp
    double giveGamma0(GaussPoint *gp);
    /*obsolete? move to updateInternalState)
     * void setTemperature(GaussPoint* gp, double T);
     */

    void setMixture(MixtureType mix);

    // returns current prestress force in gp
    double givePrestress(GaussPoint *gp);
    // returns current flow creep viscosity in gp
    double giveViscosity(GaussPoint *gp);

    //modified_zh 13.6.2004 for vol/dev split
    /**
     * Computes coefficients of viscoelasticity and viscous flow influence on moduli K and G
     * Original model assumes total stress causing creep, with moDeviatoricCreep, different
     * coefficients must be used for volumetric and deviatoric part of stress-strain relation.
     * Correct computation is now possible only in 3D mode (Poisson ratio not fixed -> zero
     * stress can't be imposed directly)
     * @param mode VM_Total:
     * @param kv  returns coefficient for total / volumetric stress
     * @param gv  returns coefficient for deviatoric stress
     *
     * @param mode VM_Incremental:
     * @param kv  returns coefficient for viscous creep (ev)
     * @param gv  returns coefficient for flow creep    (fv)
     */                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       //ev        //fv
    void giveKvCoeffs(GaussPoint *gp, TimeStep *atTime, double &kv, double &gv, ValueModeType mode);
    /// Original service for compatibility with non-3D ananlysis, uses gv from giveKvCoeffs.
    double giveKvCoeff(GaussPoint *gp, TimeStep *atTime);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *atTime);

    // zero without hydration; with hydration, zero initial value is used
    void giveShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                   GaussPoint *gp, TimeStep *atTime, ValueModeType mode);
    // strain from prescribed prestress (called by giveEigenStrainVector)
    void givePrestressStrainVector(FloatArray &answer, MatResponseForm form,
                                   GaussPoint *gp, TimeStep *atTime, ValueModeType mode);
    // sums viscous creep strains, flow creep strains and thermal dilation strains and prestress- if any
    // uses temporary creep status ???
    void giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                               GaussPoint *gp, TimeStep *atTime, ValueModeType mode);

    virtual void giveRealStressVector(FloatArray &, MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

    virtual void printOutputAt(FILE *file, TimeStep *atTime);

    // returns the time increment dt of step atTime
    double giveTimeIncrement(TimeStep *atTime);
    // returns the time of step atTime
    double giveTime(TimeStep *atTime);

    /**
     * Updates auxiliary status values (temperature, hydration degree, viscosity) for given time
     * computed at step start, time saved in auxStateAt
     */
    void initAuxStatus(GaussPoint *gp, TimeStep *atTime);

    /**
     * Initializes the status in given integration point.
     * The non-isothermal status is NOT initialized in each iteration, it is only initialized at step start in initAuxStatus.
     */
    virtual void initTempStatus(GaussPoint *gp);

    virtual void initGpForNewStep(GaussPoint *gp);

    virtual void updateYourself(GaussPoint *gp, TimeStep *atTime);

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);

    /**
     * Computes reduced strain vector in given integration point, generated by internal processes in
     * material, which are independent on loading in particular integration point.
     * Default implementation takes only into account temperature induced strains.
     * @param answer returned strain vector
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     * @param determines response mode (Total or incremental)
     */
    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *atTime, ValueModeType mode);

    // gives the base linear elastic material with current hydration state in given gp
    // for reduced stress modes without plasticity
    LinearElasticMaterial *giveLinearElasticMaterial(GaussPoint *gp, TimeStep *atTime);

    // === identification and auxiliary functions ===

    virtual const char *giveClassName() const { return "HellmichMaterial"; }
    virtual classType giveClassID() const { return HellmichMaterialClass; }

    virtual double give(int aProperty, GaussPoint *gp);

    // === Postprocessing functions ===
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const;
};

#else // #ifdef __TM_MODULE

class HellmichMaterial : public StructuralMaterial
{
public:
    /// Constructor
    HellmichMaterial(int n, Domain *d);
    virtual ~HellmichMaterial() {}

    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm f, GaussPoint *gp,
                              const FloatArray &strain, TimeStep *tstep) { answer.resize(0); }

    virtual const char *giveClassName() const { return "HellmichMaterial"; }
    virtual classType giveClassID() const { return HellmichMaterialClass; }
};
#endif
} // end namespace oofem
#endif // hellmat_h
