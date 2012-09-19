/* TODO
 *
 * could only one HydrationModel instance be used during staggered analysis for both tm and sm analysis?
 * if initialized in one of the slave problems, it can hardly be accessed from the other domain,
 * if initialized in master staggered problem, it could be assigned to both of the slave problems (setHydrationModel(*)), but it's not very clean
 * in the future, when more types of hydration models are implemented, a separate input file for the hydration model could be used.
 *
 * the staggered problem should be extended to enable
 * 1) different time step lengths - done by dtf
 * 2) alternate stepping - several time steps performed in one analysis, while only one in the other
 *    e.g. plastic evolution - short time x tm analysis without load changes
 *
 * should provide a function to compute we (water consumption) based on chemical cement composition and/or concrete mixture
 * w/c and a/c or [c/V .. result of wc,ac,rho] ratio.
 *
 * UpdateInternalState:
 *   maybe hydration model should ensure that the hydration degree is NOT recalculated in isothermal analysis for each gp,
 *  when each gp refers to the material level status. Maybe it will be necessary to save time stamp in status?
 *  but even then, there is no way to check if the state vector isn't changed.
 * ==> keep it at caller (only HellMat is now the case)
 * in case of the  planned several hydration model statuses + interpolation, internal state update must be done for each gp,
 *  and HydrationModel must keep track of the status variables extents, to select the appropriate range for interpolation.
 *  ??? / choose 'typical' integration points at beginning of analysis, compute hydration there
 \ use extreme computed values, which may be at different places in time
 *
 */

#include "hydram.h"
#include "gausspnt.h"
#include "datastream.h"
#include "mathfem.h"
#include "contextioerr.h"

#ifndef __MAKEDEPEND
 #include <cstdio>
#endif

namespace oofem {
// ======= class HydrationModelStatus implementation =======
HydrationModelStatus :: HydrationModelStatus(int n, Domain *d, GaussPoint *g) : MaterialStatus(n, d, g)
{
    hydrationDegree = 0;
    tempHydrationDegree = 0;
}

void
HydrationModelStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    fprintf(file, " ksi %.5f", hydrationDegree);
}

void
HydrationModelStatus :: initTempStatus()
// initializes temp status - resets temp variables to saved equilibrium state
// Normally, this is used to clear temp status variables during equilibrium iteration
// (e.g. plastic multiplier, active yield surface, ...)
// In HellmichMaterial, there is number of auxiliary status values that only need to be evaluated once for particular
// time, and have to recalculated only if the iteration is restarted.

// In the hydration model, this means resetting to hydration degree at start of step.
// This should only be done at step restart or restoring context (InitStepIncrements->initGpForNewStep)

// Because hydration model initTempStatus is never called directly, it should be determined in master status
// initTempStatus(), whether hydration model init is necessary.
{
    tempHydrationDegree = hydrationDegree;
}

void
HydrationModelStatus :: updateYourself(TimeStep *atTime)
// update after new equilibrium state reached - prepare status for saving context and for next time increment
{
    hydrationDegree = tempHydrationDegree;
}

//modified_zh 25.6.2004 obj = NULL is set in .h
contextIOResultType
HydrationModelStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
// saves current context(state) into stream
{
    if ( !stream->write(& hydrationDegree, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

//modified_zh 25.6.2004 obj = NULL is set in .h
contextIOResultType
HydrationModelStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
// restores current context(state) from stream
{
    if ( !stream->read(& hydrationDegree, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    tempHydrationDegree = hydrationDegree; // for oofeg
    return CIO_OK;
}


// ======= class HydrationModel implementation =======
// default constructor - Lafarge mixture
HydrationModel :: HydrationModel() : Material(0, NULL) {
    setMixture(mtLafarge);
    useFindRoot = frMixed;
}

// constructor
HydrationModel :: HydrationModel(MixtureType mix, FindRootMethod usefr) : Material(0, NULL)
{
    setMixture(mix);
    if ( ( usefr ) && ( usefr <= frMixed ) ) {
        useFindRoot = usefr;
    } else {
        _error("constructor: unknown FindRootMethod\n");
    }
}

IRResultType
HydrationModel :: initializeFrom(InputRecord *ir)
/**
 * Initializes the hydration model according to object description stored in input record.
 * Hydration options (constant value, cast time) are read in interface::initialize
 */
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                            // Required by IR_GIVE_FIELD macro
    double value;

    //hydration>0  ->  initial hydration degree
    initialHydrationDegree = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, initialHydrationDegree, IFT_HydrationModel_hydration, "hydration"); // Macro
    if ( initialHydrationDegree >= 0. ) {
        OOFEM_LOG_INFO("HydrationModel: Hydration from %.2f.", initialHydrationDegree);
    } else {
        _error("Hydration degree input incorrect, use 0..1 to set initial material hydration degree.");
    }

    if ( ir->hasField(IFT_HydrationModel_c60mix, "c60mix") ) {
        OOFEM_LOG_INFO("HydrationModel: Model parameters for Skanska C60/75 mixture.");
        setMixture(mtC60);
    }

    timeScale = 1.;
    value = -1.;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_HydrationModel_timeScale, "timescale"); // Macro
    if ( value >= 0. ) {
        timeScale = value;
        OOFEM_LOG_INFO("HydrationModel: Time scale set to %.0f", timeScale);
    }

    // Optional direct input of material parameters
    le = 0;
    value = -1.;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_HydrationModel_hheat, "hheat"); // Macro
    if ( value >= 0 ) {
        le = value;
        OOFEM_LOG_INFO("HydrationModel: Latent heat of hydration set to %.0f", le);
    }

    value = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_HydrationModel_cv, "cv"); // Macro
    if ( value >= 0 ) {
        cv = value;
        OOFEM_LOG_INFO("HydrationModel: Cement content set to %.0f kg/m3", cv);
        we = 0.23 * cv;
    }

    value = -1.;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_HydrationModel_water, "water"); // Macro
    if ( value >= 0 ) {
        we = value;
    }

    if ( cv || ( value >= 0 ) ) {
        OOFEM_LOG_INFO("HydrationModel: Water consumption for hydration set to %.0f kg/m3", we);
    }

    return IRRT_OK;
}

void
HydrationModel :: setMixture(MixtureType mix)
{
    // === Initializes material parameters for given mixture ===
    mixture = mix;
    // Normalized chemical affinity regression function coefficients
    // (Lafarge mixture)
    if ( mix == mtLafarge ) {
        aa = 7.313;
        ba = 10.46;
        ca = 169.3;
        da = 4.370;

        e0 = 0.05; // ksi_0
    }
    // Huber mixture
    else if ( mix == mtHuber ) {
        aa = 15.93;
        ba = 7.33;
        ca = 0.855e8;
        da = 26.7;

        e0 = 0.10; // ksi_0
    }
    // Skanska C60/75 mixture
    else if ( ( mix == mtC60 ) || ( mix == mtC100 ) ) { // used for C100 too
        aa = 8.5;
        ba = 5.0;
        ca = 300;
        da = 10;

        e0 = 0.05; // ksi_0
    } else {
        _error("Unknown mixture type!");
    }

    ear = 4000; // activation term [K]
    if ( !le ) {
        le = 190000;   // latent heat [kJ/m3]
    }

    // rc = 2428; // heat capacity [kJ/m3K] - set in master material input
}


// === Material functions ===
double
HydrationModel :: affinity(double ksi)
// Returns the normalized chemical affinity A~(ksi) [1/s]
{
    if ( ksi < e0 ) {
        ksi = e0;
    }

    return ( aa * ( 1 - exp(-ba * ksi) ) / ( 1 + ca * pow(ksi, da) ) );
}

double
HydrationModel :: dAdksi(double ksi)
// Returns the derivation of chemical affinity dA~/dksi(ksi)
{
    double ksinad, enaksi;

    if ( ksi < e0 ) {
        return ( 0 );
    }

    ksinad = pow(ksi, da);
    enaksi = exp(-ba * ksi);

    return ( aa * ( ca * da * ksinad * ( enaksi - 1 ) / ksi + ba * enaksi * ( ca * ksinad + 1 ) ) /
            pow(1 + ca * ksinad, 2) );
}

double
HydrationModel :: localResidual(double dks) // G(dksi) 4.12
//!!! uses auxiliary ksi, dt, T in hydration model instance
{
    return ( dks - auxdt * affinity(auxksi + dks) * exp(-ear / auxT) * ( 1. + auxh * auxh ) / 2. );
}

double
HydrationModel :: dksidT(double ksi, double T, double h, double dt)
/*
 * !!! should use state vector to enable adding state variables, but is OK on the level of dIntSource/dState
 * G = dksi - dt * aff(ksi) * exp(-ear/T) * (1+h^2)/2 = 0
 * dksi/dT = - dG/dT / dG/dksi, dksi->0                        dksi/dksi=1
 *       = - -dt*aff(ksi)*(--ear/T^2)*exp(-ear/T)*(1+h^2)/2  /  (1 - dt*daff/dksi * exp(-ear/T) * (1+h^2)/2)
 *
 * dh/dT, dh/dksi uz nejde, od toho je iterace ksi(h,T)
 */
{
    double aux;
    aux = dt * exp(-ear / T) * ( 1 + h * h ) / 2;
    return ( aux * affinity(ksi) * ear / ( T * T ) ) /
           ( 1 - aux * dAdksi(ksi) );
}

double
HydrationModel :: dksidh(double ksi, double T, double h, double dt)
/*
 * G = dks - dt * aff(ksi) * exp(-ear/T) * (1+h^2)/2 = 0
 * dksi/dh = - dG/dh / dG/dksi
 *       = - -dt*aff(ksi) *exp(-ear/T) * h  /  (1 - dt*daff/dksi * exp(-ear/T) * (1+h^2)/2)
 */
{
    double aux;
    aux = dt * exp(-ear / T);
    return ( aux * affinity(ksi) * h ) /
           ( 1 - aux * dAdksi(ksi) * ( 1 + h * h ) / 2 );
}

MaterialStatus *
HydrationModel :: giveStatus(GaussPoint *gp) const
/**
 * Returns the hydration model status obtained from gp associated material status or from model associated status in case of isothermal analysis
 * Creates the hydration model status if necessary.
 */
{
    HydrationModelStatusInterface *hmi = ( HydrationModelStatusInterface * ) this->giveStatus(gp)->giveInterface(HydrationModelStatusInterfaceType);
    HydrationModelStatus *status = NULL;
    if ( hmi ) {
        status = hmi->giveHydrationModelStatus();
        if ( status == NULL ) {
            status = ( HydrationModelStatus * ) this->CreateStatus(gp);
            hmi->setHydrationModelStatus(status);
        }
    } else {
        _error("giveStatus: Master status undefined.");
    }

    return status;
}

void
HydrationModel :: computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
/** Returns the hydration heat generated in gp during time step atTime.
 *  Can be overriden to return also water consumption as second component of the internal source vector.
 */
//!!! Works only for current TimeStep, no check done
//!!! Should suffice to return val(1) without moisture analysis.
{
    val.resize(2);
    val(0) = le * giveHydrationDegree(gp, atTime, mode); // heat SOURCE
    val(1) = -we *giveHydrationDegree(gp, atTime, mode); // water CONSUMPTION
}

double
HydrationModel :: giveCharacteristicValue(const FloatArray &vec, MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime)
/// returns coefficients for LHS contribution from internal sources (dHeat/dT, dWaterSource/dw)
// Transport status needs to be obtained from master status, it's better to pass as a parameter
// to enable usage from structural material
{
    double answer = 0;

    if ( ( rmode == IntSource ) || ( rmode == IntSource_hh ) || ( rmode == IntSource_ww ) || ( rmode == IntSource_wh ) || ( rmode == IntSource_hw ) ) {
        if ( vec.isEmpty() ) {
            _error("giveCharacteristicValue: undefined state vector.");
        }

        answer = computeIntSource(vec, gp, atTime, rmode);
    } else {
        _error("giveCharacteristicValue: wrong MatResponseMode.");
    }

    return answer;
}

double
HydrationModel :: computeHydrationDegreeIncrement(double ksi, double T, double h, double dt)
/// Computes the hydration degree increment according to given hydration degree, temperature and time increment.
//!!! sets aux values in material
{
    double result = 0.0;

    if ( ksi < 1.0 ) {
        auxksi = ksi;
        auxT = T;
        auxh = h;
        auxdt = dt;
        switch ( useFindRoot ) {
        case frRegula:  result = regulafindroot();
            break;
        case frBinTree: result = bintreefindroot();
            break;
        case frMixed:   result = mixedfindroot();
            break;
#ifdef DEBUG
        default: {
            OOFEM_ERROR2("HydrationModel :: dksi - unknown FindRootMethod %d", useFindRoot);
        }
#endif
        }

        if ( ksi + result > 1.0 ) {
#ifdef VERBOSE
            OOFEM_LOG_INFO("temp=%.12f, dksi %.15f -> %f\n", T, result, 1.0 - ksi);
#endif
            result = 1.0 - ksi;
        }
    } else {
        result = 0.;
    }

    return ( result );
}

double
HydrationModel :: computeIntSource(const FloatArray &vec, GaussPoint *gp, TimeStep *atTime, MatResponseMode rmode)
/**
 * Computes and returns the derivatives of the material-generated Internal Source with respect to the tm state vector.
 * Used in IntSourceLHS matrix for quadratic convergency of global iteration.
 * State vector must contain relative humidity, not water content
 *
 * Called from giveCharacteristicValue - IntSource_hh, IntSource_ww, IntSource_wh, IntSource_hw.
 * IntSource_hh = -lksi * dksidT (dHeat / dT)
 * IntSource_ww = wksi * dksidh (dWater / dh)
 * IntSource_hw = -lksi * dksidh (dHeat / dh)
 * IntSource_wh = -wksi * dksidT (dWater / dT)
 */
{
    // prepare ksi, T, h, dt
    double ksi, T, h, dt;
    ksi = giveHydrationDegree(gp, atTime, VM_Total);
    T = vec(0);
    h = vec(1);
    ///!!! timeScale
    dt = atTime->giveTimeIncrement() * timeScale;

    if ( ksi < 1.0 ) {
        switch ( rmode ) {
        case IntSource:
        case IntSource_hh: return -le *dksidT(ksi, T, h, dt);

        case IntSource_ww: return we * dksidh(ksi, T, h, dt);

        case IntSource_hw: return -le *dksidh(ksi, T, h, dt);

        case IntSource_wh: return -we *dksidT(ksi, T, h, dt);

        default: _error("computeIntSource: Wrong MatResponseMode.");
        }
    }

    return 0;
}

// === public services ===

double
HydrationModel :: giveHydrationDegree(GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
// returns the hydration degree in integration point gp
{
    HydrationModelStatus *status = ( HydrationModelStatus * ) giveStatus(gp);
    double ksi = status->giveTempHydrationDegree();
    if ( mode == VM_Incremental ) {
        ksi -= status->giveHydrationDegree();
    }

    return ksi;
}

void
HydrationModel :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *atTime)
/** Recalculates the hydration degree according to the given new state vector and time increment, using equilib status
 *  State vector is supposed to contain [1]->temperature, [2]->relative humidity!
 *  caller should ensure that this is called only when state vector is changed
 */
{
    double ksi, dksi, T = 0., h = 1., dt;
    // get hydration model status associated with integration point
    HydrationModelStatus *status = ( HydrationModelStatus * ) giveStatus(gp);

    if ( vec.giveSize() ) {
        T = vec(0);
        if ( vec.giveSize() > 1 ) {
            h = vec(1);
        } else {
            h = 1; // assume saturated if undefined
        }
    } else {
        _error("updateInternalState: undefined state vector.");
    }

    ksi = status->giveHydrationDegree();
    if ( !ksi && initialHydrationDegree ) {
        ksi = initialHydrationDegree;
        status->setHydrationDegree(ksi);
    }

    //!!! timeScale
    if ( ( dt = atTime->giveTimeIncrement() ) > 0. ) {
        dksi = computeHydrationDegreeIncrement(ksi, T, h, atTime->giveTimeIncrement() * timeScale);
    } else {
        dksi = 0.;
    }

    status->setTempHydrationDegree(ksi + dksi);
}


// === Auxiliary functions === (root finding)
double
HydrationModel :: regulafindroot()
{
    double x0, y0, yl, xl = 0., xr = 1.;

    do {
        yl = localResidual(xl);
        x0 = yl * ( xl - xr ) / ( localResidual(xr) - yl ) + xl;
        y0 = localResidual(x0);
        if ( y0 < 0 ) {
            xl = x0;
        } else {
            xr = x0;
        }

#ifdef VERBOSEFINDROOT
        OOFEM_LOG_INFO("regulafindroot: x=%.15f, chyba %.15f \n", x0, y0);
#endif
    } while ( fabs(y0) > ROOT_PRECISION_DKSI );

    return ( x0 );
}

double
HydrationModel :: bintreefindroot()
{
    double xl = 0., xr = 1., x0, y0;

    do {
        x0 = ( xl + xr ) / 2;
        y0 = localResidual(x0);
        if ( y0 < 0 ) {
            xl = x0;
        } else {
            xr = x0;
        }

#ifdef VERBOSEFINDROOT
        OOFEM_LOG_INFO("bintreefindroot: x=%.15f, chyba %.15f \n", x0, y0);
#endif
    } while ( fabs(y0) > ROOT_PRECISION_DKSI );

    return ( x0 );
}



double
HydrationModel :: mixedfindroot()
{
    double x0 = 0., y0, yl, xl = 0., xr = 1.;
    int jcount, done = 0;

    do {
        for ( jcount = 0; ( jcount < BINARY_TREE_STEPS ); jcount++ ) {
            x0 = ( xl + xr ) / 2;
            y0 = localResidual(x0);
            if ( fabs(y0) < ROOT_PRECISION_DKSI ) {
                done = 1;
                break;
            }

            if ( y0 < 0 ) {
                xl = x0;
            } else {
                xr = x0;
            }
        }

        if ( !done ) {
            yl = localResidual(xl);
            x0 = yl * ( xl - xr ) / ( localResidual(xr) - yl ) + xl;
            y0 = localResidual(x0);

            if ( fabs(y0) < ROOT_PRECISION_DKSI ) {
                break;
            } else
            if ( y0 < 0 ) {
                xl = x0;
            } else {
                xr = x0;
            }

#ifdef VERBOSEFINDROOT
            OOFEM_LOG_INFO("mixedfindroot: x=%.15f, chyba %.15f \n", x0, y0);
#endif
        }
    } while ( !done );

    return ( x0 );
}

MaterialStatus *
HydrationModel :: CreateStatus(GaussPoint *gp) const
{
    return new HydrationModelStatus(1, this->giveDomain(), gp);
}

// ======= HydrationModelStatusInterface implementation =======
/// Updates the equilibrium variables to temporary values.
void
HydrationModelStatusInterface :: updateYourself(TimeStep *atTime)
{
    if ( hydrationModelStatus ) {
        hydrationModelStatus->updateYourself(atTime);
    }
}
/// Outputs the status variables
void
HydrationModelStatusInterface :: printOutputAt(FILE *file, TimeStep *atTime)
{
    if ( hydrationModelStatus ) {
        hydrationModelStatus->printOutputAt(file, atTime);
    }
}

// ======= HydrationModelInterface implementation =======

IRResultType
HydrationModelInterface :: initializeFrom(InputRecord *ir)
/**
 * Creates the hydration model according to object description stored in input record
 * and invokes hydration model initialization from the input record
 */
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro
    double value;

    //!!! should use separate field, e.g. hydramname #hydramnumber
    // Hydration>0  ->  Model starting at value, hydration<0 -> Constant at given value
    value = -2.;
    constantHydrationDegree = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_HydrationModelInterface_hydration, "hydration"); // Macro
    if ( value >= 0. ) {
        OOFEM_LOG_INFO("HydratingMaterial: creating HydrationModel.");
        if ( !( hydrationModel = new HydrationModel() ) ) {
            ( ( Material * ) this )->_error("Could not create HydrationModel instance.");
        }

        hydrationModel->initializeFrom(ir);
    }
    // constant hydration degree
    else if ( value >= -1. ) {
        constantHydrationDegree = -value;
        OOFEM_LOG_INFO("HydratingMaterial: Hydration degree set to %.2f.", -value);
    } else {
        OOFEM_ERROR("HydrationModelInterface :: initializeFrom - Hydration degree input incorrect, use -1..<0 for constant hydration degree, 0..1 to set initial material hydration degree.");
    }

    // Material cast time - start of hydration
    // 11/3/2004 OK *unfinished in Hellmat, needs to be checked in hm_Interface->updateInternalState
    castAt = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, castAt, IFT_HydrationModelInterface_castAt, "castat"); // Macro
    if ( castAt >= 0. ) {
        OOFEM_LOG_INFO("HydratingMaterial: Hydration starts at time %.2g.", castAt);
    }

    return IRRT_OK;
}

void
HydrationModelInterface :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *atTime)
/**
 * Calls hydrationModel->updateInternalState, if the material is already cast
 * In case the cast time lies within the span of current timestep, the timestep increment is set to (time-castAt)
 */
{
    if ( hydrationModel ) {
        TimeStep *hydraTime = new TimeStep( ( const TimeStep ) *atTime );
        int notime = 0;
        if ( atTime->giveTargetTime() - atTime->giveTimeIncrement() < castAt ) {
            if ( atTime->giveTargetTime() >= castAt ) {
                hydraTime->setTimeIncrement(atTime->giveTargetTime() - castAt);
            } else {
                notime = 1;
            }
        }

        if ( !notime ) {
            hydrationModel->updateInternalState(vec, gp, hydraTime);
        }

        delete hydraTime;
    }
}

double
HydrationModelInterface :: giveHydrationDegree(GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
/// Calls the hydrationModel to return the hydration degree in integration point gp or returns the constant hydration degree
{
    if ( hydrationModel ) {
        return hydrationModel->giveHydrationDegree(gp, atTime, mode);
    } else {
        return constantHydrationDegree;
    }
}
} // end namespace oofem
