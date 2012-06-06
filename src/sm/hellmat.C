// (ViscoPlastic material model (short-term creep, long-term creep, microcracking - multisurface plasticity (Drucker-Prager with hardening + tension cut-off)
// direct displacement control, incremental formulation
// isothermal chemical hardening capability in external class HydrationModel (constant temperature as input material parameter)
//     could be used in thermo- analysis, reading only temperature from thermo- model, evaluating hydration degree increment for given time increment
//  ?? standalone thermochemical analysis - temperature history as input / both ksi & T input ... thermochemo a few steps ahead - temperature interpolation for given time / coordinates
// file hellmat.C

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #include <string.h>
#endif

#include "gausspnt.h"
#include "timestep.h"
#include "loadtime.h"
#include "isolinearelasticmaterial.h"
#include "structuralcrosssection.h"
#include "hellmat.h"
#include "mathfem.h"
#include "oofemtxtinputrecord.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
#ifdef __TM_MODULE

// --------- class AgingIsoLEMaterial implementation---------
AgingIsoLEMaterial :: AgingIsoLEMaterial(int n, Domain *d, double E, double nu) : IsotropicLinearElasticMaterial(n, d, E, nu) { }
void AgingIsoLEMaterial :: setE(double newE) {
    E = newE;
    G = 0.5 * E / ( 1. + nu );
}

// --------- class HellmichMaterial implementation---------
HellmichMaterial :: HellmichMaterial(int n, Domain *d) : StructuralMaterial(n, d), HydrationModelInterface()
// Constructor - does nothing, initialization is done in initializeFrom(ir);
{ }

void
HellmichMaterial :: initializeParameters()
// Initializes the auxiliary material constants for plasticity according to material options and parameters
{
    alpha = sqrt(2. / 3.) * ( kappa - 1 ) / ( 2 * kappa - 1 ); // 3.20

    tc0 = delta * fc;
    kDP0 = alpha * omega * fc;                       // 3.21

    chi1u = -( 9 * alpha * sqrt(2.) * ( ae * epscu - fc ) ) / // 3.25
            ( ae * ( 3 * alpha * sqrt(2.) - 2 * sqrt(3.) ) ); // 3.26

    auxkdp = sqrt(2. / 3.) * kappa / ( 2. * kappa - 1. );
    if ( options & moSqrtHardeningLaw ) {
        auxd = auxkdp * ( 1. - omega ) * 3. * alpha / chi1u;
    } else {
        auxd = -2. * auxkdp * ( 1. - omega ) * 3. * alpha / pow(chi1u, 2.);
    }
}


void
HellmichMaterial :: createMaterialGp()
/*
 * Creates isothermal analysis (material-level) gp that contains the auxiliary status values constant in material.
 * Uses virtual Element instance to have reference to material.
 * Sets gp number to 0 to let status services know it's the material-level gp.
 */
{
    Element *elem = new Element( 0, giveDomain() );
    IntegrationRule *ir = new IntegrationRule(0, elem);
    char eirstr [ 50 ];
    sprintf( eirstr, "mat %d crosssect 1 nodes 0", giveNumber() );
    OOFEMTXTInputRecord eir(eirstr);

    // initialize elem to have reference to this material
    elem->initializeFrom(& eir);
    // number 0, no coords, zero weight, unknown material mode
    materialGp = new GaussPoint(ir, 0, NULL, 0, _3dMat);
    if ( !materialGp ) {
        _error("Could not create the material-level gp.");
    }

    ( ( HellmichMaterialStatus * ) giveStatus(materialGp) )->setInitialTemperature(initialTemperature); // set material temperature in K
 #ifdef VERBOSE_HELLMAT
    // output to check if it worked
    printf("\n Hellmat::createMaterialGp:");
    printf("\n  Auxiliary element input record: %s", eirstr);
    printf( "\n  Created auxiliary element with material %s.", elem->giveMaterial()->giveClassName() );
    printf( "\n  Created material-level gp with material %s.", materialGp->giveMaterial()->giveClassName() );
    printf( "\n  Material temperature set to %f.", giveTemperature(NULL) );
 #endif
}

GaussPoint *
HellmichMaterial :: giveMaterialGp() {
    if ( !materialGp ) {
        createMaterialGp();
    }

    return materialGp;
}

IRResultType
HellmichMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                            // Required by IR_GIVE_FIELD macro
    double value;
    int intvalue;

    // parent classes initialization
    this->StructuralMaterial :: initializeFrom(ir); // reads material density d
    this->HydrationModelInterface :: initializeFrom(ir); // initializes the hydration model if needed

    // === init material options ===
    // default options set in hellmat.h
    options = moDefault;

    // === material constants (Hellmich original) ===
    ny = 0.2;    // Hellmich 3.2.4.2 (pg 38): Lafarge shotcrete
    kappa = 1.16; //   f_b-biaxial c/ f_c-uniaxial
    omega = 0.25; //   f_y-yield c / f_c
    delta = 0.1; //   f_t-tensile / f_c
    epscu = 0.0022; //   strain at f_c,8
    fc    = 40.1e6; //   f_c,8 [Pa] (39.6 at 90 days)  // Huber ~ 30e6
    ae    = 41.3e9; //   E,8   [Pa]   // Huber 40.3e9

    ashr  = 4.96e-4; //5e-4;   // shrinkage strains - Huber mixture!
    bshr  = -10.91e-4; //-11e-4; //  eps(ksi) = a,shr + b,shr*ksi, eps < 0;

    kshr  = 8e-4;   //8e-4;   // reversible humidity shrinkage parameter Kappa

    ksi0 = 0.05;    // percolation threshold

    // Creep
    modulusH = 1e6 / 7.; // [Pa]; Laplante 1/5 MPa, B3 model 1/7 MPa;
    ur = 2700;     // [K] U/R (activation energy / universal gas const)
    refT = 20 + 273; // [K] reference temperature for long-term creep
    jv = 24e-12;   // [Pa^(-1)] viscous creep compliance
    tw = 28. * 3600; // [s] short-term creep characteristic time - 28 hours
    c = 1e-17;     // [K-1?] drying shrinkage parameter 1e-16..1e-18?


    // Optionally change basic material properties
    IR_GIVE_OPTIONAL_FIELD(ir, ae, IFT_HellmichMaterial_E, "e"); // Macro
    printf("\nHellMat: Ultimate Young modulus E=%.4g.", ae);

    if ( ir->hasField(IFT_HellmichMaterial_linearE, "lineare") ) {
        options = options | moLinearEModulus;
        printf("\nHellMat: Forcing linear relation E(ksi).");
    }

    IR_GIVE_OPTIONAL_FIELD(ir, ny, IFT_HellmichMaterial_nu, "n"); // Macro
    printf("\nHellMat: Poisson const nu=%.3f.", ny);
    IR_GIVE_OPTIONAL_FIELD(ir, epscu, IFT_HellmichMaterial_epscu, "epscu"); // Macro // 0.0022
    IR_GIVE_OPTIONAL_FIELD(ir, fc, IFT_HellmichMaterial_fc, "fc");    // Macro // 40.1e6  88e6
    IR_GIVE_FIELD(ir, value, IFT_HellmichMaterial_tAlpha, "talpha"); // Macro
    propertyDictionary->add(tAlpha, value);

    initialTemperature = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, initialTemperature, IFT_HellmichMaterial_isoT, "isot"); // Macro
    if ( initialTemperature >= 0 ) {
        printf("\nHellMat: Isothermal analysis at %.2f K.", initialTemperature);
        options = options | moIsothermal;
        tTimeFunction = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, tTimeFunction, IFT_HellmichMaterial_Tltf, "tltf"); // Macro
        if ( tTimeFunction ) {
            printf("\nHellMat: Time function number %d used as temperature history.", tTimeFunction);
        }

        hTimeFunction = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, tTimeFunction, IFT_HellmichMaterial_hltf, "hltf"); // Macro
        if ( tTimeFunction ) {
            printf("\nHellMat: Time function number %d used as moisture history.", tTimeFunction);
        }
    } else {
        options = options & ~moIsothermal;
        IR_GIVE_OPTIONAL_FIELD(ir, initialTemperature, IFT_HellmichMaterial_iniT, "init"); // Macro
    }

    materialGpUpdateFlag = 0; // needs init
    materialGpInitAt = materialGpSaveAt = materialGpRestoreAt = 0; // before first step
    // Options for non-isothermal analysis - base of the temperature field, flattening of the temperature field, hemo material for sorption isotherm functions
    flatTemperature = 0;
    if ( !( options & moIsothermal ) ) {
        value = -1;
        // for temperature field in Celsius degrees, hydration analysis requires absolute temperatures
        IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_HellmichMaterial_baseT, "baset"); // Macro
        if ( value >= 0 ) {
            temperatureFieldBase = value;
        } else {
            temperatureFieldBase = 0; //273.15;
        }

        if ( temperatureFieldBase ) {
            printf("\nHellMat: Temperature field moved by %.2f.", temperatureFieldBase);
        }


        flatTemperature = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, flatTemperature, IFT_HellmichMaterial_flatT, "flatt"); // Macro
        if ( flatTemperature ) {
            if ( flatTemperature == -1 ) {
                printf("\nHellMat: Temperature field flattened xz->xy");
            } else {
                printf("\nHellMat: Temperature field flattened in axis %d.", ( flatTemperature / 2 ) + 1);
            }
        }

        hemoMaterial = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, hemoMaterial, IFT_HellmichMaterial_hemomat, "hemomat"); // Macro
        if ( hemoMaterial ) {
            printf("\nHellMat: Using material %d for sorption isotherms functions.", hemoMaterial);
        }
    }

    // Hydration switch; options are read in HydrationModelInterface
    // !!! not very consistent with interface concept?
    // mix setting not available in input record
    value = -1.;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_HellmichMaterial_hydration, "hydration"); // Macro
    if ( value >= 0. ) {
        options = options | moHydration;
    } else {
        options = options & ~moHydration;
    }

    if ( options & moHydration ) {
        setMixture(mtLafarge);
    }

    // Stress-strain plot output (plastic.out)
    pssIndex = -1;
    pssElement = 0;
    pssGaussPoint = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, pssIndex, IFT_HellmichMaterial_plotss, "plotss"); // Macro
    if ( pssIndex > -1 ) {
        options = options | moPlotStressStrain;
        printf("\nHellMat: Stress-strain plot for component %d.", pssIndex);
        FILE *sstream = fopen(OUTPUTFILE_SS, "w");
        if ( sstream ) {
            plotStressStrain(sstream, NULL, NULL, pssIndex, -1, 0);
            fclose(sstream);
        }
    } else {
        options = options & ~moPlotStressStrain;
    }

    if ( options & moPlotStressStrain ) {
        if ( ir->hasField(IFT_HellmichMaterial_pssiter, "pssiter") ) {
            options = options | moPlotStressStrainIter;
            printf("\nHellMat: Stress-Strain plot in each iteration.");
        } else {
            options = options & ~moPlotStressStrainIter;
        }

        IR_GIVE_OPTIONAL_FIELD(ir, pssElement, IFT_HellmichMaterial_psselem,  "psselem"); // Macro
        if ( pssElement ) {
            printf("\nHellMat: Stress-Strain plot for Element %d.", pssElement);
        }

        IR_GIVE_OPTIONAL_FIELD(ir, pssGaussPoint, IFT_HellmichMaterial_pssgp, "pssgp"); // Macro
        if ( pssGaussPoint ) {
            printf("\nHellMat: Stress-Strain plot for GaussPoint %d.", pssGaussPoint);
        }
    }

    // !!! Prestress (as initial stress)
    prestress = 0;
    prestressFrom = -1;
    prestressTo = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, prestress, IFT_HellmichMaterial_prestress, "prestress"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, prestressFrom, IFT_HellmichMaterial_prestressFrom, "prestressfrom"); // Macro
    prestressTo = prestressFrom;
    IR_GIVE_OPTIONAL_FIELD(ir, prestressTo, IFT_HellmichMaterial_prestressTo, "prestressto"); // Macro
    if ( prestress ) {
        printf("\n Hellmat: prestress %.5e applied from %.2f to %.2f", prestress, prestressFrom, prestressTo);
    }

    // Input time scale
    timeScale = 1.;
    value = -1.;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_HellmichMaterial_timeScale, "timescale"); // Macro
    if ( value > 0. ) {
        timeScale = value;
        printf("\nHellMat: Time scale set to %.0f", timeScale);
    }

    prestressFrom *= timeScale;
    prestressTo *= timeScale;

    // Select chemical shrinkage and moisture-dependent volume changes
    intvalue = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, intvalue, IFT_HellmichMaterial_shr, "shr"); // Macro
    if ( ( intvalue == 0 ) || ir->hasField(IFT_HellmichMaterial_noshr, "noshr") ) { // noshr set
        options = options & ~moShrinkage;
        options = options & ~moHumidityStrain;
        options = options & ~moDryingShrinkage;
        printf("\nHellMat: No autogenous shrinkage, no humidity shrinkage.");
    } else if ( intvalue > 0 ) {
        if ( intvalue & 1 ) {
            options = options | moShrinkage;
            printf("\nHellMat: Autogenous shrinkage enabled.");
        } else {
            options = options & ~moShrinkage;
            printf("\nHellMat: No autogenous shrinkage.");
        }

        if ( intvalue & 2 ) {
            options = options | moHumidityStrain;
            printf("\nHellMat: Humidity shrinkage enabled.");
        } else {
            options = options & ~moHumidityStrain;
            printf("\nHellMat: No humidity shrinkage.");
        }

        if ( intvalue & 4 ) {
            options = options | moDryingShrinkage;
            printf("\nHellMat: Drying shrinkage enabled.");
        } else {
            options = options & ~moDryingShrinkage;
            printf("\nHellMat: No drying shrinkage.");
        }
    }

    // if autogenous shrinkage is used, check ashr options
    if ( options & moShrinkage ) {
        if ( !( options & moHydration ) ) {
            printf("\nHellMat: No chemical shrinkage without hydration.");
        } else { // optionally change autogenous shrinkage parameters (eshr = ashr + bshr * ksi; eshr <= 0 )
            IR_GIVE_OPTIONAL_FIELD(ir, ashr, IFT_HellmichMaterial_ashr, "ashr"); // Macro //   5e-4;
            // Set bshr 0 to disable empiric hydration-dependent autogenous shrinkage
            IR_GIVE_OPTIONAL_FIELD(ir, bshr, IFT_HellmichMaterial_bshr,  "bshr"); // Macro // -11e-4;
        }
    }

    // if humidity shrinkage is enabled, check hshr parameter Kappa
    if ( options & moHumidityStrain ) {
        IR_GIVE_OPTIONAL_FIELD(ir, kshr, IFT_HellmichMaterial_kshr, "kshr"); // Macro // 8e-4;
        if ( ir->hasField(IFT_HellmichMaterial_kshr, "kshr") ) {
            printf("\nHellMat: Humidity volume change parameter Kappa set to %.3e", kshr);
        }
    }

    // if drying shrinkage is enabled, check parameter C
    if ( options & moDryingShrinkage ) {
        IR_GIVE_OPTIONAL_FIELD(ir, c, IFT_HellmichMaterial_dryingc, "dc"); // Macro // 1e-17;
        if ( ir->hasField(IFT_HellmichMaterial_dryingc, "dc") ) {
            printf("\nHellMat: Drying shrinkage coefficient set to %.3e", c);
        }
    }

    // Switch off creep
    if ( ir->hasField(IFT_HellmichMaterial_nocreep, "nocreep") ) {
        options = options & ~moCreep;
        printf("\nHellMat: No creep.");
    }

    // Use parameters for the Skanska C60/75 mixture
    if ( ir->hasField(IFT_HellmichMaterial_c60mix, "c60mix") ) {
        printf("\nHellMat: Model parameters for Skanska C60/75 mixture.");
        setMixture(mtC60);
        modulusH = 1e6 / 9.5;
        jv = 30e-12;
        fc = 88e6;
        epscu = 0.003;
    }

    // optionally change material creep parameters
    if ( options & moCreep ) {
        IR_GIVE_OPTIONAL_FIELD(ir, modulusH, IFT_HellmichMaterial_modulusH, "modulush"); // Macro // 1e6/7  1e6/9,5
        IR_GIVE_OPTIONAL_FIELD(ir, ur, IFT_HellmichMaterial_ur, "ur"); // Macro // 2700
        IR_GIVE_OPTIONAL_FIELD(ir, jv, IFT_HellmichMaterial_jv, "jv"); // Macro // 24e-12 30e-12
        IR_GIVE_OPTIONAL_FIELD(ir, tw, IFT_HellmichMaterial_tw, "tw"); // Macro // 28 days

        // deviatoric creep options
        intvalue = -1;
        IR_GIVE_OPTIONAL_FIELD(ir, intvalue, IFT_HellmichMaterial_devc, "devc"); // Macro
        if ( intvalue >= 0 ) { // devc set
            options = options & ~( moDeviatoricCreepE | moDeviatoricCreepF ); // init. turn off dev. creep
            intvalue = intvalue << 11; // 0/1/2/3 -> 0/2048/4096/6144
            options = options | intvalue;
            printf("\nHellmat: deviatoric creep");
            if ( !intvalue ) {
                printf(" disabled.");
            } else {
                printf(" enabled for");
                if ( options & moDeviatoricCreepE ) {
                    printf(" viscous creep");
                }

                if ( options & moDeviatoricCreepF ) {
                    printf(" flow creep");
                }
            }
        }
    }

    // Initialize auxiliary material constants according to material parameters
    initializeParameters();

    // === Plasticity ===
    // Switch off plasticity
    if ( ir->hasField(IFT_HellmichMaterial_noplast, "noplast") ) {
        options = options & ~moPlasticity;
        printf("\nHellMat: No plasticity.");
    }

    if ( options & moPlasticity ) {
        // Switch off volumetric dependence of the yield surface (~ J2 plasticity)
        if ( ir->hasField(IFT_HellmichMaterial_zeroalpha, "zeroalpha") ) {
            alpha = 0;
            options = options | moHardening;
            printf("\nHellMat: J2 plasticity, no hardening.");
        } else if ( ir->hasField(IFT_HellmichMaterial_nohardening, "nohardening") ) {
            options = options & ~moHardening;
            printf("\nHellMat: no plastic hardening.");
        }

        // Use numeric derivative in local plastic iteration
        if ( ( options & moPlasticity ) && ( ir->hasField(IFT_HellmichMaterial_approxnewton, "approxnewton") ) ) {
            options = options | moApproxNewton;
            printf("\nHellMat: Approximate derivation used in local plastic iteration.");
        }

        // Compute plastic multiplier directly, not using the local Newton iteration
        // !!! not sure whether works correctly in 3D corner
        if ( ( options & moPlasticity ) && ( ir->hasField(IFT_HellmichMaterial_computedl, "computedl") ) ) {
            options = options | moComputedl;
            printf("\nHellMat: Direct computation of plastic multiplier dlambda1.");
        } else {
            options = options & ~moComputedl;
            printf("\nHellMat: Newton solution of plastic multiplier dlambda1.");
        }
    }

    printf("\nHellMat: Material options %d", options);
    linMat = new AgingIsoLEMaterial(this->giveDomain()->giveNumberOfMaterialModels() + 1, this->giveDomain(), ae, ny);

    printf("\n");
    return IRRT_OK;
}

HellmichMaterial :: ~HellmichMaterial()
// Destructor
{
    if ( linMat ) {
        delete linMat;
    }
}

/* From OOFEM documentation, but I have not used it
 * Interface*
 * HellmichMaterial::giveInterface (InterfaceType type)
 * /// Returns this as HydrationModelInterface
 * {
 * if (type == HydrationModelInterfaceType) return (HydrationModelInterface*) this;
 * else return NULL;
 * }
 */

double
HellmichMaterial :: agingE(double ksi)
// Returns the Young modulus of material at given hydration degree ksi
// Hellmich A.14
{
    if ( ( options & moLinearEModulus ) || ( mixture == mtHuber ) ) {
        return ( ae * ksi );                                             // Huber mixture, input forced
    } else {
        return ( ae * sqrt(ksi) ); // Lafarge mixture
    }
}

void
HellmichMaterial :: elasticStiffness(FloatArray &stress, FloatArray &strain, GaussPoint *gp, TimeStep *atTime, MatResponseForm form, int coeff)
/*
 * Returns the isotropic linear elastic stress vector for the given strain vector.
 * Uses the attached aging linear elastic material for material modes other than 3dMat or 1dMat.
 * If coeff is true, uses kv and gv creep coefficients according to moDeviatoricCreep settings
 */
{
    int i;
    double ksi, Gv, Kv, evKG;
    MaterialMode mmode = gp->giveMaterialMode();
    ksi = giveHydrationDegree(gp, atTime, VM_Total);

    if ( mmode == _1dMat ) {
        if ( form == ReducedForm ) {
            stress.resize(1);
        } else {
            stress.resize(6);
            stress.zero();
        }

        stress(0) = agingE(ksi) * strain(0);
    } else if ( mmode == _3dMat ) {
        double kv, gv;
        if ( coeff ) {
            giveKvCoeffs(gp, atTime, kv, gv, VM_Total);
        } else {
            kv = gv = 1.;
        }

        stress.resize(6);
        Kv = agingK(ksi) / kv;
        Gv = agingG(ksi) / gv;
        // sigv * [K-2/3*G = ny*E/((1+ny)*(1-2*ny))];
        evKG = ( strain(0) + strain(1) + strain(2) ) * ( Kv - 2. / 3. * Gv );

        for ( i = 0; i < 3; i++ ) {
            stress(i) = evKG + 2 *Gv *strain(i);
        }

        for ( i = 3; i < 6; i++ ) {
            stress(i) = Gv * strain(i);
        }
    } else if ( strain.containsOnlyZeroes() ) {
        stress = strain;
    } else {                                                                // linear elastic material
        FloatMatrix d;
        LinearElasticMaterial *lMat;
        lMat = giveLinearElasticMaterial(gp, atTime);
        if ( lMat->hasMaterialModeCapability(mmode) ) {
            lMat->giveCharacteristicMatrix(d, form, ElasticStiffness, gp, atTime);
            stress.beProductOf(d, strain);
        } else {
            _error("elasticStiffness: unsupported material mode.");
        }
    }
}

void HellmichMaterial :: elasticCompliance(FloatArray &strain, FloatArray &stress, GaussPoint *gp, TimeStep *atTime, MatResponseForm form, int coeff)
// Returns the isotropic linear elastic strain vector for the given stress vector
{
    int i;
    double ksi, E, Gv, Kv, svE;
    MaterialMode mmode = gp->giveMaterialMode();
    ksi = giveHydrationDegree(gp, atTime, VM_Total);
    // adjust hydration degree to have a non-zero stiffness
    if ( ksi < HYDRATION_MINDEGREE ) {
        ksi = HYDRATION_MINDEGREE;
    }

    E = agingE(ksi);
    if ( !E ) {
        _error("Zero stiffness in elastic compliance!");
    }

    if ( mmode == _1dMat ) {
        if ( form == ReducedForm ) {
            strain.resize(1);
        } else {
            strain.resize(6);
            strain.zero();
        }

        strain(0) = stress(0) / E;
    } else if ( mmode == _3dMat ) {
        double kv, gv;
        if ( coeff ) {
            giveKvCoeffs(gp, atTime, kv, gv, VM_Total);
        } else {
            kv = gv = 1.;
        }

        Gv = agingG(ksi) / gv;
        Kv = agingK(ksi) / kv;

        strain.resize(6);
        svE = ( stress(0) + stress(1) + stress(2) ) / 3. * ( 1 / ( 3 * Kv ) - 1 / ( 2 * Gv ) );
        for ( i = 0; i < 3; i++ ) {
            strain(i) = svE + stress(i) / ( 2 * Gv );
        }

        for ( i = 3; i < 6; i++ ) {
            strain(i) = stress(i) / Gv;
        }
    } else if ( stress.containsOnlyZeroes() ) {
        strain = stress;
    } else {
        FloatMatrix d;
        LinearElasticMaterial *lMat;
        lMat = giveLinearElasticMaterial(gp, atTime);
        if ( lMat->hasMaterialModeCapability(mmode) ) {
            lMat->giveCharacteristicComplianceMatrix(d, form, ElasticStiffness, gp, atTime);
            strain.beProductOf(d, stress);
        } else {
            _error("elasticCompliance: unsupported material mode.");
        }
    }
}


double HellmichMaterial :: prestressValue(double time)
{
    if ( prestress ) {
        if ( time < prestressFrom ) {
            return 0;
        } else if ( ( prestressTo > prestressFrom ) && ( time < prestressTo ) ) {
            return prestress * ( time - prestressFrom ) / ( prestressTo - prestressFrom );
        } else {
            return prestress;
        }
    } else {
        return 0.;
    }
}


double HellmichMaterial :: autoShrinkageCoeff(double ksi)
// A.15
{
    if ( bshr && ( ksi > -ashr / bshr ) ) {
        return ( ashr + bshr * ksi );
    } else {
        return ( 0. );
    }
}


double HellmichMaterial :: dryingShrinkageCoeff(double h)
// added 17/02/2004
// shrinkage as function of relative humidity
// eshrh = K (1 - h^3), K = -0,0008 [Desky... Beton TKS 2/2002]
{
    if ( h < 1. ) {
        return kshr * ( 1 - h * h * h );
    } else {
        return 0.;
    }
}


// --------- creep ---------
inline double HellmichMaterial :: twTime(double ksi)
// Linear dependence of the short-term creep characteristic time on the hydration degree
// 3.48
{
    return ( tw * ksi );
}


double HellmichMaterial :: computeGamma0(double t, double T)
/*
 * Micropresstress force at time and temperature of material setting point (ksi0).
 * Integrated analytically for isothermal conditions
 */
// 2.80, p=2
{
    return ( 1. / ( modulusH * t * exp( -ur * ( 1 / T - 1 / refT ) ) ) );
}


double HellmichMaterial :: computeViscousSlipIncrement(GaussPoint *gp, TimeStep *atTime)
/*
 * Computes and returns the viscous slip increment for given time increment, temperature, humidity and previous prestress and slip values.
 * Positive root of quadratic equation for viscous slip increment.
 */
// 4.36
{
    double expterm, ba, ca;
    double base, gamman, dt, dT, T, dh, h;
    dt = giveTimeIncrement(atTime);
    HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) giveStatus(gp);
    // might also get directly from status - should be called only for materialGp from initAuxStatus in isothermal case
    T = giveTemperature(gp);
    dT = giveTemperatureChange(gp, VM_Incremental);
    h = giveHumidity(gp, VM_Total);
    dh = giveHumidity(gp, VM_Incremental);
    base = status->giveGamma0();
    gamman = status->giveViscousSlip();
    if ( base == 0. ) {
        return ( 0. );
    }

    expterm = dt * exp( -2. * ur * ( 1. / T - 1. / refT ) );
    ca = pow(modulusH, 2.) * expterm;
    ba = ( 1. + 2. * base * modulusH * expterm ) / ca;
    // gamman = 0 ~ 1e-14, fabs = 0 ~ 10, dg = 1e-18~1e-16 ... c = 1e-17 (guess 1e-16..1e-18)
    if ( options & moDryingShrinkage ) {
        ca = ( gamman + base * base * expterm - c * fabs(dT * log(h)  +  T * dh / h) ) / ca;
    } else {
        ca = ( gamman + base * base * expterm ) / ca;
    }

    return ( ( ba - sqrt(ba * ba - 4. * ca) ) / 2. - gamman );
}

double HellmichMaterial :: computeViscosity(GaussPoint *gp, TimeStep *atTime)
// returns the flow creep viscosity 1/ny,f for given viscous slip -> microprestress
// viscosity is saved in nonisothermal status, could be removed and computed from gamma directly when needed
{
    double base, gamma, T;
    T = giveTemperature(gp);
    HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) giveStatus(gp);
    base = status->giveGamma0();
    gamma = status->giveTempViscousSlip();
    // c'=1e-6 - unit constant for Pa
    // c' * Microprestress * exp(...)
    return ( 1e-6 * ( base - modulusH * gamma ) * exp( -2 * ur * ( 1 / T - 1 / refT ) ) );
}

// --------- chemical hardening ---------
double HellmichMaterial :: fcStrength(double ksi)
//
// 3.24
{
    if ( ksi <= ksi0 ) {
        return ( 0 );
    } else {
        return ( fc * ( ksi - ksi0 ) / ( 1 - ksi0 ) );
    }
}

double HellmichMaterial :: rcThreshold(double chi1, double ksi)
//
// 3.23
{
    if ( chi1 < chi1u ) {
        if ( options & moSqrtHardeningLaw ) {
            if ( chi1 >= 0 ) {
                return ( fcStrength(ksi) *
                        ( omega + ( 1 - omega ) * sqrt( chi1 * ( 2 * chi1u - chi1 ) ) / chi1u ) );
            } else { // symmetric negative branch needed for stability of return mapping iteration
                return ( fcStrength(ksi) *
                        ( omega - ( 1 - omega ) * sqrt( chi1 * ( -2 * chi1u - chi1 ) ) / chi1u ) );
            }
        } else { // quadratic hardening law - original
            return ( fcStrength(ksi) *
                    ( omega + ( 1 - omega ) * ( 1 - pow( ( chi1 - chi1u ) / chi1u, 2 ) ) ) );
        }
    } else {
        return ( fcStrength(ksi) );
    }
}


double HellmichMaterial :: dRcdchi1(double chi1, double ksi)
// full dRcdchi1 = fcStrength(ksi)*(1-omega) * (chi1u-chi1)/(chi1u*sqrt(chi1*(2*chi1u-chi1)))
{
    if ( chi1 < chi1u ) {
        if ( options & moSqrtHardeningLaw ) {
            if ( chi1 >= 0 ) {
                if ( chi1 < NEWTON_DERIVATIVEDX * chi1u ) {
                    chi1 = NEWTON_DERIVATIVEDX * chi1u;
                }

                return ( fcStrength(ksi) * ( chi1u - chi1 ) / sqrt( chi1 * ( 2 * chi1u - chi1 ) ) );
            } else {
                if ( fabs(chi1) < NEWTON_DERIVATIVEDX * chi1u ) {
                    chi1 = -NEWTON_DERIVATIVEDX * chi1u;
                }

                return ( fcStrength(ksi) * ( chi1u + chi1 ) / sqrt( chi1 * ( -2 * chi1u - chi1 ) ) );
            }
        } else { // quadratic hardening law - original
            return ( fcStrength(ksi) * ( chi1 - chi1u ) );
        }
    } else {
        return ( 0 );
    }
}

// --------- plasticity ---------
double HellmichMaterial :: invariantI1(FloatArray &stress)
{
    return ( stress(0) + stress(1) + stress(2) );
}

void HellmichMaterial :: deviator(FloatArray &dev, FloatArray &stress)
{
    dev.resize(6);

    int i;
    double sv;

    sv = invariantI1(stress) / 3;
    for ( i = 0; i < 3; i++ ) {
        dev(i) = stress(i) - sv;
    }

    for ( i = 3; i < 6; i++ ) {
        dev(i) = stress(i);
    }
}

double HellmichMaterial :: deviatorNorm(FloatArray &stress)
{
    int i;
    double sv, result;
    FloatArray deviator(6);

    sv = ( stress(0) + stress(1) + stress(2) ) / 3;
    for ( i = 0; i < 3; i++ ) {
        deviator(i) = stress(i) - sv;
    }

    for ( i = 3; i < 6; i++ ) {
        deviator(i) = stress(i);
    }

    result = 0;
    for ( i = 3; i < 6; i++ ) {
        result += pow(deviator(i), 2);
    }

    result *= 2.;
    for ( i = 0; i < 3; i++ ) {
        result += pow(deviator(i), 2);
    }

    return ( sqrt(result) );
}

double HellmichMaterial :: f1(double dlambda1)
// DP yield value for projection algorithm - uses auxiliary values of ksi, chi
{
    return ( tempa + tempb * dlambda1 - auxkdp * rcThreshold(auxchi1 + ( ( options & moHardening ) > 0 ) * dlambda1 * 3 * alpha, auxksi) );
}

double HellmichMaterial :: dfdx(double dlambda1)
{
    if ( ( options & moHardening ) && ( auxchi1 < chi1u ) ) {
        return ( tempb - auxd * dRcdchi1(auxchi1 + 3 * alpha * dlambda1, auxksi) );
    } else {
        return ( tempb );
    }
}

double HellmichMaterial :: computedl(double trial, double EKv)
// solves the quadratic equation for plastic multiplier dl1:
// f1 = trial - EKv*dl1 - [fc*omega + fc*(1-omega)/chi1u * sqrt((chi1+3*alfa*dl)*(2*chi1u - (chi1+3*alfa*dl)))];
// trial == tempa/auxkdp, EKv == -tempb/auxkdp;
{
    double fca = fcStrength(auxksi);
    if ( auxchi1 < chi1u ) {
        double fcol, fcom, poma, pomb, pomc;
        fcol = fca * ( 1 - omega ) / chi1u;
        fcol *= fcol;
        fcom = fca * omega - trial;
        poma = 9. * alpha * alpha + EKv * EKv / fcol;
        pomb = 6. * alpha * ( auxchi1 - chi1u ) + 2. * EKv * fcom / fcol;
        pomc = auxchi1 * ( auxchi1 - 2. * chi1u ) + fcom * fcom / fcol;

        return ( -pomb - sqrt(pomb * pomb - 4. * poma * pomc) ) / ( 2. * poma );
    } else {
        return -( trial + fca ) / EKv;
    }
}

void HellmichMaterial :: projection(ActiveSurface &active, double &dlambda1, double &dlambda2, FloatArray &trialStress)
// projects given trialStress on the given active yield surface in 3D stress space
// uses auxksi, auxchi1, auxKv set in giveRealStressVector... independent on gp,status,timestep
// returns activeSurface and computed flowIncrements
{
    double i1;

 #ifdef DEBUG
    if ( active == asNone ) {
        dlambda1 = dlambda2 = 0;
        printf("HellmichMaterial::projection: No active yield surface!\n");
        return;
    }

 #endif


    i1 = invariantI1(trialStress);

    if ( active == asTC ) {
        dlambda1 = 0;
        dlambda2 = auxKv * ( i1 - delta * fcStrength(auxksi) )
                   / ( 9 * agingK(auxksi) );
        return;
    } else if ( active == asDP ) {
        dlambda2 = 0;

        tempa = alpha * i1 + deviatorNorm(trialStress);
        tempb = -( 9 * alpha * alpha * agingK(auxksi) / auxKv + 2 * agingG(auxksi) / auxGv );
    } else if ( active == asCorner ) {
        tempa = deviatorNorm(trialStress) + alpha *delta *fcStrength(auxksi);
        tempb = -2 * agingG(auxksi) / auxGv;
    } else {
        _error("HellmichMaterial::projection: Unknown ActiveSurface value.");
    }


    if ( options & moComputedl ) {
        dlambda1 = computedl(tempa / auxkdp, -tempb / auxkdp);
    } else
    if ( options & moApproxNewton ) {
        dlambda1 = approxnewtonfindroot();
    } else {
        dlambda1 = newtonfindroot();
    }

    if ( active == asCorner ) {
        if ( dlambda1 < 0 ) {
 #ifdef VERBOSE
            //   printf("HellmichMaterial::projection: corner region - only TC surface active.\n");
 #endif
            active = asTC;
            projection(active, dlambda1, dlambda2, trialStress);
        } else {
            dlambda2 = auxKv * ( i1 - delta * fcStrength(auxksi) ) / ( 9 * agingK(auxksi) )
                       - dlambda1 * alpha;
            if ( dlambda2 < 0 ) {
 #ifdef VERBOSE
                //    printf("HellmichMaterial::projection: corner region - only DP surface active.\n");
 #endif
                active = asDP;
                projection(active, dlambda1, dlambda2, trialStress);
            }
        }
    }

 #ifdef DEBUG
    // check
    if ( ( dlambda1 < 0 ) || ( ( dlambda2 < 0 ) ) ) {
        printf("HellmichMaterial::projection: negative yield increment!\n \
  dl1=%f, dl2=%f, activeSurface=%d", dlambda1, dlambda2, active);
        _error("");
    }

 #endif
}

void HellmichMaterial :: stressReturn(FloatArray &stress, FloatArray &trialStress, GaussPoint *gp, TimeStep *atTime)
// uses auxksi, auxchi1, auxKv - is set in giveRealStressVector / in main
// sets (now material!) gp->status->activeSurface, plasticStrain, flowIncrement, activeSurface
{
    double inv, snorm, dl1, dl2;
    FloatArray auxepsp(6), depsp(6);
    HellmichMaterialStatus *status;
    ActiveSurface as;
    int i;

    status = ( HellmichMaterialStatus * ) giveStatus(gp);

    inv = invariantI1(trialStress);
    snorm = deviatorNorm(trialStress);
    deviator(depsp, trialStress);

    as = asNone;
    if ( alpha * inv + snorm - auxkdp * rcThreshold(auxchi1, auxksi) > 0 ) {
        as = ( ActiveSurface ) ( as | asDP );
    }

    if ( inv - delta * fcStrength(auxksi) > 0 ) {
        as = ( ActiveSurface ) ( as | asTC );
    }

    // elastoplastic
    if ( as != asNone ) {
        projection(as, dl1, dl2, trialStress);

        // evaluate plastic strain increment
        if ( as & asDP ) {
            depsp(0) = dl1 * ( depsp(0) / snorm + alpha );
            depsp(1) = dl1 * ( depsp(1) / snorm + alpha );
            depsp(2) = dl1 * ( depsp(2) / snorm + alpha );
            snorm /= 2 * dl1; // eng. notation 2*dl1, tensor 1*dl1
            // depsp_ij = dl1 * sigtr_ij/snorm
            depsp(3) /= snorm;
            depsp(4) /= snorm;
            depsp(5) /= snorm;
        } else {
            depsp.zero();
        }

        if ( as & asTC ) {
            depsp(0) += dl2;
            depsp(1) += dl2;
            depsp(2) += dl2;
        }

        // original stress = trial - 1/Kv * [Ce]::depsp
        // devcreep stress = trial - [Cev]::depsp
        elasticStiffness(stress, depsp, gp, atTime, FullForm, true);
        for ( i = 0; i < 6; i++ ) {
            stress(i) = trialStress(i) - stress(i);
        }

        // plastic status
        if ( options & moHardening ) {
            status->setTempHardeningVar(auxchi1 + 3 * alpha * dl1);
        }

        status->givePlasticStrainVector(auxepsp);
        auxepsp.add(depsp);

        // temporary plastic status update
        status->setTempPlasticStrainVector(auxepsp);
        status->setActiveSurface(as);
        status->setFlowIncrement(dl1);
    } else {
        // elastic
        stress = trialStress;
        status->setActiveSurface(asNone);
        status->setFlowIncrement(0.);
    }
}

void HellmichMaterial :: give1dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                       MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime)
//
// consistent tangents for 1D
// Young modulus E + Kv creep coeff + plasticity
//
{
    ActiveSurface as;
    double Eep = 0, Ev, drc, ksi;
    HellmichMaterialStatus *status;
 #ifdef DEBUG
    if ( gp->giveMaterialMode() != _1dMat ) {
        _error("give1dStiffMtrx: wrong stress-strain mode!");
    }

 #endif
    ksi = giveHydrationDegree(gp, atTime, VM_Total);
    if ( ksi < HYDRATION_MINDEGREE ) {
        ksi = HYDRATION_MINDEGREE;
    }

    Ev = agingE(ksi) / giveKvCoeff(gp, atTime);
    status = ( HellmichMaterialStatus * ) giveStatus(gp);

    as = status->giveActiveSurface();
    if ( ( as == asNone ) || ( rMode == ElasticStiffness ) ) { // elastic
        Eep = Ev;
    } else if ( as == asTC ) { // tension
        Eep = 0;
    } else if ( as == asDP ) { // compression
        drc = ( 1 - omega ) / chi1u *dRcdchi1(status->giveTempHardeningVar(), ksi);
        Eep = Ev * ( 1  -  1 / ( 1 + 3 * ( kappa - 1 ) * drc / ( Ev * kappa ) ) );
    } else {
        _error("Unknown active surface.");
    }

    if ( form == ReducedForm ) {
        answer.resize(1, 1);
    } else { // FullForm
        answer.resize(6, 6);
        answer.zero();
    }

    answer(0, 0) = Eep;
}

void HellmichMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                       MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime)
//
// consistent tangents for 3D
// isotropic linear elastic matrix + Kv creep coeff + plasticity
//
// with plasticity works correctly only for 3D material mode
//
// for deviatoric creep, G and K has different Kv coefficient.
{
    ActiveSurface as;
    double dl1, ksi, chi;
    double kv, gv, Gv, GGv, Kv, KGv, GGv2, curent;
    int i, j;
    HellmichMaterialStatus *status;
    FloatArray stress;

 #ifdef DEBUG
    if ( gp->giveMaterialMode() != _3dMat ) {
        _error("give3dStiffMtrx: wrong stress-strain mode!");
    }

 #endif
    if ( !( ( rMode == TangentStiffness ) || ( rMode == ElasticStiffness ) || ( rMode == SecantStiffness ) ) ) {
        _error("HellmichMaterial:give3dStiffness.: Unsupported material response mode.");
    }

    answer.resize(6, 6);
    answer.zero();


    ksi = giveHydrationDegree(gp, atTime, VM_Total);
    if ( ksi < HYDRATION_MINDEGREE ) {
        ksi = HYDRATION_MINDEGREE;
    }

    giveKvCoeffs(gp, atTime, kv, gv, VM_Total); // obtain visc coefficients for volumetric(kv) and deviatoric(gv) stress
    status = ( HellmichMaterialStatus * ) giveStatus(gp);
    as = status->giveActiveSurface();

    Gv = agingG(ksi) / gv; // deviatoric visc coeff included in base value of G
    GGv = 2 * Gv;
    Kv = agingK(ksi) / kv; // volumetric visc coeff included in base value of K
    KGv = Kv - 2. / 3. * Gv; // for K*1x1 - 1/3*2G*1x1

    // ======   elastic  ======
    // isotropic linear elastic stiffness, using current E, G, K; if CREEP then reduced by Kv - viscous flow coefficient
    // use also if forced elastic response
    /*
     * G = 1/2 E(ksi) / (1 + ny)
     * K = 1/3 E(ksi) / (1-2*ny)
     * orig. [De] = 1/Kv ( 2G*[Idev] + K*[1x1] )
     * devc. [De] = 2Gv*[Idev] + Kv*[1x1]  = 2Gv*[I] + (Kv-2/3Gv)*[1x1]  = GGv*[I] + KGv*[1x1]
     */
    if ( ( as == asNone ) || ( rMode == ElasticStiffness ) ) {
        for ( i = 0; i < 6; i++ ) {
            for ( j = 0; j < 6; j++ ) {
                if ( i == j ) {
                    curent = GGv;    // I
                } else {
                    curent = 0;
                }

                if ( ( i < 3 ) && ( j < 3 ) ) {
                    curent += KGv;   // 1x1
                }

                answer(i, j) = curent;
            }
        }
    }
    // ======   elastoplastic   =======
    // --- tension cut-off ---
    // [Dep] = 2Gv*[Idev] = GGv * ([I] - 1/3 [1x1])
    else if ( as == asTC ) {
        for ( i = 0; i < 6; i++ ) {
            for ( j = 0; j < 6; j++ ) { // Idev = I  -  1/3 1x1
                if ( i == j ) {
                    curent = 1;        // I
                } else {
                    curent = 0;
                }

                if ( ( i < 3 ) && ( j < 3 ) ) {
                    curent -= 1. / 3.; // 1x1
                }

                answer(i, j) = curent * GGv;
            }
        }
    } else {
        // --- Drucker-Prager active ---
        FloatArray normdev(6), flowdir(6);
        double snorm, df;
        // get reduced stress vector and plastic status
        stress = status->giveTempStressVector();
        dl1 = status->giveFlowIncrement();
        chi = status->giveTempHardeningVar();

        GGv2 = GGv * GGv; // GGv2 = (2G/gv)^2    (orig. ~= 4G^2 / Kv^2)
        df = 9. * Kv * alpha * alpha + GGv;
        if ( options & moHardening ) {
            df += auxd * dRcdchi1(chi, ksi);
        }

        deviator(normdev, stress);
        snorm = deviatorNorm(stress);
        for ( i = 0; i < 6; i++ ) {
            normdev(i) /= snorm;
        }

        status->giveTrialStressVector(stress);
        snorm = deviatorNorm(stress);

        for ( i = 0; i < 6; i++ ) {
            flowdir(i) = GGv * normdev(i) + ( i < 3 ) * 3 * Kv * alpha;
        }

        // only DP active
        // [Dep] = (1 - 2Gv*dl1/snorm)*2Gv*[Idev]  +  Kv*[1x1]  -  1/df [fxf]  +  (2Gv)^2 dl1/snorm * [nxn]
        // [Dep] = (2Gv - (2Gv)^2*dl1/snorm)*I + (Kv-2/3Gv + 1/3 (2Gv)^2*dl1/snorm)*1x1 - ...
        if ( as == asDP ) {
            for ( i = 0; i < 6; i++ ) {
                for ( j = 0; j < 6; j++ ) {
                    if ( i == j ) {
                        curent = GGv - GGv2 * dl1 / snorm;    // I
                    } else {
                        curent = 0;                                  // Idev + 1x1
                    }

                    if ( ( i < 3 ) && ( j < 3 ) ) {
                        curent += KGv  +  GGv2 / 3. * dl1 / snorm; // 1x1
                    }

                    // -1/df [fxf]  +  (2Gv)^2 dl1/snorm * [nxn]
                    curent += -flowdir(i) * flowdir(j) / df  +  GGv2 *normdev(i) * normdev(j) * dl1 / snorm;

                    answer(i, j) = curent;
                }
            }
        } else {
            // corner region
            // [Dep] = (1 - 2Gv*dl1/snorm)*2Gv*[Idev]  -  2Gv/df * [fxn]  +  (2Gv)^2 dl1/snorm * [nxn]  +  3Kv*alpha * [1xn] * 2Gv/df
            df -= 9. * Kv * alpha * alpha; // df = 2Gv [+ auxd*dRcdchi1]
            for ( i = 0; i < 6; i++ ) {
                for ( j = 0; j < 6; j++ ) {
                    if ( i == j ) {
                        curent = GGv - GGv2 * dl1 / snorm;    // I
                    } else {
                        curent = 0;                                  // Idev
                    }

                    if ( ( i < 3 ) && ( j < 3 ) ) {
                        curent += -( GGv - GGv2 * dl1 / snorm ) / 3; // 1x1
                    }

                    // -2Gv/df * [fxn]  +  (2Gv)^2 dl1/snorm * [nxn]  +  3Kv*alpha * [1xn] * 2Gv/df
                    curent += -flowdir(i) * GGv * normdev(j) / df  +  GGv2 *normdev(i) * normdev(j) * dl1 / snorm;
                    if ( i < 3 ) {
                        curent += 3 *Kv *GGv *alpha *normdev(j) / df; // 1xn
                    }

                    answer(i, j) = curent;
                }
            }
        }
    }

    // change to account for engineering strain notation
    for ( i = 0; i < 6; i++ ) {
        for ( j = 3; j < 6; j++ ) {
            answer(i, j) /= 2.;
        }
    }
}

// when solving for dlambda in projection algorithm, approximates df/dx with dx = NEWTON_DERIVATIVEDX
///@todo Check if approx. is faster/stable
double HellmichMaterial :: approxnewtonfindroot()
// derivative extrapolated from 'small' dx - useless in general case
{
    double x0 = 0., y0, dfdx, prec;

    y0 = f1(x0);
    prec = ROOT_PRECISION_DLAMBDA * fabs(y0);
    if ( prec < FINDROOT_SMALLNUM ) {
        prec = FINDROOT_SMALLNUM;
    }

    do {
        dfdx = ( f1(x0 + NEWTON_DERIVATIVEDX) - y0 ) / NEWTON_DERIVATIVEDX;
        x0 -= y0 / dfdx;
        y0 = f1(x0);
 #ifdef VERBOSEFINDROOT
        printf("approxnewtonfindroot: x=%.15g, dfdx= %.15g, chyba %.2g \n", x0, dfdx, y0);
 #endif
    } while ( fabs(y0) > prec );

    return ( x0 );
}

double HellmichMaterial :: newtonfindroot()
{
    double x0 = 0., y0, df, prec;

    y0 = f1(x0);
    prec = ROOT_PRECISION_DLAMBDA * fabs(y0);
    if ( prec < FINDROOT_SMALLNUM ) {
        prec = FINDROOT_SMALLNUM;
    }

    do {
        df = dfdx(x0);
 #ifdef VERBOSEFINDROOT
        printf("newtonfindroot: x=%.15g, dfdx= %.15g, chyba %.2g \n", x0, df, y0);
 #endif
        x0 -= y0 / df;
        y0 = f1(x0);
    } while ( fabs(y0) > prec );

 #ifdef VERBOSEFINDROOT
    printf("          root: x=%.15g, dfdx= %.15g, chyba %.2g \n", x0, df, y0);
 #endif
    return ( x0 );
}


void principalStresses(double *answer, double *s, stressStrainPrincMode mode)
// adapted from oofem1.3 StructuralMaterial::computePrincipalValues
// This function computes Principal values of streses.
// streses are stored in vector form in array s.
// Engineering notation is used.
//
// Problem size: 3D problem, array s contains:
//               {Sxx,Syy,Szz,Syz,Szx,Sxy}
// Return Values:
//
//    array sp -> principal strains or stresses
{
    double swap;
    int nonzeroFlag = 0;
    double I1 = 0.0, I2 = 0.0, I3 = 0.0, s1, s2, s3;
    int i, j;

    for ( i = 0; i < 6; i++ ) {
        if ( fabs(s [ i ]) > 1.e-20 ) {
            nonzeroFlag = 1;
        }
    }

    if ( nonzeroFlag == 0 ) {
        return;
    }

    I1 = s [ 0 ] + s [ 1 ] + s [ 2 ];
    I2 = s [ 0 ] * s [ 1 ] + s [ 1 ] * s [ 2 ] + s [ 2 ] * s [ 0 ] -
         ( s [ 3 ] * s [ 3 ] + s [ 4 ] * s [ 4 ] + s [ 5 ] * s [ 5 ] );
    I3 = s [ 0 ] * s [ 1 ] * s [ 2 ] + 2. * s [ 3 ] * s [ 4 ] * s [ 5 ] -
         ( s [ 0 ] * s [ 3 ] * s [ 3 ] + s [ 1 ] * s [ 4 ] * s [ 4 ] + s [ 2 ] * s [ 5 ] * s [ 5 ] );

    /*
     * Call cubic3r to ensure, that all three real eigenvalues will be found,
     * because we have symmetric tensor.
     * This allows to overcome various rounding errors when solving general
     * cubic equation.
     */
    cubic3r( ( double ) -1., I1, -I2, I3, & s1, & s2, & s3, & i );

    if ( i > 0 ) {
        answer [ 0 ] = s1;
    }

    if ( i > 1 ) {
        answer [ 1 ] = s2;
    }

    if ( i > 2 ) {
        answer [ 2 ] = s3;
    }

 #ifdef DEBUG
    else {
        printf("principalStresses: 3 valid roots not found!\n");
    }
 #endif

    // sort results
    for ( i = 1; i < 3; i++ ) {
        for ( j = 0; j < 2; j++ ) {
            if ( answer [ j + 1 ] > answer [ j ] ) {
                swap = answer [ j + 1 ];
                answer [ j + 1 ] = answer [ j ];
                answer [ j ] = swap;
            }
        }
    }
}


LoadTimeFunction *HellmichMaterial :: giveTTimeFunction()
// Returns the load-time function of the receiver.
{
    if ( !tTimeFunction ) {
        return NULL;
    }

    return this->giveDomain()->giveLoadTimeFunction(tTimeFunction);
}

LoadTimeFunction *HellmichMaterial :: givehTimeFunction()
// Returns the load-time function of the receiver.
{
    if ( !hTimeFunction ) {
        return NULL;
    }

    return this->giveDomain()->giveLoadTimeFunction(hTimeFunction);
}


void HellmichMaterial :: plotReturn(FILE *outputStream, GaussPoint *gp, TimeStep *atTime)
// stress return visualisation output - only correct for 3D analysis
{
    double inv, snorm, ksi, psigma [ 3 ];
    HellmichMaterialStatus *status;
    FloatArray trialStress(6), stress(6);

 #ifdef DEBUG
    if ( outputStream == NULL ) {
        _error("plotReturn: can't write into NULL stream");
    }

 #endif

    status = ( HellmichMaterialStatus * ) giveStatus(gp);
    status->giveTrialStressVector(trialStress);
    stress = status->giveTempStressVector();
    ksi = giveHydrationDegree(gp, atTime, VM_Total);

    // principal stress space visualization - MathCad
    principalStresses(psigma, trialStress.givePointer(), principal_stress);
    fprintf(outputStream, "%.15g %.15g %.15g\n", psigma [ 0 ], psigma [ 1 ], psigma [ 2 ]);
    principalStresses(psigma, stress.givePointer(), principal_stress);
    fprintf(outputStream, "%.15g %.15g %.15g\n", psigma [ 0 ], psigma [ 1 ], psigma [ 2 ]);

    // I1-|s| stress space visualization - Excel
    fprintf( outputStream, "%d\n", status->giveActiveSurface() );
    fprintf(outputStream, "%.15g\n", alpha);
    fprintf( outputStream, "%.15g\n",
            auxkdp * rcThreshold(status->giveHardeningVar(), ksi) );
    fprintf( outputStream, "%.15g\n",
            auxkdp * rcThreshold(status->giveTempHardeningVar(), ksi) );
    fprintf( outputStream, "%.15g\n", delta * fcStrength(ksi) );

    inv = invariantI1(trialStress);
    snorm = deviatorNorm(trialStress);
    fprintf(outputStream, "%.15g %.15g\n", inv, snorm); // trial stress
    inv = invariantI1(stress);
    snorm = deviatorNorm(stress);
    fprintf(outputStream, "%.15g %.15g\n", inv, snorm); // corrected stress
}

void HellmichMaterial :: plotStressStrain(FILE *outputStream, GaussPoint *gp, TimeStep *atTime, int idx, int id, double err)
// stress-strain graph output
{
    FloatArray auxVector;
    HellmichMaterialStatus *status;
    double ksi, eshr, evisc = 0, eflow = 0, epl = 0, stress, strain, chi = 0;

    if ( id == -1 ) { // print header
        fprintf(outputStream, "sigma%d eps%d epsp chi1 ksi elem gp evisc eflow time(h) eshr\n", idx, idx);
        return; // no data yet
    }

    if ( ( pssElement && ( id != pssElement ) ) || ( pssGaussPoint && ( ( int ) err != pssGaussPoint ) ) ) {
        return; // not this gp / element
    }

    ksi = giveHydrationDegree(gp, atTime, VM_Total);
    status = ( HellmichMaterialStatus * ) giveStatus(gp);
    if ( status->giveTempStressVector().giveSize() ) {
        stress = status->giveTempStressVector() (idx);
    } else {
        stress = 0;
    }

    if ( status->giveTempStrainVector().giveSize() ) {
        strain = status->giveTempStrainVector() (idx);
    } else {
        strain = 0;
    }

    if ( ( ( options & moHasShrinkage ) == moHasShrinkage ) || ( options & moHumidityStrain ) ) {
        giveShrinkageStrainVector(auxVector, ReducedForm, gp, atTime, VM_Total);
        if ( auxVector.giveSize() ) {
            eshr = auxVector(idx);
        } else {
            eshr = 0;
        }
    } else {
        eshr = 0;
    }

    if ( options & moPlasticity ) {
        status->giveTempPlasticStrainVector(auxVector);
        if ( auxVector.giveSize() > idx ) {
            epl = auxVector(idx);
        }

        chi = status->giveTempHardeningVar() / chi1u;
    }

    if ( options & moCreep ) {
        status->giveTempViscousStrainVector(auxVector);
        if ( auxVector.giveSize() > idx ) {
            evisc = auxVector(idx);
        }

        status->giveTempFlowStrainVector(auxVector);
        if ( auxVector.giveSize() > idx ) {
            eflow = auxVector(idx);
        }
    }

    fprintf(outputStream, "%.15g %.15g %.15g %.15g %.15f %d %.2g %.15g %.15g %15g %.15g\n",
            stress, strain, epl, chi, ksi, id, err, evisc, eflow, atTime->giveTargetTime(), eshr);
}

void HellmichMaterial :: plotStressPath(FILE *outputStream, GaussPoint *gp, TimeStep *atTime, int id, bool trial)
// I1-|s| stress path output
{
    HellmichMaterialStatus *status;
    FloatArray stress(6);
    double ksi;

    status = ( HellmichMaterialStatus * ) giveStatus(gp);

    if ( id == -1 ) { // print header
        ksi = giveHydrationDegree(gp, atTime, VM_Total);
        fprintf(outputStream, "%.15g\n", alpha);
        fprintf( outputStream, "%.15g\n", auxkdp * rcThreshold(status->giveHardeningVar(), ksi) );
        fprintf( outputStream, "%.15g\n", delta * fcStrength(ksi) );
    }

    if ( trial ) {
        status->giveTrialStressVector(stress);
        fprintf( outputStream, "%.15g %.15g %d\n", invariantI1(stress), deviatorNorm(stress), status->giveActiveSurface() );
    }

    stress = status->giveTempStressVector();
    fprintf(outputStream, "%.15g %.15g %d\n", invariantI1(stress), deviatorNorm(stress), id);
}

HeMoTKMaterial *
HellmichMaterial :: giveHeMoMaterial()
{
    if ( !hemoMaterial || !giveDomain()->giveMaterial(hemoMaterial) ) {
        _error("giveHeMoMaterial: undefined reference to heat and moisture transport material.");
    }

    return ( ( HeMoTKMaterial * ) giveDomain()->giveMaterial(hemoMaterial) );
}

void
HellmichMaterial :: printOutputAt(FILE *file, TimeStep *atTime)
{
    StateCounterType c = atTime->giveSolutionStateCounter();
    if ( materialGpOutputAt != c ) {
        materialGpOutputAt = c;
        giveStatus( giveMaterialGp() )->printOutputAt(file, atTime);
    }
}

double HellmichMaterial :: giveTimeIncrement(TimeStep *atTime)
// returns the time increment of step atTime
// needed to enable time scaling
{
    return atTime->giveTimeIncrement() * timeScale;
}

double HellmichMaterial :: giveTime(TimeStep *atTime)
// returns the time of step atTime
{
    return atTime->giveTargetTime() * timeScale;
}

void
HellmichMaterial :: initAuxStatus(GaussPoint *gp, TimeStep *atTime)
/*
 * Auxiliary status variables need to be calculated once for each gp in one time step. Should be called at step start.
 * Called for both gp and material-level status
 * Should be added after step restart - probably Material->InitGpForNewStep()?
 * Now called at:
 * giveRealStressVector()
 * computeStressIndependentStrainVector()
 * giveCharacteristicMatrix()
 * Checks and sets the nonisoData status update flag (auxStateUpdateFlag)
 */
// 1) temperature (only iso now)
// 2) hydration degree (const or iso from hydration model)
// 3) creep (for given temperature and hydr. degree - ok)
//    sets Gamma0 if ksi>=ksi0;
{
    // === Material-level check ===
    // Init materialGp and set the material-level status update flag if not done yet
    // 0  -> initialize, set to 1
    // 1  -> nothing
    if ( ( options & moIsothermal ) && !materialGpUpdateFlag ) {
 #ifdef VERBOSE_HELLMAT
        printf( "Hellmat: Material-level auxiliary status initialized for time=%.0f, flag 0->1\n", giveTime(atTime) );
 #endif
        materialGpUpdateFlag = 1;
        materialGpInitAt = atTime->giveSolutionStateCounter();
        initAuxStatus(giveMaterialGp(), atTime);
    }

    // === Check if aux status is uninitialized ===  works for both gp and material-level status
    HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) giveStatus(gp);
    if ( status->giveUpdateFlag() ) {
        return;
    } else {
        status->setUpdateFlag(1);
    }

    // === Proceed with aux status update ===
    double ksi, dksi, ksi0tim, ksi0T, dt, tim, w = 0., h = 1., T = 0., dT, g0 = -1., auxg;
    FloatArray gcoords, et;
    dt = giveTimeIncrement(atTime);
    tim = giveTime(atTime);
    // ---- Evaluate current temperature and moisture ----
    if ( options & moIsothermal ) { // === Isothermal ===
        // constant temperature in space
        // temperature history might be prescribed, computed from adiabatic/quasiadiabatic test sim., ....
        // if history ltf is not prescribed, only takes previous value
        // history is expected in original time scale
        if ( giveTTimeFunction() ) {
            T = ( ( LoadTimeFunction * ) giveTTimeFunction() )->evaluate(atTime, VM_Total);
        } else {
            T = giveTemperature(gp);
        }

        if ( givehTimeFunction() ) {
            h = ( ( LoadTimeFunction * ) givehTimeFunction() )->evaluate(atTime, VM_Total);
        } else {
            h = giveHumidity(gp, VM_Total);
        }
    } else {   // === Non-isothermal ===
        // Get temperature from transportProblem via Field Manager
        // Might also get hydration degree, so that it need not be computed twice
        FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
        Field *tf;
        StructuralElement *elem;
        // == Temperature ==
        if ( ( tf = fm->giveField(FT_Temperature) ) ) {
            // temperature field registered
            elem = ( StructuralElement * ) gp->giveElement();
            elem->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
            if ( flatTemperature ) {
                if ( flatTemperature == -1 ) { // temperature xy -> xz (2D cross-section for 3D girder)
                    gcoords.at(2) = gcoords.at(3);
                    gcoords.at(3) = 0.;
                } else {
                    if ( flatTemperature & 1 ) {
                        gcoords.at(1) = 0.;
                    }

                    if ( ( gcoords.giveSize() > 1 ) && ( flatTemperature & 2 ) ) {
                        gcoords.at(2) = 0.;
                    }

                    if ( ( gcoords.giveSize() > 2 ) && ( flatTemperature & 4 ) ) {
                        gcoords.at(3) = 0.;
                    }
                }
            }

            tf->evaluateAt(et, gcoords, VM_Total, atTime);
            T = et.at(1) + temperatureFieldBase;
        } else { // no temperature field
            if ( initialTemperature > 0 ) {
                T = initialTemperature;
            } else {
                _error("InitAuxStatus: Temperature field not registered, initial temperature not set. Use isoT or iniT setting.");
            }

 #ifdef VERBOSE_HELLMAT
            printf("\nInitAuxStatus: Temperature field not registered, using initial teperature %.2f.", T);
 #endif
        }

        // set initial temperature to actual temperature at time of cast
        if ( T && ( tim >= castAt ) && ( status->giveInitialTemperature() <= 0 ) ) {
            status->setInitialTemperature(T);
 #ifdef VERBOSE_HELLMAT
            printf("\nInitAuxStatus: Setting initial temperature of gp %d to %.3f", gp->giveNumber(), T);
 #endif
        } // end set initial temperature

        // == Moisture ==
        if ( tf && ( tf = fm->giveField(FT_HumidityConcentration) ) ) {
            // temperature AND humidity concentration field registered
            // using gcoords defined in temperature field
            tf->evaluateAt(et, gcoords, VM_Total, atTime);
            if ( et.giveSize() > 0 ) {
                w = et.at(1);
            } else {
                _error("\nInitAuxStatus: Moisture value couldn't be obtained.");
            }

            // compute relative humidity from moisture content
            if ( w ) {
                h = giveHeMoMaterial()->inverse_sorption_isotherm(w);
            } else {
                h = 0.;
            }
        } else { // no moisture field
            h = 1.;
 #ifdef VERBOSE_HELLMAT
            printf("\nInitAuxStatus: Moisture field not registered, using saturated conditions.");
 #endif
        } // == end Moisture ==

    } // === end non-isothermal ===

    status->setTempTemperature(T);
    status->setTempHumidity(h);
    // prepare the new state vector to contain current temperature and humidity
    et.resize(2);
    et.at(1) = T;
    et.at(2) = h;
    // compute hydration degree increment
    HydrationModelInterface :: updateInternalState(et, gp, atTime);
    //??? how to evaluate the long-term creep process start? Must be interpolated through the whole step
    // multistepping for hydration degree increment is removed from hydrationModel, because it is not consistent with the tm analysis.
    // determine the time from start of hydration in current step
    // 17/02/2004 ... is removed, might be implemented anew.
    // get hydration degree increment
    ksi = giveHydrationDegree(gp, atTime, VM_Total);
    dksi = giveHydrationDegree(gp, atTime, VM_Incremental);
    if ( options & moCreep ) {
        // read the base prestress value
        g0 = status->giveGamma0();
        if ( ( g0 <= 0 ) && ( ksi >= ksi0 ) && dksi ) { // if base prestress is not set yet and hydration
            // Evaluate the base micro-prestress value
            // compute interpolated time; Temperature? at ksi=ksi0 ...
            //!!! a little bit inconsistent in case tim-castAt<dt
            dT = giveTemperatureChange(gp, VM_Incremental);
            ksi0tim = tim - dt * ( ksi - ksi0 ) / dksi;
            ksi0T = T - dT * ( tim - ksi0tim ) / dt;
            g0 = computeGamma0(ksi0tim, ksi0T);
 #ifdef VERBOSE_HELLMAT
            printf("Starting long-term creep process at ksi=%.4f for time=%f, T=%f.\n", ksi, ksi0tim, ksi0T);
 #endif
            // save the base prestress value
            status->setGamma0(g0);
        } // end base prestress evaluation

        // update creep status stored in gp
        // uses temperature and hydration degree saved in gp above
        if ( g0 > 0. ) {
            auxg = status->giveViscousSlip();
            auxg += computeViscousSlipIncrement(gp, atTime);
            status->setTempViscousSlip(auxg);
            status->setViscosity( computeViscosity(gp, atTime) );
        } // end creep status update

    } // end creep

}
double
HellmichMaterial :: giveHydrationDegree(GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
/*
 * Returns the hydration degree in given integration point. Overrides the giveHydrationDegree from HydrationModelInterface
 * to enable isothermal/noniso switch.
 */
{
    if ( options & moHydration ) {
        if ( options & moIsothermal ) {
            gp = giveMaterialGp();
        }

        return hydrationModel->giveHydrationDegree(gp, atTime, mode);
    } else {
        return ( mode == VM_Total ) ? constantHydrationDegree : 0.;
    }
}

double HellmichMaterial :: giveTemperature(GaussPoint *gp)
/*
 * Returns the temperature in given integration point.
 */
{
    if ( options & moIsothermal ) {
        gp = giveMaterialGp();
    }

    return ( ( HellmichMaterialStatus * ) giveStatus(gp) )->giveTempTemperature();
}

double HellmichMaterial :: giveTemperatureChange(GaussPoint *gp, ValueModeType mode)
{
    if ( options & moIsothermal ) {
        gp = giveMaterialGp();
    }

    HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) giveStatus(gp);
    if ( mode == VM_Total ) {
        return ( status->giveTempTemperature() - status->giveInitialTemperature() );
    } else if ( mode == VM_Incremental ) {
        return ( status->giveTempTemperature() - status->giveTemperature() );
    }

    return 0.;
}

double HellmichMaterial :: giveHumidity(GaussPoint *gp, ValueModeType mode)
{
    if ( options & moIsothermal ) {
        gp = giveMaterialGp();
    }

    HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) giveStatus(gp);
    if ( mode == VM_Total ) {
        return status->giveTempHumidity();
    } else if ( mode == VM_Incremental ) {
        return status->giveTempHumidity() - status->giveHumidity();
    } else if ( mode == VM_Velocity ) { // VM_Previous
        return status->giveHumidity();
    }

    return 1.;
}

//10.6.2004 - unused
double HellmichMaterial :: giveViscousSlip(GaussPoint *gp)
/*
 * Returns the viscous slip in given integration point.
 */
{
    if ( options & moIsothermal ) {
        gp = giveMaterialGp();
    }

    return ( ( HellmichMaterialStatus * ) giveStatus(gp) )->giveViscousSlip();
}

double HellmichMaterial :: giveGamma0(GaussPoint *gp)
/*
 * Returns the base microprestress in given integration point.
 */
{
    if ( options & moIsothermal ) {
        gp = giveMaterialGp();
    }

    return ( ( HellmichMaterialStatus * ) giveStatus(gp) )->giveGamma0();
}

// === Creep ===
double HellmichMaterial :: givePrestress(GaussPoint *gp)
// returns the microprestress force Gamma
{
    if ( options & moCreep ) {
        if ( options & moIsothermal ) {
            gp = giveMaterialGp();
        }

        HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) giveStatus(gp);
        double g0;
        g0 = status->giveGamma0();
        if ( g0 <= 0 ) {
            return ( 0. ); // gamma0 not set yet
        } else {
            return ( g0 - modulusH * status->giveTempViscousSlip() );
        }
    } else {
        return 0.; // no creep
    }
}

double HellmichMaterial :: giveViscosity(GaussPoint *gp)
// returns the flow creep viscosity 1/ny,f for given integration point
// 4.34
{
    if ( options & moIsothermal ) {
        gp = giveMaterialGp();
    }

    HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) giveStatus(gp);
    if ( status->giveGamma0() <= 0. ) {
        return ( 0. );                     // gamma0 not set yet
    } else {
        return ( status->giveViscosity() );
    }
}

double HellmichMaterial :: giveKvCoeff(GaussPoint *gp, TimeStep *atTime)
// coefficient accounting for the influence of visco-elasticity and viscous flow
// on the modulus of elasticity of the stress-strain relation
// 4.39
{
    if ( options & ( moDeviatoricCreepE | moDeviatoricCreepF ) ) {
        _error("giveKvCoeff: can't compute total ceep coefficient, different for volumetric/deviatoric stress.");
    }

    double kv, gv;
    giveKvCoeffs(gp, atTime, kv, gv, VM_Total);
    return kv;
}
//ev        //fv
void HellmichMaterial :: giveKvCoeffs(GaussPoint *gp, TimeStep *atTime, double &kv, double &gv, ValueModeType mode)
/*
 * Coefficient accounting for the influence of visco-elasticity and viscous flow
 * on the modulus of elasticity of the stress-strain relation
 * zh 13.6.2004:
 * kv - coefficient for volumetric stress
 * gv - coefficient for deviatoric stress
 * With moDeviatoricCreepE/F option, viscous/flow creep is evaluated only based on deviatoric part of stress tensor,
 * so the kv coefficient doesn't contain the flow creep part.
 */

{
    if ( !( ( mode == VM_Total ) || ( mode == VM_Incremental ) ) ) {
        _error2("giveKvCoeffs: unknown mode %d (use VM_Total for kv/gv, VM_Incremental for ev/fv)", mode);
    }

    if ( options & moCreep ) {
        if ( options & moIsothermal ) {
            gp = giveMaterialGp();
        }

        double dt, ksi, ev, fv;
        dt = giveTimeIncrement(atTime);
        if ( dt > 0 ) {
            ksi = giveHydrationDegree(gp, atTime, VM_Total);

            // viscous creep coefficient ev
            ev = 1. / ( 1. + twTime(ksi) / dt ) * agingE(ksi) * jv;

            // compute flow creep coefficient fv
            if ( ( ( HellmichMaterialStatus * ) giveStatus(gp) )->giveGamma0() > 0 ) {
                fv = dt * agingE(ksi) * giveViscosity(gp);
            } else {
                fv = 0.;
            }

            if ( mode == VM_Total ) {
                gv = 1. + ev + fv; // everything is always included in deviatoric coefficient
                kv = 1.;
                if ( !( options & moDeviatoricCreepE ) ) {
                    kv += ev;                       // add visc creep to volumetric coeff
                }

                if ( !( options & moDeviatoricCreepF ) ) {
                    kv += fv;                       // add flow creep to volumetric coeff
                }
            } else {
                kv = ev;
                gv = fv;
            }

            return;
        }
    }

    if ( mode == VM_Total ) {
        kv = gv = 1.;
    } else {
        kv = gv = 0.;
    }
}

// === Temperature strains ===
void HellmichMaterial :: giveThermalDilatationVector(FloatArray &answer,
                                                     GaussPoint *gp, TimeStep *atTime)
//
// returns a FloatArray(6) of initial strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
{
    answer.resize(6);
    answer.zero();
    double aux = this->give(tAlpha, gp);
    answer(0) = aux;
    answer(1) = aux;
    answer(2) = aux;
}

void HellmichMaterial :: giveShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                                   GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
{
    int size = 6;
    answer.resize(0);
    MaterialMode mmode = gp->giveMaterialMode();

    if ( ( ( options & moHasShrinkage ) == moHasShrinkage ) || ( options & moHumidityStrain ) ) {
        double eshr = 0., eshrh, ksi;
        int i;

        if ( ( mmode == _3dShell ) || ( mmode ==  _3dBeam ) || ( mmode == _2dPlate ) || ( mmode == _2dBeam ) ) {
            answer.resize(12);
            size = 12;
        } else {
            answer.resize(6);
            size = 6;
        }

        // shrinkage as function of hydration degree
        if ( ( options & moHasShrinkage ) == moHasShrinkage ) {
            ksi = giveHydrationDegree(gp, atTime, VM_Total);
            eshr = autoShrinkageCoeff(ksi);
            if ( mode == VM_Incremental ) {
                if ( giveTimeIncrement(atTime) > 0 ) {
                    eshr -= autoShrinkageCoeff( ksi - giveHydrationDegree(gp, atTime, VM_Incremental) );
                } else {
                    eshr = 0;
                }
            }
        }

        // shrinkage as function of relative humidity
        if ( options & moHumidityStrain ) {
            // eshrh = K (1 - h^3), K = -0,0008 [Desky... Beton TKS 2/2002]
            eshrh = dryingShrinkageCoeff( giveHumidity(gp, VM_Total) );
            if ( mode == VM_Incremental ) {
                eshrh -= dryingShrinkageCoeff( giveHumidity(gp, VM_Velocity) );              // Previous
            }

            if ( eshrh ) {
                eshr += eshrh;
            }
        }

        for ( i = 0; i < 3; i++ ) {
            answer(i) = eshr;
        }

        for ( i = 3; i < size; i++ ) {
            answer(i) = 0;
        }

        if ( form == ReducedForm ) {
            FloatArray auxvec = answer;
            ( ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection() )->giveReducedCharacteristicVector(answer, gp, auxvec);
        }
    }
}

void HellmichMaterial :: givePrestressStrainVector(FloatArray &answer, MatResponseForm form,
                                                   GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
{
    answer.resize(0);
    if ( prestress ) {
        double ksi, prestr, time = giveTime(atTime);
        prestr = prestressValue(time);
        if ( mode == VM_Incremental ) {
            if ( giveTimeIncrement(atTime) > 0 ) {
                prestr -= prestressValue( time - giveTimeIncrement(atTime) );
            } else {
                prestr = 0;
            }
        }

        if ( prestr ) {
            ksi = giveHydrationDegree(gp, atTime, VM_Total);
            if ( ksi < HYDRATION_MINDEGREE ) {
                _error("Can't apply prestress: zero hydration degree!");
            }

            answer.resize(6);
            answer.zero();
            answer(0) = prestr / agingE(ksi);
        }
    }

    if ( form == ReducedForm ) {
        FloatArray auxvec = answer;
        ( ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection() )->giveReducedCharacteristicVector(answer, gp, auxvec);
    }
}

void HellmichMaterial :: giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                               GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
// returns temperature strains based on current temperature, not on prescribed temperature loads
// + creep strains: (Kv-1)C:sig_n - 1/(1+tau/dt)*eps_v,n
{
    FloatArray stress, auxvec, redvec, et;
    double dT, dt, ev, fv;
    StructuralElement *elem = ( StructuralElement * ) gp->giveElement();
    StructuralCrossSection *crossSection =  dynamic_cast< StructuralCrossSection * >( gp->giveCrossSection() );
    MaterialMode mmode = gp->giveMaterialMode();

    // do not generate any temperature
    if ( atTime->giveTargetTime() < castAt ) {
        if ( ( mmode == _3dShell ) || ( mmode ==  _3dBeam ) || ( mmode == _2dPlate ) || ( mmode == _2dBeam ) ) {
            answer.resize(12);
        } else {
            answer.resize(6);
        }

        answer.zero();
        return;
    }


    // StructuralTemperatureLoad
    elem->computeResultingIPTemperatureAt(et, atTime, gp, mode);

    // thermal dilatation from intrinsic temperature (independent on StructuralTemperatureLoad)
    answer.resize(0);
    dT = giveTemperatureChange(gp, mode);


    if ( et.isNotEmpty() ) {
        dT += et.at(1);
    }

    if ( dT ) {
        giveThermalDilatationVector(auxvec, gp, atTime); // {tAlpha, tAlpha, tAlpha, 0, 0, 0}
        if ( ( mmode == _2dBeam ) || ( mmode == _3dBeam ) || ( mmode == _3dShell ) || ( mmode == _2dPlate ) ) {
            double thick = crossSection->give(CS_Thickness);
            if ( mmode == _2dBeam ) {
                answer.resize(12);
                answer.zero();
                answer.at(1) = auxvec.at(1) * ( dT );
                if ( et.giveSize() > 1 ) {
                    answer.at(8) = auxvec.at(1) * et.at(2) / thick; // kappa_x
                }
            } else if ( mmode == _3dBeam ) {
                answer.resize(12);
                answer.zero();
                double width = crossSection->give(CS_Width);
                answer.at(1) = auxvec.at(1) * ( dT );
                if ( et.giveSize() > 1 ) {
                    answer.at(8) = auxvec.at(1) * et.at(2) / thick; // kappa_y
                    if ( et.giveSize() > 2 ) {
                        answer.at(9) = auxvec.at(1) * et.at(3) / width; // kappa_z
                    }
                }
            } else if ( mmode == _2dPlate ) {
                if ( et.giveSize() > 1 ) {
                    answer.resize(12);
                    answer.zero();

                    if ( et.giveSize() > 1 ) {
                        answer.at(7) = auxvec.at(1) * et.at(2) / thick; // kappa_x
                        answer.at(8) = auxvec.at(2) * et.at(2) / thick; // kappa_y
                    }
                }
            } else if ( mmode == _3dShell ) {
                answer.resize(12);
                answer.zero();

                answer.at(1) = auxvec.at(1) * ( dT );
                answer.at(2) = auxvec.at(2) * ( dT );
                if ( et.giveSize() > 1 ) {
                    answer.at(7) = auxvec.at(1) * et.at(2) / thick; // kappa_x
                    answer.at(8) = auxvec.at(2) * et.at(2) / thick; // kappa_y
                }
            }
        } else {
            if ( auxvec.giveSize() ) {
                answer = auxvec;
                answer.times(dT);
            }
        }
    }


    if ( prestress ) {
        givePrestressStrainVector(auxvec, FullForm, gp, atTime, mode); // {prestress/E, 0, 0, 0, 0, 0}
        if ( auxvec.giveSize() ) {
            if ( !answer.giveSize() ) {
                answer.resize(6);
                answer.zero();
            }

            answer.add(auxvec);
        }
    }

    dt = giveTimeIncrement(atTime);
    if ( ( options & moCreep ) && ( dt > 0 ) ) {
        // reduced stress vector
        HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) giveStatus(gp);
        stress = status->giveStressVector();
        // compute creep coefficients
        giveKvCoeffs(gp, atTime, ev, fv, VM_Incremental);

        // compute redvec tensor for visc creep - stress or stress deviator; deviatoric works only for 3D
        if ( options & moDeviatoricCreepE ) {
            deviator(redvec, stress);
        } else {
            redvec = stress;
        }

        redvec.times(ev); // viscous creep
        // compute auxvec tensor for flow creep
        if ( options & moDeviatoricCreepF ) {
            deviator(auxvec, stress);
        } else {
            auxvec = stress;
        }

        auxvec.times(fv); // flow creep

        redvec.add(auxvec); // Sum into redvec

        // Compute auxvec creep strains, don't use creep coeffs
        elasticCompliance(auxvec, redvec, gp, atTime, ReducedForm, false);
        redvec = auxvec;

        // subtract previous viscous creep strains dev_n(auxvec)
        status->giveViscousStrainVector(auxvec);
        auxvec.times( 1. / ( 1. + twTime( giveHydrationDegree(gp, atTime, VM_Total) ) / dt ) );
        redvec.subtract(auxvec); // can't use ev, here isn't JvE

        // convert to FullForm and add to answer
        StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();
        crossSection->giveFullCharacteristicVector(auxvec, gp, redvec);
        if ( !answer.giveSize() ) {
            if ( ( mmode == _3dShell ) || ( mmode ==  _3dBeam ) || ( mmode == _2dPlate ) || ( mmode == _2dBeam ) ) {
                answer.resize(12);
            } else {
                answer.resize(6);
            }

            //answer.resize(6); answer.zero();
        }

        answer.add(auxvec);
    }
}

void HellmichMaterial :: giveRealStressVector(FloatArray &answer,
                                              MatResponseForm form, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime)
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current strain increment
// calls stress return method - updates temporary status variables
//
// formulated in reduced stress-strain space
{
    FloatArray auxStrain, strainIncrement, elasticStrain;
    FloatArray auxStress, trialStressVector, fullStressVector(6), redStressVector;
    HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();
    ValueModeType mode = VM_Incremental;
    MaterialMode mmode = gp->giveMaterialMode();
    double E, dt;
    // clear temporary status variables - maybe unnecessary
    initTempStatus(gp);
    // if needed, update auxiliary status values in gp and material
    initAuxStatus(gp, atTime);
    // === evaluate the trial stress ===
    if ( mode == VM_Incremental ) {
        // trial stress = sigma_n  +  1/Kv * (Ce:(depsep = deps - deT - deshr + evaux*dev_n) - (evaux*Jv + dt*viscosity)*E*sigma_n)
        // !!! auxiliary status - stored in material for stress return
        // e.g. parallel - goes to hell (several gp's giveRealStressVector use the same material auxiliary variables)
        auxksi = giveHydrationDegree(gp, atTime, VM_Total);
        auxchi1 = status->giveHardeningVar();
        giveKvCoeffs(gp, atTime, auxKv, auxGv, VM_Total);
        E = agingE(auxksi);
        // --- Strain increment ---
        auxStrain = totalStrain;
        auxStrain.subtract( status->giveStrainVector() );
        // subtract eigenstrains and shrinkage strains (eps_T, eps_v, eps_f, eps_shr)
        // OOFEM: calls StructuralElement->computeStressIndependentStrainVector (reduced)
        //  there sums computeTemperatureStrainVectorAt + CrossSection->computeStressIndependentStrainVector
        //  ...TemperatureEps multiplies mat->giveThermalDilatationVector * computeResultingIPTemperature...sum of thermal loads
        //  ...CrossSect sums mat->giveShrinkageStrainVector...(eps_shr) + mat->giveEigenStrainVector...(thermal + creep)
        giveStressDependentPartOfStrainVector(strainIncrement, gp, auxStrain, atTime, VM_Incremental);
        // --- Trial Stress ---
        // original: sigtr = sig_n + [De]/Kv : de;
        // devcreep: sigtr = sig_n + [Dev] : de;
        //  1D: dsig = E/Kv*de, 3D: dsig = 1/Kv*[(1x1)K + 2G Idev]:de, jinak 1/Kv * (linMat->De):de
        elasticStiffness(trialStressVector, strainIncrement, gp, atTime, ReducedForm, true);
        trialStressVector.add( status->giveStressVector() );
    }
    // ------ ( Total Mode ) -------
    else if ( mode == VM_Total ) {
        if ( options & moCreep ) {
            _error("Can't use VM_Total with creep.");
        }

 #ifdef DEBUG
        _warning1("HellmichMaterial::giveRealStress - VM_Total!\n");
 #endif

        auxKv = auxGv = 1.; // no creep
        auxStrain = totalStrain;
        giveStressDependentPartOfStrainVector(elasticStrain, gp, auxStrain, atTime, VM_Total);
        if ( options & moPlasticity ) {
            status->givePlasticStrainVector(auxStrain);
            elasticStrain.subtract(auxStrain);
        }

        // trial stress = Ce * eps_n+1; no creep - can use Kv and Gv, but doesn't
        elasticStiffness(trialStressVector, elasticStrain, gp, atTime, ReducedForm, false);
    }
    // ------ ( Total Konec ) ------
    else {
        _error("giveRealStressVector: Unknown ValueModeType!");
    }

    // save reduced Trial Stress vector
    if ( options & moPlasticity ) {
        status->setTrialStressVector(trialStressVector);
    }

    // ======  1D  ======
    // uniaxial stress - elastic limit in tension (ft) and compression (Rc with hardening)
    // formulated in reduced stress and strain
    if ( mmode == _1dMat ) {
        // auxiliary variables for 1D
        double depsp, sig, sigtr, tc, rc, fca, dl;
        ActiveSurface as;

        sigtr = trialStressVector(0);
        sig = sigtr;

        if ( ( options & moPlasticity ) && ( auxksi > ksi0 ) ) {
            fca = fcStrength(auxksi);
            tc = delta * fca;
            rc = rcThreshold(auxchi1, auxksi);
            if ( sig > tc ) { // Tension - cut-off
                as = asTC;
                dl = 0;
                depsp = ( sig - tc ) * auxKv / ( 9 * agingK(auxksi) );
                sig = tc;
            } else if ( ( fca > 0 ) && ( -sig > rc ) ) { // Compression - DP
                as = asDP;
                double EKv = auxkdp * E / auxKv;
                dl = computedl(-sig, EKv);
                depsp = -auxkdp * dl;
                sig += EKv * dl;
            } else { // Elastic
                as = asNone;
                dl = 0;
                depsp = 0;
            }

            // save Plastic status
            if ( options & moHardening ) {
                status->setTempHardeningVar(auxchi1 + 3 * alpha * dl);
            }

            status->givePlasticStrainVector(auxStrain);
            auxStrain(0) += depsp;
            status->setTempPlasticStrainVector(auxStrain);

            // temporary plastic status update
            status->setActiveSurface(as);
            status->setFlowIncrement(dl);
        }

        redStressVector.resize(1);
        redStressVector(0) = sig;
        fullStressVector.zero();
        fullStressVector(0) = sig;
    }
    // ====== 1D Konec ======
    else if ( mmode == _PlaneStress ) {
        // ====== 2D Plane Stress - not working! ======
        // maybe should work 3D projection + stiffness matrix inversion->reduction?
        // auxiliary variables for 2D
        // !! CORNER !!
        _error("2D Plane stress plasticity not implemented!");
        double sigx, sigy, tau, tc, rc, fca, dl, fDP, Itrial, snorm;
        FloatArray depsp(3);
        ActiveSurface as;

        if ( ( options & moPlasticity ) && ( auxksi > ksi0 ) ) {
            sigx = trialStressVector(0);
            sigy = trialStressVector(1);
            tau = trialStressVector(2);
            fca = fcStrength(auxksi);
            tc = delta * fca;
            if ( sigx + sigy > tc ) { // Tension - cut-off ... ok
                as = asTC;
                dl = 0;
                depsp(0) = ( sigx + sigy - tc ) * auxGv / ( 4 * agingG(auxksi) ) * ( 1 - ny ) / ( 1 + ny );
                depsp(1) = depsp(0);
                depsp(2) = 0;
                sigx -= 0.5 * ( sigx + sigy - tc );
                sigy -= 0.5 * ( sigx + sigy - tc );
            } else {
                rc = rcThreshold(auxchi1, auxksi);
                Itrial = sigx + sigy;
                snorm = sqrt(2. / 3. * ( sigx * sigx + sigx * sigy + sigy * sigy ) + 2 * tau * tau);
                fDP = alpha * Itrial + snorm - auxkdp * rc;
                if ( ( fca > 0 ) && ( -fDP > 0 ) ) { // Compression - DP ... ??
                    as = asDP;

                    tempb = -1 / auxKv * ( 9 * alpha * alpha * agingK(auxksi) + 2 * agingG(auxksi) );
                    dl = computedl( ( alpha * Itrial + snorm ) / auxkdp, tempb / auxkdp );

                    depsp(0) = dl * ( alpha + ( 2. / 3. * sigx - 1. / 3. * sigy ) / snorm );
                    depsp(1) = dl * ( alpha + ( 2. / 3. * sigy - 1. / 3. * sigx ) / snorm );
                    depsp(2) = 2 * dl * tau / snorm;

                    // devcreep - use [Dev(Kv, Gv)] instead of [De]/Kv
                    elasticStiffness(redStressVector, depsp, gp, atTime, ReducedForm, true);
                    redStressVector.add(trialStressVector);
                } else { // Elastic
                    as = asNone;
                    dl = 0;
                    depsp.zero();
                }
            }

            // save Plastic status
            if ( options & moHardening ) {
                status->setTempHardeningVar(auxchi1 + 3 * alpha * dl);
            }

            status->givePlasticStrainVector(auxStrain);
            auxStrain.add(depsp);
            status->setTempPlasticStrainVector(auxStrain);

            // temporary plastic status update
            status->setActiveSurface(as);
            status->setFlowIncrement(dl);
            if ( as != asDP ) {
                redStressVector.resize(3);
                redStressVector(0) = sigx;
                redStressVector(1) = sigy;
                redStressVector(2) = tau;
            }
        } else {
            redStressVector = trialStressVector;
        }

        fullStressVector.zero();
        fullStressVector(0) = redStressVector(0);
        fullStressVector(1) = redStressVector(1);
        fullStressVector(5) = redStressVector(2);
    }
    // ====== 2D Konec ======
    else {
        if ( ( options & moPlasticity ) && ( auxksi > ksi0 ) ) {
            if ( mmode == _3dMat ) {
                // stress return for 3D; sets temp plastic status, overwrites fullStressVector
                stressReturn(redStressVector, trialStressVector, gp, atTime);
                fullStressVector = redStressVector; // 3D - full=reduced
            } else {
                _error("giveRealStress: unsupported material mode");
            }
        } else { // elastic
            redStressVector = trialStressVector;
            crossSection->giveFullCharacteristicVector(fullStressVector, gp, redStressVector);
        }
    }

    // update creep strains - reduced form
    dt = giveTimeIncrement(atTime);
    if ( ( options & moCreep ) && ( dt > 0 ) ) {              // moDeviatoricCreepE
        // viscous creep - ev = ev,n + 1/(1+tw/dt)*(Jv*E*compl([sigma | deviator(sigma)])-ev,n)
        // compute auxStrain stress tensor for visc creep - stress or stress deviator
        if ( options & moDeviatoricCreepE ) {
            deviator(auxStress, redStressVector);
        } else {
            auxStress = redStressVector;
        }

        elasticCompliance(elasticStrain, auxStress, gp, atTime, ReducedForm, false); // using no Kv, Gv!
        elasticStrain.times(jv * E);
        status->giveViscousStrainVector(auxStrain);
        elasticStrain.subtract(auxStrain);
        elasticStrain.times( 1. / ( 1. + twTime(auxksi) / dt ) );

        elasticStrain.add(auxStrain);
        status->setTempViscousStrainVector(elasticStrain);
        // moDeviatoricCreepF
        // flow creep - ef = ef,n + dt*viscosity*E*compl([sigma | deviator(sigma)])
        // compute auxvec tensor for flow creep
        if ( options & moDeviatoricCreepF ) {
            deviator(auxStress, redStressVector);
        } else {
            auxStress = redStressVector;
        }

        elasticCompliance(elasticStrain, auxStress, gp, atTime, ReducedForm, false); // using no Kv, Gv!
        elasticStrain.times(dt * giveViscosity(gp) * E);
        status->giveFlowStrainVector(auxStrain);
        elasticStrain.add(auxStrain);
        status->setTempFlowStrainVector(elasticStrain);
    } // end creep

    // return full/reduced stress vector
    if ( form == FullForm ) {
        answer = fullStressVector;
    } else {
        answer = redStressVector;
    }

    // update temp status stress and strain vector (reduced)
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    // Stress-strain plot output
    if ( options & moPlotStressStrainIter ) {
        FILE *sstream;
        sstream = fopen(OUTPUTFILE_SS, "a");
        if ( sstream ) {
            plotStressStrain( sstream, gp, atTime, pssIndex, gp->giveElement()->giveNumber(), ( double ) ( gp->giveNumber() ) );
            fclose(sstream);
        }
    }
}

void
HellmichMaterial :: computeStressIndependentStrainVector(FloatArray &answer,
                                                         GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
{
    FloatArray fullAnswer, et;

    // if needed, update auxiliary status values in gp and material
    initAuxStatus(gp, atTime);

    giveShrinkageStrainVector(fullAnswer, FullForm, gp, atTime, mode); // Autogenous shrinkage and humidity volume changes
    giveEigenStrainVector(et, FullForm, gp, atTime, mode); // Material temperature and temperature loads, creep
    if ( et.giveSize() ) {
        if ( fullAnswer.giveSize() ) {
            fullAnswer.add(et);
        } else {
            fullAnswer = et;
        }
    }

    // export reduced answer
    if ( fullAnswer.giveSize() ) {
        this->giveReducedCharacteristicVector(answer, gp, fullAnswer);
        return;
    }

    answer.resize(0);
}

contextIOResultType HellmichMaterial :: saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp)
// saves full status for this material, also invokes saving
// for sub-objects of this.
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex: can't write into NULL stream");
    }

    // save material gp status if necessary
    if ( gp && options & moIsothermal ) {
        StateCounterType c = gp ->giveElement()->giveDomain()->giveEngngModel()->giveCurrentStep()->giveSolutionStateCounter();
        if ( materialGpSaveAt != c ) {
            materialGpSaveAt = c;
            if ( ( iores = saveIPContext( stream, mode, giveMaterialGp() ) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        }
    }

    // write parent data
    if ( ( iores = StructuralMaterial :: saveIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // save hydration model data - maybe should check moHydration?
    // needs to save only in case nonisodata is present = moIsothermal option is not set in gp options
    if ( gp && !( ( ( HellmichMaterialStatus * ) giveStatus( gp ) )->giveMaterialOptions() & moIsothermal ) ) {
      if ( ( iores = HydrationModelInterface :: saveContext(stream, mode, (void*) gp) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}
contextIOResultType HellmichMaterial :: restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp)
// restores full status for this material, also invokes restoring for sub-objects of this.
// !!! gp is passed in obj
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("restoreContext: can't read from NULL stream");
    }

    // read material gp status if necessary
    if ( gp && options & moIsothermal ) {
        StateCounterType c = gp ->giveElement()->giveDomain()->giveEngngModel()->giveCurrentStep()->giveSolutionStateCounter();
        if ( materialGpRestoreAt != c ) {
            materialGpRestoreAt = c;
            if ( ( iores = restoreIPContext( stream, mode, giveMaterialGp() ) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        }
    }

    // read parent data
    if ( ( iores = StructuralMaterial :: restoreIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read hydration model data - maybe should check moHydration?
    if ( gp && !( ( ( HellmichMaterialStatus * ) giveStatus( gp ) )->giveMaterialOptions() & moIsothermal ) ) {
      if ( ( iores = HydrationModelInterface :: restoreContext(stream, mode, gp) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}

int HellmichMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
// should be: without plasticity, return giveElasticMaterial()->hasMaterialModeCapability(mode);
// with plasticity, currently only 3D material response is supported, should like to implement 1D plasticity
// (may be independent on 3D functions, or redefining them to check mode).
//
{
    if ( options & ( moDeviatoricCreepE | moDeviatoricCreepF ) ) { // deviatoric creep supported only in 3D analysis
        if ( mode == _3dMat ) {
            return 1;
        } else {
            return 0;
        }
    } else if ( options & moPlasticity ) { // 3D and 1D plasticity fully implemented
        if ( ( mode == _3dMat ) || ( mode == _1dMat ) ) {
            return 1;
        } else {
            return 0;
        }
    } else {
        return ( linMat->hasMaterialModeCapability(mode) );
    }

    /* LinearElastic:
     *
     * if ((mode == _3dMat) || (mode == _PlaneStress) ||
     * (mode == _PlaneStrain) || (mode == _1dMat) ||
     * (mode == _2dPlateLayer) || (mode == _2dBeamLayer) ||
     * (mode == _3dShellLayer) || (mode == _2dPlate) ||
     * (mode == _2dBeam) || (mode == _3dShell) ||
     * (mode == _3dBeam) || (mode == _PlaneStressRot) ||
     * (mode == _1dFiber)) return 1;
     * return 0;
     */
}

void
HellmichMaterial :: initGpForNewStep(GaussPoint *gp)
/*
 * Step restart initialization
 * Standard function - calls initTempStatus(gp)
 * Furthermore, calls initAuxStatus(gp, currentStep) for each gp, in any case of status->nonisodata->updateflag
 * also calls initAuxStatus(materialGp), but only if materialGpInitAt is not uptodate.
 */
{
    // get current time step - needed for initAuxStatus
    TimeStep *atTime = gp->giveElement()->giveDomain()->giveEngngModel()->giveCurrentStep();

    // initialize the material-level status
    if ( ( options & moIsothermal ) && ( materialGpInitAt < atTime->giveSolutionStateCounter() ) ) {
        // avoid multiple calls
        materialGpInitAt = atTime->giveSolutionStateCounter();
        // avoid calling material gp init again from initAuxStatus(materialGp)
        materialGpUpdateFlag = 1;
        // call init for materialGp
        initGpForNewStep( giveMaterialGp() );
    }

    // call parent method - calls initTempStatus(gp)
    StructuralMaterial :: initGpForNewStep(gp);
    // Force initialization of auxiliary status
    ( ( HellmichMaterialStatus * ) giveStatus(gp) )->setUpdateFlag(0);
    initAuxStatus(gp, atTime);
}

void
HellmichMaterial :: initTempStatus(GaussPoint *gp)
{
    HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) this->giveStatus(gp);
    status->initTempStatus();
}

void
HellmichMaterial :: updateYourself(GaussPoint *gp, TimeStep *atTime)
{
    // update the material-level status if necessary (once per material)
    if ( ( options & moIsothermal ) && materialGpUpdateFlag ) {
 #ifdef VERBOSE_HELLMAT
        printf( "Hellmat: Material-level status updated for time=%.0f, flag 1->0.\n", giveTime(atTime) );
 #endif
        materialGpUpdateFlag = 0;
        updateYourself(giveMaterialGp(), atTime);
    }

    HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) this->giveStatus(gp);
    if ( status ) {
        status->updateYourself(atTime);
    }
}

void
HellmichMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                             MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime)
{
    MaterialMode mMode = gp->giveMaterialMode();

    // if needed, update auxiliary status values in gp and material
    initAuxStatus(gp, atTime);

    // if (giveTempHydrationDegree(gp)<HYDRATION_MINDEGREE)
    //  _warning("HellmichMaterial::giveCharMtrx: No stiffness - zero hydration degree.",10);

    if ( mMode == _1dMat ) {
        give1dMaterialStiffnessMatrix(answer, form, rMode, gp, atTime);
    } else if ( mMode == _3dMat ) {
        give3dMaterialStiffnessMatrix(answer, form, rMode, gp, atTime);
    } else if ( !( options & moPlasticity ) ) {
        giveLinearElasticMaterial(gp, atTime)->giveCharacteristicMatrix(answer, form, rMode, gp, atTime);
        answer.times( 1 / giveKvCoeff(gp, atTime) );
    } else {
        _error("giveCharMtrx: unsupported stress-strain mode!");
    }

 #ifdef VERBOSE_TANGENT
    if ( options & moPlotStressStrainIter ) {
        HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) giveStatus(gp);
        if ( mMode == _1dMat ) {
            printf( "Time: %.2f, sig: %.10e, Eep: %.10e\n", giveTime(atTime),
                   ( status->giveTempStressVector().giveSize() ) ? ( status->giveTempStressVector() )(0) : 0, answer(0, 0) );
        } else {
            answer.printYourself();
        }
    }

 #endif
}

void HellmichMaterial :: setMixture(MixtureType mix)
{
    mixture = mix;
    // Set also hydration model to use given mixture
    /* !!! Ensure that hydrationModel is available at time of calling setMixture
     * - hydrationModel is set up at beginning of initializeFrom if applicable */
    if ( hydrationModel ) {
        hydrationModel->setMixture(mix);
    }
}

LinearElasticMaterial *
HellmichMaterial :: giveLinearElasticMaterial(GaussPoint *gp, TimeStep *atTime)
{
    if ( !linMat ) {
        _error("giveLinearElasticMaterial: base linear elastic material is not initialized.");
    }

    double ksi = giveHydrationDegree(gp, atTime, VM_Total);
    if ( ksi < HYDRATION_MINDEGREE ) {
        ksi = HYDRATION_MINDEGREE;
    }

    linMat->setE( agingE(ksi) );
    return linMat;
}

int
HellmichMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
// returns the equilibrated IP value. should be called after step update or context restore - temp=equlib
{
    HellmichMaterialStatus *status = ( HellmichMaterialStatus * ) this->giveStatus(aGaussPoint);
    if ( type == IST_PlasticStrainTensor ) {
        if ( options & moPlasticity ) {
            status->givePlasticStrainVector(answer);
        } else {
            answer.resize( this->giveSizeOfReducedStressStrainVector( aGaussPoint->giveMaterialMode() ) );
            answer.zero();
        }

        return 1;
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        int indx;
        FloatArray st(6), s;
        if ( options & moPlasticity ) {
            status->givePlasticStrainVector(s);
            for ( int i = 1; i <= s.giveSize(); i++ ) {
                indx = this->giveStressStrainComponentIndOf(ReducedForm, aGaussPoint->giveMaterialMode(), i);
                if ( indx ) {
                    st.at(indx) = s.at(i);
                }
            }

            this->computePrincipalValues(answer, st, principal_strain);
        } else {
            answer.resize(3);
            answer.zero();
        }

        return 1;
    } else if ( type == IST_DamageTensor ) { // return hardening internal variable as first component of damage tensor
        answer.resize(1);
        answer.at(1) = status->giveHardeningVar() / chi1u;
        return 1;
    } else if ( type == IST_HydrationDegree ) {
        answer.resize(1);
        answer.at(1) = giveHydrationDegree(aGaussPoint, atTime, VM_Total);
        return 1;
    } else if ( type == IST_Temperature ) {
        answer.resize(1);
        answer.at(1) = giveTemperature(aGaussPoint);
        return 1;
    } else if ( type == IST_Humidity ) {
        answer.resize(1);
        answer.at(1) = giveHumidity(aGaussPoint, VM_Total);
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}
MaterialStatus *
HellmichMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new HellmichMaterialStatus(1, this->giveDomain(), gp);
}

double HellmichMaterial :: give(int aProperty, GaussPoint *gp)
/*
 * Returns the value of the property aProperty.
 * Only the final values of the Young and shear modulus andPoisson ratio are returned.
 */
{
    if ( ( aProperty == 'G' ) || ( aProperty == Gyz ) || ( aProperty == Gxz ) ||
        ( aProperty == Gxy ) ) {
        return agingG(1.0);
    }

    if ( ( aProperty == 'E' ) || ( aProperty == Ex ) || ( aProperty == Ey ) ||
        ( aProperty == Ez ) ) {
        return agingE(1.0);
    }

    if ( ( aProperty == NYxy ) || ( aProperty == NYxz ) || ( aProperty == NYyz ) ) {
        return ny;
    }

    if ( ( aProperty == 'n' ) || ( aProperty == NYzx ) || ( aProperty == NYzy ) ||
        ( aProperty == NYyx ) ) {
        return ny;
    } else {
        return this->Material :: give(aProperty, gp);
    }
}

// === Postprocessing ===
InternalStateValueType
HellmichMaterial :: giveIPValueType(InternalStateType type)
{
    // strains components packed in enginnering notation
    if ( ( type == IST_PlasticStrainTensor ) || ( type == IST_PrincipalPlasticStrainTensor ) ) {
        return ISVT_TENSOR_S3E;
    } else if ( type == IST_DamageTensor ) {
        return ISVT_TENSOR_G;
    } else if ( type == IST_HydrationDegree ) {
        return ISVT_SCALAR;
    } else if ( type == IST_Temperature ) {
        return ISVT_SCALAR;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}

int
HellmichMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_PlasticStrainTensor ) {
        this->giveStressStrainMask(answer, FullForm, mmode);
        return 1;
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        answer.resize(6);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        return 1;
    } else if ( type == IST_DamageTensor ) { // used for internal variables
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( ( type == IST_HydrationDegree ) || ( type == IST_Temperature ) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}

int
HellmichMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_PlasticStrainTensor ) {
        return this->giveSizeOfReducedStressStrainVector( aGaussPoint->giveMaterialMode() );
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        return 3;
    } else if ( ( type == IST_DamageTensor ) || ( type == IST_HydrationDegree ) || ( type == IST_Temperature ) ) {
        return 1;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}
// end of HellmichMaterial implemantation

// ========================  STATUS implementation ===================================
PlastData :: PlastData() : plasticStrainVector(), tempPlasticStrainVector(), trialStressVector()
{
    hardeningVar = 0;
    tempHardeningVar = 0;
    activeSurface = asNone;
    flowIncrement = 0;
}


NonisoData :: NonisoData()
/*
 * Initializes the non-isothermal status variables.
 * Also used as material level status for isothermal analysis. Temperature is initialized via setInitialTemperature.
 */
{
    // creep
    gamma0 = -1.;
    viscousSlip = 0;
    tempViscousSlip = 0;
    viscosity = 0;
    humidity = tempHumidity = 1.;
    // Set status update flag to zero - means auxiliary non-isothermal status values
    //  will be initialized when initAuxStatus is called.
    auxStatusUpdateFlag = 0;
}

HellmichMaterialStatus :: HellmichMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g), HydrationModelStatusInterface()
{
    // Get material options; gp number 0 => material level status
    MaterialOptions options = giveMaterialOptions();

    // Non-isothermal
    if ( !( options & moIsothermal ) ) {
        nonisoData = new NonisoData();
        // uninitialized
        setInitialTemperature(-1.);
    } else {
        nonisoData = NULL;
    }

    // Plasticity
    if ( options & moPlasticity ) {
        plastData = new PlastData();
    } else {
        plastData = NULL;
    }

    // Creep
    if ( options & moCreep ) {
        creepData = new CreepData();
    } else {
        creepData = NULL;
    }
}

HellmichMaterialStatus :: ~HellmichMaterialStatus()
{
    delete plastData;
    delete creepData;
    delete nonisoData;
}

MaterialOptions
HellmichMaterialStatus :: giveMaterialOptions()
// Returns the material options mask.
{
    if ( gp->giveNumber() ) {
        return ( ( HellmichMaterial * ) gp->giveMaterial() )->giveOptions();
    } else {
        return moMaterialLevel;
    }
}

void
HellmichMaterialStatus :: setInitialTemperature(double v)
{
    if ( nonisoData ) {
        nonisoData->temperature = v;
        nonisoData->tempTemperature = v;
        nonisoData->initialTemperature = v;
    }
}

// necessary for proper cast to interface, can't be done from outside
Interface *
HellmichMaterialStatus :: giveInterface(InterfaceType type)
{
    if ( type == HydrationModelStatusInterfaceType ) {
        return ( HydrationModelStatusInterface * ) this;
    } else {
        return NULL;
    }
}

// === Temporary status initialization ===
void PlastData :: initTempStatus(GaussPoint *gp)
{
    // initialization
    if ( !plasticStrainVector.giveSize() ) {
        int n = ( ( HellmichMaterial * ) ( gp->giveMaterial() ) )->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() );

        plasticStrainVector.resize(n);
    }

    // initTemp
    tempPlasticStrainVector = plasticStrainVector;
    tempHardeningVar = hardeningVar;
    // clear temp plastic status
    activeSurface = asNone;
    flowIncrement = 0;
}

void CreepData :: initTempStatus(GaussPoint *gp)
{
    // initialization
    if ( !viscousStrainVector.giveSize() ) {
        int n = ( ( HellmichMaterial * ) ( gp->giveMaterial() ) )->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() );
        viscousStrainVector.resize(n);
        flowStrainVector.resize(n);
    }

    // initTemp
    tempViscousStrainVector = viscousStrainVector;
    tempFlowStrainVector = flowStrainVector;
}

void HellmichMaterialStatus :: initTempStatus()
/*
 * Sets temp variables to last equilibrium state. Called at start of each iteration.
 * Allocates the status vectors if necessary.
 */
{
    // StructuralMaterialStatus initialized even case of material level status, but that shouldn't matter
    StructuralMaterialStatus :: initTempStatus();
    if ( plastData ) {
        plastData->initTempStatus(gp);
    }

    if ( creepData ) {
        creepData->initTempStatus(gp);
    }

    // nonisoData is not initialized in each iteration, it's initialized in initAuxStatus at step start
}

// === Status update ===
void PlastData :: updateYourself()
{
    plasticStrainVector = tempPlasticStrainVector;
    hardeningVar = tempHardeningVar;
}
void CreepData :: updateYourself()
{
    viscousStrainVector = tempViscousStrainVector;
    flowStrainVector = tempFlowStrainVector;
}
void NonisoData :: updateYourself()
{
    temperature = tempTemperature;
    humidity = tempHumidity;

    viscousSlip = tempViscousSlip;

    // Set status update flag to zero - means auxiliary non-isothermal status values
    //  will be initialized when initAuxStatus is called.
    auxStatusUpdateFlag = 0;
}
void HellmichMaterialStatus :: updateYourself(TimeStep *atTime)
// Update after equilibrium has been reached: temp -> equilib.
{
    HydrationModelStatusInterface :: updateYourself(atTime);
    StructuralMaterialStatus :: updateYourself(atTime);

    if ( nonisoData ) {
        nonisoData->updateYourself();
    }

    if ( plastData ) {
        plastData->updateYourself();
    }

    if ( creepData ) {
        creepData->updateYourself();
    }
}

void HellmichMaterialStatus :: printOutputAt(FILE *stream, TimeStep *atTime)
// OOFEM output routine
{
    FloatArray helpVec, fullHelpVec;
    int i, n;
    HellmichMaterial *mat = ( HellmichMaterial * ) ( gp->giveMaterial() );
    MaterialOptions options = giveMaterialOptions();
    ActiveSurface as;
    StructuralCrossSection *cs = ( StructuralCrossSection * ) ( gp->giveCrossSection() );
    // output material gp status if necessary
    if ( options & moIsothermal ) {
        mat->printOutputAt(stream, atTime);
    }

    if ( ( options & moPlotStressStrain ) && !( options & moPlotStressStrainIter ) ) {
        FILE *sstream;
        sstream = fopen(OUTPUTFILE_SS, "a");
        if ( sstream ) {
            mat->plotStressStrain( sstream, gp, atTime, mat->givePssIndex(), gp->giveElement()->giveNumber(), ( double ) ( gp->giveNumber() ) );
            fclose(sstream);
        }
    }

    StructuralMaterialStatus :: printOutputAt(stream, atTime);

    fprintf(stream, "   status HellMat");
    if ( nonisoData ) {
        fprintf( stream, " Temp %.3f", giveTemperature() );
        if ( ( ( options & moHasShrinkage ) == moHasShrinkage ) || ( options & moHumidityStrain ) ) {
            mat->giveShrinkageStrainVector(helpVec, ReducedForm, gp, atTime, VM_Total);
            if ( helpVec.isEmpty() ) {
                helpVec.resize(1);
                helpVec.zero();
            }

            fprintf( stream, " eshr %.5e", helpVec.at(1) );
        }
    }

    if ( plastData ) {
        as = giveActiveSurface();
        if ( !as ) {
            fprintf(stream, " Elastic");
        } else {
            fprintf(stream, " Plast");
            if ( as == asDP ) {
                fprintf(stream, "DP");
            } else if ( as == asTC ) {
                fprintf(stream, "TC");
            } else if ( as == asCorner ) {
                fprintf(stream, "CO");
            }

            fprintf( stream, " chi1 %.5e eps_pl ", giveHardeningVar() );
            givePlasticStrainVector(helpVec);
            cs->giveFullCharacteristicVector(fullHelpVec, gp, helpVec);
            n = fullHelpVec.giveSize();
            for ( i = 0; i < n; i++ ) {
                fprintf( stream, " %.4e", fullHelpVec(i) );
            }
        }
    }

    if ( creepData ) {
        fprintf(stream, "   Creep");
        if ( nonisoData ) {
            fprintf( stream, " slip %.5e", giveViscousSlip() );
        }

        fprintf(stream, " eps_visc");
        cs->giveFullCharacteristicVector(helpVec, gp, creepData->viscousStrainVector);
        n = helpVec.giveSize();
        for ( i = 0; i < n; i++ ) {
            fprintf( stream, " %.4e", helpVec(i) );
        }

        fprintf(stream, " eps_flow ");
        cs->giveFullCharacteristicVector(helpVec, gp, creepData->flowStrainVector);
        for ( i = 0; i < n; i++ ) {
            fprintf( stream, " %.4e", helpVec(i) );
        }
    }

    // output hydration model status, if present
    HydrationModelStatusInterface :: printOutputAt(stream, atTime);

    fprintf(stream, "\n");
}

contextIOResultType
PlastData :: saveContext(DataStream *stream, ContextMode mode)
// Saves the plasticity status data
{
    contextIOResultType iores;
    if ( ( iores = plasticStrainVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->write(& hardeningVar, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    int _as = activeSurface;
    if ( !stream->write(& _as, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
contextIOResultType
CreepData :: saveContext(DataStream *stream, ContextMode mode)
// Saves the plasticity status data
{
    contextIOResultType iores;
    if ( ( iores = viscousStrainVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = flowStrainVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
contextIOResultType
NonisoData :: saveContext(DataStream *stream, ContextMode mode)
// Saves the non-isothermal/mat.level status data
{
    //contextIOResultType iores;
    if ( !stream->write(& initialTemperature, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& temperature, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& humidity, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& gamma0, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& viscousSlip, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& viscosity, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
contextIOResultType
HellmichMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
// saves full information stored in this Status
// no temp variables stored
{
    contextIOResultType iores;
    // Get material options
    MaterialOptions options = giveMaterialOptions();
    // save flag for each module to enable consistency check at restoreContext
    int dataFlag;
    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // Plasticity
    if ( ( ( options & moPlasticity ) == 0 ) ^ ( plastData == NULL ) ) {
        _error("saveContext: plastic status inconsistent with material options.");
    }

    dataFlag = ( plastData ) ? 1 : 0;
    if ( !stream->write(& dataFlag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( dataFlag && ( ( iores = plastData->saveContext(stream, mode) ) != CIO_OK ) ) {
        THROW_CIOERR(iores);
    }

    // Creep
    if ( ( ( options & moCreep ) == 0 ) ^ ( creepData == NULL ) ) {
        _error("saveContext: creep status inconsistent with material options.");
    }

    dataFlag = ( creepData ) ? 1 : 0;
    if ( !stream->write(& dataFlag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( dataFlag && ( ( iores = creepData->saveContext(stream, mode) ) != CIO_OK ) ) {
        THROW_CIOERR(iores);
    }

    // Non-isothermal or material-level auxiliary status
    if ( ( ( options & moIsothermal ) != 0 ) ^ ( nonisoData == NULL ) ) {
        _error("saveContext: non-isothermal status inconsistent with material options.");
    }

    dataFlag = ( nonisoData ) ? 1 : 0;
    if ( !stream->write(& dataFlag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( dataFlag && ( ( iores = nonisoData->saveContext(stream, mode) ) != CIO_OK ) ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

contextIOResultType
PlastData :: restoreContext(DataStream *stream, ContextMode mode)
// Restores the plasticity status data
{
    contextIOResultType iores;
    if ( ( iores = plasticStrainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->read(& hardeningVar, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    int _as;
    if ( !stream->read(& _as, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    activeSurface = ( ActiveSurface ) _as;
    return CIO_OK;
}
contextIOResultType
CreepData :: restoreContext(DataStream *stream, ContextMode mode)
// Restores the plasticity status data
{
    contextIOResultType iores;
    if ( ( iores = viscousStrainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = flowStrainVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
contextIOResultType
NonisoData :: restoreContext(DataStream *stream, ContextMode mode)
// Restores the non-isothermal/mat.level status data
{
    if ( !stream->read(& initialTemperature, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& temperature, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& humidity, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& gamma0, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& viscousSlip, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& viscosity, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
contextIOResultType
HellmichMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
// Restores full information stored in stream to this Status
{
    contextIOResultType iores;
    // read flag for each module to enable consistency check
    int dataFlag;
    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // Plasticity
    if ( !stream->read(& dataFlag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( dataFlag == 0 ) ^ ( plastData == NULL ) ) {
        _error("restoreContext: creep status inconsistent with material options.");
    }

    if ( dataFlag && ( ( iores = plastData->restoreContext(stream, mode) ) != CIO_OK ) ) {
        THROW_CIOERR(iores);
    }

    // Creep
    if ( !stream->read(& dataFlag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( dataFlag == 0 ) ^ ( creepData == NULL ) ) {
        _error("restoreContext: creep status inconsistent with material options.");
    }

    if ( dataFlag && ( ( iores = creepData->restoreContext(stream, mode) ) != CIO_OK ) ) {
        THROW_CIOERR(iores);
    }

    // Non-isothermal or material-level auxiliary status
    if ( !stream->read(& dataFlag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( dataFlag == 0 ) ^ ( nonisoData == NULL ) ) {
        _error("restoreContext: non-isothermal status inconsistent with material options.");
    }

    if ( dataFlag && ( ( iores = nonisoData->restoreContext(stream, mode) ) != CIO_OK ) ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

#else // #ifdef __TM_MODULE
HellmichMaterial :: HellmichMaterial(int n, Domain *d) : StructuralMaterial(n, d) {
    _error("Can't create instance of this class, TM module required");
}
#endif
} // end namespace oofem
