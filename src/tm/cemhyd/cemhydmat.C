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

/* Majority of this file was developed at the National Institute of
 * Standards and Technology by employees of the Federal
 * Government in the course of their official duties. Pursuant
 * to title 17 Section 105 of the United States Code this
 * software is not subject to copyright protection and is in
 * the public domain. CEMHYD3D is an experimental system. NIST
 * assumes no responsibility whatsoever for its use by other
 * parties, and makes no guarantees, expressed or implied,
 * about its quality, reliabi#define OUTFILESlity, or any other characteristic.
 * We would appreciate acknowledgement if the software is used.
 * This software can be redistributed and/or modified freely
 * provided that any derivative works bear some notice that
 * they are derived from it, and any modified versions bear
 * some notice that they have been modified.
 */

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "cemhydmat.h"
#include "homogenize.h"

#ifdef __TM_MODULE //OOFEM transport module
 #include "domain.h"
 #include "flotmtrx.h"
 #include "gausspnt.h"
 #include "oofem_limits.h"
 #ifndef __MAKEDEPEND
  #include <cstdlib>
 #endif
#endif

namespace oofem {
/* This software was developed at the National Institute of */
/* Standards and Technology by employees of the Federal */
/* Government in the course of their official duties. Pursuant */
/* to title 17 Section 105 of the United States Code this */
/* software is not subject to copyright protection and is in */
/* the public domain. CEMHYD3D is an experimental system. NIST */
/* assumes no responsibility whatsoever for its use by other */
/* parties, and makes no guarantees, expressed or implied, */
/* about its quality, reliability, or any other characteristic. */
/* We would appreciate acknowledgement if the software is used. */
/* This software can be redistributed and/or modified freely */
/* provided that any derivative works bear some notice that */
/* they are derived from it, and any modified versions bear */
/* some notice that they have been modified. */
/* Modified 3/97 to allow placement of pozzolanic, inert and fly ash particles */
/* Modified 9/98 to allow placement of various forms of gypsum */
/* Documented version produced 1/00 */
/* Modified by smilauer@cml.fsv.cvut.cz to include 1 voxel particles 16.6.2005
 * Dynamical allocation of memory arrays (possible in input file)
 */

#define OUTFILES //if defined, output files are generated
#define IMAGEFILES //if defined, output percolated and unpercolated images in each cycle (directories perc/ and unperc/)
#define PRINTF //if defined, printf results simultaneously on screen

#ifdef __TM_MODULE //OOFEM transport module
 #undef OUTFILES
 #undef IMAGEFILES
 #undef PRINTF
CemhydMat :: CemhydMat(int n, Domain *d) : IsotropicHeatTransferMaterial(n, d)
{
    MasterCemhydMatStatus = NULL;
}

CemhydMat :: ~CemhydMat()
{
}

//returns hydration power [W/m3 of concrete]
void
CemhydMat :: computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
{
    double averageTemperature;
    CemhydMatStatus *ms = ( CemhydMatStatus * ) this->giveStatus(gp);
    val.resize(1);

    if ( eachGP || ms == MasterCemhydMatStatus ) {
        averageTemperature = ms->giveAverageTemperature();
        if ( mode == VM_Total ) {
            //for nonlinear solver, return the last value even no time has elapsed
            if ( atTime->giveTargetTime() != ms->LastCallTime ) {
                val.at(1) = ms->GivePower( averageTemperature, atTime->giveTargetTime() );
            } else {
                val.at(1) = ms->PartHeat;
            }
        } else {
            OOFEM_ERROR2( "Undefined mode %s\n", __ValueModeTypeToString(mode) );
        }
    } else { //return released heat from the master
        if ( mode == VM_Total ) {
            val.at(1) = MasterCemhydMatStatus->PartHeat;
        } else {
            OOFEM_ERROR2( "Undefined mode %s\n", __ValueModeTypeToString(mode) );
        }
    }

    //val.at(1) = 1500;//constant source
}

void
CemhydMat :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *atTime)
{
    CemhydMatStatus *ms = ( CemhydMatStatus * ) this->giveStatus(gp);
    ms->letTempStateVectorBe(vec);
}


int CemhydMat :: giveCycleNumber(GaussPoint *gp)
{
    CemhydMatStatus *ms = ( CemhydMatStatus * ) this->giveStatus(gp);
    if ( MasterCemhydMatStatus ) {
        ms = MasterCemhydMatStatus;
    }

    return ms->GiveCycNum();
}

double CemhydMat :: giveTimeOfCycle(GaussPoint *gp)
{
    CemhydMatStatus *ms = ( CemhydMatStatus * ) this->giveStatus(gp);
    if ( MasterCemhydMatStatus ) {
        ms = MasterCemhydMatStatus;
    }

    return ms->GiveCycTime();
}



double CemhydMat :: giveDoHActual(GaussPoint *gp)
{
    CemhydMatStatus *ms = ( CemhydMatStatus * ) this->giveStatus(gp);
    if ( MasterCemhydMatStatus ) {
        ms = MasterCemhydMatStatus;
    }

    return ms->GiveDoHActual();
}

//standard units are [Wm-1K-1]
double CemhydMat :: giveConcreteConductivity(GaussPoint *gp)
{
    CemhydMatStatus *ms = ( CemhydMatStatus * ) this->giveStatus(gp);
    double conduct;

    if ( MasterCemhydMatStatus ) {
        ms = MasterCemhydMatStatus;
    }

    if ( conductivityType == 0 ) { //given from OOFEM input file
        conduct = IsotropicHeatTransferMaterial :: give('k', gp);
    } else if ( conductivityType == 1 ) { //compute according to Ruiz, Schindler, Rasmussen. Kim, Chang: Concrete temperature modeling and strength prediction using maturity concepts in the FHWA HIPERPAV software, 7th international conference on concrete pavements, Orlando (FL), USA, 2001
        conduct = IsotropicHeatTransferMaterial :: give('k', gp) * ( 1.33 - 0.33 * ms->GiveDoHActual() );
    } else {
        OOFEM_ERROR2("Unknown conductivityType %d\n", conductivityType);
    }


    //Parallel Voigt model, 20 W/m/K for steel
    conduct = conduct * ( 1. - this->reinforcementDegree ) + 20. * this->reinforcementDegree;

    if ( !this->nowarnings.at(2) && ( conduct < 0.3 || conduct > 5 ) ) {
        OOFEM_WARNING2("Weird concrete thermal conductivity %f W/m/K\n", conduct);
    }

    conduct *= this->scaling.at(2);

    return conduct;
}

//normally it returns J/kg/K of concrete
double CemhydMat :: giveConcreteCapacity(GaussPoint *gp)
{
    CemhydMatStatus *ms = ( CemhydMatStatus * ) this->giveStatus(gp);
    double capacityConcrete;

    if ( MasterCemhydMatStatus ) {
        ms = MasterCemhydMatStatus;
    }

    if ( capacityType == 0 ) { //given from OOFEM input file
        capacityConcrete = IsotropicHeatTransferMaterial :: give('c', gp);
    } else if ( capacityType == 1 ) { //compute from CEMHYD3D according to Bentz
        capacityConcrete = ms->computeConcreteCapacityBentz();
    } else if ( capacityType == 2 ) { //compute from CEMHYD3D directly
        capacityConcrete = 1000 * ms->GiveCp();
    } else {
        OOFEM_ERROR2("Unknown capacityType %d\n", capacityType);
    }

    //Parallel Voigt model, 500 J/kg/K for steel
    capacityConcrete = capacityConcrete * ( 1. - this->reinforcementDegree ) + 500. * this->reinforcementDegree;

    if ( !this->nowarnings.at(3) && ( capacityConcrete < 500 || capacityConcrete > 2000 ) ) {
        OOFEM_WARNING2("Weird concrete heat capacity %f J/kg/K\n", capacityConcrete);
    }

    capacityConcrete *= this->scaling.at(3);

    return capacityConcrete;
}

double CemhydMat :: giveConcreteDensity(GaussPoint *gp)
{
    CemhydMatStatus *ms = ( CemhydMatStatus * ) this->giveStatus(gp);
    double concreteBulkDensity;
    if ( MasterCemhydMatStatus ) {
        ms = MasterCemhydMatStatus;
    }

    if ( densityType == 0 ) { //get from OOFEM input file
        concreteBulkDensity = IsotropicHeatTransferMaterial :: give('d', gp);
    } else if ( densityType == 1 ) { //get from XML input file
        concreteBulkDensity = ms->GiveDensity();
    } else {
        OOFEM_ERROR2("Unknown densityType %d\n", densityType);
    }

    //Parallel Voigt model, 7850 kg/m3 for steel
    concreteBulkDensity = concreteBulkDensity * ( 1. - this->reinforcementDegree ) + 7850. * this->reinforcementDegree;

    if ( !this->nowarnings.at(1) && ( concreteBulkDensity < 1000 || concreteBulkDensity > 4000 ) ) {
        OOFEM_WARNING2("Weird concrete density %f kg/m3\n", concreteBulkDensity);
    }

    concreteBulkDensity *= this->scaling.at(1);

    return concreteBulkDensity;
}

void
CemhydMat :: giveCharacteristicMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MaterialMode mMode = gp->giveMaterialMode();
    double conduct = giveConcreteConductivity(gp);

    switch  ( mMode ) {
    case _1dHeat:
        answer.resize(1, 1);
        answer.at(1, 1) = conduct;
    case _2dHeat:
        answer.resize(2, 2);
        answer.at(1, 1) = conduct;
        answer.at(2, 2) = conduct;
        return;

    case _3dHeat:
        answer.resize(3, 3);
        answer.at(1, 1) = conduct;
        answer.at(2, 2) = conduct;
        answer.at(3, 3) = conduct;
        return;

    default:
        OOFEM_ERROR2( "giveCharacteristicMatrix : unknown mode (%s)\n", __MaterialModeToString(mMode) );
    }
}

double
CemhydMat :: giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    if ( mode == Capacity ) {
        return ( giveConcreteCapacity(gp) * giveConcreteDensity(gp) );
    } else if ( mode == IntSource ) { //for nonlinear solver, return dHeat/dTemperature
        //it suffices to compute derivative of Arrhenius equation
        //double actualTemperature = ((TransportMaterialStatus *) giveStatus(gp))->giveTempStateVector().at(1);
        double lastEquilibratedTemperature = ( ( TransportMaterialStatus * ) giveStatus(gp) )->giveStateVector().at(1);
        //double dt = atTime->giveTimeIncrement();
        double krate, EaOverR, val;
        CemhydMatStatus *ms = ( CemhydMatStatus * ) this->giveStatus(gp);
        if ( MasterCemhydMatStatus ) {
            ms = MasterCemhydMatStatus;
        }

        EaOverR = 1000. * ms->E_act / 8.314;

        if ( ms->icyc > 1 ) {
            krate = exp( -EaOverR * ( 1. / ( ms->temp_cur + 273.15 ) -  1. / ( lastEquilibratedTemperature + 273.15 ) ) );
            //use PartHeat from the last cycle as a corrector tangent, at least one cycle has elapsed
            //         if( fabs(3600*ms->time_cur - ms->PrevHydrTime) > 1.e-3 ){
            //             power = ms->heat_new-ms->PrevCycHeat / (3600*ms->time_cur - ms->PrevHydrTime);//[J/s = W] per gram of cement
            //             power *= 1000 * ms->Mass_cement_concrete;//W/m3/s
            //         } else {
            //             power = 1.e-6;
            //         }
        } else {
            krate = 1.;
        }

        val = EaOverR * krate / ( ms->temp_cur + 273.15 ) / ( ms->temp_cur + 273.15 );

        return val;
    } else {
        OOFEM_ERROR2( "giveCharacteristicValue : unknown mode (%s)\n", __MatResponseModeToString(mode) );
    }

    return 0.;
}


int
CemhydMat :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    CemhydMatStatus *ms;
    if ( MasterCemhydMatStatus ) {
        ms = MasterCemhydMatStatus;
    } else {
        ms = ( CemhydMatStatus * ) this->giveStatus(aGaussPoint);
    }

    if ( type == IST_HydrationDegree ) {
        answer.resize(1);
        answer.at(1) = this->giveDoHActual(aGaussPoint);
        return 1;
    } else if ( type == IST_Density ) {
        answer.resize(1);
        answer.at(1) = this->giveConcreteDensity(aGaussPoint);
        return 1;
    } else if ( type == IST_ThermalConductivityIsotropic ) {
        answer.resize(1);
        answer.at(1) = this->giveConcreteConductivity(aGaussPoint);
        return 1;
    } else if ( type == IST_HeatCapacity ) {
        answer.resize(1);
        answer.at(1) = this->giveConcreteCapacity(aGaussPoint);
        return 1;
    } else if ( type == IST_AverageTemperature ) {
        answer.resize(1);
        answer.at(1) = ms->giveAverageTemperature();
        return 1;
    } else if ( type == IST_YoungModulusVirginPaste ) {
        answer.resize(1);
        answer.at(1) = ms->last_values [ 2 ];
        return 1;
    } else if ( type == IST_PoissonRatioVirginPaste ) {
        answer.resize(1);
        answer.at(1) = ms->last_values [ 3 ];
        return 1;
    } else if ( type == IST_YoungModulusConcrete ) {
        answer.resize(1);
        answer.at(1) = ms->last_values [ 4 ];
        return 1;
    } else if ( type == IST_PoissonRatioConcrete ) {
        answer.resize(1);
        answer.at(1) = ms->last_values [ 5 ];
        return 1;
    } else {
        return TransportMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

int
CemhydMat :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_HydrationDegree ) {
        return 1;
    } else {
        return TransportMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}


InternalStateValueType
CemhydMat :: giveIPValueType(InternalStateType type)
{
    return TransportMaterial :: giveIPValueType(type);
}

int
CemhydMat :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( ( type == IST_HydrationDegree ) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return TransportMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}

int CemhydMat :: initMaterial(Element *element) {
    IntegrationRule *iRule;
    GaussPoint *gp;
    CemhydMatStatus *ms;
    int i;


    iRule = element->giveDefaultIntegrationRulePtr();
    for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        if ( !MasterCemhydMatStatus && !eachGP ) {
            ms = new CemhydMatStatus(1, domain, gp, NULL, this, 1);
            MasterCemhydMatStatus = ms;
        } else if ( eachGP ) {
            ms = new CemhydMatStatus(1, domain, gp, MasterCemhydMatStatus, this, 1);
        } else {
            ms = new CemhydMatStatus(1, domain, gp, NULL, this, 0);
        }

        gp->setMaterialStatus(ms);
    }

    return 1;
}

void CemhydMat :: clearWeightTemperatureProductVolume(Element *element)
{
    IntegrationRule *iRule;
    GaussPoint *gp;
    CemhydMatStatus *ms;
    int i;
    iRule = element->giveDefaultIntegrationRulePtr();

    for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        ms = ( CemhydMatStatus * ) this->giveStatus(gp);
        ms->setAverageTemperatureVolume(0.0, 0.0);
    }
}

void CemhydMat :: storeWeightTemperatureProductVolume(Element *element, TimeStep *tStep)
{
    IntegrationRule *iRule;
    GaussPoint *gp;
    CemhydMatStatus *ms;
    int i;
    double dV;
    FloatArray vecTemperature;
    iRule = element->giveDefaultIntegrationRulePtr();

    if ( !eachGP ) {
        for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule->getIntegrationPoint(i);
            ms = ( CemhydMatStatus * ) this->giveStatus(gp);
            //when more GPs are lumped to a master GP
            dV  = element->computeVolumeAround(gp);
            element->giveIPValue(vecTemperature, gp, IST_Temperature, tStep);
            MasterCemhydMatStatus->setAverageTemperatureVolume(MasterCemhydMatStatus->giveAverageTemperature() + dV * vecTemperature.at(1), MasterCemhydMatStatus->giveTotalVolume() + dV);
        }
    }
}

void CemhydMat :: averageTemperature(void) {
    if ( !eachGP ) {
        MasterCemhydMatStatus->setAverageTemperatureVolume( MasterCemhydMatStatus->giveAverageTemperature() / MasterCemhydMatStatus->giveTotalVolume(), MasterCemhydMatStatus->giveTotalVolume() );
    }
}

IRResultType CemhydMat :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro
    castingTime = 0.;
    int i;

    this->IsotropicHeatTransferMaterial :: initializeFrom(ir); //read d,k,c
    conductivityType = 0;
    capacityType = 0;
    densityType = 0;
    eachGP = 0;
    nowarnings.resize(4);
    nowarnings.zero();
    scaling.resize(3);
    for ( i = 1; i <= scaling.giveSize(); i++ ) {
        scaling.at(i) = 1.;
    }

    reinforcementDegree = 0.;
    //if you want computation of material properties directly from CEMHYD3D, sum up 1 for density, 2 for conductivity, 4 for capacity
    IR_GIVE_OPTIONAL_FIELD(ir, conductivityType, IFT_CemhydMat_conductivitytype, "conductivitytype"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, capacityType, IFT_CemhydMat_capacitytype, "capacitytype"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, densityType, IFT_CemhydMat_densitytype, "densitytype"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, eachGP, IFT_CemhydMat_eachgp, "eachgp"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, nowarnings, IFT_CemhydMat_nowarnings, "nowarnings"); // Macro
    if ( nowarnings.giveSize() != 4 ) {
        OOFEM_ERROR2( "Incorrect size %d of nowarnings", nowarnings.giveSize() );
    }

    IR_GIVE_OPTIONAL_FIELD(ir, scaling, IFT_CemhydMat_scaling, "scaling"); // Macro
    if ( scaling.giveSize() != 3 ) {
        OOFEM_ERROR2( "Incorrect size %d of scaling", nowarnings.giveSize() );
    }

    IR_GIVE_OPTIONAL_FIELD(ir, reinforcementDegree, IFT_CemhydMat_reinforcementDegree, "reinforcementdegree"); // Macro
    IR_GIVE_FIELD(ir, XMLfileName, IFT_CemhydMat_inputFileName, "file");

    return IRRT_OK;
}


MaterialStatus *
CemhydMat :: CreateStatus(GaussPoint *gp) const
{
    OOFEM_ERROR("Use function CemhydMat :: initMaterial instead");

    return NULL;
}


//constructor allowing to copy a microstructure from another CemhydMatStatus
//particular instance of CemhydMat in an integration point
CemhydMatStatus :: CemhydMatStatus(int n, Domain *d, GaussPoint *gp, CemhydMatStatus *CemStat, CemhydMat *cemhydmat, bool withMicrostructure) : TransportMaterialStatus(n, d, gp)
{
    int i, j, k;
    PartHeat = 0.;
    //to be sure, set all pointers to NULL
    mic = NULL;
    mic_CSH = NULL;
    micorig = NULL;
    micpart = NULL;
    mask = NULL;
    ArrPerc = NULL;
    ConnNumbers = NULL;
    cshage = NULL;
    faces = NULL;
    PhaseFrac = NULL;
    last_values = NULL;
    phase = NULL;
    
    adiafile = NULL;
    thfile = NULL;
    elasfile = NULL;
    heatfile = NULL;
    chsfile = NULL;
    ptmpfile = NULL;
    movfile = NULL;
    pHfile = NULL;
    micfile = NULL;
    fileperc = NULL;
    percfile = NULL;
    disprobfile = NULL;
    phasfile = NULL;
    perc_phases = NULL;
    CSHfile = NULL;
    infoperc = NULL;
    infoUnperc = NULL;

    mic = NULL;
    micorig = NULL;
    micpart = NULL;
    cement = NULL;
    cemreal = NULL;
    clust = NULL;
    mask = NULL;
    curvature = NULL;
    mic_CSH = NULL;
    ArrPerc = NULL;
    ConnNumbers = NULL;
    faces = NULL;

#ifdef TINYXML
    xmlFile = NULL;
#endif

    headant = NULL;
    tailant = NULL;
    soluble = NULL;
    creates = NULL;
    CSH_vicinity = NULL;
    molarvcsh = NULL;
    watercsh = NULL;
    xsph = NULL;
    ysph = NULL;
    zsph = NULL;
    iv = NULL;
    discount = NULL;
    count = NULL;
    disprob = NULL;
    disbase = NULL;
    specgrav = NULL;
    molarv = NULL;
    heatf = NULL;
    waterc = NULL;
    pHeffect = NULL;
    soluble = NULL;
    creates = NULL;

    
    //set common variables in constructor
 #ifdef PRINTF
    printf("Constructor of CemhydMatStatus on GP %p, withMicrostructure %d, copy from CemhydMatStatus %p\n", gp, withMicrostructure, CemStat);
    fflush(stdout);
 #endif
    this->gp = gp;
    if ( withMicrostructure ) {
        this->initializeMicrostructure();
        if ( !CemStat ) {
            this->readInputFileAndInitialize(cemhydmat->XMLfileName.c_str(), 1);
        } else { //copy 3D microstructure
            this->readInputFileAndInitialize(cemhydmat->XMLfileName.c_str(), 0); //read input but do not reconstruct 3D microstructure
            for ( k = 0; k < SYSIZE; k++ ) {
                for ( j = 0; j < SYSIZE; j++ ) {
                    for ( i = 0; i < SYSIZE; i++ ) {
                        micpart [ i ] [ j ] [ k ] = CemStat->micpart [ i ] [ j ] [ k ];
                        micorig [ i ] [ j ] [ k ] = CemStat->micorig [ i ] [ j ] [ k ];
                        mic [ i ] [ j ] [ k ] = micorig [ i ] [ j ] [ k ];
                    }
                }
            }
        }
    }
}
#endif //__TM_MODULE


void CemhydMatStatus :: initializeMicrostructure()
{
    icyc = 1; //set the cycle counter
    time_cur = 0.; //hydration time [h]
    heat_new = 0.;
    heat_cf = 0.;
    Cp_now = 0.900; //initial guess [J/g of concrete], which is later updated

    averageTemperature = 0.;
    IPVolume = 0.;
    init_material_time = 0.;

    xoff [ 0 ] = 1;
    xoff [ 1 ] = 0;
    xoff [ 2 ] = 0;
    xoff [ 3 ] = -1;
    xoff [ 4 ] = 0;
    xoff [ 5 ] = 0;
    xoff [ 6 ] = 1;
    xoff [ 7 ] = 1;
    xoff [ 8 ] = -1;
    xoff [ 9 ] = -1;
    xoff [ 10 ] = 0;
    xoff [ 11 ] = 0;
    xoff [ 12 ] = 0;
    xoff [ 13 ] = 0;
    xoff [ 14 ] = 1;
    xoff [ 15 ] = 1;
    xoff [ 16 ] = -1;
    xoff [ 17 ] = -1;
    xoff [ 18 ] = 1;
    xoff [ 19 ] = 1;
    xoff [ 20 ] = 1;
    xoff [ 21 ] = 1;
    xoff [ 22 ] = -1;
    xoff [ 23 ] = -1;
    xoff [ 24 ] = -1;
    xoff [ 25 ] = -1;
    xoff [ 26 ] = 0;

    yoff [ 0 ] = 0;
    yoff [ 1 ] = 1;
    yoff [ 2 ] = 0;
    yoff [ 3 ] = 0;
    yoff [ 4 ] = -1;
    yoff [ 5 ] = 0;
    yoff [ 6 ] = 1;
    yoff [ 7 ] = -1;
    yoff [ 8 ] = 1;
    yoff [ 9 ] = -1;
    yoff [ 10 ] = 1;
    yoff [ 11 ] = -1;
    yoff [ 12 ] = 1;
    yoff [ 13 ] = -1;
    yoff [ 14 ] = 0;
    yoff [ 15 ] = 0;
    yoff [ 16 ] = 0;
    yoff [ 17 ] = 0;
    yoff [ 18 ] = 1;
    yoff [ 19 ] = -1;
    yoff [ 20 ] = 1;
    yoff [ 21 ] = -1;
    yoff [ 22 ] = 1;
    yoff [ 23 ] = 1;
    yoff [ 24 ] = -1;
    yoff [ 25 ] = -1;
    yoff [ 26 ] = 0;

    zoff [ 0 ] = 0;
    zoff [ 1 ] = 0;
    zoff [ 2 ] = 1;
    zoff [ 3 ] = 0;
    zoff [ 4 ] = 0;
    zoff [ 5 ] = -1;
    zoff [ 6 ] = 0;
    zoff [ 7 ] = 0;
    zoff [ 8 ] = 0;
    zoff [ 9 ] = 0;
    zoff [ 10 ] = 1;
    zoff [ 11 ] = 1;
    zoff [ 12 ] = -1;
    zoff [ 13 ] = -1;
    zoff [ 14 ] = 1;
    zoff [ 15 ] = -1;
    zoff [ 16 ] = 1;
    zoff [ 17 ] = -1;
    zoff [ 18 ] = 1;
    zoff [ 19 ] = 1;
    zoff [ 20 ] = -1;
    zoff [ 21 ] = -1;
    zoff [ 22 ] = 1;
    zoff [ 23 ] = -1;
    zoff [ 24 ] = 1;
    zoff [ 25 ] = -1;
    zoff [ 26 ] = 0;

    iy = 0; //random generator ran1()
    LastCallTime = -1.e6; //start from begining (the LastCall in the beginning is -1.e6)
    alpha_cur = 0.0;
    alpha_last = 0.; //last degree of hydration
    LastTargTime = 0.; //stores last call of TargTime
    TargDoHelas = 0.; //stores DoH at which to perform analytic homogenization
    //definitions originally from CemhydMat.h

    //for generation of microstructures
    CEM = 100;       /* and greater */
    CEMID = 1;       /* phase identifier for cement */
    C2SID = 2;       /* phase identified for C2S cement */
    GYPID = 5;       /* phase identifier for gypsum */
    HEMIHYDRATE = 6; /* phase identifier for hemihydrate */
    POZZID = 8;      /* phase identifier for pozzolanic material - REACTIVE */
    INERTID = 9;     /* phase identifier for inert material - UNREACTIVE */
    SLAGID = 10;     /* phase identifier for slag - REACTIVE */
    AGG = 28;        /* phase identifier for flat aggregate - UNREACTIVE */
    FLYASH = 30;     /* phase identifier for all fly ash components - UNREACTIVE*/

    //for hydration part
    POROSITY = 0;
    C3S = 1;
    C2S = 2;
    C3A = 3;
    C4AF = 4;
    GYPSUM = 5;
    HEMIHYD = 6;
    ANHYDRITE = 7;
    POZZ = 8;
    INERT = 9;
    SLAG = 10;
    ASG = 11; /* aluminosilicate glass */
    CAS2 = 12;
    CH = 13;
    CSH = 14;
    C3AH6 = 15;
    ETTR = 16;
    ETTRC4AF = 17; /* Iron-rich stable ettringite phase */
    AFM = 18;
    FH3 = 19;
    POZZCSH = 20;
    SLAGCSH = 21; /* Slag gel-hydration product */
    CACL2 = 22;
    FREIDEL = 23; /* Freidel's salt */
    STRAT = 24; /* stratlingite (C2ASH8) */
    GYPSUMS = 25; /* Gypsum formed from hemihydrate and anhydrite */
    CACO3 = 26;
    AFMC = 27;
    INERTAGG = 28;
    ABSGYP = 29;
    DIFFCSH = 30;
    DIFFCH = 31;
    DIFFGYP = 32;
    DIFFC3A = 33;
    DIFFC4A = 34;
    DIFFFH3 = 35;
    DIFFETTR = 36;
    DIFFCACO3 = 37;
    DIFFAS = 38;
    DIFFANH = 39;
    DIFFHEM = 40;
    DIFFCAS2 = 41;
    DIFFCACL2 = 42;
    EMPTYP = 45; /*Empty porosity due to self desiccation*/
    HDCSH = 46;
    OFFSET = 50; /*Offset for highlighted potentially soluble pixel*/

    NEIGHBORS = 26; /* number of neighbors to consider (6, 18, or 26) in dissolution */
    BoxSize = 1; /*int describing vicinity of CSH*/
    SolidLimit = 27; /*how many solid phase voxels must be in a box (max. <=(2*BoxSize+1)^3)*/
    MAXTRIES = 150000; /* maximum number of random tries for sphere placement */
    MAXCYC_SEAL = 30000; /* Maximum number of cycles of sealed hydration (originally MAXCYC in disrealnew.c */
    NUMSIZES = 100;   /* maximum number of different particle sizes */
    MAXSPH = 10000; /* maximum number of elements in a spherical template */

    Cp_pozz = 0.75;
    Cp_CH = 0.75;
    Cp_h2o = 4.18; /* Cp for free water */
    Cp_bh2o = 2.2; /* Cp for bound water */
    WN = 0.23;    /* water bound per gram of cement during hydration */
    WCHSH = 0.06; /* water imbibed per gram of cement during chemical shrinkage (estimate) */
    CUBEMIN = 3;    /* Minimum cube size for checking pore size */

    DISBIAS = 30.0; /* Dissolution bias- to change all dissolution rates */
    DISMIN = 0.001; /* Minimum dissolution for C3S dissolution */
    DISMIN2 = 0.00025; /* Minimum dissolution for C2S dissolution */
    DISMINSLAG = 0.0001; /* Minimum dissolution for SLAG dissolution */
    DISMINASG = 0.0005; /* Minimum dissolution for ASG dissolution */
    DISMINCAS2 = 0.0005; /* Minimum dissolution for CAS2 dissolution */
    DISMIN_C3A_0 = 0.002; /* Minimum dissolution for C3A dissolution */
    DISMIN_C4AF_0 = 0.0005; /* Minimum dissolution for C4AF dissolution */

    C3AH6GROW = 0.01; /* Probability for C3AH6 growth */
    CHGROW = 1.0;    /* Probability for CH growth */
    CHGROWAGG = 1.0;  /* Probability for CH growth on aggregate surface */
    ETTRGROW = 0.002; /* Probability for ettringite growth */
    C3AETTR = 0.001;  /* Probability for reaction of diffusing C3A with ettringite  */
    C3AGYP = 0.001;    /*  Probability for diffusing C3A to react with diffusing gypsum */
    SOLIDC3AGYP = 0.5;  /* Probability of solid C3A to react with diffusing sulfate */
    SOLIDC4AFGYP = 0.1; /* Probability of solid C4AF to react with diffusing sulfate */
    PPOZZ = 0.05;    /* base probability for pozzolanic reaction */
    PCSH2CSH = 0.002; /* probability for CSH dissolution */
    A0_CHSOL = 1.325; /* Parameters for variation of CH solubility with */
    A1_CHSOL = 0.008162; /* temperature (data from Taylor- Cement Chemistry) */

    BURNT = 70;    /* label for a burnt pixel <255 (in char type arrays) */
    SIZE2D = 49000; /* size of matrices for holding burning locations */

    SIZESET = 100000;
    AGRATE = 0.25;      /* Probability of gypsum absorption by CSH */
    VOLFACTOR = 0.00001; /* dm per pixel Note- dm*dm*dm = Liters */
    MASSFACTOR = 0.0001; /* cm per pixel - specific gravities in g/cm^3 */
    MMNa = 22.9898;
    MMK = 39.102;
    MMNa2O = 61.979;
    MMK2O = 94.203;
    BNa = 0.00031;   /* From Taylor paper in liters (31 mL/1000/ 100 g) */
    BK = 0.00020;   /* From Taylor paper in liters (20 mL/1000/ 100 g) */
    BprimeNa = 0.0030; /* From Taylor paper in liters (3 mL/1000/ 1 g POZZ) */
    BprimeK = 0.0033; /* From Taylor paper in liters (3.3 mL/1000/ 1 g POZZ) */
    KspCH25C = 0.00000646;
    KspGypsum = 0.0000263;
    KspSyngenite = 0.00000010;
    SpecgravSyngenite = 2.607; /* Source Taylor, H.F.W., Cement Chemistry */
    KperSyn = 2.0; /* moles of K+ per mole of syngenite */
    activeA0 = 0.0366;   /* A at 295 K (from Ken) */
    activeB0 = 0.01035; /* B at 295 K (from Ken) */
    zCa = 2.;
    zSO4 = 2.;
    zOH = 1.;
    zNa = 1.;
    zK = 1.;
    aK = 1.33;
    aCa = 1.;
    aOH = 3.;
    aNa = 3.;
    aSO4 = 4.5;   /* Estimate as S ionic radii + O ionic diameter */
    lambdaOH_0 = 198.0; /* Units: S cm-cm eq.^(-1) */
    lambdaNa_0 = 50.1;
    lambdaK_0 = 73.5;
    lambdaSO4_0 = 39.5;
    lambdaCa_0 = 29.5;   /* Note that CRC has 60./2 for this */
    GOH = 0.353; /* Units: (eq.^2 mol/L)^(-0.5) */
    GK = 0.548;
    GNa = 0.733;
    GCa = 0.771;
    GSO4 = 0.877;
    cm2perL2m = 0.1; /* Conversion from cm2/Liter to 1/m */
    EPSS = 6.e-8;
    MAXIT = 100;
    EPSP = 2.0e-6;
    MAXM = 100;

    /*random generator*/
    IA = 16807;
    IM = 2147483647;
    IQ = 127773;
    IR = 2836;
    NTAB = 32;
    EPS = 1.2E-07;

    PhaseFrac = new double [ 34 ];
    last_values = new double [ 6 ];
    last_values [ 2 ] = 0.001; //Young's modulus of virgin paste
    last_values [ 3 ] = 0.499924; //Poisson's ratio of virgin paste
    last_values [ 4 ] = 0.001; //Young's modulus of concrete
    last_values [ 5 ] = 0.499924; //Poisson's ratio of concrete
    phase = new long int [ 51 ];

    CSH_vicinity = new unsigned int [ ( 2 * BoxSize + 1 ) * ( 2 * BoxSize + 1 ) * ( 2 * BoxSize + 1 ) + 1 ];
    molarvcsh = new float [ MAXCYC_SEAL ];
    watercsh = new float [ MAXCYC_SEAL ];
    xsph = new int [ MAXSPH ];
    ysph = new int [ MAXSPH ];
    zsph = new int [ MAXSPH ];

    iv = new int [ NTAB ];

    discount = new long int [ EMPTYP + 1 ];
    count = new long int [ HDCSH + 1 ];
    disprob = new float [  HDCSH + 1 ];
    disbase = new float [ EMPTYP + 1 ];
    specgrav = new float [ EMPTYP + 1 ];
    molarv = new float [ EMPTYP + 1 ];
    heatf = new float [ EMPTYP + 1 ];
    waterc = new float [ EMPTYP + 1 ];
    pHeffect = new float [ EMPTYP + 1 ];

    //set zero to arrays
    for ( int i = 0; i <= EMPTYP; i++ ) {
        heatf [ i ] = 0.;
        waterc [ i ] = 0.;
    }

    for ( int i = 0; i <= HDCSH; i++ ) {
        disprob [ i ] = 0.;
    }

    soluble = new int [ EMPTYP + 1 ];
    creates = new int [ EMPTYP + 1 ];
}



///destructor
CemhydMatStatus :: ~CemhydMatStatus(void) {
#ifdef OUTFILES
    if ( fileperc != NULL ) {
        fclose(fileperc);
    }

    if ( disprobfile != NULL ) {
        fclose(disprobfile);
    }

    //fclose(percfile);
    if ( heatfile != NULL ) {
        fclose(heatfile);
    }

    if ( pHfile != NULL ) {
        fclose(pHfile);
    }

    if ( phasfile != NULL ) {
        fclose(phasfile);
    }

    if ( perc_phases != NULL ) {
        fclose(perc_phases);
    }

    if ( adiafile != NULL ) {
        fclose(adiafile);
    }

    if ( elasfile != NULL ) {
        fclose(elasfile);
    }

    if ( CSHfile != NULL ) {
        fclose(CSHfile);
    }

    if ( infoperc != NULL ) {
        fclose(infoperc);
    }

    if ( infoUnperc != NULL ) {
        fclose(infoUnperc);
    }

#endif
    delete [] PhaseFrac;
    delete [] last_values;
    delete [] phase;

#ifdef CMLFILE
    delete F; //delete cmlfile
#endif
#ifdef TINYXML
    if ( xmlFile != NULL ) {
        delete xmlFile;
    }
#endif

    delete [] CSH_vicinity;
    delete [] molarvcsh;
    delete [] watercsh;

    delete [] xsph;
    delete [] ysph;
    delete [] zsph;
    delete [] iv;
    delete [] discount;
    delete [] count;
    delete [] disprob;
    delete [] disbase;
    delete [] specgrav;
    delete [] molarv;
    delete [] heatf;
    delete [] waterc;
    delete [] pHeffect;
    delete [] soluble;
    delete [] creates;

    /* Eliminate the whole list */
    struct ants *curant;

    if ( headant != NULL ) { //if hydration did not start sucessfully
        while ( headant->nextant != NULL ) {
            curant = headant->nextant;
            free(headant);
            headant = curant;
#ifdef PRINTF
            printf("Deallocating headant\n");
#endif
        }
    }

    if ( tailant != NULL ) { //if hydration did not start sucessfully
        free(tailant);
    }

#ifdef PRINTF
    printf("Deallocating arrays\n");
    fflush(stdout);
#endif

    dealloc_char_3D(mic, SYSIZE);
    dealloc_int_3D(mic_CSH, SYSIZE);
    dealloc_char_3D(micorig, SYSIZE);
    dealloc_long_3D(micpart, SYSIZE);
    dealloc_int_3D(mask, SYSIZE + 1);
    dealloc_int_3D(ArrPerc, SYSIZE);
    dealloc_int_3D(ConnNumbers, SYSIZE);
    dealloc_shortint_3D(cshage, SYSIZE);
    dealloc_shortint_3D(faces, SYSIZE);
}

void CemhydMatStatus :: alloc_char_3D(char ***( & mic ), long SYSIZE) {
    mic = new char ** [ SYSIZE ];
    for ( int x = 0; x < SYSIZE; x++ ) {
        mic [ x ] = new char * [ SYSIZE ];
        for ( int y = 0; y < SYSIZE; y++ ) {
            mic [ x ] [ y ] = new char [ SYSIZE ];
            if ( mic [ x ] [ y ] == NULL ) {
                printf("Cannot allocate memory (file %s, line %d)\n", __FILE__, __LINE__);
            }
        }
    }
}

void CemhydMatStatus :: dealloc_char_3D(char ***( & mic ), long SYSIZE) {
    if ( mic != NULL ) {
        for ( int x = 0; x < SYSIZE; x++ ) {
            for ( int y = 0; y < SYSIZE; y++ ) {
                delete [] mic [ x ] [ y ];
            }

            delete [] mic [ x ];
        }

        delete [] mic;
    }
}

void CemhydMatStatus :: alloc_long_3D(long ***( & mic ), long SYSIZE) {
    mic = new long ** [ SYSIZE ];
    for ( int x = 0; x < SYSIZE; x++ ) {
        mic [ x ] = new long * [ SYSIZE ];
        for ( int y = 0; y < SYSIZE; y++ ) {
            mic [ x ] [ y ] = new long [ SYSIZE ];
            if ( mic [ x ] [ y ] == NULL ) {
                printf("Cannot allocate memory (file %s, line %d)\n", __FILE__, __LINE__);
            }
        }
    }
}


void CemhydMatStatus :: dealloc_long_3D(long ***( & mic ), long SYSIZE) {
    if ( mic != NULL ) {
        for ( int x = 0; x < SYSIZE; x++ ) {
            for ( int y = 0; y < SYSIZE; y++ ) {
                delete [] mic [ x ] [ y ];
            }

            delete [] mic [ x ];
        }

        delete [] mic;
    }
}

void CemhydMatStatus :: alloc_int_3D(int ***( & mic ), long SYSIZE) {
    mic = new int ** [ SYSIZE ];
    for ( int x = 0; x < SYSIZE; x++ ) {
        mic [ x ] = new int * [ SYSIZE ];
        for ( int y = 0; y < SYSIZE; y++ ) {
            mic [ x ] [ y ] = new int [ SYSIZE ];
            if ( mic [ x ] [ y ] == NULL ) {
                printf("Cannot allocate memory (file %s, line %d)\n", __FILE__, __LINE__);
            }
        }
    }
}


void CemhydMatStatus :: dealloc_int_3D(int ***( & mic ), long SYSIZE) {
    if ( mic != NULL ) {
        for ( int x = 0; x < SYSIZE; x++ ) {
            for ( int y = 0; y < SYSIZE; y++ ) {
                delete [] mic [ x ] [ y ];
            }

            delete [] mic [ x ];
        }

        delete [] mic;
    }
}

void CemhydMatStatus :: alloc_shortint_3D(short int ***( & mic ), long SYSIZE) {
    mic = new short int ** [ SYSIZE ];
    for ( int x = 0; x < SYSIZE; x++ ) {
        mic [ x ] = new short int * [ SYSIZE ];
        for ( int y = 0; y < SYSIZE; y++ ) {
            mic [ x ] [ y ] = new short int [ SYSIZE ];
            if ( mic [ x ] [ y ] == NULL ) {
                printf("Cannot allocate memory (file %s, line %d)\n", __FILE__, __LINE__);
            }
        }
    }
}


void CemhydMatStatus :: dealloc_shortint_3D(short int ***( & mic ), long SYSIZE) {
    if ( mic != NULL ) {
        for ( int x = 0; x < SYSIZE; x++ ) {
            for ( int y = 0; y < SYSIZE; y++ ) {
                delete [] mic [ x ] [ y ];
            }

            delete [] mic [ x ];
        }

        delete [] mic;
    }
}

void CemhydMatStatus :: alloc_double_3D(double ***( & mic ), long SYSIZE) {
    mic = new double ** [ SYSIZE ];
    for ( int x = 0; x < SYSIZE; x++ ) {
        mic [ x ] = new double * [ SYSIZE ];
        for ( int y = 0; y < SYSIZE; y++ ) {
            mic [ x ] [ y ] = new double [ SYSIZE ];
            if ( mic [ x ] [ y ] == NULL ) {
                printf("Cannot allocate memory (file %s, line %d)\n", __FILE__, __LINE__);
            }
        }
    }
}


void CemhydMatStatus :: dealloc_double_3D(double ***( & mic ), long SYSIZE) {
    if ( mic != NULL ) {
        for ( int x = 0; x < SYSIZE; x++ ) {
            for ( int y = 0; y < SYSIZE; y++ ) {
                delete [] mic [ x ] [ y ];
            }

            delete [] mic [ x ];
        }

        delete [] mic;
    }
}

#ifdef TINYXML
//functions to read int, double and string with error checking
void CemhydMatStatus :: QueryNumAttributeExt(TiXmlDocument *xmlFile, const char *elementName, int position, int &val) {
    int success;
    char key [ 256 ];
    TiXmlHandle docHandle = TiXmlHandle(xmlFile);
    TiXmlElement *elemSelected = docHandle.FirstChild(elementName).ToElement();
    if ( elemSelected == NULL ) {
        printf("Cannot find entry %s, terminating, file %s, line %d\n", elementName, __FILE__, __LINE__);
        exit(0);
    }

    sprintf(key, "key%d", position);
    success = elemSelected->QueryIntAttribute(key, & val);
    if ( success != TIXML_SUCCESS ) {
        printf("Cannot read int value or attribute %s from the entry %s, terminating, file %s, line %d\n", key, elementName, __FILE__, __LINE__);
        exit(0);
    }
}

void CemhydMatStatus :: QueryNumAttributeExt(TiXmlDocument *xmlFile, const char *elementName, int position, long int &val) {
    int temp;
    QueryNumAttributeExt(xmlFile, elementName, position, temp);
    val = static_cast< long int >(temp);
}

void CemhydMatStatus :: QueryNumAttributeExt(TiXmlDocument *xmlFile, const char *elementName, const char *key, int &val) {
    int success;
    TiXmlHandle docHandle = TiXmlHandle(xmlFile);
    TiXmlElement *elemSelected = docHandle.FirstChild(elementName).ToElement();
    if ( elemSelected == NULL ) {
        printf("Cannot find entry %s, terminating, file %s, line %d\n", elementName, __FILE__, __LINE__);
        exit(0);
    }

    success = elemSelected->QueryIntAttribute(key, & val);
    if ( success != TIXML_SUCCESS ) {
        printf("Cannot read int value or attribute %s from the entry %s, terminating, file %s, line %d\n", key, elementName, __FILE__, __LINE__);
        exit(0);
    }
}


void CemhydMatStatus :: QueryNumAttributeExt(TiXmlDocument *xmlFile, const char *elementName, int position, double &val) {
    int success;
    char key [ 256 ];
    TiXmlHandle docHandle = TiXmlHandle(xmlFile);
    TiXmlElement *elemSelected = docHandle.FirstChild(elementName).ToElement();
    if ( elemSelected == NULL ) {
        printf("Cannot find entry %s, terminating, file %s, line %d\n", elementName, __FILE__, __LINE__);
        exit(0);
    }

    sprintf(key, "key%d", position);
    success = elemSelected->QueryDoubleAttribute(key, & val);
    if ( success != TIXML_SUCCESS ) {
        printf("Cannot read double value or attribute %s from the entry %s, terminating, file %s, line %d\n", key, elementName, __FILE__, __LINE__);
        exit(0);
    }
}

void CemhydMatStatus :: QueryNumAttributeExt(TiXmlDocument *xmlFile, const char *elementName, const char *key, double &val) {
    int success;
    TiXmlHandle docHandle = TiXmlHandle(xmlFile);
    TiXmlElement *elemSelected = docHandle.FirstChild(elementName).ToElement();
    if ( elemSelected == NULL ) {
        printf("Cannot find entry %s, terminating, file %s, line %d\n", elementName, __FILE__, __LINE__);
        exit(0);
    }

    success = elemSelected->QueryDoubleAttribute(key, & val);
    if ( success != TIXML_SUCCESS ) {
        printf("Cannot read double value or attribute %s from the entry %s, terminating, file %s, line %d\n", key, elementName, __FILE__, __LINE__);
        exit(0);
    }
}

void CemhydMatStatus :: QueryStringAttributeExt(TiXmlDocument *xmlFile, const char *elementName, int position, char *chars) {
    int success;
    char key [ 256 ];
    TiXmlHandle docHandle = TiXmlHandle(xmlFile);
    std :: string str1;
    TiXmlElement *elemSelected = docHandle.FirstChild(elementName).ToElement();
    if ( elemSelected == NULL ) {
        printf("Cannot find entry %s, terminating, file %s, line %d\n", elementName, __FILE__, __LINE__);
        exit(0);
    }

    sprintf(key, "key%d", position);
    //success = elemSelected->QueryStringAttribute(key, & str1);
    // Since ubuntu/debian is still stuck at 2.5.3, lacking QueryStringAttribute.
    // Change with above whenever packages are updated.
    const char *cstr = elemSelected->Attribute(key);
    if (cstr) {
        success = TIXML_SUCCESS;
    } else {
        success = TIXML_NO_ATTRIBUTE;
    }
    if ( success != TIXML_SUCCESS ) {
        printf("Cannot read string value or key %s from the entry %s, terminating, file %s, line %d\n", key, elementName, __FILE__, __LINE__);
        exit(0);
    }
    str1 = std :: string(cstr);

    strcpy( chars, str1.c_str() );
}
#endif

/* read input parameters in file, use XML or cmlfile construction
 * allocate necessary arrays (especially those dependent on SYSIZE)
 * returns (1) if generation of particles in the Representative Volume Element (RVE) was unsuccessful
 */
int CemhydMatStatus :: readInputFileAndInitialize(const char *inp, bool generateMicrostructure)
{
    int read_micr;
#ifdef CMLFILE
    F = new cmlfile(inp);
    // set number of keywords
    F->set_labels(54);
    F->set_label_string(0, "Rand_seed_num");
    F->set_label_string(1, "Input_img_file");
    F->set_label_string(2, "Input_id_file");
    F->set_label_string(3, "Saturated_sealed");
    F->set_label_string(4, "Induction_time");
    F->set_label_string(5, "Ea_cement");
    F->set_label_string(6, "Ea_pozz");
    F->set_label_string(7, "Ea_slag");
    F->set_label_string(8, "Beta");
    F->set_label_string(9, "Mass_SCM_FA_CA_inert_frac");
    F->set_label_string(10, "Mass_cem");
    F->set_label_string(11, "Cp_SCM_FA_CA_inert");
    F->set_label_string(12, "Cp_cem");
    F->set_label_string(13, "Given_microstructure");
    F->set_label_string(14, "Output_initial_microstructure");
    F->set_label_string(15, "Output_initial_microstructure_img_file");
    F->set_label_string(16, "Output_initial_microstructure_id_file");
    F->set_label_string(17, "Cycle_freq_perc_pore");
    F->set_label_string(18, "Cycle_freq_perc_sol");
    F->set_label_string(19, "Total_sodium");
    F->set_label_string(20, "Total_potassium");
    F->set_label_string(21, "Readily_soluble_sodium");
    F->set_label_string(22, "Readily_soluble_potassium");
    F->set_label_string(23, "Diffusion_steps_per_cycle");
    F->set_label_string(24, "CH_nucleation_probability");
    F->set_label_string(25, "CH_scale_factor");
    F->set_label_string(26, "Gypsum_nucleation_probability");
    F->set_label_string(27, "Gypsum_scale_factor");
    F->set_label_string(28, "C3AH6_nucleation_probability");
    F->set_label_string(29, "C3AH6_scale_factor");
    F->set_label_string(30, "FH3_nucleation_probability");
    F->set_label_string(31, "FH3_scale_factor");
    F->set_label_string(32, "Microstructure_size");
    F->set_label_string(33, "Adiabatic_conditions");
    F->set_label_string(34, "Vol_cement_clinker_gypsum");
    F->set_label_string(35, "Vol_cement_SCM");
    F->set_label_string(36, "Vol_water");
    F->set_label_string(37, "Vol_FA");
    F->set_label_string(38, "Vol_CA");
    F->set_label_string(39, "Vol_inert_filler");
    F->set_label_string(40, "Vol_entrained_entrapped_air");
    F->set_label_string(41, "Grain_average_FA");
    F->set_label_string(42, "Grain_average_CA");
    F->set_label_string(43, "ITZ_thickness");
    F->set_label_string(44, "ITZ_Young_red");
    F->set_label_string(45, "Young_SCM");
    F->set_label_string(46, "Poisson_SCM");
    F->set_label_string(47, "Young_FA");
    F->set_label_string(48, "Poisson_FA");
    F->set_label_string(49, "Young_CA");
    F->set_label_string(50, "Poisson_CA");
    F->set_label_string(51, "Young_inert");
    F->set_label_string(52, "Poisson_inert");
    F->set_label_string(53, "Calculate_elastic_homogenization");

    // these keywords with #id will be required
    F->require(0);
    //F->require( 1 ) ;
    //F->require( 2 ) ;
    F->require(3);
    F->require(4);
    F->require(5);
    F->require(6);
    F->require(7);
    F->require(8);
    F->require(9);
    F->require(10);
    F->require(11);
    F->require(12);
    F->require(13);
    F->require(14);
    F->require(17);
    F->require(18);
    F->require(19);
    F->require(20);
    F->require(21);
    F->require(22);
    F->require(23);
    F->require(24);
    F->require(27);
    F->require(28);
    F->require(29);
    F->require(30);
    F->require(31);
    F->require(32);
    F->require(33);
    F->require(34);
    F->require(35);
    F->require(36);
    F->require(37);
    F->require(38);
    F->require(39);
    F->require(40);
    F->require(41);
    F->require(42);
    F->require(43);
    F->require(44);
    F->require(45);
    F->require(46);
    F->require(47);
    F->require(48);
    F->require(49);
    F->require(50);
    F->require(51);
    F->require(52);
    F->require(53);
    // set number and names of sections
    F->set_sections(1);
    F->set_section_string(0, "CEMHYD_generate_particles");

    F->check_requirements();

    if ( F->error_in_requirements() ) {
        printf("Cemhyd input file %s is not complete (file %s, line %d)\n", inp, __FILE__, __LINE__);
        exit(0);
    }

    F->get_value(0, ( long & )iseed);
#endif
#ifdef TINYXML
    xmlFile = new TiXmlDocument(inp);
    countKey = 0;
    if ( !xmlFile->LoadFile() ) {
        printf("\nError reading XML file %s or nonletter symbols used, such as =, (file %s, line %d)\n", inp, __FILE__, __LINE__);
        exit(0);
    }

    QueryNumAttributeExt(xmlFile, "Rand_seed_num", 0, iseed);
    QueryNumAttributeExt(xmlFile, "Microstructure_size", 0, SYSSIZE);
    QueryNumAttributeExt(xmlFile, "Given_microstructure", 0, read_micr);
#endif

    nseed = iseed;
    seed = ( & nseed );
    //printf("iseed = %d", iseed);

    //read SYSSIZE of microstructure and allocate arrays
#ifdef CMLFILE
    F->get_value(32, ( long & )SYSSIZE);
#endif
    if ( SYSSIZE < 10 ) {
        printf("Can not run small microstructure %d (< 10 voxels a side), file %s, line %d)\n", SYSSIZE, __FILE__, __LINE__);
        exit(0);
    }

    SYSIZE = SYSSIZE;
    SYSIZEM1 = ( SYSIZE - 1 ); /* System size -1 */
    SYSIZE_POW3 = ( SYSIZE * SYSIZE * SYSIZE );
    NPARTC = ( long ) ( 700000 * ( double ) SYSIZE_POW3 / 1000000. );
    BURNTG = ( NPARTC > 100 ? NPARTC : 100 );
    CUBEMAX = ( SYSIZE > 6 ? 7 : 3 );
    DETTRMAX = ( 1200. * SYSIZE_POW3 / 1000000. ); /* Maximum allowed # of ettringite diffusing species */
    DGYPMAX = ( 2000. * SYSIZE_POW3 / 1000000. ); /*  Maximum allowed # of gypsum diffusing species */
    DCACO3MAX = ( 1000 * SYSIZE_POW3 / 1000000. ); /*  Maximum allowed # of CaCO3 diffusing species */
    DCACL2MAX = ( 2000 * SYSIZE_POW3 / 1000000. ); /* Maximum allowed # of CaCl2 diffusing species */
    DCAS2MAX = ( 2000 * SYSIZE_POW3 / 1000000. ); /* Maximum allowed # of CAS2 diffusing species */
    CHCRIT = ( 50. * SYSIZE_POW3 / 1000000. ); /* Scale parameter to adjust CH dissolution probability */
    C3AH6CRIT = ( 10. * SYSIZE_POW3 / 1000000. ); /* Scale par. to adjust C3AH6 dissolution prob. */
    CSHSCALE = ( 70000. * SYSIZE_POW3 / 1000000. ); /*scale factor for CSH controlling induction */
    C3AH6_SCALE = ( 2000. * SYSIZE_POW3 / 1000000. ); /*scale factor for C3AH6 controlling induction of aluminates */

    alloc_char_3D(mic, SYSIZE);
    alloc_int_3D(mic_CSH, SYSIZE);
    alloc_char_3D(micorig, SYSIZE);
    alloc_long_3D(micpart, SYSIZE);
    alloc_int_3D(mask, SYSIZE + 1);
    alloc_int_3D(ArrPerc, SYSIZE);
    alloc_int_3D(ConnNumbers, SYSIZE);
    alloc_shortint_3D(cshage, SYSIZE);
    alloc_shortint_3D(faces, SYSIZE);

#ifdef CMLFILE
    F->get_value(13, ( long & )read_micr);
#endif

    if ( !read_micr && generateMicrostructure == 1 ) { //generate new microstructure
        if ( genpartnew() == 1 ) { //read input file for RVE generation, if unsuccessful microstructure generation, return
            return 1;
        }

#ifdef PRINTF
        printf("MONOPHASE microstructure created\n");
#endif
        distrib3d(); //read autocorrelation functions and distribute clinker phases in RVE
    }

    readhydrparam(); //read hydration parameters

    return 0;
}


/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */
double CemhydMatStatus :: ran1(int *idum)
/* Calls: no routines */
/* Called by: gsphere,makefloc */
{
    int j, k;
    NDIV = 1.0 / ( 1.0 + ( IM - 1.0 ) / NTAB );
    RNMX = ( 1.0 - EPS );
    AM = ( 1.0 / IM );

    if ( ( * idum <= 0 ) || ( iy == 0 ) ) {
        * idum = ( -* idum > * idum ) ? -* idum : * idum; //MAX(-*idum>*idum);
        for ( j = NTAB + 7; j >= 0; j-- ) {
            k = * idum / IQ;
            * idum = IA * ( * idum - k * IQ ) - IR * k;
            if ( * idum < 0 ) {
                * idum += IM;
            }

            if ( j < NTAB ) {
                iv [ j ] = * idum;
            }
        }

        iy = iv [ 0 ];
    }

    k = * idum / IQ;
    * idum = IA * ( * idum - k * IQ ) - IR * k;
    if ( * idum < 0 ) {
        * idum += IM;
    }

    j = ( int ) ( iy * NDIV );
    iy = iv [ j ];
    iv [ j ] = * idum;
    return AM * iy < RNMX ? AM * iy : RNMX; //MIN(AM*iy,RNMX);
}


/* routine to add a flat plate aggregate in the microstructure */
void CemhydMatStatus :: addagg(void)
/* Calls: no other routines */
/* Called by: main program */
{
    int ix, iy, iz;
    int agglo, agghi;

    /* Be sure aggregate size is an even integer */
    do {
        printf("Enter thickness of aggregate to place (an even integer) \n");
#ifdef CMLFILE
        F->get_next_line_in_section(0, ( long & )aggsize);
#endif
#ifdef TINYXML
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, aggsize);
#endif

        //fscanf(in, "%d",&aggsize);
#ifdef PRINTF
        printf("%d\n", aggsize);
#endif
    } while ( ( ( aggsize % 2 ) != 0 ) || ( aggsize > ( SYSSIZE - 2 ) ) );

    if ( aggsize != 0 ) {
        agglo = ( SYSSIZE / 2 ) - ( ( aggsize - 2 ) / 2 );
        agghi = ( SYSSIZE / 2 ) + ( aggsize / 2 );

        /* Aggregate is placed in yz plane */
        for ( ix = agglo; ix <= agghi; ix++ ) {
            for ( iy = 1; iy <= SYSSIZE; iy++ ) {
                for ( iz = 1; iz <= SYSSIZE; iz++ ) {
                    /* Mark aggregate in both particle and microstructure images */
                    cement [ ix ] [ iy ] [ iz ] = AGG;
                    cemreal [ ix ] [ iy ] [ iz ] = AGG;
                }
            }
        }
    }
}



/* routine to check or perform placement of sphere of ID phasein */
/* centered at location (xin,yin,zin) of radius radd */
/* wflg=1 check for fit of sphere */
/* wflg=2 place the sphere */
/* phasein and phase2 are phases to assign to cement and cemreal images resp. */
int CemhydMatStatus :: chksph(int xin, int yin, int zin, int radd, int wflg, int phasein, int phase2)
/* Calls: no other routines */
/* Called by: gsphere */
{
    int nofits, xp, yp, zp, i, j, k;
    float dist, xdist, ydist, zdist, ftmp;

    nofits = 0;   /* Flag indicating if placement is possible */
    /* Check all pixels within the digitized sphere volume */
    for ( i = xin - radd; ( ( i <= xin + radd ) && ( nofits == 0 ) ); i++ ) {
        xp = i;
        /* use periodic boundary conditions for sphere placement */
        if ( xp < 1 ) {
            xp += SYSSIZE;
        } else if ( xp > SYSSIZE ) {
            xp -= SYSSIZE;
        }

        ftmp = ( float ) ( i - xin );
        xdist = ftmp * ftmp;
        for ( j = yin - radd; ( ( j <= yin + radd ) && ( nofits == 0 ) ); j++ ) {
            yp = j;
            /* use periodic boundary conditions for sphere placement */
            if ( yp < 1 ) {
                yp += SYSSIZE;
            } else if ( yp > SYSSIZE ) {
                yp -= SYSSIZE;
            }

            ftmp = ( float ) ( j - yin );
            ydist = ftmp * ftmp;
            for ( k = zin - radd; ( ( k <= zin + radd ) && ( nofits == 0 ) ); k++ ) {
                zp = k;
                /* use periodic boundary conditions for sphere placement */
                if ( zp < 1 ) {
                    zp += SYSSIZE;
                } else if ( zp > SYSSIZE ) {
                    zp -= SYSSIZE;
                }

                ftmp = ( float ) ( k - zin );
                zdist = ftmp * ftmp;

                /* Compute distance from center of sphere to this pixel */
                dist = sqrt(xdist + ydist + zdist);
                if ( ( dist - 0.5 ) <= ( float ) radd ) {
                    /* Perform placement */
                    if ( wflg == 2 ) {
                        cement [ xp ] [ yp ] [ zp ] = phasein;
                        cemreal [ xp ] [ yp ] [ zp ] = phase2;
                    }
                    /* or check placement */
                    else if ( ( wflg == 1 ) && ( cement [ xp ] [ yp ] [ zp ] != POROSITY ) ) {
                        nofits = 1;
                    }
                }

                /* Check for overlap with aggregate */
                if ( ( wflg == 1 ) && ( ( fabs( xp - ( ( float ) ( SYSSIZE + 1 ) / 2.0 ) ) ) < ( ( float ) aggsize / 2.0 ) ) ) {
                    nofits = 1;
                }
            }
        }
    }

    /* return flag indicating if sphere will fit */
    return ( nofits );
}



/* routine to place spheres of various sizes and phases at random */
/* locations in 3-D microstructure */
/* numgen is number of different size spheres to place */
/* numeach holds the number of each size class */
/* sizeeach holds the radius of each size class */
/* pheach holds the phase of each size class */
int CemhydMatStatus :: gsphere(int numgen, long int *numeach, int *sizeeach, int *pheach) {
    /* Calls: makesph, ran1 */
    /* Called by: create */
    int count, x, y, z, radius, ig, tries, phnow;
    long int jg, i;
    float testgyp, typegyp;

    /* Generate spheres of each size class in turn (largest first) */
    for ( ig = 0; ig < numgen; ig++ ) {
        phnow = pheach [ ig ]; /* phase for this class */
        radius = sizeeach [ ig ]; /* radius for this class */
        /* loop for each sphere in this size class */
        for ( jg = 1; jg <= numeach [ ig ]; jg++ ) {
            tries = 0;
            /* Stop after MAXTRIES random tries */
            do {
                tries += 1;
                /* generate a random center location for the sphere */
                x = ( int ) ( ( float ) SYSSIZE * ran1(seed) ) + 1;
                y = ( int ) ( ( float ) SYSSIZE * ran1(seed) ) + 1;
                z = ( int ) ( ( float ) SYSSIZE * ran1(seed) ) + 1;
                /* See if the sphere will fit at x,y,z */
                /* Include dispersion distance when checking */
                /* to insure requested separation between spheres */
                count = chksph(x, y, z, radius + dispdist, 1, npart + CEM, 0);
                if ( ( tries > MAXTRIES ) && ( dispdist == 2 ) ) {
                    tries = 0;
                    dispdist += 1;
                }

                if ( tries > MAXTRIES ) {
                    printf("Could not place sphere %d after %ld random attempts \n", npart, MAXTRIES);
                    printf("Skipping this microstructure parameters\n");
                    for ( i = 1; i <= npart; i++ ) {
                        free(clust [ i ]);
                    }

                    return ( 1 );
                }
            } while ( count != 0 );

            /* place the sphere at x,y,z */
            npart += 1;
            if ( npart >= NPARTC ) {
                printf("Too many spheres being generated \n");
                printf("User needs to increase value of NPARTC at top of C-code\n");
                printf("Skipping this microstructure parameters\n");
                return ( 1 );
            }

            /* Allocate space for new particle info */
            clust [ npart ] = ( struct cluster * ) malloc( sizeof( struct cluster ) );
            clust [ npart ]->partid = npart;
            clust [ npart ]->clustid = npart;
            /* Default to cement placement */
            clust [ npart ]->partphase = CEMID;
            clust [ npart ]->x = x;
            clust [ npart ]->y = y;
            clust [ npart ]->z = z;
            clust [ npart ]->r = radius;
            clusleft += 1;
            if ( phnow == 1 ) {
                testgyp = ran1(seed);
                /* Do not use dispersion distance when placing particle */
                if ( ( ( testgyp > probgyp ) && ( ( target_sulfate - n_sulfate ) < ( target_total - n_total ) ) ) || ( n_sulfate > target_sulfate ) || ( volpart [ radius ] > ( target_sulfate - n_sulfate ) ) || ( numeach [ ig ] <= 2 ) ) {
                    count = chksph(x, y, z, radius, 2, npart + CEM - 1, CEMID);
                    n_total += volpart [ radius ];
                } else {
                    /* Place particle as gypsum */
                    typegyp = ran1(seed);
                    n_total += volpart [ radius ];
                    n_sulfate += volpart [ radius ];
                    if ( ( probanh >= 1.0 ) || ( ( typegyp < probanh ) && ( n_anhydrite < target_anhydrite ) && ( volpart [ radius ] <= ( target_anhydrite - n_anhydrite ) ) ) ) {
                        /* Place particle as anhydrite */
                        n_anhydrite += volpart [ radius ];
                        count = chksph(x, y, z, radius, 2, npart + CEM - 1, ANHYDRITE);
                        clust [ npart ]->partphase = ANHYDRITE;
                    } else if ( ( ( probanh + probhem ) >= 1.0 ) || ( ( typegyp < ( probanh + probhem ) ) && ( n_hemi < target_hemi ) && ( volpart [ radius ] <= ( target_hemi - n_hemi ) ) ) ) {
                        /* Place particle as hemihydrate */
                        n_hemi += volpart [ radius ];
                        count = chksph(x, y, z, radius, 2, npart + CEM - 1, HEMIHYDRATE);
                        clust [ npart ]->partphase = HEMIHYDRATE;
                    } else {
                        count = chksph(x, y, z, radius, 2, npart + CEM - 1, GYPID);
                        /* Correct phase ID of particle */
                        clust [ npart ]->partphase = GYPID;
                    }
                }
            }
            /* place as inert, CaCO3, C2S, slag, or pozzolanic material */
            else {
                count = chksph(x, y, z, radius, 2, npart + CEM - 1, phnow);
                /* Correct phase ID of particle */
                clust [ npart ]->partphase = phnow;
            }

            clust [ npart ]->nextpart = NULL;
        }
    }

    //deallocate
    for ( i = 1; i <= npart; i++ ) {
        free(clust [ i ]);
        //printf("Dealloc clust %ld of %d", i,npart);
    }

    return ( 0 );
}


/* routine to obtain user input and create a starting microstructure */
int CemhydMatStatus :: create(void) {
    /* Calls: gsphere */
    /* Called by: main program */
    int numsize;
    int *sphrad, *sphase;
    long int *sphnum;
    long int inval1;
    int isph, inval;

    sphrad = new int [ NUMSIZES ];
    sphase = new int [ NUMSIZES ];
    sphnum = new long int [ NUMSIZES ];


    do {
#ifdef PRINTF
        printf("Enter number of different size spheres to use(max. is %d) \n", NUMSIZES);
#endif
#ifdef CMLFILE
        F->get_next_line_in_section(0, ( long & )numsize);
#endif
#ifdef TINYXML
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, numsize);
#endif
        //fscanf(in, "%d",&numsize);
#ifdef PRINTF
        printf("%d \n", numsize);
#endif
    } while ( ( numsize > NUMSIZES ) || ( numsize < 0 ) );

    do {
#ifdef PRINTF
        printf("Enter dispersion factor (separation distance in pixels) for spheres (0-2) \n");
        printf("0 corresponds to totally random placement \n");
#endif
#ifdef CMLFILE
        F->get_next_line_in_section(0, ( long & )dispdist);
#endif
#ifdef TINYXML
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, dispdist);
#endif
        //fscanf(in, "%d",&dispdist);
#ifdef PRINTF
        printf("%d \n", dispdist);
#endif
    } while ( ( dispdist < 0 ) || ( dispdist > 2 ) );

    do {
#ifdef PRINTF
        printf("Enter probability for gypsum particles on a random particle basis (0.0-1.0) \n");
#endif
#ifdef CMLFILE
        F->get_next_line_in_section(0, probgyp);
#endif
#ifdef TINYXML
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", "dihydrate", probgyp);
#endif
        //fscanf(in, "%f",&probgyp);
#ifdef PRINTF
        printf("%f \n", probgyp);
#endif
    } while ( ( probgyp < 0.0 ) || ( probgyp > 1.0 ) );

    do {
#ifdef PRINTF
        printf("Enter probabilities for hemihydrate and anhydrite forms of gypsum (0.0-1.0) \n");
#endif
#ifdef CMLFILE
        F->get_next_line_in_section(0, probhem);
        F->get_next_line_in_section(0, probanh);
#endif
#ifdef TINYXML
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", "hemihydrate", probhem);
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", "anhydrite", probanh);
#endif
        //fscanf(in, "%f %f",&probhem,&probanh);
#ifdef PRINTF
        printf("%f %f\n", probhem, probanh);
#endif
    } while ( ( probhem < 0.0 ) || ( probhem > 1.0 ) || ( probanh < 0.0 ) || ( probanh > 1.0 ) || ( ( probanh + probhem ) > 1.001 ) );

    if ( ( numsize > 0 ) && ( numsize < ( NUMSIZES + 1 ) ) ) {
#ifdef PRINTF
        printf("Enter number, radius, and phase ID for each sphere class (largest radius 1st) \n");
        printf("Phases are %d- Cement and (random) calcium sulfate, %d- C2S, %d- Gypsum, %d- hemihydrate %d- anhydrite %d- Pozzolanic, %d- Inert, %d- Slag, %d- CaCO3 %d- Fly Ash \n", CEMID, C2SID, GYPID, HEMIHYDRATE, ANHYDRITE, POZZID, INERTID, SLAGID, CACO3, FLYASH);
#endif
        /* Obtain input for each size class of spheres */
        for ( isph = 0; isph < numsize; isph++ ) {
#ifdef PRINTF
            printf("Enter number of spheres of class %d \n", isph + 1);
#endif
#ifdef CMLFILE
            F->get_next_line_in_section(0, inval1);
#endif
#ifdef TINYXML
            QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, inval1);
            //inval1 = static_cast<long int>(inval);
#endif
            //fscanf(in, "%ld",&inval1);
#ifdef PRINTF
            printf("%ld \n", inval1);
#endif
            sphnum [ isph ] = inval1;

            //      do{
#ifdef PRINTF
            printf("Enter radius of spheres of class %d \n", isph + 1);
            printf("(Integer <=%d please) \n", SYSSIZE / 3);
#endif
#ifdef CMLFILE
            F->get_next_line_in_section(0, ( long & )inval);
#endif
#ifdef TINYXML
            QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, inval);
#endif
            //fscanf(in, "%d",&inval);
            if ( inval > ( SYSSIZE / 3 ) ) {
                printf("Given radius %d exceeded maximum radius of %d, terminating\n", inval, SYSSIZE / 3);
                exit(0);
            }

#ifdef PRINTF
            printf("%d \n", inval);
#endif
            //} while ((inval<0)||(inval>(SYSSIZE/3)));

            sphrad [ isph ] = inval;
            do {
#ifdef PRINTF
                printf("Enter phase of spheres of class %d \n", isph + 1);
#endif
#ifdef CMLFILE
                F->get_next_line_in_section(0, ( long & )inval);
#endif
#ifdef TINYXML
                QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, inval);
#endif
                //fscanf(in, "%d",&inval);
#ifdef PRINTF
                printf("%d \n", inval);
#endif
            } while ( ( inval != CEMID ) && ( inval != C2SID ) && ( inval != GYPID ) && ( inval != HEMIHYDRATE ) && ( inval != ANHYDRITE ) && ( inval != POZZID ) && ( inval != INERTID ) && ( inval != SLAGID ) && ( inval != FLYASH ) && ( inval != CACO3 ) && ( inval != AGG ) && ( inval != ASG ) );

            sphase [ isph ] = inval;
            if ( inval == CEMID ) {
                target_total += sphnum [ isph ] * volpart [ sphrad [ isph ] ];
            }
        }

        /* Determine target pixel counts for calcium sulfate forms */
        target_sulfate = ( int ) ( ( float ) target_total * probgyp );
        target_anhydrite = ( int ) ( ( float ) target_total * probgyp * probanh );
        target_hemi = ( int ) ( ( float ) target_total * probgyp * probhem );
        if ( gsphere(numsize, sphnum, sphrad, sphase) == 1 ) { //unsuccessful generation of microstructure due to excessive amount of particles or too dense
            delete [] sphrad;
            delete [] sphase;
            delete [] sphnum;
            return ( 1 );
        }
    }

    delete [] sphrad;
    delete [] sphase;
    delete [] sphnum;


    return ( 0 );
}


/* Routine to draw a particle during flocculation routine */
/* See routine chksph for definition of parameters */
void CemhydMatStatus :: drawfloc(int xin, int yin, int zin, int radd, int phasein, int phase2) {
    /* Calls: no other routines */
    /* Called by: makefloc */
    int xp, yp, zp, i, j, k;
    float dist, xdist, ydist, zdist, ftmp;

    /* Check all pixels within the digitized sphere volume */
    for ( i = xin - radd; ( i <= xin + radd ); i++ ) {
        xp = i;
        /* use periodic boundary conditions for sphere placement */
        if ( xp < 1 ) {
            xp += SYSSIZE;
        } else if ( xp > SYSSIZE ) {
            xp -= SYSSIZE;
        }

        ftmp = ( float ) ( i - xin );
        xdist = ftmp * ftmp;
        for ( j = yin - radd; ( j <= yin + radd ); j++ ) {
            yp = j;
            /* use periodic boundary conditions for sphere placement */
            if ( yp < 1 ) {
                yp += SYSSIZE;
            } else if ( yp > SYSSIZE ) {
                yp -= SYSSIZE;
            }

            ftmp = ( float ) ( j - yin );
            ydist = ftmp * ftmp;
            for ( k = zin - radd; ( k <= zin + radd ); k++ ) {
                zp = k;
                /* use periodic boundary conditions for sphere placement */
                if ( zp < 1 ) {
                    zp += SYSSIZE;
                } else if ( zp > SYSSIZE ) {
                    zp -= SYSSIZE;
                }

                ftmp = ( float ) ( k - zin );
                zdist = ftmp * ftmp;

                /* Compute distance from center of sphere to this pixel */
                dist = sqrt(xdist + ydist + zdist);
                if ( ( dist - 0.5 ) <= ( float ) radd ) {
                    /* Update both cement and cemreal images */
                    cement [ xp ] [ yp ] [ zp ] = phasein;
                    cemreal [ xp ] [ yp ] [ zp ] = phase2;
                }
            }
        }
    }
}


/* Routine to check particle placement during flocculation */
/* for particle of size radd centered at (xin,yin,zin) */
/* Returns flag indicating if placement is possible */
int CemhydMatStatus :: chkfloc(int xin, int yin, int zin, int radd) {
    /* Calls: no other routines */
    /* Called by: makefloc */
    int nofits, xp, yp, zp, i, j, k;
    float dist, xdist, ydist, zdist, ftmp;

    nofits = 0;   /* Flag indicating if placement is possible */

    /* Check all pixels within the digitized sphere volume */
    for ( i = xin - radd; ( ( i <= xin + radd ) && ( nofits == 0 ) ); i++ ) {
        xp = i;
        /* use periodic boundary conditions for sphere placement */
        if ( xp < 1 ) {
            xp += SYSSIZE;
        } else if ( xp > SYSSIZE ) {
            xp -= SYSSIZE;
        }

        ftmp = ( float ) ( i - xin );
        xdist = ftmp * ftmp;
        for ( j = yin - radd; ( ( j <= yin + radd ) && ( nofits == 0 ) ); j++ ) {
            yp = j;
            /* use periodic boundary conditions for sphere placement */
            if ( yp < 1 ) {
                yp += SYSSIZE;
            } else if ( yp > SYSSIZE ) {
                yp -= SYSSIZE;
            }

            ftmp = ( float ) ( j - yin );
            ydist = ftmp * ftmp;
            for ( k = zin - radd; ( ( k <= zin + radd ) && ( nofits == 0 ) ); k++ ) {
                zp = k;
                /* use periodic boundary conditions for sphere placement */
                if ( zp < 1 ) {
                    zp += SYSSIZE;
                } else if ( zp > SYSSIZE ) {
                    zp -= SYSSIZE;
                }

                ftmp = ( float ) ( k - zin );
                zdist = ftmp * ftmp;

                /* Compute distance from center of sphere to this pixel */
                dist = sqrt(xdist + ydist + zdist);
                if ( ( dist - 0.5 ) <= ( float ) radd ) {
                    if ( ( cement [ xp ] [ yp ] [ zp ] != POROSITY ) ) {
                        /* Record ID of particle hit */
                        nofits = cement [ xp ] [ yp ] [ zp ];
                    }
                }

                /* Check for overlap with aggregate */
                if ( ( fabs( xp - ( ( float ) ( SYSSIZE + 1 ) / 2.0 ) ) ) < ( ( float ) aggsize / 2.0 ) ) {
                    nofits = AGG;
                }
            }
        }
    }

    /* return flag indicating if sphere will fit */
    return ( nofits );
}


/* routine to perform flocculation of particles */
void CemhydMatStatus :: makefloc(void)
/* Calls: drawfloc, chkfloc, ran1 */
/* Called by: main program */
{
    int partdo, numfloc;
    int nstart;
    int nleft = 0, ckall;
    int xm, ym, zm, moveran;
    int xp, yp, zp, rp, clushit, valkeep;
    int iclus;
    int *cluspart;
    struct cluster *parttmp, *partpoint, *partkeep = NULL;

    cluspart = new int [ NPARTC ];


    nstart = npart; /* Counter of number of flocs remaining */
    for ( iclus = 1; iclus <= npart; iclus++ ) {
        cluspart [ iclus ] = iclus;
    }

    do {
#ifdef PRINTF
        printf("Enter number of flocs desired at end of routine (>0) \n");
#endif
#ifdef CMLFILE
        F->get_next_line_in_section(0, ( long & )numfloc);
#endif
#ifdef TINYXML
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, numfloc);
#endif
        //fscanf(in, "%d",&numfloc);
#ifdef PRINTF
        printf("%d\n", numfloc);
#endif
    } while ( numfloc <= 0 );

    while ( nstart > numfloc ) {
        nleft = 0;

        /* Try to move each cluster in turn */
        for ( iclus = 1; iclus <= npart; iclus++ ) {
            if ( clust [ iclus ] == NULL ) {
                nleft += 1;
            } else {
                xm = ym = zm = 0;
                /* Generate a random move in one of 6 principal directions */
                moveran = ( int ) ( 6. * ran1(seed) );
                switch ( moveran ) {
                case 0:
                    xm = 1;
                    break;
                case 1:
                    xm = ( -1 );
                    break;
                case 2:
                    ym = 1;
                    break;
                case 3:
                    ym = ( -1 );
                    break;
                case 4:
                    zm = 1;
                    break;
                case 5:
                    zm = ( -1 );
                    break;
                default:
                    break;
                }

                /* First erase all particles in cluster */
                partpoint = clust [ iclus ];
                while ( partpoint != NULL ) {
                    xp = partpoint->x;
                    yp = partpoint->y;
                    zp = partpoint->z;
                    rp = partpoint->r;
                    drawfloc(xp, yp, zp, rp, 0, 0);
                    partpoint = partpoint->nextpart;
                }

                ckall = 0;
                /* Now try to draw cluster at new location */
                partpoint = clust [ iclus ];
                while ( ( partpoint != NULL ) && ( ckall == 0 ) ) {
                    xp = partpoint->x + xm;
                    yp = partpoint->y + ym;
                    zp = partpoint->z + zm;
                    rp = partpoint->r;
                    ckall = chkfloc(xp, yp, zp, rp);
                    partpoint = partpoint->nextpart;
                }

                if ( ckall == 0 ) {
                    /* Place cluster particles at new location */
                    partpoint = clust [ iclus ];
                    while ( partpoint != NULL ) {
                        xp = partpoint->x + xm;
                        yp = partpoint->y + ym;
                        zp = partpoint->z + zm;
                        rp = partpoint->r;
                        valkeep = partpoint->partphase;
                        partdo = partpoint->partid;
                        drawfloc(xp, yp, zp, rp, partdo + CEM - 1, valkeep);
                        /* Update particle location */
                        partpoint->x = xp;
                        partpoint->y = yp;
                        partpoint->z = zp;
                        partpoint = partpoint->nextpart;
                    }
                } else {
                    /* A cluster or aggregate was hit */
                    /* Draw particles at old location */
                    partpoint = clust [ iclus ];
                    /* partkeep stores pointer to last particle in list */
                    while ( partpoint != NULL ) {
                        xp = partpoint->x;
                        yp = partpoint->y;
                        zp = partpoint->z;
                        rp = partpoint->r;
                        valkeep = partpoint->partphase;
                        partdo = partpoint->partid;
                        drawfloc(xp, yp, zp, rp, partdo + CEM - 1, valkeep);
                        partkeep = partpoint;
                        partpoint = partpoint->nextpart;
                    }

                    /* Determine the cluster hit */
                    if ( ckall != AGG ) {
                        clushit = cluspart [ ckall - CEM + 1 ];
                        /* Move all of the particles from cluster clushit to cluster iclus */
                        parttmp = clust [ clushit ];
                        /* Attach new cluster to old one */
                        partkeep->nextpart = parttmp;
                        while ( parttmp != NULL ) {
                            cluspart [ parttmp->partid ] = iclus;
                            /* Relabel all particles added to this cluster */
                            parttmp->clustid = iclus;
                            parttmp = parttmp->nextpart;
                        }

                        /* Disengage the cluster that was hit */
                        clust [ clushit ] = NULL;
                        nstart -= 1;
                    }
                }
            }
        }

        printf("Number left was %d but number of clusters is %d \n", nleft, nstart);
    }

    /* end of while loop */
    clusleft = nleft;

    delete [] cluspart;
}


/* routine to assess global phase fractions present in 3-D system */
void CemhydMatStatus :: measure(void)
/* Calls: no other routines */
/* Called by: main program */
{
    long int npor, nc2s, ngyp, ncem, nagg, npozz, ninert, nflyash, nanh, nhem, ncaco3, nslag;
    int i, j, k, valph;

    /* counters for the various phase fractions */
    npor = 0;
    ngyp = 0;
    ncem = 0;
    nagg = 0;
    ninert = 0;
    nslag = 0;
    nc2s = 0;
    npozz = 0;
    nflyash = 0;
    nanh = 0;
    nhem = 0;
    ncaco3 = 0;

    /* Check all pixels in 3-D microstructure */
    for ( i = 1; i <= SYSSIZE; i++ ) {
        for ( j = 1; j <= SYSSIZE; j++ ) {
            for ( k = 1; k <= SYSSIZE; k++ ) {
                valph = cemreal [ i ] [ j ] [ k ];
                if ( valph == POROSITY ) {
                    npor += 1;
                } else if ( valph == CEMID ) {
                    ncem += 1;
                } else if ( valph == C2SID ) {
                    nc2s += 1;
                } else if ( valph == GYPID ) {
                    ngyp += 1;
                } else if ( valph == ANHYDRITE ) {
                    nanh += 1;
                } else if ( valph == HEMIHYDRATE ) {
                    nhem += 1;
                } else if ( valph == AGG ) {
                    nagg += 1;
                } else if ( valph == POZZID ) {
                    npozz += 1;
                } else if ( valph == SLAGID ) {
                    nslag += 1;
                } else if ( valph == INERTID ) {
                    ninert += 1;
                } else if ( valph == FLYASH ) {
                    nflyash += 1;
                } else if ( valph == CACO3 ) {
                    ncaco3 += 1;
                }
            }
        }
    }

    /* Output results */
    printf("\n Phase counts are: \n");
    printf("Porosity= %ld \n", npor);
    printf("Cement= %ld \n", ncem);
    printf("C2S= %ld \n", nc2s);
    printf("Gypsum= %ld \n", ngyp);
    printf("Anhydrite= %ld \n", nanh);
    printf("Hemihydrate= %ld \n", nhem);
    printf("Pozzolan= %ld \n", npozz);
    printf("Inert= %ld \n", ninert);
    printf("Slag= %ld \n", nslag);
    printf("CaCO3= %ld \n", ncaco3);
    printf("Fly Ash= %ld \n", nflyash);
    printf("Aggregate= %ld \n", nagg);
}


/* Routine to measure phase fractions as a function of distance from */
/* aggregate surface*/
void CemhydMatStatus :: measagg(void)
/* Calls: no other routines */
/* Called by: main program */
{
    int phase [ 40 ], ptot;
    int icnt, ixlo, ixhi, iy, iz, phid, idist;
    FILE *aggfile;

    /* By default, results are sent to output file called agglist.out */
    aggfile = fopen("agglist.out", "w");
    printf("Distance  Porosity  Cement C2S  Gypsum  Anhydrite Hemihydrate Pozzolan  Inert   Slag CaCO3  Fly Ash\n");
    fprintf(aggfile, "Distance  Porosity  Cement  C2S Gypsum  Anhydrite Hemihydrate Pozzolan  Inert  Slag  CaCO3   Fly Ash\n");

    /* Increase distance from aggregate in increments of one */
    for ( idist = 1; idist <= ( SYSSIZE - aggsize ) / 2; idist++ ) {
        /* Pixel left of aggregate surface */
        ixlo = ( ( SYSSIZE - aggsize + 2 ) / 2 ) - idist;
        /* Pixel right of aggregate surface */
        ixhi = ( ( SYSSIZE + aggsize ) / 2 ) + idist;

        /* Initialize phase counts for this distance */
        for ( icnt = 0; icnt < 39; icnt++ ) {
            phase [ icnt ] = 0;
        }

        ptot = 0;

        /* Check all pixels which are this distance from aggregate surface */
        for ( iy = 1; iy <= SYSSIZE; iy++ ) {
            for ( iz = 1; iz <= SYSSIZE; iz++ ) {
                phid = cemreal [ ixlo ] [ iy ] [ iz ];
                ptot += 1;
                if ( phid <= FLYASH ) {
                    phase [ phid ] += 1;
                }

                phid = cemreal [ ixhi ] [ iy ] [ iz ];
                ptot += 1;
                if ( phid <= FLYASH ) {
                    phase [ phid ] += 1;
                }
            }
        }

        /* Output results for this distance from surface */
        printf("%d   %d   %d   %d  %d  %d %d %d  %d %d %d %d\n", idist, phase [ 0 ], phase [ CEMID ], phase [ C2SID ], phase [ GYPID ], phase [ ANHYDRITE ], phase [ HEMIHYDRATE ], phase [ POZZID ], phase [ INERTID ], phase [ SLAGID ], phase [ CACO3 ], phase [ FLYASH ]);
        fprintf(aggfile, "%d   %d   %d  %d  %d %d %d  %d %d %d  %d %d\n", idist, phase [ 0 ], phase [ CEMID ], phase [ C2SID ], phase [ GYPID ], phase [ ANHYDRITE ], phase [ HEMIHYDRATE ], phase [ POZZID ], phase [ INERTID ], phase [ SLAGID ], phase [ CACO3 ], phase [ FLYASH ]);
    }

    fclose(aggfile);
}


/* routine to assess the connectivity (percolation) of a single phase */
/* Two matrices are used here: one for the current burnt locations */
/* the other to store the newly found burnt locations */
void CemhydMatStatus :: connect(void)
/* Calls: no other routines */
/* Called by: main program */
{
    long int inew, ntop, nthrough, ncur, nnew, ntot;
    int i, j, k, nmatx [ 29000 ], nmaty [ 29000 ], nmatz [ 29000 ];
    int xcn, ycn, zcn, npix, x1, y1, z1, igood, nnewx [ 29000 ], nnewy [ 29000 ], nnewz [ 29000 ];
    int jnew, icur;

    do {
        printf("Enter phase to analyze 0) pores 1) Cement  \n");
#ifdef CMLFILE
        F->get_next_line_in_section(0, ( long & )npix);
#endif
#ifdef TINYXML
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, npix);
#endif
        //fscanf(in, "%d",&npix);
        printf("%d \n", npix);
    } while ( ( npix != 0 ) && ( npix != 1 ) );

    /* counters for number of pixels of phase accessible from top surface */
    /* and number which are part of a percolated pathway */
    ntop = 0;
    nthrough = 0;

    /* percolation is assessed from top to bottom only */
    /* and burning algorithm is periodic in x and y directions */
    k = 1;
    for ( i = 1; i <= SYSSIZE; i++ ) {
        for ( j = 1; j <= SYSSIZE; j++ ) {
            ncur = 0;
            ntot = 0;
            igood = 0; /* Indicates if bottom has been reached */
            if ( ( ( cement [ i ] [ j ] [ k ] == npix ) && ( ( cement [ i ] [ j ] [ SYSSIZE ] == npix ) ||
                                                            ( cement [ i ] [ j ] [ SYSSIZE ] == ( npix + BURNTG ) ) ) ) ||
                ( ( cement [ i ] [ j ] [ SYSSIZE ] >= CEM ) &&
                 ( cement [ i ] [ j ] [ k ] >= CEM ) && ( cement [ i ] [ j ] [ k ] < BURNTG ) && ( npix == 1 ) ) ) {
                /* Start a burn front */
                cement [ i ] [ j ] [ k ] += BURNTG;
                ntot += 1;
                ncur += 1;
                /* burn front is stored in matrices nmat* */
                /* and nnew* */
                nmatx [ ncur ] = i;
                nmaty [ ncur ] = j;
                nmatz [ ncur ] = 1;
                /* Burn as long as new (fuel) pixels are found */
                do {
                    nnew = 0;
                    for ( inew = 1; inew <= ncur; inew++ ) {
                        xcn = nmatx [ inew ];
                        ycn = nmaty [ inew ];
                        zcn = nmatz [ inew ];

                        /* Check all six neighbors */
                        for ( jnew = 1; jnew <= 6; jnew++ ) {
                            x1 = xcn;
                            y1 = ycn;
                            z1 = zcn;
                            if ( jnew == 1 ) {
                                x1 -= 1;
                                if ( x1 < 1 ) {
                                    x1 += SYSSIZE;
                                }
                            } else if ( jnew == 2 ) {
                                x1 += 1;
                                if ( x1 > SYSSIZE ) {
                                    x1 -= SYSSIZE;
                                }
                            } else if ( jnew == 3 ) {
                                y1 -= 1;
                                if ( y1 < 1 ) {
                                    y1 += SYSSIZE;
                                }
                            } else if ( jnew == 4 ) {
                                y1 += 1;
                                if ( y1 > SYSSIZE ) {
                                    y1 -= SYSSIZE;
                                }
                            } else if ( jnew == 5 ) {
                                z1 -= 1;
                            } else if ( jnew == 6 ) {
                                z1 += 1;
                            }

                            /* Nonperiodic in z direction so be sure to remain in the 3-D box */
                            if ( ( z1 >= 1 ) && ( z1 <= SYSSIZE ) ) {
                                if ( ( cement [ x1 ] [ y1 ] [ z1 ] == npix ) || ( ( cement [ x1 ] [ y1 ] [ z1 ] >= CEM ) &&
                                                                                 ( cement [ x1 ] [ y1 ] [ z1 ] < BURNTG ) && ( npix == 1 ) ) ) {
                                    ntot += 1;
                                    cement [ x1 ] [ y1 ] [ z1 ] += BURNTG;
                                    nnew += 1;
                                    if ( nnew >= 29000 ) {
                                        printf("error in size of nnew \n");
                                    }

                                    nnewx [ nnew ] = x1;
                                    nnewy [ nnew ] = y1;
                                    nnewz [ nnew ] = z1;
                                    /* See if bottom of system has been reached */
                                    if ( z1 == SYSSIZE ) {
                                        igood = 1;
                                    }
                                }
                            }
                        }
                    }

                    if ( nnew > 0 ) {
                        ncur = nnew;
                        /* update the burn front matrices */
                        for ( icur = 1; icur <= ncur; icur++ ) {
                            nmatx [ icur ] = nnewx [ icur ];
                            nmaty [ icur ] = nnewy [ icur ];
                            nmatz [ icur ] = nnewz [ icur ];
                        }
                    }
                } while ( nnew > 0 );

                ntop += ntot;
                if ( igood == 1 ) {
                    nthrough += ntot;
                }
            }
        }
    }

    printf("Phase ID= %d \n", npix);
    printf("Number accessible from top= %ld \n", ntop);
    printf("Number contained in through pathways= %ld \n", nthrough);

    /* return the burnt sites to their original phase values */
    for ( i = 1; i <= SYSSIZE; i++ ) {
        for ( j = 1; j <= SYSSIZE; j++ ) {
            for ( k = 1; k <= SYSSIZE; k++ ) {
                if ( cement [ i ] [ j ] [ k ] >= BURNTG ) {
                    cement [ i ] [ j ] [ k ] -= BURNTG;
                }
            }
        }
    }
}


/* Routine to output final microstructure to file */
void CemhydMatStatus :: outmic(void)
/* Calls: no other routines */
/* Called by: main program */
{
    FILE *outfile, *partfile;
    char filen [ 80 ], filepart [ 80 ];
    int ix, iy, iz, valout;

#ifdef PRINTF
    printf("Enter name of file to save microstructure to \n");
#endif
#ifdef CMLFILE
    F->get_next_line_in_section(0, filen);
#endif
#ifdef TINYXML
    QueryStringAttributeExt(xmlFile, "Generate_microstructure", countKey++, filen);
#endif
    //fscanf(in, "%s",filen);
    printf("%s\n", filen);

    outfile = fopen(filen, "w");

#ifdef PRINTF
    printf("Enter name of file to save particle IDs to \n");
#endif
#ifdef CMLFILE
    F->get_next_line_in_section(0, filepart);
#endif
#ifdef TINYXML
    QueryStringAttributeExt(xmlFile, "Generate_microstructure", countKey++, filepart);
#endif
    //fscanf(in, "%s",filepart);
#ifdef PRINTF
    printf("%s\n", filepart);
#endif

    partfile = fopen(filepart, "w");

    for ( iz = 1; iz <= SYSSIZE; iz++ ) {
        for ( iy = 1; iy <= SYSSIZE; iy++ ) {
            for ( ix = 1; ix <= SYSSIZE; ix++ ) {
                valout = cemreal [ ix ] [ iy ] [ iz ];
                fprintf(outfile, "%1d\n", valout); //img file
                valout = cement [ ix ] [ iy ] [ iz ];
                if ( valout < 0 ) {
                    valout = 0;
                }

                fprintf(partfile, "%d\n", valout); //id file
            }
        }
    }

    fclose(outfile);
    fclose(partfile);
}


int CemhydMatStatus :: genpartnew(void) {
    int userc;    /* User choice from menu */
    int ig, jg, kg;

    n_sulfate = 0;
    target_sulfate = 0;
    n_total = 0;
    target_total = 0;
    n_anhydrite = 0;
    target_anhydrite = 0;
    n_hemi = 0;
    target_hemi = 0;

    alloc_long_3D(cement, SYSIZE + 1);
    alloc_long_3D(cemreal, SYSIZE + 1);

    clust = new cluster * [ NPARTC ];

    /* Initialize volume array */
    volpart [ 0 ] = 1;
    volpart [ 1 ] = 19;
    volpart [ 2 ] = 81;
    volpart [ 3 ] = 179;
    volpart [ 4 ] = 389;
    volpart [ 5 ] = 739;
    volpart [ 6 ] = 1189;
    volpart [ 7 ] = 1791;
    volpart [ 8 ] = 2553;
    volpart [ 9 ] = 3695;
    volpart [ 10 ] = 4945;
    volpart [ 11 ] = 6403;
    volpart [ 12 ] = 8217;
    volpart [ 13 ] = 10395;
    volpart [ 14 ] = 12893;
    volpart [ 15 ] = 15515;
    volpart [ 16 ] = 18853;
    volpart [ 17 ] = 22575;
    volpart [ 18 ] = 26745;
    volpart [ 19 ] = 31103;
    volpart [ 20 ] = 36137;
    volpart [ 21 ] = 41851;
    volpart [ 22 ] = 47833;
    volpart [ 23 ] = 54435;
    volpart [ 24 ] = 61565;
    volpart [ 25 ] = 69599;
    volpart [ 26 ] = 78205;
    volpart [ 27 ] = 87271;
    volpart [ 28 ] = 97233;
    volpart [ 29 ] = 107783;
    volpart [ 30 ] = 119009;
    volpart [ 31 ] = 131155;
    volpart [ 32 ] = 143761;
    volpart [ 33 ] = 157563;
    volpart [ 34 ] = 172317;
    volpart [ 35 ] = 187511;
    volpart [ 36 ] = 203965;

    //printf("Enter random number seed value (a negative integer) \n");
    //fscanf(in, "%d",&iseed);
    //printf("%d \n",iseed);
    nseed = iseed;
    seed = & nseed;

    /* Initialize counters and system parameters */
    npart = 0;
    aggsize = 0;
    clusleft = 0;

    /* clear the 3-D system to all porosity to start */
    for ( ig = 1; ig <= SYSSIZE; ig++ ) {
        for ( jg = 1; jg <= SYSSIZE; jg++ ) {
            for ( kg = 1; kg <= SYSSIZE; kg++ ) {
                cement [ ig ] [ jg ] [ kg ] = POROSITY; //particle ID
                cemreal [ ig ] [ jg ] [ kg ] = POROSITY; //microstructure
            }
        }
    }

    /* present menu and execute user choice */
    do {
#ifdef PRINTF
        printf(" \n Input User Choice \n");
        printf("1) Exit \n");
        printf("2) Add spherical particles (cement,gypsum, pozzolans, etc.) to microstructure \n");
        printf("3) Flocculate system by reducing number of particle clusters \n");
        printf("4) Measure global phase fractions \n");
        printf("5) Add an aggregate to the microstructure \n");
        printf("6) Measure single phase connectivity (pores or solids) \n");
        printf("7) Measure phase fractions vs. distance from aggregate surface \n");
        printf("8) Output current microstructure to file \n");
#endif
#ifdef CMLFILE
        F->get_next_line_in_section(0, ( long & )userc);
#endif
#ifdef TINYXML
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, userc);
#endif
        //fscanf(in, "%d",&userc);
        //printf("%d \n",userc);
        fflush(stdout);

        switch ( userc ) {
        case 2:
            if ( create() == 1 ) { //unsuccessful generation of microstructure due to excessive amount of particles or too dense
                delete [] clust;
                dealloc_long_3D(cement, SYSIZE + 1);
                dealloc_long_3D(cemreal, SYSIZE + 1);
                return ( 1 );
            }

            break;
        case 3:
            makefloc();
            break;
        case 4:
            measure();
            break;
        case 5:
            addagg();
            break;
        case 6:
            connect();
            break;
        case 7:
            if ( aggsize != 0 ) {
                measagg();
            } else {
                printf("No aggregate present. \n");
            }

            break;
        case 8:
            outmic();
            break;
        default:
            break;
        }
    } while ( userc != 1 );

    //store ID in an array
    for ( ig = 0; ig < SYSSIZE; ig++ ) {
        for ( jg = 0; jg < SYSSIZE; jg++ ) {
            for ( kg = 0; kg < SYSSIZE; kg++ ) {
                micpart [ ig ] [ jg ] [ kg ] = cement [ ig + 1 ] [ jg + 1 ] [ kg + 1 ];
            }
        }
    }

    dealloc_long_3D(cement, SYSIZE + 1);
    delete [] clust;
    return ( 0 );
}


/*************************************************************************/
/****************************DISTRIB3D************************************/
/*************************************************************************/

/* routine to create a template for the sphere of interest of radius size */
/* to be used in curvature evaluation */
/* Called by: runsint */
/* Calls no other routines */
int CemhydMatStatus :: maketemp(int size)
{
    int icirc, xval, yval, zval;
    float xtmp, ytmp;
    float dist;

    /* determine and store the locations of all pixels in the 3-D sphere */
    icirc = 0;
    for ( xval = ( -size ); xval <= size; xval++ ) {
        xtmp = ( float ) ( xval * xval );
        for ( yval = ( -size ); yval <= size; yval++ ) {
            ytmp = ( float ) ( yval * yval );
            for ( zval = ( -size ); zval <= size; zval++ ) {
                dist = sqrt( xtmp + ytmp + ( float ) ( zval * zval ) );
                if ( dist <= ( ( float ) size + 0.5 ) ) {
                    icirc += 1;
                    if ( icirc >= MAXSPH ) {
                        printf("Too many elements in sphere \n");
                        printf("Must change value of MAXSPH parameter \n");
                        printf("Currently set at %ld \n", MAXSPH);
                        exit(1);
                    }

                    xsph [ icirc ] = xval;
                    ysph [ icirc ] = yval;
                    zsph [ icirc ] = zval;
                }
            }
        }
    }

    /* return the number of pixels contained in sphere of radius (size+0.5) */
    return ( icirc );
}

/* routine to count phase fractions (porosity and solids) */
/* Called by main routine */
/* Calls no other routines */
void CemhydMatStatus :: phcount(void) {
    long int npore, nsolid [ 37 ];
    int ix, iy, iz;

    npore = 0;
    for ( ix = 1; ix < 37; ix++ ) {
        nsolid [ ix ] = 0;
    }

    /* check all pixels in the 3-D system */
    for ( ix = 1; ix <= SYSSIZE; ix++ ) {
        for ( iy = 1; iy <= SYSSIZE; iy++ ) {
            for ( iz = 1; iz <= SYSSIZE; iz++ ) {
                if ( mask [ ix ] [ iy ] [ iz ] == 0 ) {
                    npore += 1;
                } else {
                    nsolid [ mask [ ix ] [ iy ] [ iz ] ] += 1;
                }
            }
        }
    }

    printf("Pores are: %ld \n", npore);
    printf("Solids are: %ld %ld %ld %ld %ld %ld\n", nsolid [ 1 ], nsolid [ 2 ],
           nsolid [ 3 ], nsolid [ 4 ], nsolid [ 5 ], nsolid [ 6 ]);
}

/* routine to return number of surface faces exposed to porosity */
/* for pixel located at (xin,yin,zin) */
/* Called by rhcalc */
/* Calls no other routines */
int CemhydMatStatus :: surfpix(int xin, int yin, int zin) {
    int npix, ix1, iy1, iz1;

    npix = 0;

    /* check each of the six immediate neighbors */
    /* using periodic boundary conditions */
    ix1 = xin - 1;
    if ( ix1 < 1 ) {
        ix1 += SYSSIZE;
    }

    if ( mask [ ix1 ] [ yin ] [ zin ] == 0 ) {
        npix += 1;
    }

    ix1 = xin + 1;
    if ( ix1 > SYSSIZE ) {
        ix1 -= SYSSIZE;
    }

    if ( mask [ ix1 ] [ yin ] [ zin ] == 0 ) {
        npix += 1;
    }

    iy1 = yin - 1;
    if ( iy1 < 1 ) {
        iy1 += SYSSIZE;
    }

    if ( mask [ xin ] [ iy1 ] [ zin ] == 0 ) {
        npix += 1;
    }

    iy1 = yin + 1;
    if ( iy1 > SYSSIZE ) {
        iy1 -= SYSSIZE;
    }

    if ( mask [ xin ] [ iy1 ] [ zin ] == 0 ) {
        npix += 1;
    }

    iz1 = zin - 1;
    if ( iz1 < 1 ) {
        iz1 += SYSSIZE;
    }

    if ( mask [ xin ] [ yin ] [ iz1 ] == 0 ) {
        npix += 1;
    }

    iz1 = zin + 1;
    if ( iz1 > SYSSIZE ) {
        iz1 -= SYSSIZE;
    }

    if ( mask [ xin ] [ yin ] [ iz1 ] == 0 ) {
        npix += 1;
    }

    return ( npix );
}

/* routine to return the current hydraulic radius for phase phin */
/* Calls surfpix */
/* Called by runsint */
float CemhydMatStatus :: rhcalc(int phin) {
    int ix, iy, iz;
    long int porc, surfc;
    float rhval;

    porc = surfc = 0;

    /* Check all pixels in the 3-D volume */
    for ( ix = 1; ix <= SYSSIZE; ix++ ) {
        for ( iy = 1; iy <= SYSSIZE; iy++ ) {
            for ( iz = 1; iz <= SYSSIZE; iz++ ) {
                if ( mask [ ix ] [ iy ] [ iz ] == phin ) {
                    porc += 1;
                    surfc += surfpix(ix, iy, iz);
                }
            }
        }
    }

    printf("Phase area count is %ld \n", porc);
    printf("Phase surface count is %ld \n", surfc);
    rhval = ( float ) porc * 6. / ( 4. * ( float ) surfc );
    printf("Hydraulic radius is %f \n", rhval);
    return ( rhval );
}

/* routine to return count of pixels in a spherical template which are phase */
/* phin or porosity (phase=0) */
/* Calls no other routines */
/* Called by sysinit */
int CemhydMatStatus :: countem(int xp, int yp, int zp, int phin) {
    int xc, yc, zc;
    int cumnum, ic;

    cumnum = 0;
    for ( ic = 1; ic <= nsph; ic++ ) {
        xc = xp + xsph [ ic ];
        yc = yp + ysph [ ic ];
        zc = zp + zsph [ ic ];
        /* Use periodic boundaries */
        if ( xc < 1 ) {
            xc += SYSSIZE;
        } else if ( xc > SYSSIZE ) {
            xc -= SYSSIZE;
        }

        if ( yc < 1 ) {
            yc += SYSSIZE;
        } else if ( yc > SYSSIZE ) {
            yc -= SYSSIZE;
        }

        if ( zc < 1 ) {
            zc += SYSSIZE;
        } else if ( zc > SYSSIZE ) {
            zc -= SYSSIZE;
        }

        if ( ( xc != xp ) || ( yc != yp ) || ( zc != zp ) ) {
            if ( ( mask [ xc ] [ yc ] [ zc ] == phin ) || ( mask [ xc ] [ yc ] [ zc ] == 0 ) ) {
                cumnum += 1;
            }
        }
    }

    return ( cumnum );
}

/* routine to initialize system by determining local curvature */
/* of all phase 1 and phase 2 pixels */
/* Calls countem */
/* Called by runsint */
void CemhydMatStatus :: sysinit(int ph1, int ph2) {
    int count, xl, yl, zl;

    count = 0;
    /* process all pixels in the 3-D box */
    for ( xl = 1; xl <= SYSSIZE; xl++ ) {
        for ( yl = 1; yl <= SYSSIZE; yl++ ) {
            for ( zl = 1; zl <= SYSSIZE; zl++ ) {
                /* determine local curvature */
                /* For phase 1 want to determine number of porosity pixels */
                /* (phase=0) in immediate neighborhood */
                if ( mask [ xl ] [ yl ] [ zl ] == ph1 ) {
                    count = countem(xl, yl, zl, 0);
                }

                /* For phase 2 want to determine number of porosity or phase */
                /* 2 pixels in immediate neighborhood */
                if ( mask [ xl ] [ yl ] [ zl ] == ph2 ) {
                    count = countem(xl, yl, zl, ph2);
                }

                if ( ( count < 0 ) || ( count >= nsph ) ) {
                    printf("Error count is %d \n", count);
                    printf("xl %d  yl %d  zl %d \n", xl, yl, zl);
                }

                /* case where we have a phase 1 surface pixel */
                /* with non-zero local curvature */
                if ( ( count >= 0 ) && ( mask [ xl ] [ yl ] [ zl ] == ph1 ) ) {
                    curvature [ xl ] [ yl ] [ zl ] = count;
                    /* update solid curvature histogram */
                    nsolid [ count ] += 1;
                }

                /* case where we have a phase 2 surface pixel */
                if ( ( count >= 0 ) && ( mask [ xl ] [ yl ] [ zl ] == ph2 ) ) {
                    curvature [ xl ] [ yl ] [ zl ] = count;
                    /* update air curvature histogram */
                    nair [ count ] += 1;
                }
            }
        }
    }

    /* end of xl loop */
}

/* routine to scan system and determine nsolid (ph2) and nair (ph1) */
/* histograms based on values in phase and curvature arrays */
/* Calls no other routines */
/* Called by runsint */
void CemhydMatStatus :: sysscan(int ph1, int ph2) {
    int xd, yd, zd, curvval;

    /* Scan all pixels in 3-D system */
    for ( xd = 1; xd <= SYSSIZE; xd++ ) {
        for ( yd = 1; yd <= SYSSIZE; yd++ ) {
            for ( zd = 1; zd <= SYSSIZE; zd++ ) {
                curvval = curvature [ xd ] [ yd ] [ zd ];

                if ( mask [ xd ] [ yd ] [ zd ] == ph2 ) {
                    nair [ curvval ] += 1;
                } else if ( mask [ xd ] [ yd ] [ zd ] == ph1 ) {
                    nsolid [ curvval ] += 1;
                }
            }
        }
    }
}

/* routine to return how many cells of solid curvature histogram to use */
/* to accomodate nsearch pixels moving */
/* want to use highest values first */
/* Calls no other routines */
/* Called by movepix */
int CemhydMatStatus :: procsol(int nsearch) {
    int valfound, i, stop;
    long int nsofar;

    /* search histogram from top down until cumulative count */
    /* exceeds nsearch */
    valfound = nsph - 1;
    nsofar = 0;
    stop = 0;
    for ( i = ( nsph - 1 ); ( ( i >= 0 ) && ( stop == 0 ) ); i-- ) {
        nsofar += nsolid [ i ];
        if ( nsofar > nsearch ) {
            valfound = i;
            stop = 1;
        }
    }

    return ( valfound );
}

/* routine to determine how many cells of air curvature histogram to use */
/* to accomodate nsearch moving pixels */
/* want to use lowest values first */
/* Calls no other routines */
/* Called by movepix */

int CemhydMatStatus :: procair(int nsearch) {
    int valfound, i, stop;
    long int nsofar;

    /* search histogram from bottom up until cumulative count */
    /* exceeds nsearch */
    valfound = 0;
    nsofar = 0;
    stop = 0;
    for ( i = 0; ( ( i < nsph ) && ( stop == 0 ) ); i++ ) {
        nsofar += nair [ i ];
        if ( nsofar > nsearch ) {
            valfound = i;
            stop = 1;
        }
    }

    return ( valfound );
}

/* routine to move requested number of pixels (ntomove) from highest */
/* curvature phase 1 (ph1) sites to lowest curvature phase 2 (ph2) sites */
/* Calls procsol and procair */
/* Called by runsint */
int CemhydMatStatus :: movepix(int ntomove, int ph1, int ph2) {
    int xloc [ 2100 ], yloc [ 2100 ], zloc [ 2100 ];
    int count1, count2, ntot, countc, i, xp, yp, zp;
    int cmin, cmax, cfg;
    int alldone;
    long int nsolc, nairc, nsum, nsolm, nairm, nst1, nst2, next1, next2;
    float pck, plsol, plair;

    alldone = 0;
    /* determine critical values for removal and placement */
    count1 = procsol(ntomove);
    nsum = 0;
    cfg = 0;
    cmax = count1;
    for ( i = nsph; i > count1; i-- ) {
        if ( ( nsolid [ i ] > 0 ) && ( cfg == 0 ) ) {
            cfg = 1;
            cmax = i;
        }

        nsum += nsolid [ i ];
    }

    /* Determine movement probability for last cell */
    plsol = ( float ) ( ntomove - nsum ) / ( float ) nsolid [ count1 ];
    next1 = ntomove - nsum;
    nst1 = nsolid [ count1 ];

    count2 = procair(ntomove);
    nsum = 0;
    cmin = count2;
    cfg = 0;
    for ( i = 0; i < count2; i++ ) {
        if ( ( nair [ i ] > 0 ) && ( cfg == 0 ) ) {
            cfg = 1;
            cmin = i;
        }

        nsum += nair [ i ];
    }

    /* Determine movement probability for last cell */
    plair = ( float ) ( ntomove - nsum ) / ( float ) nair [ count2 ];
    next2 = ntomove - nsum;
    nst2 = nair [ count2 ];

    /* Check to see if equilibrium has been reached --- */
    /* no further increase in hydraulic radius is possible */
    if ( cmin >= cmax ) {
        alldone = 1;
        printf("Stopping - at equilibrium \n");
        printf("cmin- %d  cmax- %d \n", cmin, cmax);
        return ( alldone );
    }

    /* initialize counters for performing sintering */
    ntot = 0;
    nsolc = 0;
    nairc = 0;
    nsolm = 0;
    nairm = 0;

    /* Now process each pixel in turn */
    for ( xp = 1; xp <= SYSSIZE; xp++ ) {
        for ( yp = 1; yp <= SYSSIZE; yp++ ) {
            for ( zp = 1; zp <= SYSSIZE; zp++ ) {
                countc = curvature [ xp ] [ yp ] [ zp ];
                /* handle phase 1 case first */
                if ( mask [ xp ] [ yp ] [ zp ] == ph1 ) {
                    if ( countc > count1 ) {
                        /* convert from phase 1 to phase 2 */
                        mask [ xp ] [ yp ] [ zp ] = ph2;

                        /* update appropriate histogram cells */
                        nsolid [ countc ] -= 1;
                        nair [ countc ] += 1;
                        /* store the location of the modified pixel */
                        ntot += 1;
                        xloc [ ntot ] = xp;
                        yloc [ ntot ] = yp;
                        zloc [ ntot ] = zp;
                    }

                    if ( countc == count1 ) {
                        nsolm += 1;
                        /* generate probability for pixel being removed */
                        pck = ran1(seed);
                        if ( ( pck < 0 ) || ( pck > 1.0 ) ) {
                            pck = 1.0;
                        }

                        if ( ( ( pck < plsol ) && ( nsolc < next1 ) ) || ( ( nst1 - nsolm ) < ( next1 - nsolc ) ) ) {
                            nsolc += 1;
                            /* convert phase 1 pixel to phase 2 */
                            mask [ xp ] [ yp ] [ zp ] = ph2;

                            /* update appropriate histogram cells */
                            nsolid [ count1 ] -= 1;
                            nair [ count1 ] += 1;
                            /* store the location of the modified pixel */
                            ntot += 1;
                            xloc [ ntot ] = xp;
                            yloc [ ntot ] = yp;
                            zloc [ ntot ] = zp;
                        }
                    }
                }
                /* handle phase 2 case here */
                else if ( mask [ xp ] [ yp ] [ zp ] == ph2 ) {
                    if ( countc < count2 ) {
                        /* convert phase 2 pixel to phase 1 */
                        mask [ xp ] [ yp ] [ zp ] = ph1;

                        nsolid [ countc ] += 1;
                        nair [ countc ] -= 1;
                        ntot += 1;
                        xloc [ ntot ] = xp;
                        yloc [ ntot ] = yp;
                        zloc [ ntot ] = zp;
                    }

                    if ( countc == count2 ) {
                        nairm += 1;
                        pck = ran1(seed);
                        if ( ( pck < 0 ) || ( pck > 1.0 ) ) {
                            pck = 1.0;
                        }

                        if ( ( ( pck < plair ) && ( nairc < next2 ) ) || ( ( nst2 - nairm ) < ( next2 - nairc ) ) ) {
                            nairc += 1;
                            /* convert phase 2 to phase 1 */
                            mask [ xp ] [ yp ] [ zp ] = ph1;

                            nsolid [ count2 ] += 1;
                            nair [ count2 ] -= 1;
                            ntot += 1;
                            xloc [ ntot ] = xp;
                            yloc [ ntot ] = yp;
                            zloc [ ntot ] = zp;
                        }
                    }
                }
            }

            /* end of zp loop */
        }

        /* end of yp loop */
    }

    /* end of xloop */
    printf("ntot is %d \n", ntot);
    return ( alldone );
}

/* routine to execute user input number of cycles of sintering algorithm */
/* Calls maketemp, rhcalc, sysinit, sysscan, and movepix */
/* Called by main routine */
void CemhydMatStatus :: sinter3d(int ph1id, int ph2id, float rhtarget) {
    int natonce, i, rade, j, rflag;
    int keepgo;
    long int curvsum1, curvsum2, pixsum1, pixsum2;
    float rhnow, avecurv1, avecurv2;

    /* initialize the solid and air count histograms */
    for ( i = 0; i <= 499; i++ ) {
        nsolid [ i ] = 0;
        nair [ i ] = 0;
    }

    /* Obtain needed user input */
    natonce = 200;
    rade = 3;
    rflag = 0; /* always initialize system */

    nsph = maketemp(rade);
    printf("nsph is %d \n", nsph);
    if ( rflag == 0 ) {
        sysinit(ph1id, ph2id);
    } else {
        sysscan(ph1id, ph2id);
    }

    i = 0;
    rhnow = rhcalc(ph1id);
    while ( ( rhnow < rhtarget ) && ( i < MAXCYC_SEAL ) ) {
        printf("Now: %f  Target: %f \n", rhnow, rhtarget);
        i += 1;
#ifdef PRINTF
        printf("Cycle: %d \n", i);
#endif
        keepgo = movepix(natonce, ph1id, ph2id);
        /* If equilibrium is reached, then return to calling routine */
        if ( keepgo == 1 ) {
            return;
        }

        curvsum1 = 0;
        curvsum2 = 0;
        pixsum1 = 0;
        pixsum2 = 0;
        /* Determine average curvatures for phases 1 and 2 */
        for ( j = 0; j <= nsph; j++ ) {
            pixsum1 += nsolid [ j ];
            curvsum1 += ( j * nsolid [ j ] );
            pixsum2 += nair [ j ];
            curvsum2 += ( j * nair [ j ] );
        }

        avecurv1 = ( float ) curvsum1 / ( float ) pixsum1;
        avecurv2 = ( float ) curvsum2 / ( float ) pixsum2;
        printf("Ave. solid curvature: %f \n", avecurv1);
        printf("Ave. air curvature: %f \n", avecurv2);
        rhnow = rhcalc(ph1id);
    }
}

void CemhydMatStatus :: stat3d(void) {
    int valin, ix, iy, iz;
    int ix1, iy1, iz1, k;
    long int voltot, surftot;

    for ( ix = 0; ix <= 42; ix++ ) {
        volume [ ix ] = surface [ ix ] = 0;
    }

    /* Read in image and accumulate volume totals */
    for ( iz = 1; iz <= SYSIZE; iz++ ) {
        for ( iy = 1; iy <= SYSIZE; iy++ ) {
            for ( ix = 1; ix <= SYSIZE; ix++ ) {
                valin = mask [ ix ] [ iy ] [ iz ];
                volume [ valin ] += 1;
            }
        }
    }


    for ( iz = 1; iz <= SYSIZE; iz++ ) {
        for ( iy = 1; iy <= SYSIZE; iy++ ) {
            for ( ix = 1; ix <= SYSIZE; ix++ ) {
                if ( mask [ ix ] [ iy ] [ iz ] != 0 ) {
                    valin = mask [ ix ] [ iy ] [ iz ];
                    /* Check six neighboring pixels for porosity */
                    for ( k = 1; k <= 6; k++ ) {
                        switch ( k ) {
                        case 1:
                            ix1 = ix - 1;
                            if ( ix1 < 1 ) {
                                ix1 += SYSIZE;
                            }

                            iy1 = iy;
                            iz1 = iz;
                            break;
                        case 2:
                            ix1 = ix + 1;
                            if ( ix1 > SYSIZE ) {
                                ix1 -= SYSIZE;
                            }

                            iy1 = iy;
                            iz1 = iz;
                            break;
                        case 3:
                            iy1 = iy - 1;
                            if ( iy1 < 1 ) {
                                iy1 += SYSIZE;
                            }

                            ix1 = ix;
                            iz1 = iz;
                            break;
                        case 4:
                            iy1 = iy + 1;
                            if ( iy1 > SYSIZE ) {
                                iy1 -= SYSIZE;
                            }

                            ix1 = ix;
                            iz1 = iz;
                            break;
                        case 5:
                            iz1 = iz - 1;
                            if ( iz1 < 1 ) {
                                iz1 += SYSIZE;
                            }

                            iy1 = iy;
                            ix1 = ix;
                            break;
                        case 6:
                            iz1 = iz + 1;
                            if ( iz1 > SYSIZE ) {
                                iz1 -= SYSIZE;
                            }

                            iy1 = iy;
                            ix1 = ix;
                            break;
                        default:
                            break;
                        }

                        if ( ( ix1 < 1 ) || ( iy1 < 1 ) || ( iz1 < 1 ) || ( ix1 > SYSIZE ) || ( iy1 > SYSIZE ) || ( iz1 > SYSIZE ) ) {
                            printf("%d %d %d \n", ix1, iy1, iz1);
                            exit(1);
                        }

                        if ( mask [ ix1 ] [ iy1 ] [ iz1 ] == 0 ) {
                            surface [ valin ] += 1;
                        }
                    }
                }
            }
        }
    }

#ifdef PRINTF
    printf("Phase    Volume      Surface     Volume    Surface \n");
    printf(" ID      count        count      fraction  fraction \n");
#endif
    /* Only include clinker phases in surface area fraction calculation */
    surftot = surface [ 1 ] + surface [ 2 ] + surface [ 3 ] + surface [ 4 ];
    voltot = volume [ 1 ] + volume [ 2 ] + volume [ 3 ] + volume [ 4 ];
    k = 0;
#ifdef PRINTF
    printf("  %d    %8ld     %8ld  \n", k, volume [ 0 ], surface [ 0 ]);

    for ( k = 1; k <= 4; k++ ) {
        printf("  %d    %8ld     %8ld     %.5f   %.5f\n", k, volume [ k ], surface [ k ],
               ( float ) volume [ k ] / ( float ) voltot, ( float ) surface [ k ] / ( float ) surftot);
    }

    printf("Total  %8ld     %8ld\n\n\n", voltot, surftot);

    for ( k = 5; k <= 11; k++ ) {
        printf("  %d    %8ld     %8ld\n", k, volume [ k ], surface [ k ]);
    }

    printf(" 20    %8ld     %8ld\n", volume [ 20 ], surface [ 20 ]);

    for ( k = 24; k <= 27; k++ ) {
        printf(" %d    %8ld     %8ld\n", k, volume [ k ], surface [ k ]);
    }

    printf(" 28    %8ld     %8ld\n", volume [ 28 ], surface [ 28 ]);
#endif
}

void CemhydMatStatus :: rand3d(int phasein, int phaseout, float xpt) {
    int ires;
    float s2, ss, sdiff, xtmp, ytmp;
    //static float normm[SYSIZE+1][SYSIZE+1][SYSIZE+1];
    //static float res[SYSIZE+1][SYSIZE+1][SYSIZE+1];
    double ***normm, ***res;
    //static float filter [32][32][32];
    double ***filter;
    int done, r [ 61 ];
    //static float s[61],xr[61],sum[502];
    float *s, *xr, *sum;
    double val2;
    double t1, t2, x1, x2, u1, u2, xrad, resmax, resmin;
    float xtot, filval, radius, sect, sumtot, vcrit;
    int valin, r1, r2, i1, i2, i3, i, j, k, j1, k1;
    int ido, iii, jjj, ix, iy, iz, index;
#ifdef CMLFILE
    char tempstr [ 256 ];
#endif

    alloc_double_3D(normm, SYSIZE + 1);
    alloc_double_3D(res, SYSIZE + 1);
    alloc_double_3D(filter, 32);

    s = new float [ 61 ];
    xr = new float [ 61 ];
    sum = new float [ 502 ];

    /* Create the Gaussian noise image */
    i1 = i2 = i3 = 1;
    for ( i = 1; i <= ( ( SYSIZE * SYSIZE * SYSIZE ) / 2 ); i++ ) {
        u1 = ran1(seed);
        u2 = ran1(seed);
        t1 = 2. * M_PI * u2;
        t2 = sqrt( -2. * log(u1) );
        x1 = cos(t1) * t2;
        x2 = sin(t1) * t2;
        normm [ i1 ] [ i2 ] [ i3 ] = x1;
        i1 += 1;
        if ( i1 > SYSIZE ) {
            i1 = 1;
            i2 += 1;
            if ( i2 > SYSIZE ) {
                i2 = 1;
                i3 += 1;
            }
        }

        normm [ i1 ] [ i2 ] [ i3 ] = x2;
        i1 += 1;
        if ( i1 > SYSIZE ) {
            i1 = 1;
            i2 += 1;
            if ( i2 > SYSIZE ) {
                i2 = 1;
                i3 += 1;
            }
        }
    }

    /* Now perform the convolution */
#ifdef CMLFILE
    F->get_next_line_in_section(0, ( long & )ido);
#endif
    //fscanf(in,"%d",&ido);
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, ido);
#endif

#ifdef PRINTF
    printf("Number of points in correlation file is %d \n", ido);
#endif

    for ( i = 1; i <= ido; i++ ) {
#ifdef CMLFILE
        F->get_next_line_in_section(0, tempstr);
        sscanf(tempstr, "%d %f", & valin, & val2);
#endif
#ifdef TINYXML
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, valin);
        QueryNumAttributeExt(xmlFile, "Generate_microstructure", countKey++, val2);
#endif
        //fscanf(in,"%d %f",&valin,&val2);
        r [ i ] = valin;
        s [ i ] = ( float ) val2;
        xr [ i ] = ( float ) r [ i ];
    }

    ss = s [ 1 ];
    s2 = ss * ss;
    /* Load up the convolution matrix */
    sdiff = ss - s2;
    for ( i = 0; i < 31; i++ ) {
        iii = i * i;
        for ( j = 0; j < 31; j++ ) {
            jjj = j * j;
            for ( k = 0; k < 31; k++ ) {
                xtmp = ( float ) ( iii + jjj + k * k );
                radius = sqrt(xtmp);
                r1 = ( int ) radius + 1;
                r2 = r1 + 1;
                if ( s [ r1 ] < 0.0 ) {
                    printf("%d and %d %f and %f with xtmp of %f\n", r1, r2, s [ r1 ], s [ r2 ], xtmp);
                    fflush(stdout);
                    exit(1);
                }

                xrad = radius + 1 - r1;
                filval = s [ r1 ] + ( s [ r2 ] - s [ r1 ] ) * xrad;
                filter [ i + 1 ] [ j + 1 ] [ k + 1 ] = ( filval - s2 ) / sdiff;
            }
        }
    }

    /* Now filter the image maintaining periodic boundaries */
    /*fixed for periodic boundaries (small microstructures) - smilauer 4.12.2006*/
    resmax = 0.0;
    resmin = 1.0;
    for ( i = 1; i <= SYSIZE; i++ ) {
        for ( j = 1; j <= SYSIZE; j++ ) {
            for ( k = 1; k <= SYSIZE; k++ ) {
                res [ i ] [ j ] [ k ] = 0.0;
                if ( ( float ) mask [ i ] [ j ] [ k ] == phasein ) {
                    for ( ix = 1; ix <= 31; ix++ ) {
                        i1 = i + ix - 1;
                        while ( i1 < 1 ) { //if(i1<1){i1+=SYSIZE;}
                            i1 += SYSIZE;
                        }

                        while ( i1 > SYSIZE ) { //else if(i1>SYSIZE){
                            i1 -= SYSIZE;
                        }

                        for ( iy = 1; iy <= 31; iy++ ) {
                            j1 = j + iy - 1;
                            while ( j1 < 1 ) { //if(j1<1){j1+=SYSIZE;}
                                j1 += SYSIZE;
                            }

                            while ( j1 > SYSIZE ) { //else if(j1>SYSIZE){
                                j1 -= SYSIZE;
                            }

                            for ( iz = 1; iz <= 31; iz++ ) {
                                k1 = k + iz - 1;
                                while ( k1 < 1 ) { //if(k1<1){k1+=SYSIZE;}
                                    k1 += SYSIZE;
                                }

                                while ( k1 > SYSIZE ) { //else if(k1>SYSIZE){
                                    k1 -= SYSIZE;
                                }

                                res [ i ] [ j ] [ k ] += normm [ i1 ] [ j1 ] [ k1 ] * filter [ ix ] [ iy ] [ iz ];
                            }
                        }
                    }

                    if ( res [ i ] [ j ] [ k ] > resmax ) {
                        resmax = res [ i ] [ j ] [ k ];
                    }

                    if ( res [ i ] [ j ] [ k ] < resmin ) {
                        resmin = res [ i ] [ j ] [ k ];
                    }
                }
            }

#ifdef PRINTF
            printf(".");
#endif
        }

#ifdef PRINTF
        printf("%d out of %d\n", i, SYSIZE);
#endif
    }

#ifdef PRINTF
    printf("\n");
#endif
    /* Now threshold the image */
    sect = ( resmax - resmin ) / 500.;
    for ( i = 1; i <= 500; i++ ) {
        sum [ i ] = 0.0;
    }

    xtot = 0.0;
    for ( i = 1; i <= SYSIZE; i++ ) {
        for ( j = 1; j <= SYSIZE; j++ ) {
            for ( k = 1; k <= SYSIZE; k++ ) {
                if ( ( float ) mask [ i ] [ j ] [ k ] == phasein ) {
                    xtot += 1.0;
                    index = 1 + ( int ) ( ( res [ i ] [ j ] [ k ] - resmin ) / sect );
                    if ( index > 500 ) {
                        index = 500;
                    }

                    sum [ index ] += 1.0;
                }
            }
        }
    }

    /* Determine which bin to choose for correct thresholding */
    sumtot = vcrit = 0.0;
    done = 0;
    for ( i = 1; ( ( i <= 500 ) && ( done == 0 ) ); i++ ) {
        sumtot += sum [ i ] / xtot;
        if ( sumtot > xpt ) {
            ytmp = ( float ) i;
            vcrit = resmin + ( resmax - resmin ) * ( ytmp - 0.5 ) / 500.;
            done = 1;
        }
    }

#ifdef PRINTF
    printf("Critical volume fraction is %f\n", vcrit);
#endif
    ires = 0;

    for ( k = 1; k <= SYSIZE; k++ ) {
        for ( j = 1; j <= SYSIZE; j++ ) {
            for ( i = 1; i <= SYSIZE; i++ ) {
                if ( ( float ) mask [ i ] [ j ] [ k ] == phasein ) {
                    if ( res [ i ] [ j ] [ k ] > vcrit ) {
                        mask [ i ] [ j ] [ k ] = phaseout;
                    }
                }
            }
        }
    }

    dealloc_double_3D(normm, SYSIZE + 1);
    dealloc_double_3D(res, SYSIZE + 1);
    dealloc_double_3D(filter, 32);

    delete [] s;
    delete [] xr;
    delete [] sum;
}


/*disabled sintering*/
void CemhydMatStatus :: distrib3d(void)
{
    int i, j, k, alumflag, alumval, alum2, valin;
    int output_img;
    double volin, volf [ 5 ], surff [ 5 ], rhtest, rdesire;
    char filen [ 80 ];
    //char fileout[80],filecem[80]
    FILE *infile;
    FILE *outfile_img = NULL, *outfile_id = NULL;

    alloc_int_3D(curvature, SYSIZE + 1);

    /* Seed the random number generator */
    //printf("Enter random number seed (negative integer) \n");
    //fscanf(in, "%d",&nseed);
    //printf("%d\n",nseed);
    nseed = iseed;
#ifdef PRINTF
    printf("%d\n", * seed);
#endif

    /* Read in the parameters to use */
    //printf("Enter name of cement microstructure image file\n");
    //fscanf(in, "%s",filen);
    //printf("%s\n",filen);


    /* Set up the correlation filenames
     * printf("Enter name of sil correlation files\n");
     * fscanf(in, "%s",filesil);
     * printf("%s\n",filesil);
     * printf("Enter name of c3s correlation files\n");
     * fscanf(in, "%s",filec3s);
     * printf("%s\n",filec3s);
     * printf("Enter name of c4af correlation files\n");
     * fscanf(in, "%s",filealum);
     * printf("%s\n",filealum);
     */

    alumflag = 1;
    alumval = 4;

    //assume always C4AF
    /*
     * testfile=fopen(filealum,"r");
     * if(testfile==NULL){
     * alumflag=0;
     * sprintf(filealum,"%s",filecem);
     * strcat(filealum,".c3a");
     * alumval=3;
     * }
     * else{
     * fclose(testfile);
     * }
     */
    //printf("Enter name of new cement microstructure image file\n");
    //fscanf(in, "%s",fileout);
    // printf("%s\n",fileout);
    for ( i = 1; i <= 4; i++ ) {
#ifdef CMLFILE
        F->get_next_line_in_section(0, volin);
#endif
#ifdef TINYXML
        if ( i == 1 ) {
            QueryNumAttributeExt(xmlFile, "Generate_microstructure", "C3S_unit_frac", volin);
        }

        if ( i == 2 ) {
            QueryNumAttributeExt(xmlFile, "Generate_microstructure", "C2S_unit_frac", volin);
        }

        if ( i == 3 ) {
            QueryNumAttributeExt(xmlFile, "Generate_microstructure", "C3A_unit_frac", volin);
        }

        if ( i == 4 ) {
            QueryNumAttributeExt(xmlFile, "Generate_microstructure", "C4AF_unit_frac", volin);
        }

#endif
        //fscanf(in, "%f",&volin);
        volf [ i ] = volin;
#ifdef PRINTF
        printf("Volume %f\n", volf [ i ]);
#endif
        //fscanf(in, "%f",&volin);
        surff [ i ] = volin;
#ifdef PRINTF
        printf("Surface %f\n", surff [ i ]);
#endif
    }

    /* Read in the original microstructure image file */
    if ( ( infile = fopen(filen, "r") ) != NULL ) {
        for ( k = 1; k <= SYSIZE; k++ ) {
            for ( j = 1; j <= SYSIZE; j++ ) {
                for ( i = 1; i <= SYSIZE; i++ ) {
                    if ( fscanf(infile, "%d", & valin) != 1 ) {
                        OOFEM_ERROR("CemhydMatStatus :: distrib3d reading error");
                    }

                    mask [ i ] [ j ] [ k ] = valin;
                    curvature [ i ] [ j ] [ k ] = 0;
                }
            }
        }

        fclose(infile);
    } else {
        for ( k = 1; k <= SYSIZE; k++ ) {
            for ( j = 1; j <= SYSIZE; j++ ) {
                for ( i = 1; i <= SYSIZE; i++ ) {
                    mask [ i ] [ j ] [ k ] = cemreal [ i ] [ j ] [ k ];
                }
            }
        }
    }

    stat3d();

    /* First filtering */
    volin = volf [ 1 ] + volf [ 2 ];
    if ( volin < 1.0 ) {
        rand3d(1, alumval, volin);

        /* First sintering */
        stat3d();
        rdesire = ( surff [ 1 ] + surff [ 2 ] ) * ( float ) ( surface [ 1 ] + surface [ alumval ] );
        if ( rdesire != 0.0 ) {
            if ( ( int ) rdesire < surface [ 1 ] ) {
                rhtest = ( 6. / 4. ) * ( float ) ( volume [ 1 ] ) / rdesire;
                //sinter3d(1,alumval,rhtest);
            } else {
                rdesire = ( surff [ 3 ] + surff [ 4 ] ) * ( float ) ( surface [ 1 ] + surface [ alumval ] );
                if ( rdesire != 0.0 ) {
                    rhtest = ( 6. / 4. ) * ( float ) ( volume [ alumval ] ) / rdesire;
                    //sinter3d(alumval,1,rhtest);
                }
            }
        }
    }

    /* Second filtering */
    if ( ( volf [ 1 ] + volf [ 2 ] ) > 0.0 ) {
        volin = volf [ 1 ] / ( volf [ 1 ] + volf [ 2 ] );
        if ( volin < 1.0 ) {
            rand3d(1, 2, volin);

            /* Second sintering */
            stat3d();
            rdesire = ( surff [ 1 ] / ( surff [ 1 ] + surff [ 2 ] ) ) * ( float ) ( surface [ 1 ] + surface [ 2 ] );
            if ( rdesire != 0.0 ) {
                if ( ( int ) rdesire < surface [ 1 ] ) {
                    rhtest = ( 6. / 4. ) * ( float ) ( volume [ 1 ] ) / rdesire;
                    //sinter3d(1,2,rhtest);
                } else {
                    rdesire = ( surff [ 2 ] / ( surff [ 1 ] + surff [ 2 ] ) ) * ( float ) ( surface [ 1 ] + surface [ 2 ] );
                    if ( rdesire != 0.0 ) {
                        rhtest = ( 6. / 4. ) * ( float ) ( volume [ 2 ] ) / rdesire;
                        //sinter3d(2,1,rhtest);
                    }
                }
            }
        }
    }

    /* Third (final) filtering */
    if ( alumval == 4 ) {
        volin = volf [ 4 ] / ( volf [ 4 ] + volf [ 3 ] );
        alum2 = 3;
    } else {
        volin = volf [ 3 ] / ( volf [ 4 ] + volf [ 3 ] );
        alum2 = 4;
    }

    if ( volin < 1.0 ) {
        rand3d(alumval, alum2, volin);

        /* Third (final) sintering */
        stat3d();
        if ( alumval == 4 ) {
            rdesire = ( surff [ 4 ] / ( surff [ 3 ] + surff [ 4 ] ) ) * ( float ) ( surface [ 3 ] + surface [ 4 ] );
            if ( rdesire != 0.0 ) {
                if ( ( int ) rdesire < surface [ 4 ] ) {
                    rhtest = ( 6. / 4. ) * ( float ) ( volume [ 4 ] ) / rdesire;
                    //sinter3d(alumval,alum2,rhtest);
                } else {
                    rdesire = ( surff [ 3 ] / ( surff [ 3 ] + surff [ 4 ] ) ) * ( float ) ( surface [ 3 ] + surface [ 4 ] );
                    if ( rdesire != 0.0 ) {
                        rhtest = ( 6. / 4. ) * ( float ) ( volume [ 3 ] ) / rdesire;
                        //sinter3d(alum2,alumval,rhtest);
                    }
                }
            }
        } else {
            rdesire = ( surff [ 3 ] / ( surff [ 3 ] + surff [ 4 ] ) ) * ( float ) ( surface [ 3 ] + surface [ 4 ] );
            if ( rdesire != 0.0 ) {
                if ( ( int ) rdesire < surface [ 3 ] ) {
                    rhtest = ( 6. / 4. ) * ( float ) ( volume [ 3 ] ) / rdesire;
                    //sinter3d(alumval,alum2,rhtest);
                } else {
                    rdesire = ( surff [ 4 ] / ( surff [ 3 ] + surff [ 4 ] ) ) * ( float ) ( surface [ 3 ] + surface [ 4 ] );
                    if ( rdesire != 0.0 ) {
                        rhtest = ( 6. / 4. ) * ( float ) ( volume [ 4 ] ) / rdesire;
                        //sinter3d(alum2,alumval,rhtest);
                    }
                }
            }
        }
    }

    /* Output final microstructure */

    //check whether to save it into img and id files
#ifdef CMLFILE
    F->get_value(14, ( long & )output_img);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Output_initial_microstructure", 0, output_img);
#endif


    if ( output_img ) {
#ifdef CMLFILE
        F->get_value(15, filen);
#endif
#ifdef TINYXML
        QueryStringAttributeExt(xmlFile, "Output_initial_microstructure_img_file", 0, filen);
#endif
        if ( ( outfile_img = fopen(filen, "w") ) == NULL ) {
            printf("Output img file %s can not be created (file %s, line %d)\n", filen, __FILE__, __LINE__);
            exit(1);
        }

#ifdef CMLFILE
        F->get_value(16, filen);
#endif
#ifdef TINYXML
        QueryStringAttributeExt(xmlFile, "Output_initial_microstructure_id_file", 0, filen);
#endif
        if ( ( outfile_id = fopen(filen, "w") ) == NULL ) {
            printf("Output id file %s can not be created (file %s, line %d)\n", filen, __FILE__, __LINE__);
            exit(1);
        }
    }

    for ( k = 0; k < SYSIZE; k++ ) {
        for ( j = 0; j < SYSIZE; j++ ) {
            for ( i = 0; i < SYSIZE; i++ ) {
                mic [ i ] [ j ] [ k ] = mask [ i + 1 ] [ j + 1 ] [ k + 1 ];
                micorig [ i ] [ j ] [ k ] = mic [ i ] [ j ] [ k ];
                if ( output_img ) {
                    fprintf(outfile_img, "%d\n", mic [ i ] [ j ] [ k ]);
                    fprintf(outfile_id, "%ld\n", micpart [ i ] [ j ] [ k ]);
                }
            }
        }
    }

    if ( output_img ) {
        fclose(outfile_img);
        fclose(outfile_id);
    }

    //deallocate curvature
    dealloc_int_3D(curvature, SYSIZE + 1);

    //deallocate cemreal[][][]
    dealloc_long_3D(cemreal, SYSIZE + 1);
}


/**********************************************************/
/*************DISREALNEW***********************************/
/**********************************************************/

void CemhydMatStatus :: init(void)
{
    int i;
    double slagin, CHperslag;
    FILE *slagfile;

    for ( i = 0; i <= EMPTYP; i++ ) {
        creates [ i ] = 0;
        soluble [ i ] = 0;
        disprob [ i ] = 0.0;
        disbase [ i ] = 0.0;
        pHeffect [ i ] = 0.0;
    }

    /* soluble [x] - flag indicating if phase x is soluble */
    /* disprob [x] - probability of dissolution (relative diss. rate) */
    gypabsprob = 0.01;  /* One sulfate absorbed per 100 CSH units */
    /* Source is from Taylor's Cement Chemistry  for K and Na absorption */
    /* Note that for the first cycle, of the clinker phases only the */
    /* aluminates and gypsum are soluble */
    soluble [ C4AF ] = 1;
    disprob [ C4AF ] = disbase [ C4AF ] = 0.067 / DISBIAS;
    creates [ C4AF ] = POROSITY;
    pHeffect [ C4AF ] = 1.0;
    soluble [ C3S ] = 0;
    disprob [ C3S ] = disbase [ C3S ] = 0.7 / DISBIAS;
    creates [ C3S ] = DIFFCSH;
    pHeffect [ C3S ] = 1.0;
    soluble [ C2S ] = 0;
    disprob [ C2S ] = disbase [ C2S ] = 0.1 / DISBIAS;
    creates [ C2S ] = DIFFCSH;
    pHeffect [ C2S ] = 1.0;
    soluble [ C3A ] = 1;
    /* increased back to 0.4 from 0.25 7/8/99 */
    disprob [ C3A ] = disbase [ C3A ] = 0.4 / DISBIAS;
    creates [ C3A ] = POROSITY;
    pHeffect [ C3A ] = 1.0;
    soluble [ GYPSUM ] = 1;
    /* Changed from 0.05 to 0.015  9/29/98 */
    /* Changed to 0.040 10/15/98 */
    /* back to 0.05 from 0.10 7/8/99 */
    /* from 0.05 to 0.02 4/4/00 */
    /* from 0.02 to 0.025 8/13/01 */
    disprob [ GYPSUM ] = disbase [ GYPSUM ] = 0.025; /*geaendert am 04.04.00, urspr. 0.05*/
    /* dissolved gypsum distributed at random throughout microstructure */
    creates [ GYPSUM ] = POROSITY;
    soluble [ GYPSUMS ] = 1;
    pHeffect [ GYPSUM ] = 0.0;
    /* Changed from 0.05 to 0.015  9/29/98 */
    /* Changed to 0.020 10/15/98 */
    /* and also changed all sulfate based dissolution rates */
    disprob [ GYPSUMS ] = disbase [ GYPSUMS ] = disprob [ GYPSUM ];
    creates [ GYPSUMS ] = POROSITY;
    pHeffect [ GYPSUMS ] = 0.0;
    soluble [ ANHYDRITE ] = 1;
    /* set anhydrite dissolution at 4/5ths of that of gypsum */
    /* Source: Uchikawa et al., CCR, 1984 */
    disprob [ ANHYDRITE ] = disbase [ ANHYDRITE ] = 0.8 * disprob [ GYPSUM ];
    /* dissolved anhydrite distributed at random throughout microstructure */
    creates [ ANHYDRITE ] = POROSITY;
    pHeffect [ ANHYDRITE ] = 0.0;
    soluble [ HEMIHYD ] = 1;
    /* set hemihydrate dissolution at 3 times that of gypsum */
    /* Source: Uchikawa et al., CCR, 1984 */
    /* Changed to 1.5 times that of gypsum 6/1/00 */
    disprob [ HEMIHYD ] = disbase [ HEMIHYD ] = 1.5 * disprob [ GYPSUM ]; /* ge�ndert am 01.06.00, urspr. 3.0 */
    /* dissolved hemihydrate distributed at random throughout microstructure */
    creates [ HEMIHYD ] = POROSITY;
    pHeffect [ HEMIHYD ] = 0.0;
    /* CH soluble to allow for Ostwald ripening of crystals */
    soluble [ CH ] = 1;
    disprob [ CH ] = disbase [ CH ] = 0.5 / DISBIAS;
    creates [ CH ] = DIFFCH;
    pHeffect [ CH ] = 0.0;
    /* CaCO3 is only mildly soluble */
    soluble [ CACO3 ] = 1;
    disprob [ CACO3 ] = disbase [ CACO3 ] = 0.10 / DISBIAS;
    creates [ CACO3 ] = DIFFCACO3;
    pHeffect [ CACO3 ] = 0.0;
    /* Slag is not truly soluble, but use its dissolution probability for reaction probability */
    soluble [ SLAG ] = 0;
    disprob [ SLAG ] = disbase [ SLAG ] = 0.005 / DISBIAS;
    creates [ SLAG ] = 0;
    pHeffect [ SLAG ] = 1.0;
    soluble [ C3AH6 ] = 1;
    disprob [ C3AH6 ] = disbase [ C3AH6 ] = 0.01 / DISBIAS; /* changed from 0.5 to 0.01 06.09.00 */
    creates [ C3AH6 ] = POROSITY;
    pHeffect [ C3AH6 ] = 0.0;
    /* Ettringite is initially insoluble */
    soluble [ ETTR ] = 0;
    /* Changed to 0.008 from 0.020  3/11/99 */
    disprob [ ETTR ] = disbase [ ETTR ] = 0.008 / DISBIAS;
    creates [ ETTR ] = DIFFETTR;
    pHeffect [ ETTR ] = 0.0;
    /* Iron-rich ettringite is always insoluble */
    soluble [ ETTRC4AF ] = 0;
    disprob [ ETTRC4AF ] = disbase [ ETTRC4AF ] = 0.0;
    creates [ ETTRC4AF ] = ETTRC4AF;
    pHeffect [ ETTRC4AF ] = 0.0;
    /* calcium chloride is soluble */
    soluble [ CACL2 ] = 1;
    disprob [ CACL2 ] = disbase [ CACL2 ] = 0.1 / DISBIAS;
    creates [ CACL2 ] = DIFFCACL2;
    pHeffect [ CACL2 ] = 0.0;
    /* aluminosilicate glass is soluble */
    soluble [ ASG ] = 1;
    disprob [ ASG ] = disbase [ ASG ] = 0.2 / DISBIAS;
    creates [ ASG ] = DIFFAS;
    pHeffect [ ASG ] = 1.0;
    /* calcium aluminodisilicate is soluble */
    soluble [ CAS2 ] = 1;
    disprob [ CAS2 ] = disbase [ CAS2 ] = 0.2 / DISBIAS;
    creates [ CAS2 ] = DIFFCAS2;
    pHeffect [ CAS2 ] = 1.0;

    /* establish molar volumes and heats of formation */
    /* molar volumes are in cm^3/mole */
    /* heats of formation are in kJ/mole */
    /* See paper by Fukuhara et al., Cem. & Conc. Res., 11, 407-14, 1981. */
    /* values for Porosity are those of water */
    molarv [ POROSITY ] = 18.068;
    heatf [ POROSITY ] = ( -285.83 );
    waterc [ POROSITY ] = 1.0;
    specgrav [ POROSITY ] = 0.99707;

    molarv [ C3S ] = 71.129;
    heatf [ C3S ] = ( -2927.82 );
    waterc [ C3S ] = 0.0;
    specgrav [ C3S ] = 3.21;
    /* For improvement in chemical shrinkage correspondence */
    /* Changed molar volume of C(1.7)-S-H(4.0) to 108 5/24/95 */
    molarv [ CSH ] = 108.;
    heatf [ CSH ] = ( -3283.0 );
    waterc [ CSH ] = 4.0;
    specgrav [ CSH ] = 2.11;

    molarv [ CH ] = 33.1;
    heatf [ CH ] = ( -986.1 );
    waterc [ CH ] = 1.0;
    specgrav [ CH ] = 2.24;

    //uncomment to have the same values as in Pignat's thesis
    /*
     * molarv[C3S]=72.381;//changed 18.5.05
     * specgrav[C3S]=3.15;//changed 18.5.05
     * molarv[CSH]=113.5;//changed 18.5.05
     * specgrav[CSH]=2.0;//changed 18.5.05
     * molarv[CH]=33.036;//changed 18.5.05
     */

    /* Assume that calcium carbonate is in the calcite form */
    molarv [ CACO3 ] = 36.93;
    waterc [ CACO3 ] = 0.0;
    specgrav [ CACO3 ] = 2.71;
    heatf [ CACO3 ] = ( -1206.92 );

    molarv [ AFMC ] = 261.91;
    waterc [ AFMC ] = 11.0;
    specgrav [ AFMC ] = 2.17;
    /* Need to fill in heat of formation at a later date */
    heatf [ AFMC ] = ( 0.0 );

    molarv [ C2S ] = 52.513;
    heatf [ C2S ] = ( -2311.6 );
    waterc [ C2S ] = 0.0;
    specgrav [ C2S ] = 3.28;

    molarv [ C3A ] = 88.94;
    heatf [ C3A ] = ( -3587.8 );
    waterc [ C3A ] = 0.0;
    specgrav [ C3A ] = 3.038;

    molarv [ GYPSUM ] = 74.21;
    heatf [ GYPSUM ] = ( -2022.6 );
    waterc [ GYPSUM ] = 0.0;
    specgrav [ GYPSUM ] = 2.32;
    molarv [ GYPSUMS ] = 74.21;
    heatf [ GYPSUMS ] = ( -2022.6 );
    waterc [ GYPSUMS ] = 0.0;
    specgrav [ GYPSUMS ] = 2.32;

    molarv [ ANHYDRITE ] = 52.16;
    heatf [ ANHYDRITE ] = ( -1424.6 );
    waterc [ ANHYDRITE ] = 0.0;
    specgrav [ ANHYDRITE ] = 2.61;

    molarv [ HEMIHYD ] = 52.973;
    heatf [ HEMIHYD ] = ( -1574.65 );
    waterc [ HEMIHYD ] = 0.0;
    specgrav [ HEMIHYD ] = 2.74;

    molarv [ C4AF ] = 130.29;
    heatf [ C4AF ] = ( -5090.3 );
    waterc [ C4AF ] = 0.0;
    specgrav [ C4AF ] = 3.73;

    molarv [ C3AH6 ] = 150.12;
    heatf [ C3AH6 ] = ( -5548. );
    waterc [ C3AH6 ] = 6.0;
    specgrav [ C3AH6 ] = 2.52;

    /* Changed molar volume of FH3 to 69.8 (specific gravity of 3.06) 5/23/95 */
    molarv [ FH3 ] = 69.803;
    heatf [ FH3 ] = ( -823.9 );
    waterc [ FH3 ] = 3.0;
    specgrav [ FH3 ] = 3.062;

    molarv [ ETTRC4AF ] = 735.01;
    heatf [ ETTRC4AF ] = ( -17539.0 );
    waterc [ ETTRC4AF ] = 26.0;
    specgrav [ ETTRC4AF ] = 1.7076;
    /* Changed molar volue of ettringite to 735 (spec. gr.=1.7076)  5/24/95 */
    molarv [ ETTR ] = 735.01;
    heatf [ ETTR ] = ( -17539.0 );
    waterc [ ETTR ] = 26.0;
    specgrav [ ETTR ] = 1.7076;

    molarv [ AFM ] = 312.82;
    heatf [ AFM ] = ( -8778.0 );
    /* Each mole of AFM which forms requires 12 moles of water, */
    /* two of which are supplied by gypsum in forming ETTR  */
    /* leaving 10 moles to be incorporated from free water */
    waterc [ AFM ] = 10.0;
    specgrav [ AFM ] = 1.99;

    molarv [ CACL2 ] = 51.62;
    heatf [ CACL2 ] = ( -795.8 );
    waterc [ CACL2 ] = 0;
    specgrav [ CACL2 ] = 2.15;

    molarv [ FREIDEL ] = 296.662;
    /* No data available for heat of formation */
    heatf [ FREIDEL ] = ( 0.0 );
    /* 10 moles of H2O per mole of Freidel's salt */
    waterc [ FREIDEL ] = 10.0;
    specgrav [ FREIDEL ] = 1.892;

    /* Basic reaction is 2CH + ASG + 6H --> C2ASH8 (Stratlingite) */
    molarv [ ASG ] = 49.9;
    /* No data available for heat of formation */
    heatf [ ASG ] = 0.0;
    waterc [ ASG ] = 0.0;
    specgrav [ ASG ] = 3.247;

    molarv [ CAS2 ] = 100.62;
    /* No data available for heat of formation */
    heatf [ CAS2 ] = 0.0;
    waterc [ CAS2 ] = 0.0;
    specgrav [ CAS2 ] = 2.77;

    molarv [ STRAT ] = 215.63;
    /* No data available for heat of formation */
    heatf [ STRAT ] = 0.0;
    /* 8 moles of water per mole of stratlingite */
    waterc [ STRAT ] = 8.0;
    specgrav [ STRAT ] = 1.94;

    molarv [ POZZ ] = 27.0;
    /* Use heat of formation of SiO2 (quartz) for unreacted pozzolan */
    /* Source- Handbook of Chemistry and Physics */
    heatf [ POZZ ] = -907.5;
    waterc [ POZZ ] = 0.0;
    specgrav [ POZZ ] = 2.22;

    /* Data for Pozzolanic CSH based on work of Atlassi, DeLarrard, */
    /* and Jensen */
    /* gives a chemical shrinkage of 0.2 g H2O/g CSF */
    /* heat of formation estimated based on heat release of */
    /* 780 J/g Condensed Silica Fume */
    /* Changed stoichiometry to be C(1.1)SH(3.9) to see effect on */
    /* results 1/22/97 */
    /* MW is 191.8 g/mole */
    /* Changed molar volume to 101.81  3/10/97 */
    molarv [ POZZCSH ] = 101.81;
    waterc [ POZZCSH ] = 3.9;
    specgrav [ POZZCSH ] = 1.884;
    heatf [ POZZCSH ] = ( -2.1 );

    /* Assume inert filler has same specific gravity and molar volume as SiO2 */
    molarv [ INERT ] = 27.0;
    heatf [ INERT ] = 0.0;
    waterc [ INERT ] = 0.0;
    specgrav [ INERT ] = 2.2;

    molarv [ ABSGYP ] = 74.21;
    heatf [ ABSGYP ] = ( -2022.6 );
    waterc [ ABSGYP ] = 0.0;
    specgrav [ ABSGYP ] = 2.32;

    molarv [ EMPTYP ] = 18.068;
    heatf [ EMPTYP ] = ( -285.83 );
    waterc [ EMPTYP ] = 0.0;
    specgrav [ EMPTYP ] = 0.99707;

    /* Read in values for alkali characteristics and */
    /* convert them to fractions from percentages */
    // First define non-zero values due to divisions
    //totsodium=0.2/100;
    //totpotassium=0.5/100;
    //rssodium=0.03/100;
    //rspotassium=0.25/100;
    specgrav [ SLAG ] = 2.87;
    specgrav [ SLAGCSH ] = 2.35;
    molarv [ SLAG ] = 945.72;
    molarv [ SLAGCSH ] = 1717.53;
    slagcasi = 1.433;
    slaghydcasi = 1.35;
    siperslag = 15.0;
    waterc [ SLAGCSH ] = 5.533 * siperslag;
    slagc3a = 1.0;
    slagreact = 1.0;

#ifdef CMLFILE
    F->get_value(19, totsodium);
    F->get_value(20, totpotassium);
    F->get_value(21, rssodium);
    F->get_value(22, rspotassium);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Total_sodium", 0, totsodium);
    QueryNumAttributeExt(xmlFile, "Total_potassium", 0, totpotassium);
    QueryNumAttributeExt(xmlFile, "Readily_soluble_sodium", 0, rssodium);
    QueryNumAttributeExt(xmlFile, "Readily_soluble_potassium", 0, rspotassium);
#endif

    /*fscanf(alkalifile,"%f",&totsodium);
     * fscanf(alkalifile,"%f",&totpotassium);
     * fscanf(alkalifile,"%f",&rssodium);
     * fscanf(alkalifile,"%f",&rspotassium);
     */

    totsodium /= 100.;
    totpotassium /= 100.;
    rssodium /= 100.;
    rspotassium /= 100.;


    if ( pHactive ) {
        /* Read in values for slag characteristics */

        if ( ( slagfile = fopen("slagchar.dat", "r") ) == NULL ) {
            printf("Slag file can not be opened\n");
            exit(1);
        }


        if ( fscanf(slagfile, "%lf", & slagin) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        if ( fscanf(slagfile, "%lf", & slagin) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        if ( fscanf(slagfile, "%lf", & slagin) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        specgrav [ SLAG ] = slagin;
        if ( fscanf(slagfile, "%lf", & slagin) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        specgrav [ SLAGCSH ] = slagin;
        if ( fscanf(slagfile, "%lf", & slagin) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        molarv [ SLAG ] = slagin;
        if ( fscanf(slagfile, "%lf", & slagin) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        molarv [ SLAGCSH ] = slagin;
        if ( fscanf(slagfile, "%lf", & slagcasi) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        if ( fscanf(slagfile, "%lf", & slaghydcasi) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        if ( fscanf(slagfile, "%lf", & siperslag) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        if ( fscanf(slagfile, "%lf", & slagin) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        waterc [ SLAGCSH ] = slagin * siperslag;
        if ( fscanf(slagfile, "%lf", & slagc3a) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        if ( fscanf(slagfile, "%lf", & slagreact) != 1 ) {
            OOFEM_ERROR("CemhydMatStatus::init: slagfile reading error");
        }

        waterc [ SLAG ] = 0.0;
        heatf [ SLAG ] = 0.0;
        heatf [ SLAGCSH ] = 0.0;
        fclose(slagfile);
    }

    //use default values
    specgrav [ SLAG ] = 2.87;
    specgrav [ SLAGCSH ] = 2.35;
    molarv [ SLAG ] = 945.72;
    molarv [ SLAGCSH ] = 1717.53;
    slagcasi = 1.4333;
    slaghydcasi = 1.35;
    siperslag = 15.0;
    waterc [ SLAGCSH ] = 5.533 * siperslag;
    slagc3a = 1.0;
    slagreact = 1.0;

    /* Compute slag probabilities as defined above */
    CHperslag = siperslag * ( slaghydcasi - slagcasi ) + 3. * slagc3a;
    if ( CHperslag < 0.0 ) {
        CHperslag = 0.0;
    }

    p2slag = ( molarv [ SLAG ] + molarv [ CH ] * CHperslag + molarv [ POROSITY ] * ( waterc [ SLAGCSH ] - CHperslag + waterc [ C3AH6 ] * slagc3a ) - molarv [ SLAGCSH ] - molarv [ C3AH6 ] * slagc3a ) / molarv [ SLAG ];
    p1slag = 1.0 - p2slag;
    p3slag = ( molarv [ SLAGCSH ] / molarv [ SLAG ] ) - p1slag;
    p4slag = CHperslag * molarv [ CH ] / molarv [ SLAG ];
    p5slag = slagc3a * molarv [ C3A ] / molarv [ SLAG ];
    if ( p5slag > 1.0 ) {
        p5slag = 1.0;
        printf("Error in range of C3A/slag value...  reset to 1.0 \n");
    }
}

/* routine to check if a pixel located at (xck,yck,zck) is on an edge */
/* (in contact with pore space) in 3-D system */
/* Called by passone */
/* Calls no other routines */
int CemhydMatStatus :: chckedge(int xck, int yck, int zck) {
    int edgeback, x2, y2, z2;
    int ip;

    edgeback = 0;

    /* Check all six neighboring pixels */
    /* with periodic boundary conditions */
    for ( ip = 0; ( ( ip < NEIGHBORS ) && ( edgeback == 0 ) ); ip++ ) {
        x2 = xck + xoff [ ip ];
        y2 = yck + yoff [ ip ];
        z2 = zck + zoff [ ip ];
        if ( x2 >= SYSIZE ) {
            x2 = 0;
        }

        if ( x2 < 0 ) {
            x2 = SYSIZEM1;
        }

        if ( y2 >= SYSIZE ) {
            y2 = 0;
        }

        if ( y2 < 0 ) {
            y2 = SYSIZEM1;
        }

        if ( z2 >= SYSIZE ) {
            z2 = 0;
        }

        if ( z2 < 0 ) {
            z2 = SYSIZEM1;
        }

        if ( mic [ x2 ] [ y2 ] [ z2 ] == POROSITY ) {
            edgeback = 1;
        }
    }

    return ( edgeback );
}

/* routine for first pass through microstructure during dissolution */
/* low and high indicate phase ID range to check for surface sites */
/* Called by dissolve */
/* Calls chckedge */
void CemhydMatStatus :: passone(int low, int high, int cycid, int cshexflag) {
    int i, xid, yid, zid, phid, edgef, phread, cshcyc;

    /* gypready used to determine if any soluble gypsum remains */
    if ( ( low <= GYPSUM ) && ( GYPSUM <= high ) ) {
        gypready = 0;
    }

    /* Zero out count for the relevant phases */
    for ( i = low; i <= high; i++ ) {
        count [ i ] = 0;
    }

    /* Scan the entire 3-D microstructure */
    for ( xid = 0; xid < SYSIZE; xid++ ) {
        for ( yid = 0; yid < SYSIZE; yid++ ) {
            for ( zid = 0; zid < SYSIZE; zid++ ) {
                phread = mic [ xid ] [ yid ] [ zid ];
                /* Update heat data and water consumed for solid CSH */
                if ( ( cshexflag == 1 ) && ( phread == CSH ) ) {
                    cshcyc = cshage [ xid ] [ yid ] [ zid ];
                    heatsum += heatf [ CSH ] / molarvcsh [ cshcyc ];
                    molesh2o += watercsh [ cshcyc ] / molarvcsh [ cshcyc ];
                }

                /* Identify phase and update count */
                phid = 60;
                for ( i = low; ( ( i <= high ) && ( phid == 60 ) ); i++ ) {
                    if ( mic [ xid ] [ yid ] [ zid ] == i ) {
                        phid = i;
                        /* Update count for this phase */
                        count [ i ] += 1;
                        if ( ( i == GYPSUM ) || ( i == GYPSUMS ) ) {
                            gypready += 1;
                        }

                        /* If first cycle, then accumulate initial counts */
                        if ( cycid == 1 ) { //fixed (ncyc cancelled)
                            if ( i == POROSITY ) {
                                porinit += 1;
                            }
                            /* Ordered in terms of likely volume fractions */
                            /* (largest to smallest) to speed execution */
                            else if ( i == C3S ) {
                                c3sinit += 1;
                            } else if ( i == C2S ) {
                                c2sinit += 1;
                            } else if ( i == C3A ) {
                                c3ainit += 1;
                            } else if ( i == C4AF ) {
                                c4afinit += 1;
                            } else if ( i == GYPSUM ) {
                                ncsbar += 1;
                            } else if ( i == GYPSUMS ) {
                                ncsbar += 1;
                            } else if ( i == ANHYDRITE ) {
                                anhinit += 1;
                            } else if ( i == HEMIHYD ) {
                                heminit += 1;
                            } else if ( i == POZZ ) {
                                nfill += 1;
                            } else if ( i == SLAG ) {
                                slaginit += 1;
                            } else if ( i == ETTR ) {
                                netbar += 1;
                            } else if ( i == ETTRC4AF ) {
                                netbar += 1;
                            }
                        }
                    }
                }

                if ( phid != 60 ) {
                    /* If phase is soluble, see if it is in contact with porosity */
                    if ( ( cycid != 0 ) && ( soluble [ phid ] == 1 ) ) {
                        edgef = chckedge(xid, yid, zid);
                        if ( edgef == 1 ) {
                            /* Surface eligible species has an ID OFFSET greater than its original value */
                            mic [ xid ] [ yid ] [ zid ] += OFFSET;
                        }
                    }
                }
            }

            /* end of zid */
        }

        /* end of yid */
    }

    /* end of xid */
}

/* routine to locate a diffusing CSH species near dissolution source */
/* at (xcur,ycur,zcur) */
/* Called by dissolve */
/* Calls no other routines */
int CemhydMatStatus :: loccsh(int xcur, int ycur, int zcur, int extent) {
    int effort, tries, xmod, ymod, zmod;
    struct ants *antnew;

    effort = 0; /* effort indicates if appropriate location found */
    tries = 0;
    /* 500 tries in immediate vicinity */
    while ( ( effort == 0 ) && ( tries < 500 ) ) {
        tries += 1;
        xmod = ( -extent ) + ( int ) ( ( 2. * ( float ) extent + 1. ) * ran1(seed) );
        ymod = ( -extent ) + ( int ) ( ( 2. * ( float ) extent + 1. ) * ran1(seed) );
        zmod = ( -extent ) + ( int ) ( ( 2. * ( float ) extent + 1. ) * ran1(seed) );
        if ( xmod > extent ) {
            xmod = extent;
        }

        if ( ymod > extent ) {
            ymod = extent;
        }

        if ( zmod > extent ) {
            zmod = extent;
        }

        xmod += xcur;
        ymod += ycur;
        zmod += zcur;
        /* Periodic boundaries */
        if ( xmod >= SYSIZE ) {
            xmod -= SYSIZE;
        } else if ( xmod < 0 ) {
            xmod += SYSIZE;
        }

        if ( zmod >= SYSIZE ) {
            zmod -= SYSIZE;
        } else if ( zmod < 0 ) {
            zmod += SYSIZE;
        }

        if ( ymod < 0 ) {
            ymod += SYSIZE;
        } else if ( ymod >= SYSIZE ) {
            ymod -= SYSIZE;
        }

        if ( mic [ xmod ] [ ymod ] [ zmod ] == POROSITY ) {
            effort = 1;
            mic [ xmod ] [ ymod ] [ zmod ] = DIFFCSH;
            nmade += 1;
            ngoing += 1;
            /* Add CSH diffusing species to the linked list */
            antnew = ( struct ants * ) malloc( sizeof( struct ants ) );
            antnew->x = xmod;
            antnew->y = ymod;
            antnew->z = zmod;
            antnew->id = DIFFCSH;
            antnew->cycbirth = cyccnt;
            /* Now connect this ant structure to end of linked list */
            antnew->prevant = tailant;
            tailant->nextant = antnew;
            antnew->nextant = NULL;
            tailant = antnew;
        }
    }

    return ( effort );
}

/* routine to count number of pore pixels in a cube of size boxsize */
/* centered at (qx,qy,qz) */
/* Called by makeinert */
/* Calls no other routines */
int CemhydMatStatus :: countbox(int boxsize, int qx, int qy, int qz) {
    int nfound, ix, iy, iz, qxlo, qxhi, qylo, qyhi, qzlo, qzhi;
    int hx, hy, hz, boxhalf;

    boxhalf = boxsize / 2;
    nfound = 0;
    qxlo = qx - boxhalf;
    qxhi = qx + boxhalf;
    qylo = qy - boxhalf;
    qyhi = qy + boxhalf;
    qzlo = qz - boxhalf;
    qzhi = qz + boxhalf;
    /* Count the number of requisite pixels in the 3-D cube box */
    /* using periodic boundaries */
    for ( ix = qxlo; ix <= qxhi; ix++ ) {
        hx = ix;
        if ( hx < 0 ) {
            hx += SYSIZE;
        } else if ( hx >= SYSIZE ) {
            hx -= SYSIZE;
        }

        for ( iy = qylo; iy <= qyhi; iy++ ) {
            hy = iy;
            if ( hy < 0 ) {
                hy += SYSIZE;
            } else if ( hy >= SYSIZE ) {
                hy -= SYSIZE;
            }

            for ( iz = qzlo; iz <= qzhi; iz++ ) {
                hz = iz;
                if ( hz < 0 ) {
                    hz += SYSIZE;
                } else if ( hz >= SYSIZE ) {
                    hz -= SYSIZE;
                }

                /* Count if porosity, diffusing species, or empty porosity */
                if ( ( mic [ hx ] [ hy ] [ hz ] < C3S ) || ( mic [ hx ] [ hy ] [ hz ] > ABSGYP ) ) {
                    nfound += 1;
                }
            }
        }
    }

    return ( nfound );
}

/* routine to count number of special pixels in a cube of size boxsize */
/* centered at (qx,qy,qz) */
/* special pixels are those not belonging to one of the cement clinker, */
/* calcium sulfate, or pozzolanic mineral admixture phases */
/* Called by addrand */
/* Calls no other routines */
int CemhydMatStatus :: countboxc(int boxsize, int qx, int qy, int qz) {
    int nfound, ix, iy, iz, qxlo, qxhi, qylo, qyhi, qzlo, qzhi;
    int hx, hy, hz, boxhalf;

    boxhalf = boxsize / 2;
    nfound = 0;
    qxlo = qx - boxhalf;
    qxhi = qx + boxhalf;
    qylo = qy - boxhalf;
    qyhi = qy + boxhalf;
    qzlo = qz - boxhalf;
    qzhi = qz + boxhalf;
    /* Count the number of requisite pixels in the 3-D cube box */
    /* using periodic boundaries */
    for ( ix = qxlo; ix <= qxhi; ix++ ) {
        hx = ix;
        if ( hx < 0 ) {
            hx += SYSIZE;
        } else if ( hx >= SYSIZE ) {
            hx -= SYSIZE;
        }

        for ( iy = qylo; iy <= qyhi; iy++ ) {
            hy = iy;
            if ( hy < 0 ) {
                hy += SYSIZE;
            } else if ( hy >= SYSIZE ) {
                hy -= SYSIZE;
            }

            for ( iz = qzlo; iz <= qzhi; iz++ ) {
                hz = iz;
                if ( hz < 0 ) {
                    hz += SYSIZE;
                } else if ( hz >= SYSIZE ) {
                    hz -= SYSIZE;
                }

                /* Count if not cement clinker */
                if ( ( mic [ hx ] [ hy ] [ hz ] < C3S ) || ( mic [ hx ] [ hy ] [ hz ] > POZZ ) ) {
                    nfound += 1;
                }
            }
        }
    }

    return ( nfound );
}

/* routine to create ndesire pixels of empty pore space to simulate */
/* self-desiccation, or external drying */
/* Called by dissolve */
/* Calls countbox */
/* each ....togo structure contains information only at x,y,z*/
/* headtogo is the first voxel where the water is taken from*/


void CemhydMatStatus :: makeinert(long int ndesire) {
    struct togo *headtogo, *tailtogo, *newtogo, *lasttogo, *onetogo;
    long int idesire;
    int px, py, pz, placed, cntpore, cntmax;

    /* First allocate the first element of the linked list */
    headtogo = ( struct togo * ) malloc( sizeof( struct togo ) );
    headtogo->x = headtogo->y = headtogo->z = ( -1 );
    headtogo->npore = 0;
    headtogo->nexttogo = NULL;
    headtogo->prevtogo = NULL;
    tailtogo = headtogo;
    cntmax = 0;

#ifdef PRINTF
    printf("In makeinert with %ld needed elements \n", ndesire);
#endif

    fflush(stdout);
    /* Add needed number of elements to the end of the list */
    for ( idesire = 2; idesire <= ndesire; idesire++ ) {
        newtogo = ( struct togo * ) malloc( sizeof( struct togo ) );
        newtogo->npore = 0;
        newtogo->x = newtogo->y = newtogo->z = ( -1 );
        tailtogo->nexttogo = newtogo;
        newtogo->prevtogo = tailtogo;
        tailtogo = newtogo;
    }

    /* Now scan the microstructure and rank the sites */
    for ( px = 0; px < SYSIZE; px++ ) {
        for ( py = 0; py < SYSIZE; py++ ) {
            for ( pz = 0; pz < SYSIZE; pz++ ) {
                if ( mic [ px ] [ py ] [ pz ] == POROSITY ) {
                    cntpore = countbox(cubesize, px, py, pz);
                    if ( cntpore > cntmax ) {
                        cntmax = cntpore;
                    }

                    /* Store this site value at appropriate place in */
                    /* sorted linked list */
                    if ( cntpore > ( tailtogo->npore ) ) {
                        placed = 0;
                        lasttogo = tailtogo;
                        while ( placed == 0 ) {
                            newtogo = lasttogo->prevtogo;
                            if ( newtogo == NULL ) { //if at the beginning of the list (tailtogo->prevtogo==NULL)
                                placed = 2;
                            } else {
                                if ( cntpore <= ( newtogo->npore ) ) {
                                    placed = 1;
                                }
                            }

                            if ( placed == 0 ) {
                                lasttogo = newtogo;
                            }
                        }

                        onetogo = ( struct togo * ) malloc( sizeof( struct togo ) );
                        onetogo->x = px;
                        onetogo->y = py;
                        onetogo->z = pz;
                        onetogo->npore = cntpore;
                        /* Insertion at the head of the list */
                        if ( placed == 2 ) {
                            onetogo->prevtogo = NULL;
                            onetogo->nexttogo = headtogo;
                            headtogo->prevtogo = onetogo;
                            headtogo = onetogo;
                        }

                        if ( placed == 1 ) {
                            onetogo->nexttogo = lasttogo;
                            onetogo->prevtogo = newtogo;
                            lasttogo->prevtogo = onetogo;
                            newtogo->nexttogo = onetogo;
                        }

                        /* Eliminate the last element */
                        lasttogo = tailtogo;
                        tailtogo = tailtogo->prevtogo;
                        tailtogo->nexttogo = NULL;
                        //printf("DEMEM lasttogo1 addr %p\n", lasttogo);
                        free(lasttogo);
                    }
                }
            }
        }
    }

    /* Now remove the sites */
    /* starting at the head of the list */
    /* and deallocate all of the used memory */
    for ( idesire = 1; idesire <= ndesire; idesire++ ) {
        px = headtogo->x;
        py = headtogo->y;
        pz = headtogo->z;
        if ( px != ( -1 ) ) {
            mic [ px ] [ py ] [ pz ] = EMPTYP;
            count [ POROSITY ] -= 1;
            count [ EMPTYP ] += 1;
        }

        lasttogo = headtogo;
        headtogo = headtogo->nexttogo;
        //printf("DEMEM lasttogo2 addr %p\n", (const void *)lasttogo);
        free(lasttogo);
    }

    /* If only small cubes of porosity were found, then adjust */
    /* cubesize to have a more efficient search in the future */
    if ( cubesize > CUBEMIN ) {
        if ( ( 2 * cntmax ) < ( cubesize * cubesize * cubesize ) ) {
            cubesize -= 2;
        }
    }
}


/* routine to add extra SLAG CSH when SLAG reacts */
/* SLAG located at (xpres,ypres,zpres) */
/* Called by dissolve */
/* Calls moveone and edgecnt */
void CemhydMatStatus :: extslagcsh(int xpres, int ypres, int zpres) {
    int check, sump, xchr, ychr, zchr, fchr, i1, action, numnear;
    long int tries;
    int mstest = 0, mstest2 = 0;

    /* first try 6 neighboring locations until      */
    /*    a) successful                */
    /*    b) all 6 sites are tried or             */
    /*    c) 100 tries are made           */
    /* try to grow slag C-S-H as plates */
    fchr = 0;
    sump = 1;
    for ( i1 = 1; ( ( i1 <= 100 ) && ( fchr == 0 ) && ( sump != 30030 ) ); i1++ ) {
        /* determine location of neighbor (using periodic boundaries) */
        xchr = xpres;
        ychr = ypres;
        zchr = zpres;
        action = 0;
        sump *= moveone(& xchr, & ychr, & zchr, & action, sump);
        if ( action == 0 ) {
            printf("Error in value of action in extpozz \n");
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];
        /* Determine the direction of the neighbor selected and */
        /* the plates possible for growth */
        if ( xchr != xpres ) {
            mstest = 1;
            mstest2 = 2;
        }

        if ( ychr != ypres ) {
            mstest = 2;
            mstest2 = 3;
        }

        if ( zchr != zpres ) {
            mstest = 3;
            mstest2 = 1;
        }

        /* if neighbor is porosity, locate the SLAG CSH there */
        if ( check == POROSITY ) {
            if ( ( faces [ xpres ] [ ypres ] [ zpres ] == 0 ) || ( mstest == faces [ xpres ] [ ypres ] [ zpres ] ) || ( mstest2 == faces [ xpres ] [ ypres ] [ zpres ] ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = SLAGCSH;
                faces [ xchr ] [ ychr ] [ zchr ] = faces [ xpres ] [ ypres ] [ zpres ];
                count [ SLAGCSH ] += 1;
                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }

    /* if no neighbor available, locate SLAG CSH at random location */
    /* in pore space */
    tries = 0;
    while ( fchr == 0 ) {
        tries += 1;
        /* generate a random location in the 3-D system */
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];
        /* if location is porosity, locate the extra SLAG CSH there */
        if ( check == POROSITY ) {
            numnear = edgecnt(xchr, ychr, zchr, SLAG, CSH, SLAGCSH);
            /* Be sure that one neighboring species is CSH or */
            /* SLAG material */
            if ( ( tries > 5000 ) || ( numnear < 26 ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = SLAGCSH;
                count [ SLAGCSH ] += 1;
                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }
}

/* routine to implement a cycle of dissolution */
/* Called by main program */
/* Calls passone, loccsh, and makeinert */
void CemhydMatStatus :: dissolve(int cycle) {
    int nc3aext, ncshext, nchext, ngypext, nanhext, plok;
    int nsum5, nsum4, nsum3, nsum2, nhemext, nsum6, nc4aext;
    int phid, phnew, plnew, cread;
    int i, xloop, yloop, zloop, xc, yc;
    int zc, cycnew;
    long int ctest;
    long int tries;
    int placed, cshrand, maxsulfate, msface;
    long int ncshgo, nsurf, suminit;
    long int xext, nhgd, npchext, nslagc3a = 0;
    float pdis, plfh3, fchext, fc3aext, fanhext, mass_now, mass_fa_now, tot_mass, heatfill;
    float dfact, dfact1, molesdh2o, h2oinit, heat4, fhemext, fc4aext;
    float pconvert, pc3scsh, pc2scsh, calcx, calcy, calcz, tdisfact;
    float frafm, frettr, frhyg, frtot, mc3ar, mc4ar, p3init;
    struct ants *antadd;

    /* Initialize variables */
    nmade = 0;
    npchext = ncshgo = cshrand = 0; /* counter for number of csh diffusing species to */
    /* be located at random locations in microstructure */
    heat_old = heat_new; /* new and old values for heat released */

    /* Initialize dissolution and phase counters */
    nsurf = 0;
    for ( i = 0; i <= EMPTYP; i++ ) {
        discount [ i ] = 0;
        count [ i ] = 0;
    }

    /* Pass one- highlight all edge points which are soluble */
    soluble [ C3AH6 ] = 0;
    heatsum = 0.0;
    molesh2o = 0.0;
    passone(0, EMPTYP, cycle, 1);
#ifdef PRINTF
    printf("Returned from passone \n");
#endif
    fflush(stdout);

    sulf_solid = count [ GYPSUM ] + count [ GYPSUMS ] + count [ HEMIHYD ] + count [ ANHYDRITE ];
    /* If first cycle, then determine all mixture proportions based */
    /* on user input and original microstructure */
    if ( cycle == 1 ) {
        /* Mass of cement in system */
        cemmass = ( specgrav [ C3S ] * ( float ) count [ C3S ] + specgrav [ C2S ] *
                    ( float ) count [ C2S ] + specgrav [ C3A ] * ( float ) count [ C3A ] +
                    specgrav [ C4AF ] * ( float ) count [ C4AF ] );
        /*    +specgrav[GYPSUM]*
         * (float)count[GYPSUM]+specgrav[ANHYDRITE]*(float)
         * count[ANHYDRITE]+specgrav[HEMIHYD]*(float)count[HEMIHYD]); */
        cemmasswgyp = ( specgrav [ C3S ] * ( float ) count [ C3S ] + specgrav [ C2S ] *
                        ( float ) count [ C2S ] + specgrav [ C3A ] * ( float ) count [ C3A ] +
                        specgrav [ C4AF ] * ( float ) count [ C4AF ] + specgrav [ GYPSUM ] *
                        ( float ) count [ GYPSUM ] + specgrav [ ANHYDRITE ] * ( float )
                        count [ ANHYDRITE ] + specgrav [ HEMIHYD ] * ( float ) count [ HEMIHYD ] );
        flyashmass = ( specgrav [ ASG ] * ( float ) count [ ASG ] + specgrav [ CAS2 ] *
                       ( float ) count [ CAS2 ] + specgrav [ POZZ ] * ( float ) count [ POZZ ] );
        CH_mass = specgrav [ CH ] * ( float ) count [ CH ];
        /* Total mass in system neglecting single aggregate */
        tot_mass = cemmass + ( float ) count [ POROSITY ] + specgrav [ INERT ] *
                   ( float ) count [ INERT ] + specgrav [ CACL2 ] * ( float ) count [ CACL2 ] +
                   specgrav [ ASG ] * ( float ) count [ ASG ] +
                   specgrav [ SLAG ] * ( float ) count [ SLAG ] +
                   specgrav [ HEMIHYD ] * ( float ) count [ HEMIHYD ] +
                   specgrav [ ANHYDRITE ] * ( float ) count [ ANHYDRITE ] +
                   specgrav [ CAS2 ] * ( float ) count [ CAS2 ] +
                   specgrav [ CSH ] * ( float ) count [ CSH ] +
                   specgrav [ GYPSUM ] * ( float ) count [ GYPSUM ] +
                   specgrav [ ANHYDRITE ] * ( float ) count [ ANHYDRITE ] +
                   specgrav [ HEMIHYD ] * ( float ) count [ HEMIHYD ] +
                   specgrav [ GYPSUMS ] * ( float ) count [ GYPSUMS ] +
                   specgrav [ POZZ ] * ( float ) count [ POZZ ] + CH_mass;
        /* water-to-cement ratio */
        if ( cemmass != 0.0 ) {
            w_to_c = ( float ) count [ POROSITY ] / ( cemmass +
                                                      specgrav [ GYPSUM ] * ( float ) count [ GYPSUM ] +
                                                      specgrav [ ANHYDRITE ] * ( float ) count [ ANHYDRITE ] +
                                                      specgrav [ HEMIHYD ] * ( float ) count [ HEMIHYD ] );
        } else {
            w_to_c = 0.0;
        }

        /* totfract is the total cement volume count including calcium sulfates */
        /* fractwithfill is the total count of cement and solid fillers and */
        /*    mineral admixtures   DPB- 10/04 */
        totfract = count [ C3S ] + count [ C2S ];
        totfract += ( count [ C3A ] + count [ C4AF ] );
        totfract += ( count [ GYPSUM ] + count [ ANHYDRITE ] + count [ HEMIHYD ] );
        fractwithfill = totfract + count [ CACO3 ] + count [ SLAG ] + count [ INERT ];
        fractwithfill += count [ POZZ ] + count [ CAS2 ] + count [ ASG ] + count [ CACL2 ];
        totfract /= ( float ) SYSIZE;
        totfract /= ( float ) SYSIZE;
        totfract /= ( float ) SYSIZE;
        fractwithfill /= ( float ) SYSIZE;
        fractwithfill /= ( float ) SYSIZE;
        fractwithfill /= ( float ) SYSIZE;
        /* Adjust masses for presence of aggregates in concrete */
        mass_water = ( 1. - mass_agg ) * ( float ) count [ POROSITY ] / tot_mass;
        mass_CH = ( 1. - mass_agg ) * CH_mass / tot_mass;
        /* pozzolan-to-cement ratio */
        if ( cemmass != 0.0 ) {
            s_to_c = ( float ) ( count [ INERT ] * specgrav [ INERT ] +
                                 count [ CACL2 ] * specgrav [ CACL2 ] + count [ ASG ] * specgrav [ ASG ] +
                                 count [ CAS2 ] * specgrav [ CAS2 ] +
                                 count [ SLAG ] * specgrav [ SLAG ] +
                                 count [ POZZ ] * specgrav [ POZZ ] ) / cemmass;
        } else {
            s_to_c = 0.0;
        }

        /* Conversion factor to kJ/kg for heat produced */
        if ( cemmass != 0.0 ) {
            heatfill = ( float ) ( count [ INERT ] + count [ SLAG ] + count [ POZZ ] + count [ CACL2 ] + count [ ASG ] + count [ CAS2 ] ) / cemmass;       //units pixels^3/g, inverse to the density
        } else {
            heatfill = 0.0;
        }

        if ( w_to_c > 0.01 ) {
            heat_cf = 1000. / SYSIZE_POW3 * ( 0.3125 + w_to_c + heatfill ); //heat conversion factor,1/ro_c=0.3125, fixed
        } else {
            /* Need volume per 1 gram of silica fume */
            heat_cf = 1000. / SYSIZE_POW3 * ( ( 1. / specgrav [ POZZ ] ) + ( float ) ( count [ POROSITY ] + count [ CH ] + count [ INERT ] ) / ( specgrav [ POZZ ] * ( float ) count [ POZZ ] ) );       //fixed
        }

        mass_fill_pozz = ( 1. - mass_agg ) * ( float ) ( count [ POZZ ] * specgrav [ POZZ ] ) / tot_mass;
        mass_fill = ( 1. - mass_agg ) * ( float ) ( count [ INERT ] * specgrav [ INERT ] +
                                                    count [ ASG ] * specgrav [ ASG ] + count [ SLAG ] * specgrav [ SLAG ] +
                                                    count [ CAS2 ] * specgrav [ CAS2 ] + count [ POZZ ] * specgrav [ POZZ ] +
                                                    count [ CACL2 ] * specgrav [ CACL2 ] ) / tot_mass;
#ifdef PRINTF
        printf("Calculated w/c is %.4f\n", w_to_c);
        printf("Calculated s/c is %.4f \n", s_to_c);
        printf("Calculated heat conversion factor is %f \n", heat_cf);
        printf("Calculated mass fractions of water and filler are %.4f  and %.4f \n",
               mass_water, mass_fill);
#endif
    }

    molesdh2o = 0.0;
    alpha = 0.0;  /* degree of hydration of clinker minerals*/
    /* heat4 contains measured heat release for C4AF hydration from  */
    /* Fukuhara et al., Cem. and Conc. Res. article */
    heat4 = 0.0;
    mass_now = 0.0; /* total cement mass corrected for hydration */
    suminit = c3sinit + c2sinit + c3ainit + c4afinit;
    /*        suminit=c3sinit+c2sinit+c3ainit+c4afinit+ncsbar+anhinit+heminit; */
    /* ctest is number of gypsum likely to form ettringite */
    /* 1 unit of C3A can react with 2.5 units of Gypsum */
    ctest = count [ DIFFGYP ];
#ifdef PRINTF
    printf("ctest is %ld\n", ctest);
#endif
    fflush(stdout);
    if ( ( float ) ctest > ( 2.5 * ( float ) ( count [ DIFFC3A ] + count [ DIFFC4A ] ) ) ) {
        ctest = ( long int ) ( 2.5 * ( float ) ( count [ DIFFC3A ] + count [ DIFFC4A ] ) );
    }

    for ( i = 0; i <= EMPTYP; i++ ) {
        if ( ( i != 0 ) && ( i <= ABSGYP ) && ( i != INERTAGG ) && ( i != CSH ) ) {
            heatsum += ( float ) count [ i ] * heatf [ i ] / molarv [ i ];
            /* Tabulate moles of H2O consumed by reactions so far */
            molesh2o += ( float ) count [ i ] * waterc [ i ] / molarv [ i ];
        }

        /* assume that all C3A which can, does form ettringite */
        if ( i == DIFFC3A ) {
            heatsum += ( ( float ) count [ DIFFC3A ] - ( float ) ctest / 2.5 ) * heatf [ C3A ] / molarv [ C3A ];
        }

        /* assume that all C4A which can, does form ettringite */
        if ( i == DIFFC4A ) {
            heatsum += ( ( float ) count [ DIFFC4A ] - ( float ) ctest / 2.5 ) * heatf [ C4AF ] / molarv [ C4AF ];
        }

        /* assume all gypsum which can, does form ettringite */
        /* rest will remain as gypsum */
        if ( i == DIFFGYP ) {
            heatsum += ( float ) ( count [ DIFFGYP ] - ctest ) * heatf [ GYPSUM ] / molarv [ GYPSUM ];
            /* 3.3 is the molar expansion from GYPSUM to ETTR */
            heatsum += ( float ) ctest * 3.30 * heatf [ ETTR ] / molarv [ ETTR ];
            molesdh2o += ( float ) ctest * 3.30 * waterc [ ETTR ] / molarv [ ETTR ];
        } else if ( i == DIFFCH ) {
            heatsum += ( float ) count [ DIFFCH ] * heatf [ CH ] / molarv [ CH ];
            molesdh2o += ( float ) count [ DIFFCH ] * waterc [ CH ] / molarv [ CH ];
        } else if ( i == DIFFFH3 ) {
            heatsum += ( float ) count [ DIFFFH3 ] * heatf [ FH3 ] / molarv [ FH3 ];
            molesdh2o += ( float ) count [ DIFFFH3 ] * waterc [ FH3 ] / molarv [ FH3 ];
        } else if ( i == DIFFCSH ) {
            /* use current CSH properties */
            heatsum += ( float ) count [ DIFFCSH ] * heatf [ CSH ] / molarvcsh [ cycle ];
            molesdh2o += ( float ) count [ DIFFCSH ] * watercsh [ cycle ] / molarvcsh [ cycle ];
        } else if ( i == DIFFETTR ) {
            heatsum += ( float ) count [ DIFFETTR ] * heatf [ ETTR ] / molarv [ ETTR ];
            molesdh2o += ( float ) count [ DIFFETTR ] * waterc [ ETTR ] / molarv [ ETTR ];
        } else if ( i == DIFFCACL2 ) {
            heatsum += ( float ) count [ DIFFCACL2 ] * heatf [ CACL2 ] / molarv [ CACL2 ];
            molesdh2o += ( float ) count [ DIFFCACL2 ] * waterc [ CACL2 ] / molarv [ CACL2 ];
        } else if ( i == DIFFAS ) {
            heatsum += ( float ) count [ DIFFAS ] * heatf [ ASG ] / molarv [ ASG ];
            molesdh2o += ( float ) count [ DIFFAS ] * waterc [ ASG ] / molarv [ ASG ];
        } else if ( i == DIFFCAS2 ) {
            heatsum += ( float ) count [ DIFFCAS2 ] * heatf [ CAS2 ] / molarv [ CAS2 ];
            molesdh2o += ( float ) count [ DIFFCAS2 ] * waterc [ CAS2 ] / molarv [ CAS2 ];
        }
        /* assume that all diffusing anhydrite leads to gypsum formation */
        else if ( i == DIFFANH ) {
            heatsum += ( float ) count [ DIFFANH ] * heatf [ GYPSUMS ] / molarv [ GYPSUMS ];
            /* 2 moles of water per mole of gypsum formed */
            molesdh2o += ( float ) count [ DIFFANH ] * 2.0 / molarv [ GYPSUMS ];
        }
        /* assume that all diffusing hemihydrate leads to gypsum formation */
        else if ( i == DIFFHEM ) {
            heatsum += ( float ) count [ DIFFHEM ] * heatf [ GYPSUMS ] / molarv [ GYPSUMS ];
            /* 1.5 moles of water per mole of gypsum formed */
            molesdh2o += ( float ) count [ DIFFHEM ] * 1.5 / molarv [ GYPSUMS ];
        } else if ( i == C3S ) {
            alpha += ( float ) ( c3sinit - count [ C3S ] );
            mass_now += specgrav [ C3S ] * ( float ) count [ C3S ];
            heat4 += .517 * ( float ) ( c3sinit - count [ C3S ] ) * specgrav [ C3S ];
        } else if ( i == C2S ) {
            alpha += ( float ) ( c2sinit - count [ C2S ] );
            mass_now += specgrav [ C2S ] * ( float ) count [ C2S ];
            heat4 += .262 * ( float ) ( c2sinit - count [ C2S ] ) * specgrav [ C2S ];
        } else if ( i == C3A ) {
            alpha += ( float ) ( c3ainit - count [ C3A ] );
            mass_now += specgrav [ C3A ] * ( float ) count [ C3A ];
            mc3ar = ( c3ainit - ( float ) count [ C3A ] ) / molarv [ C3A ];
            mc4ar = ( c4afinit - ( float ) count [ C4AF ] ) / molarv [ C4AF ];
            if ( ( mc3ar + mc4ar ) > 0.0 ) {
                frhyg = ( mc3ar / ( mc3ar + mc4ar ) ) * ( float ) count [ C3AH6 ] / molarv [ C3AH6 ];
            } else {
                frhyg = 0.0;
            }

            frettr = ( float ) count [ ETTR ] / molarv [ ETTR ];
            frafm = 3 * ( float ) count [ AFM ] / molarv [ AFM ];
            frtot = frafm + frettr + frhyg;
            if ( frtot > 0.0 ) {
                frettr /= frtot;
                frafm /= frtot;
                frhyg /= frtot;
                heat4 += frafm * 1.144 * ( float ) ( c3ainit - count [ C3A ] ) * specgrav [ C3A ];
                heat4 += frhyg * 0.908 * ( float ) ( c3ainit - count [ C3A ] ) * specgrav [ C3A ];
                heat4 += frettr * 1.672 * ( float ) ( c3ainit - count [ C3A ] ) * specgrav [ C3A ];
            }
        } else if ( i == C4AF ) {
            alpha += ( float ) ( c4afinit - count [ C4AF ] );
            mass_now += specgrav [ C4AF ] * ( float ) count [ C4AF ];
            mc3ar = ( c3ainit - ( float ) count [ C3A ] ) / molarv [ C3A ];
            mc4ar = ( c4afinit - ( float ) count [ C4AF ] ) / molarv [ C4AF ];
            if ( ( mc3ar + mc4ar ) > 0.0 ) {
                frhyg = ( mc4ar / ( mc3ar + mc4ar ) ) * ( float ) count [ C3AH6 ] / molarv [ C3AH6 ];
            } else {
                frhyg = 0.0;
            }

            frettr = ( float ) count [ ETTRC4AF ] / molarv [ ETTRC4AF ];
            frtot = frettr + frhyg;
            if ( frtot > 0.0 ) {
                frettr /= frtot;
                frhyg /= frtot;
                heat4 += frhyg * .418 * ( float ) ( c4afinit - count [ C4AF ] ) * specgrav [ C4AF ];
                heat4 += frettr * .725 * ( float ) ( c4afinit - count [ C4AF ] ) * specgrav [ C4AF ];
            }
        }
        /*                  else if(i==GYPSUM){
         *               alpha+=(float)(ncsbar-count[GYPSUM]);
         *               mass_now+=specgrav[GYPSUM]*(float)count[GYPSUM];
         *               }   */
        /* 0.187 kJ/g anhydrite for anhydrite --> gypsum conversion */
        else if ( i == ANHYDRITE ) {
            /*  alpha+=(float)(anhinit-count[ANHYDRITE]);
             *  mass_now+=specgrav[ANHYDRITE]*(float)count[ANHYDRITE]; */
            heat4 += .187 * ( float ) ( anhinit - count [ ANHYDRITE ] ) * specgrav [ ANHYDRITE ];
            /* 2 moles of water consumed per mole of anhydrite reacted */
            molesh2o += ( float ) ( anhinit - count [ ANHYDRITE ] ) * 2.0 / molarv [ ANHYDRITE ];
        }
        /* 0.132 kJ/g hemihydrate for hemihydrate-->gypsum conversion */
        else if ( i == HEMIHYD ) {
            /*                     alpha+=(float)(heminit-count[HEMIHYD]);
             *                     mass_now+=specgrav[HEMIHYD]*(float)count[HEMIHYD]; */
            heat4 += .132 * ( float ) ( heminit - count [ HEMIHYD ] ) * specgrav [ HEMIHYD ];
            /* 1.5 moles of water consumed per mole of anhydrite reacted */
            molesh2o += ( float ) ( heminit - count [ HEMIHYD ] ) * 1.5 / molarv [ HEMIHYD ];
        }
    }

    mass_fa_now = specgrav [ ASG ] * ( float ) count [ ASG ];
    mass_fa_now += specgrav [ CAS2 ] * ( float ) count [ CAS2 ];
    mass_fa_now += specgrav [ POZZ ] * ( float ) count [ POZZ ];
    if ( suminit != 0 ) {
        alpha = alpha / ( float ) suminit;
    } else {
        alpha = 0.0;
    }

    /* Current degree of hydration on a mass basis */
    if ( cemmass != 0.0 ) {
        alpha_cur = 1.0 - ( mass_now / cemmass );
    } else {
        alpha_cur = 0.0;
    }

    if ( flyashmass != 0.0 ) {
        alpha_fa_cur = 1.0 - ( mass_fa_now / flyashmass );
    } else {
        alpha_fa_cur = 0.0;
    }

    h2oinit = ( float ) porinit / molarv [ POROSITY ];

    /* Assume 780 J/g S for pozzolanic reaction */
    /* Each unit of silica fume consumes 1.35 units of CH, */
    /* so divide npr by 1.35 to get silca fume which has reacted */
    heat4 += 0.78 * ( ( float ) npr / 1.35 ) * specgrav [ POZZ ];

    /* Assume 800 J/g S for slag reaction */
    /* Seems reasonable with measurements of Biernacki and Richardson */
    heat4 += 0.8 * ( ( float ) nslagr ) * specgrav [ SLAG ];

    /* Assume 800 J/g AS for stratlingite formation (DeLarrard) */
    /* Each unit of AS consumes 1.3267 units of CH, */
    /* so divide nasr by 1.3267 to get ASG which has reacted */
    heat4 += 0.80 * ( ( float ) nasr / 1.3267 ) * specgrav [ ASG ];

    /* Should be additional code here for heat release due to CAS2 to */
    /* stratlingite conversion, but data unavailable at this time */

    /* Adjust heat sum for water left in system */
    water_left = ( long int ) ( ( h2oinit - molesh2o ) * molarv [ POROSITY ] + 0.5 );
    countkeep = count [ POROSITY ];
    heatsum += ( h2oinit - molesh2o - molesdh2o ) * heatf [ POROSITY ];
    if ( cyccnt == 0 ) {
#ifdef OUTFILES
        fprintf(heatfile, "Cycle time(h) alpha_vol alpha_mass heat4(kJ/kg_solid) Gsratio2 G-s_ratio\n");
#endif
    }

    heat_new = heat4; /* use heat4 for all adiabatic calculations */
    /* due to best agreement with calorimetry data */

    if ( cyccnt == 0 ) {
        //fprintf(chsfile,"Cycle  time(h) alpha_mass  Chemical shrinkage (ml/g cement)\n");
    }

    chs_new = ( ( float ) ( count [ EMPTYP ] + count [ POROSITY ] - water_left ) * heat_cf / 1000. );
    /*  if((molesh2o>h2oinit)&&(sealed==1)){  */
    if ( ( ( water_left + water_off ) < 0 ) && ( sealed == 1 ) ) {
#ifdef PRINTF
        printf("All water consumed at cycle %d \n", cyccnt);
        fflush(stdout);
#endif
        //exit(1);the simulation can be stopped
    }

    //remove defined amount of water (outer drying)
    //makeinert(10);

    /* Attempt to create empty porosity to account for self-desiccation */
    if ( ( sealed == 1 ) && ( ( count [ POROSITY ] - water_left ) > 0 ) ) {
        poretodo = ( count [ POROSITY ] - pore_off ) - ( water_left - water_off );
        poretodo -= slagemptyp;
        if ( poretodo > 0 ) {
            makeinert(poretodo);
            poregone += poretodo;
        }
    }

    /* Output phase counts */
    /* phasfile for reactant and product phases */
#ifdef OUTFILES
    if ( cyccnt == 0 ) {
        fprintf(phasfile, "#Cycle DoH Time Porosity C3S C2S C3A C4AF GYPSUM HEMIHYD ANHYDRITE POZZ INERT SLAG ASG CAS2 CH CSH C3AH6 ETTR ETTRC4AF AFM FH3 POZZCSH SLAGCSH CACL2 FREIDEL STRAT GYPSUMS CACO3 AFMC AGG ABSGYP EMPTYP HDCSH water_left \n");
        fprintf(perc_phases, "#Cyc Time[h] DoH SOL_per SOL_tot| Porosity C3S C2S C3A C4AF GYPSUM HEMIHYD ANHYDRITE POZZ INERT SLAG ASG CAS2 CH CSH C3AH6 ETTR ETTRC4AF AFM FH3 POZZCSH SLAGCSH CACL2 FREIDEL STRAT GYPSUMS CACO3 AFMC AGG ABSGYP EMPTYP HDCSH\n");
    }

    fprintf(phasfile, "%d %.4f %.4f ", cyccnt, alpha_cur, time_cur);
    fprintf(disprobfile, "%d %.4f %.4f ", cyccnt, alpha_cur, time_cur);


    CSHbox(CSH_vicinity);


    for ( i = 0; i <= HDCSH; i++ ) {
        if ( ( i < DIFFCSH ) || ( i >= EMPTYP ) ) {
            if ( i == CSH ) {
                fprintf(phasfile, "%ld ", count [ CSH ] - count [ HDCSH ]);
            } else {
                fprintf(phasfile, "%ld ", count [ i ]);
            }

            fprintf(disprobfile, "%f ", disprob [ i ]);
        }

        // printf("%ld ",count[i]);
    }

    fprintf(disprobfile, "\n");
    fprintf(phasfile, "%ld\n", water_left);
    fflush(phasfile);
    fflush(disprobfile);
#endif

#ifdef PRINTF
    printf("\n");
#endif

    if ( cycle == 0 ) {
        return;
    }

#ifdef PRINTF
 #ifdef __TM_MODULE
    printf("Cycle cyccnt %d, DoH %f, time %f, GP %p\n", cyccnt, alpha_cur, time_cur, this->gp);
 #else
    printf("Cycle cyccnt %d, DoH %f, time %f\n", cyccnt, alpha_cur, time_cur);
 #endif
#endif

    cyccnt += 1;
    fflush(stdout);
    /* Update current volume count for CH */
    chold = chnew;
    chnew = count [ CH ];

    /* See if ettringite is soluble yet */
    /* Gypsum 80% consumed, changed 06.09.00 from 90% to 80% */
    /* Gypsum 75% consumed, changed 09.09.01 from 80% to 75% */
    /* or system temperature exceeds 70 C */
    if ( ( ( ncsbar + anhinit + heminit ) != 0.0 ) || ( temp_cur >= 70.0 ) ) {
        /* Account for all sulfate sources and forms */
        if ( ( soluble [ ETTR ] == 0 ) && ( ( temp_cur >= 70.0 ) || ( count [ AFM ] != 0 ) ||
                                           ( ( ( ( float ) count [ GYPSUM ] + 1.42 * ( float ) count [ ANHYDRITE ] + 1.4 *
                                                ( float ) count [ HEMIHYD ] + ( float ) count [ GYPSUMS ] ) / ( ( float ) ncsbar +
                                                                                                               1.42 * ( float ) anhinit + 1.4 * ( float ) heminit + ( ( float ) netbar / 3.30 ) ) ) < 0.25 ) ) ) {
            soluble [ ETTR ] = 1;
#ifdef PRINTF
            printf("Ettringite is soluble beginning at cycle %d \n", cycle);
#endif
            /* identify all new soluble ettringite */
            passone(ETTR, ETTR, 2, 0);
        }
    }

    /* end of soluble ettringite test */

    /* Adjust ettringite solubility */
    /* if too many ettringites already in solution */
    if ( count [ DIFFETTR ] > DETTRMAX ) {
        disprob [ ETTR ] = 0.0;
    } else {
        disprob [ ETTR ] = disbase [ ETTR ];
    }

    /* Adjust CaCl2 solubility */
    /* if too many CaCl2 already in solution */
    if ( count [ DIFFCACL2 ] > DCACL2MAX ) {
        disprob [ CACL2 ] = 0.0;
    } else {
        disprob [ CACL2 ] = disbase [ CACL2 ];
    }

    /* Adjust CaCO3 solubility */
    /* if too many CaCO3 already in solution */
    if ( ( count [ DIFFCACO3 ] > DCACO3MAX ) && ( soluble [ ETTR ] == 0 ) ) {
        disprob [ CACO3 ] = 0.0;
    } else if ( count [ DIFFCACO3 ] > ( 4 * DCACO3MAX ) ) {
        disprob [ CACO3 ] = 0.0;
    } else {
        disprob [ CACO3 ] = disbase [ CACO3 ];
    }

    /* Adjust solubility of CH */
    /* based on amount of CH currently diffusing */
    /* Note that CH is always soluble to allow some */
    /* Ostwald ripening of the CH crystals */
    if ( ( float ) count [ DIFFCH ] >= CHCRIT ) {
        disprob [ CH ] = disbase [ CH ] * CHCRIT / ( float ) count [ DIFFCH ];
    } else {
        disprob [ CH ] = disbase [ CH ];
    }

    /* Adjust solubility of CH for temperature */
    /* Fit to data provided in Taylor, Cement Chemistry */
    /* Scale to a reference temperature of 25 C */
    /* and adjust based on availability of pozzolan */

#ifdef PRINTF
    printf("CH dissolution probability goes from %f ", disprob [ CH ]);
#endif
    disprob [ CH ] *= ( ( A0_CHSOL - A1_CHSOL * temp_cur ) / ( A0_CHSOL - A1_CHSOL * 25.0 ) );
    if ( ( ppozz > 0.0 ) && ( nfill > 0 ) ) {
        disprob [ CH ] *= ppozz / PPOZZ;
    }

#ifdef PRINTF
    printf("to %f \n", disprob [ CH ]);
#endif

    /* Adjust solubility of ASG and CAS2 phases */
    /* based on pH rise during hydration */
    /* To be added at a later date */
    disprob [ ASG ] = disbase [ ASG ];
    disprob [ CAS2 ] = disbase [ CAS2 ];
    /* Address solubility of C3AH6 */
    /* If lots of gypsum or reactive ettringite, allow C3AH6 to dissolve */
    /* to generate diffusing C3A species */
    if ( ( ( count [ GYPSUM ] + count [ GYPSUMS ] ) > ( int ) ( ( ( float ) ncsbar +
                                                                 1.42 * ( float ) anhinit + 1.4 * ( float ) heminit ) * 0.05 ) ) ||
        ( count [ ETTR ] > 500 * SYSIZE_POW3 / 1000000. ) ) {
        soluble [ C3AH6 ] = 1;
        passone(C3AH6, C3AH6, 2, 0);
        /* Base C3AH6 solubility on maximum sulfate in solution */
        /* from gypsum or ettringite available for dissolution */
        /* The more the sulfate, the higher this solubility should be */
        maxsulfate = count [ DIFFGYP ];
        if ( ( maxsulfate < count [ DIFFETTR ] ) && ( soluble [ ETTR ] == 1 ) ) {
            maxsulfate = count [ DIFFETTR ];
        }

        /* Adjust C3AH6 solubility based on potential gypsum which will dissolve */
        /*gypready and maxsulfate is size dependent */
        if ( maxsulfate < ( int ) ( ( float ) gypready * disprob [ GYPSUM ] * ( float ) count [ POROSITY ] / SYSIZE_POW3 ) ) {   //fixed
            maxsulfate = ( int ) ( ( float ) gypready * disprob [ GYPSUM ] * ( float ) count [ POROSITY ] / SYSIZE_POW3 );   //fixed
        }

        if ( maxsulfate > 0 ) {
            disprob [ C3AH6 ] = disbase [ C3AH6 ] * ( float ) ( maxsulfate / C3AH6CRIT );
            if ( disprob [ C3AH6 ] > 0.5 ) {
                disprob [ C3AH6 ] = 0.5;
            }
        } else {
            disprob [ C3AH6 ] = disbase [ C3AH6 ];
        }
    } else {
        soluble [ C3AH6 ] = 0;
    }

    /* See if silicates are soluble yet */
    if ( ( soluble [ C3S ] == 0 ) && ( ( cycle > 1 ) || ( count [ ETTR ] > 0 ) || ( count [ AFM ] > 0 ) || ( count [ ETTRC4AF ] > 0 ) ) ) {
        soluble [ C2S ] = 1;
        soluble [ C3S ] = 1;
        /* identify all new soluble silicates */
        passone(C3S, C2S, 2, 0);
    }

    /* end of soluble silicate test */
    /* Adjust solubility of C3S and C2S with CSH concentration */
    /* for simulation of induction period */

    tdisfact = A0_CHSOL - temp_cur * A1_CHSOL;
    /*   printf("tdisfact is %f\n",tdisfact);
     *   fflush(stdout); */

    /* Calculation of cs_acc; acceleration of C3S and C2S reaction by CaSO4 */
    /* Calculation of ca_acc; acceleration of C3A and C4AF reaction by CaSO4 */
    /* November 2004 --- modified to be on a sulfate per unit initial cement */
    /*    per unit porosity basis */
    /*  Try using current porosity count to see if this helps in initial */
    /*   hydration rates of sealed vs. saturated  */
    pfract = ( float ) count [ POROSITY ] / ( float ) ( SYSIZE * SYSIZE * SYSIZE );
    if ( ( ncsbar + anhinit + heminit ) == 0.0 ) {
        cs_acc = 1.0;
        ca_acc = 1.0;
        dismin_c3a = 5.0 * DISMIN_C3A_0;
        dismin_c4af = 5.0 * DISMIN_C4AF_0;
    } else {
        sulf_conc = sulf_cur * tfractw05 * pfractw05 / totfract / pfract;
        sulf_conc *= ( 1000000. / SYSIZE_POW3 ); //fixed
        if ( sulf_conc < 10.0 ) {
            cs_acc = 1.0;
            ca_acc = 1.0;
            dismin_c3a = DISMIN_C3A_0;
            dismin_c4af = DISMIN_C4AF_0;
        } else if ( sulf_conc < 20.0 ) {
            cs_acc = 1.0 + ( ( sulf_conc ) - 10.0 ) / 10.0;
            ca_acc = 1.0;
            dismin_c3a = DISMIN_C3A_0;
            dismin_c4af = DISMIN_C4AF_0;
        } else {
            cs_acc = 1.0 + ( double ) log10(sulf_conc - 10.0);
            ca_acc = 1.0;
            dismin_c3a = ( 6.0 - ( double ) log10(sulf_conc) ) * DISMIN_C3A_0;
            dismin_c4af = ( 6.0 - ( double ) log10(sulf_conc) ) * DISMIN_C4AF_0;
            if ( dismin_c3a < DISMIN_C3A_0 ) {
                dismin_c3a = DISMIN_C3A_0;
            }

            if ( dismin_c4af < DISMIN_C4AF_0 ) {
                dismin_c4af = DISMIN_C4AF_0;
            }
        }
    }

    /* Suggest change WCSCALE/w_to_c to (0.3125+WCSCALE)/(0.3125+w_to_c) */
    /* to have the induction period scaling depend on volume of CSH */
    /* produced per volume (not mass) of cement */
    /*    dfact=tdisfact*((float)count[CSH]/((float)CSHSCALE*(0.3125+WCSCALE)/(w_to_c+0.3125)))*((float)count[CSH]/((float)CSHSCALE*(0.3125+WCSCALE)/(w_to_c+0.3125)))*cs_acc;*/
    /* October 2004 --- changed to truly scale with volume of cement in */
    /* system for both plain portland cements and filled systems */

    dfact = tdisfact * ( ( float ) count [ CSH ] / ( ( float ) CSHSCALE * surffract * totfract / tfractw04 ) ) * ( ( float ) count [ CSH ] / ( ( float ) CSHSCALE * surffract * totfract / tfractw04 ) ) * cs_acc;

    disprob [ C3S ] = DISMIN + dfact * disbase [ C3S ];
    disprob [ C2S ] = DISMIN2 + dfact * disbase [ C2S ];
    if ( disprob [ C3S ] > ( 1. * disbase [ C3S ] ) ) {
        disprob [ C3S ] = ( 1. * disbase [ C3S ] );
    }

    if ( disprob [ C2S ] > ( 1. * disbase [ C2S ] ) ) {
        disprob [ C2S ] = ( 1. * disbase [ C2S ] );
    }

    /* Also adjust slag and fly ash dissolution rates here */
    /* Really slow down initial slag and fly ash dissolutions */
    /* Ultimately should be linked to pH of pore solution, most likely */
    disprob [ SLAG ] = slagreact * ( DISMINSLAG + dfact * disbase [ SLAG ] ) / 10.0;
    if ( disprob [ SLAG ] > ( slagreact * disbase [ SLAG ] ) ) {
        disprob [ SLAG ] = ( slagreact * disbase [ SLAG ] );
    }

    if ( disprob [ C3S ] == disbase [ C3S ] ) {
        disprob [ SLAG ] = slagreact * disbase [ SLAG ];
    }

    disprob [ ASG ] = DISMINASG + dfact * disbase [ ASG ] / 5.0;
    if ( disprob [ ASG ] > ( 1. * disbase [ ASG ] ) ) {
        disprob [ ASG ] = ( 1. * disbase [ ASG ] );
    }

    if ( disprob [ C3S ] == disbase [ C3S ] ) {
        disprob [ ASG ] = disbase [ ASG ];
    }

    disprob [ CAS2 ] = DISMINCAS2 + dfact * disbase [ CAS2 ] / 5.0;
    if ( disprob [ CAS2 ] > ( 1. * disbase [ CAS2 ] ) ) {
        disprob [ CAS2 ] = ( 1. * disbase [ CAS2 ] );
    }

    if ( disprob [ C3S ] == disbase [ C3S ] ) {
        disprob [ CAS2 ] = disbase [ CAS2 ];
    }

    /* Adjust CAS2 solubility */
    /* if too many CAS2 already in solution */
    if ( count [ DIFFCAS2 ] > DCAS2MAX ) {
        disprob [ CAS2 ] = 0.0;
    }

#ifdef PRINTF
    printf("Silicate probabilities: %f %f\n", disprob [ C3S ], disprob [ C2S ]);
#endif
    fflush(stdout);
    /* Assume that aluminate dissolution controlled by formation */
    /* of impermeable layer proportional to CSH concentration */
    /* if sulfates are present in the system */

    /*      dfact1=tdisfact*((float)count[CSH]/CSHSCALE)*((float)count[CSH]/CSHSCALE)*ca_acc;  */
    if ( ( ncsbar + heminit + anhinit ) > 1000 * SYSIZE_POW3 / 1000000. ) { //fixed
        /*    dfact1=tdisfact*((float)count[CSH]/((float)CSHSCALE*(0.3125+WCSCALE)/(0.3125+w_to_c)))*((float)count[CSH]/((float)CSHSCALE*(0.3125+WCSCALE)/(0.3125+w_to_c)))*ca_acc; */
        /* October 2004 --- changed to truly scale with volume of cement in */
        /* system for both plain portland cements and filled systems */
        dfact1 = tdisfact * ( ( float ) count [ CSH ] / ( ( float ) CSHSCALE * surffract * totfract / tfractw04 ) ) * ( ( float ) count [ CSH ] / ( ( float ) CSHSCALE * surffract * totfract / tfractw04 ) ) * ca_acc;
        disprob [ C3A ] = dismin_c3a + dfact1 * disbase [ C3A ];
        disprob [ C4AF ] = dismin_c4af + dfact1 * disbase [ C4AF ];
        if ( disprob [ C3A ] > ( 1. * disbase [ C3A ] ) ) {
            disprob [ C3A ] = ( 1. * disbase [ C3A ] );
        }

        if ( disprob [ C4AF ] > ( 1. * disbase [ C4AF ] ) ) {
            disprob [ C4AF ] = ( 1. * disbase [ C4AF ] );
        }

        /* Location to add in dissolution reduction in calcium sulfate phases */
        /* if needed */
        disprob [ GYPSUM ] = ( disbase [ GYPSUM ] / 15. ) + dfact1 * disbase [ GYPSUM ];
        if ( disprob [ GYPSUM ] > ( disbase [ GYPSUM ] ) ) {
            disprob [ GYPSUM ] = ( disbase [ GYPSUM ] );
        }

        disprob [ GYPSUMS ] = ( disbase [ GYPSUMS ] / 15. ) + dfact1 * disbase [ GYPSUMS ];
        if ( disprob [ GYPSUMS ] > ( disbase [ GYPSUMS ] ) ) {
            disprob [ GYPSUMS ] = ( disbase [ GYPSUMS ] );
        }

        /* Adjust gypsum solubility */
        /* if too many diffusing gypsums already in solution */
        if ( count [ DIFFGYP ] > DGYPMAX ) {
            disprob [ GYPSUM ] = disprob [ GYPSUMS ] = 0.0;
        }

        disprob [ HEMIHYD ] = ( disbase [ HEMIHYD ] / 15. ) + dfact1 * disbase [ HEMIHYD ];
        if ( disprob [ HEMIHYD ] > ( disbase [ HEMIHYD ] ) ) {
            disprob [ HEMIHYD ] = ( disbase [ HEMIHYD ] );
        }

        disprob [ ANHYDRITE ] = ( disbase [ ANHYDRITE ] / 15. ) + dfact1 * disbase [ ANHYDRITE ];
        if ( disprob [ ANHYDRITE ] > ( disbase [ ANHYDRITE ] ) ) {
            disprob [ ANHYDRITE ] = ( disbase [ ANHYDRITE ] );
        }
    } else {
        /* Cause flash set by increasing dissolution rates of C3A and C4AF */
        /* each by a factor of four */
        disprob [ C3A ] = 4. * disbase [ C3A ];
        disprob [ C4AF ] = 4. * disbase [ C4AF ];
        disprob [ GYPSUM ] = disbase [ GYPSUM ];
        disprob [ HEMIHYD ] = disbase [ HEMIHYD ];
        disprob [ ANHYDRITE ] = disbase [ ANHYDRITE ];
    }

    /* Reduce dissolution probabilities based on saturation of system */
    if ( ( count [ EMPTYP ] > 0 ) && ( ( count [ POROSITY ] + count [ EMPTYP ] ) < 220000. * ( double ) ( SYSIZE_POW3 / 1000000. ) ) ) {       //fixed
        if ( countpore == 0 ) {
            countpore = count [ EMPTYP ];
        }

        saturation = ( float ) ( count [ POROSITY ] ) / ( float ) ( count [ POROSITY ] + ( count [ EMPTYP ] - countpore ) );
        /* Roughly according to results of Jensen, powers for RH
         * sensitivity are:
         * C3S-19
         * C2S-29
         * C3A, C4AF-6    */
        /* Adjust fly ash silicates (ASG and CAS2) and pozzolanic reactivity */
        /* by same factor as C3S (also CH) */
        disprob [ C3S ] *= ( saturation * saturation );
        disprob [ C3S ] *= ( saturation * saturation );
        disprob [ C3S ] *= ( saturation * saturation );
        disprob [ C3S ] *= ( saturation * saturation );
        disprob [ C3S ] *= ( saturation * saturation );
        disprob [ C3S ] *= ( saturation * saturation );
        disprob [ C3S ] *= ( saturation * saturation );
        disprob [ C3S ] *= ( saturation * saturation );
        disprob [ C3S ] *= ( saturation * saturation );
        disprob [ C3S ] *= ( saturation );
        disprob [ SLAG ] *= ( saturation * saturation );
        disprob [ SLAG ] *= ( saturation * saturation );
        disprob [ SLAG ] *= ( saturation * saturation );
        disprob [ SLAG ] *= ( saturation * saturation );
        disprob [ SLAG ] *= ( saturation * saturation );
        disprob [ SLAG ] *= ( saturation * saturation );
        disprob [ SLAG ] *= ( saturation * saturation );
        disprob [ SLAG ] *= ( saturation * saturation );
        disprob [ SLAG ] *= ( saturation * saturation );
        disprob [ SLAG ] *= ( saturation );
        disprob [ CH ] *= ( saturation * saturation );
        disprob [ CH ] *= ( saturation * saturation );
        disprob [ CH ] *= ( saturation * saturation );
        disprob [ CH ] *= ( saturation * saturation );
        disprob [ CH ] *= ( saturation * saturation );
        disprob [ CH ] *= ( saturation * saturation );
        disprob [ CH ] *= ( saturation * saturation );
        disprob [ CH ] *= ( saturation * saturation );
        disprob [ CH ] *= ( saturation * saturation );
        disprob [ CH ] *= ( saturation );
        disprob [ ASG ] *= ( saturation * saturation );
        disprob [ ASG ] *= ( saturation * saturation );
        disprob [ ASG ] *= ( saturation * saturation );
        disprob [ ASG ] *= ( saturation * saturation );
        disprob [ ASG ] *= ( saturation * saturation );
        disprob [ ASG ] *= ( saturation * saturation );
        disprob [ ASG ] *= ( saturation * saturation );
        disprob [ ASG ] *= ( saturation * saturation );
        disprob [ ASG ] *= ( saturation * saturation );
        disprob [ ASG ] *= ( saturation );
        disprob [ CAS2 ] *= ( saturation * saturation );
        disprob [ CAS2 ] *= ( saturation * saturation );
        disprob [ CAS2 ] *= ( saturation * saturation );
        disprob [ CAS2 ] *= ( saturation * saturation );
        disprob [ CAS2 ] *= ( saturation * saturation );
        disprob [ CAS2 ] *= ( saturation * saturation );
        disprob [ CAS2 ] *= ( saturation * saturation );
        disprob [ CAS2 ] *= ( saturation * saturation );
        disprob [ CAS2 ] *= ( saturation * saturation );
        disprob [ CAS2 ] *= ( saturation );
        ppozz *= ( saturation * saturation );
        ppozz *= ( saturation * saturation );
        ppozz *= ( saturation * saturation );
        ppozz *= ( saturation * saturation );
        ppozz *= ( saturation * saturation );
        ppozz *= ( saturation * saturation );
        ppozz *= ( saturation * saturation );
        ppozz *= ( saturation * saturation );
        ppozz *= ( saturation * saturation );
        ppozz *= ( saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );

        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation * saturation );
        disprob [ C2S ] *= ( saturation );
        disprob [ C3A ] *= ( saturation * saturation );
        disprob [ C3A ] *= ( saturation * saturation );
        disprob [ C3A ] *= ( saturation * saturation );
        disprob [ C4AF ] *= ( saturation * saturation );
        disprob [ C4AF ] *= ( saturation * saturation );
        disprob [ C4AF ] *= ( saturation * saturation );
    }

#ifdef PRINTF
    printf("Silicate and aluminate probabilities: %f %f %f %f\n", disprob [ C3S ], disprob [ C2S ], disprob [ C3A ], disprob [ C4AF ]);
    printf("cs_acc is %f  and ca_acc is %f sulf_cur is %ld\n", cs_acc, ca_acc, sulf_cur);
#endif
    fflush(stdout);
    /* Pass two- perform the dissolution of species */
    /* Determine the pH factor to use */
    pHfactor = 0.0;
    if ( ( pHactive == 1 ) && ( count [ CSH ] > ( ( CSHSCALE * surffract * surffract * totfract * totfract / tfractw04 / tfractw04 ) / 8.0 ) ) ) {
        pHfactor = 1.5;
        if ( pH_cur > 12.5 ) {
            pHfactor = 1.0;
        }

        /* 2/02*/
        /*          pHfactor=0.30*(((14.1-pH_cur)/0.3)-1.0); */
        /* 3/02*/
        /*  pHfactor=0.60*(((14.1-pH_cur)/0.5)-1.0); */
        if ( pH_cur > 12.75 ) {
            pHfactor = 0.667;
        }

        if ( pH_cur > 13.00 ) {
            pHfactor = 0.333;
        }

        if ( pH_cur > 13.25 ) {
            pHfactor = ( 0.0 );
        }

        if ( pH_cur > 13.75 ) {
            pHfactor = ( -0.25 );
        }

        pHfactor += concsulfate; /* influence of sulfate on reactivity */
    }

    nhgd = 0;
    /* Update molar volume ratios for CSH formation */
    pc3scsh = molarvcsh [ cyccnt ] / molarv [ C3S ] - 1.0;
    pc2scsh = molarvcsh [ cyccnt ] / molarv [ C2S ] - 1.0;
    /* Once again, scan all pixels in microstructure */
    slagemptyp = 0;
    for ( xloop = 0; xloop < SYSIZE; xloop++ ) {
        for ( yloop = 0; yloop < SYSIZE; yloop++ ) {
            for ( zloop = 0; zloop < SYSIZE; zloop++ ) {
                if ( mic [ xloop ] [ yloop ] [ zloop ] > OFFSET ) {
                    phid = mic [ xloop ] [ yloop ] [ zloop ] - OFFSET;
                    /* attempt a one-step random walk to dissolve */
                    plnew = ( int ) ( ( float ) NEIGHBORS * ran1(seed) );
                    if ( ( plnew < 0 ) || ( plnew >= NEIGHBORS ) ) {
                        plnew = NEIGHBORS - 1;
                    }

                    xc = xloop + xoff [ plnew ];
                    yc = yloop + yoff [ plnew ];
                    zc = zloop + zoff [ plnew ];
                    if ( xc < 0 ) {
                        xc = ( SYSIZEM1 );
                    }

                    if ( yc < 0 ) {
                        yc = ( SYSIZEM1 );
                    }

                    if ( xc >= SYSIZE ) {
                        xc = 0;
                    }

                    if ( yc >= SYSIZE ) {
                        yc = 0;
                    }

                    if ( zc < 0 ) {
                        zc = ( SYSIZEM1 );
                    }

                    if ( zc >= SYSIZE ) {
                        zc = 0;
                    }

                    /* Generate probability for dissolution */
                    pdis = ran1(seed);
                    /* Bias dissolution for one pixel particles as */
                    /* indicated by a pixel value of zero in the */
                    /* particle microstructure image */
                    if ( ( ( pdis <= ( disprob [ phid ] / ( 1. + pHfactor * pHeffect [ phid ] ) ) ) || ( ( pdis <= ( onepixelbias * disprob [ phid ] / ( 1. + pHfactor * pHeffect [ phid ] ) ) ) && ( micpart [ xloop ] [ yloop ] [ zloop ] == 0 ) ) ) && ( mic [ xc ] [ yc ] [ zc ] == POROSITY ) ) { //not fixed yet for one-voxel particles
                        discount [ phid ] += 1;
                        cread = creates [ phid ];
                        count [ phid ] -= 1;
                        mic [ xloop ] [ yloop ] [ zloop ] = POROSITY;
                        if ( phid == C3AH6 ) {
                            nhgd += 1;
                        }

                        /* Special dissolution for C4AF */
                        if ( phid == C4AF ) {
                            plfh3 = ran1(seed);
                            if ( ( plfh3 < 0.0 ) || ( plfh3 > 1.0 ) ) {
                                plfh3 = 1.0;
                            }

                            /* For every C4AF that dissolves, 0.5453 */
                            /* diffusing FH3 species should be created */
                            if ( plfh3 <= 0.5453 ) {
                                cread = DIFFFH3;
                            }
                        }

                        if ( cread == POROSITY ) {
                            count [ POROSITY ] += 1;
                        }

                        if ( cread != POROSITY ) {
                            nmade += 1;
                            ngoing += 1;
                            phnew = cread;
                            count [ phnew ] += 1;
                            mic [ xc ] [ yc ] [ zc ] = phnew;
                            antadd = ( struct ants * ) malloc( sizeof( struct ants ) );
                            antadd->x = xc;
                            antadd->y = yc;
                            antadd->z = zc;
                            antadd->id = phnew;
                            antadd->cycbirth = cyccnt;
                            /* Now connect this ant structure to end of linked list */
                            antadd->prevant = tailant;
                            tailant->nextant = antadd;
                            antadd->nextant = NULL;
                            tailant = antadd;
                        }

                        /* Extra CSH diffusing species based on current temperature */
                        if ( ( phid == C3S ) || ( phid == C2S ) ) {
                            plfh3 = ran1(seed);
                            if ( ( ( phid == C2S ) && ( plfh3 <= pc2scsh ) ) || ( plfh3 <= pc3scsh ) ) {
                                cshboxsize = ( int ) ( 3. + 5. * ( 40. - temp_cur ) / 20. );
                                if ( cshboxsize < 1 ) {
                                    cshboxsize = 1;
                                }

                                if ( cshboxsize >= SYSIZE ) {
                                    cshboxsize = SYSIZE;
                                }

                                //11.12.2006 smilauer
                                placed = loccsh(xc, yc, zc, cshboxsize);
                                if ( placed != 0 ) {
                                    count [ DIFFCSH ] += 1;
                                    count [ POROSITY ] -= 1;
                                } else {
                                    cshrand += 1;
                                }
                            }
                        }

                        if ( ( phid == C2S ) && ( pc2scsh > 1.0 ) ) {
                            plfh3 = ran1(seed);
                            if ( plfh3 <= ( pc2scsh - 1.0 ) ) {
                                cshboxsize = ( int ) ( 3. + 5. * ( 40. - temp_cur ) / 20. );
                                if ( cshboxsize < 1 ) {
                                    cshboxsize = 1;
                                }

                                if ( cshboxsize >= SYSIZE ) {
                                    cshboxsize = SYSIZE;
                                }

                                //11.12.2006 smilauer
                                placed = loccsh(xc, yc, zc, cshboxsize);
                                if ( placed != 0 ) {
                                    count [ DIFFCSH ] += 1;
                                    count [ POROSITY ] -= 1;
                                } else {
                                    cshrand += 1;
                                }
                            }
                        }
                    } else {
                        mic [ xloop ] [ yloop ] [ zloop ] -= OFFSET;
                    }
                }

                /* end of if edge loop */
                /* Now check if CSH to pozzolanic CSH conversion is possible */
                /* Only if CH is less than 15% in volume */
                /* Only if CSH is in contact with at least one porosity */
                /* and user wishes to use this option */
                if ( ( count [ POZZ ] >= 13000 ) && ( chnew < ( 0.15 * SYSIZE * SYSIZE * SYSIZE ) ) && ( csh2flag == 1 ) ) {
                    if ( mic [ xloop ] [ yloop ] [ zloop ] == CSH ) {
                        if ( ( countbox(3, xloop, yloop, zloop) ) >= 1 ) {
                            pconvert = ran1(seed);
                            if ( pconvert < PCSH2CSH ) {
                                count [ CSH ] -= 1;
                                plfh3 = ran1(seed);
                                /* molarvcsh units of C1.7SHx goes to */
                                /* 101.81 units of C1.1SH3.9 */
                                /* with 19.86 units of CH */
                                /* so p=calcy */
                                calcz = 0.0;
                                cycnew = cshage [ xloop ] [ yloop ] [ zloop ];
                                calcy = molarv [ POZZCSH ] / molarvcsh [ cycnew ];
                                if ( calcy > 1.0 ) {
                                    calcz = calcy - 1.0;
                                    calcy = 1.0;
                                    printf("Problem of not creating enough pozzolanic CSH during CSH conversion \n");
                                    printf("Current temperature is %f C\n", temp_cur);
                                }

                                if ( plfh3 <= calcy ) {
                                    mic [ xloop ] [ yloop ] [ zloop ] = POZZCSH;
                                    count [ POZZCSH ] += 1;
                                } else {
                                    mic [ xloop ] [ yloop ] [ zloop ] = DIFFCH;
                                    nmade += 1;
                                    ncshgo += 1;
                                    ngoing += 1;
                                    count [ DIFFCH ] += 1;
                                    antadd = ( struct ants * ) malloc( sizeof( struct ants ) );
                                    antadd->x = xloop;
                                    antadd->y = yloop;
                                    antadd->z = zloop;
                                    antadd->id = DIFFCH;
                                    antadd->cycbirth = cyccnt;
                                    /* Now connect this ant structure to end of linked list */
                                    antadd->prevant = tailant;
                                    tailant->nextant = antadd;
                                    antadd->nextant = NULL;
                                    tailant = antadd;
                                }

                                /* Possibly need even more pozzolanic CSH */
                                /* Would need a diffusing pozzolanic
                                 * CSH species??? */
                                /*                    if(calcz>0.0){
                                 * plfh3=ran1(seed);
                                 * if(plfh3<=calcz){
                                 * cshrand+=1;
                                 * }
                                 * }    */


                                plfh3 = ran1(seed);
                                calcx = ( 19.86 / molarvcsh [ cycnew ] ) - ( 1. - calcy );
                                /* Ex. 0.12658=(19.86/108.)-(1.-0.94269) */
                                if ( plfh3 < calcx ) {
                                    npchext += 1;
                                }
                            }
                        }
                    }
                }

                /* See if slag can react --- in contact with at least one porosity */
                if ( mic [ xloop ] [ yloop ] [ zloop ] == SLAG ) {
                    if ( ( countbox(3, xloop, yloop, zloop) ) >= 1 ) {
                        pconvert = ran1(seed);
                        if ( pconvert < ( disprob [ SLAG ] / ( 1. + pHfactor * pHeffect [ SLAG ] ) ) ) {
                            nslagr += 1;
                            count [ SLAG ] -= 1;
                            discount [ SLAG ] += 1;
                            /* Check on extra C3A generation */
                            plfh3 = ran1(seed);
                            if ( plfh3 < p5slag ) {
                                nslagc3a += 1;
                            }

                            /* Convert slag to reaction products */
                            plfh3 = ran1(seed);
                            if ( plfh3 < p1slag ) {
                                mic [ xloop ] [ yloop ] [ zloop ] = SLAGCSH;
                                /* Assign a plate axes identifier to this slag C-S-H voxel */
                                msface = ( int ) ( 3. * ran1(seed) + 1. );
                                if ( msface > 3 ) {
                                    msface = 1;
                                }

                                faces [ xloop ] [ yloop ] [ zloop ] = msface;
                                count [ SLAGCSH ] += 1;
                            } else {
                                if ( sealed == 1 ) {
                                    /* Create empty porosity at slag site */
                                    slagemptyp += 1;
                                    mic [ xloop ] [ yloop ] [ zloop ] = EMPTYP;
                                    count [ EMPTYP ] += 1;
                                } else {
                                    mic [ xloop ] [ yloop ] [ zloop ] = POROSITY;
                                    count [ POROSITY ] += 1;
                                }
                            }

                            /* Add in extra SLAGCSH as needed */
                            p3init = p3slag;
                            while ( p3init > 1.0 ) {
                                extslagcsh(xloop, yloop, zloop);
                                p3init -= 1.0;
                            }

                            plfh3 = ran1(seed);
                            if ( plfh3 < p3init ) {
                                extslagcsh(xloop, yloop, zloop);
                            }
                        }
                    }
                }
            }

            /* end of zloop */
        }

        /* end of yloop */
    }

    /* end of xloop */

    if ( ncshgo != 0 ) {
        printf("CSH dissolved is %ld \n", ncshgo);
    }

    if ( npchext > 0 ) {
        printf("npchext is %ld at cycle %d \n", npchext, cycle);
    }

    /* Now add in the extra diffusing species for dissolution */
    /* Expansion factors from Young and Hansen and */
    /* Mindess and Young (Concrete) */
    ncshext = cshrand;
    if ( cshrand != 0 ) {
        printf("cshrand is %d \n", cshrand);
    }

    /* CH, Gypsum, and diffusing C3A are added at totally random
     * locations as opposed to at the dissolution site
     * molar_weight[CH]=74 g/mol
     * molar_weight[C3S]=228 g/mol
     * only for C3S hydration - reason is the Pignat's thesis
     * C3S+5.3 H2O -> C1.7SH4 + 1.3 CH
     */

    fchext = 1.3 * 74 * specgrav [ C3S ] / ( specgrav [ CH ] * 228 ) * ( float ) discount [ C3S ] + 0.191 * ( float ) discount [ C2S ] +
             0.2584 * ( float ) discount [ C4AF ]; //originally 0.61

    nchext = ( int ) fchext;
    if ( fchext > ( float ) nchext ) {
        pdis = ran1(seed);
        if ( ( fchext - ( float ) nchext ) > pdis ) {
            nchext += 1;
        }
    }

    nchext += npchext;
    /* Adjust CH addition for slag consumption and maintain deficit as needed */
    slagcum += discount [ SLAG ];
    chgone = ( int ) ( p4slag * ( float ) slagcum );
    nchext -= chgone;
    slagcum -= ( int ) ( ( float ) chgone / p4slag );
    nchext -= DIFFCHdeficit;
    DIFFCHdeficit = 0;
    if ( nchext < 0 ) {
        DIFFCHdeficit -= nchext;
        nchext = 0;
    }

    fc3aext = discount [ C3A ] + 0.5917 * ( float ) discount [ C3AH6 ];
    nc3aext = ( int ) ( fc3aext + nslagc3a );
    if ( fc3aext > ( float ) nc3aext ) {
        pdis = ran1(seed);
        if ( ( fc3aext - ( float ) nc3aext ) > pdis ) {
            nc3aext += 1;
        }
    }

    fc4aext = 0.696 * ( float ) discount [ C4AF ];
    nc4aext = ( int ) fc4aext;
    if ( fc4aext > ( float ) nc4aext ) {
        pdis = ran1(seed);
        if ( ( fc4aext - ( float ) nc4aext ) > pdis ) {
            nc4aext += 1;
        }
    }

    /* both forms of GYPSUM form same DIFFGYP species */
    ngypext = discount [ GYPSUM ] + discount [ GYPSUMS ];
    /* Convert to diffusing anhydrite at volume necessary for final */
    /* gypsum formation  (1 anhydrite --> 1.423 gypsum) */
    /* Since hemihydrate can now react with C3A, etc., can't */
    /* do expansion here any longer  7/99 */
    /*    fanhext=1.423*(float)discount[ANHYDRITE]; */
    fanhext = ( float ) discount [ ANHYDRITE ];
    nanhext = ( int ) fanhext;
    if ( fanhext > ( float ) nanhext ) {
        pdis = ran1(seed);
        if ( ( fanhext - ( float ) nanhext ) > pdis ) {
            nanhext += 1;
        }
    }

    /* Convert to diffusing hemiydrate at volume necessary for final */
    /* gypsum formation  (1 hemihydrate --> 1.4 gypsum) */
    /* Since hemihydrate can now react with C3A, etc., can't */
    /* do expansion here any longer  7/99 */
    fhemext = ( float ) discount [ HEMIHYD ];
    /*    fhemext=1.3955*(float)discount[HEMIHYD];  */

    nhemext = ( int ) fhemext;
    if ( fhemext > ( float ) nhemext ) {
        pdis = ran1(seed);
        if ( ( fhemext - ( float ) nhemext ) > pdis ) {
            nhemext += 1;
        }
    }

    count [ DIFFGYP ] += ngypext;
    count [ DIFFANH ] += nanhext;
    count [ DIFFHEM ] += nhemext;
    count [ DIFFCH ] += nchext;
    count [ DIFFCSH ] += ncshext;
    count [ DIFFC3A ] += nc3aext;
    count [ DIFFC4A ] += nc4aext;

    nsum2 = nchext + ncshext;
    nsum3 = nsum2 + nc3aext;
    nsum4 = nsum3 + nc4aext;
    nsum5 = nsum4 + ngypext;
    nsum6 = nsum5 + nhemext;
    fflush(stdout);

    //add extra diffusing species - check for infinite looping if porosity is not available - smilauer 01/16/2008
    for ( xext = 1; xext <= ( nsum6 + nanhext ); xext++ ) {
        plok = 0;
        tries = 0;
        do {
            xc = ( int ) ( ( float ) SYSIZE * ran1(seed) );   //locate random place in the microstructure
            yc = ( int ) ( ( float ) SYSIZE * ran1(seed) );
            zc = ( int ) ( ( float ) SYSIZE * ran1(seed) );
            if ( xc >= SYSIZE ) {
                xc = 0;
            }

            if ( yc >= SYSIZE ) {
                yc = 0;
            }

            if ( zc >= SYSIZE ) {
                zc = 0;
            }

            tries++;

            if ( mic [ xc ] [ yc ] [ zc ] == POROSITY ) {
                plok = 1;
                phid = DIFFCH;
                count [ POROSITY ] -= 1;
                if ( xext > nsum6 ) {
                    phid = DIFFANH;
                } else if ( xext > nsum5 ) {
                    phid = DIFFHEM;
                } else if ( xext > nsum4 ) {
                    phid = DIFFGYP;
                } else if ( xext > nsum3 ) {
                    phid = DIFFC4A;
                } else if ( xext > nsum2 ) {
                    phid = DIFFC3A;
                } else if ( xext > nchext ) {
                    phid = DIFFCSH;
                }

                mic [ xc ] [ yc ] [ zc ] = phid;
                nmade += 1;
                ngoing += 1;
                antadd = ( struct ants * ) malloc( sizeof( struct ants ) );
                //printf("MEM addr %p\n", (const void *)antadd);
                antadd->x = xc;
                antadd->y = yc;
                antadd->z = zc;
                antadd->id = phid;
                antadd->cycbirth = cyccnt;
                /* Now connect this ant structure to end of linked list */
                antadd->prevant = tailant;
                tailant->nextant = antadd;
                antadd->nextant = NULL;
                tailant = antadd;
            }

            if ( tries > 10 * SYSIZE_POW3 ) {
#ifdef PRINTF
                printf("Could not place extra diffusing species (too few POROSITY voxels), continuing (line %d),\n", __LINE__);
#endif
                continue;
            }
        } while ( plok == 0 );
    }

    /* end of xext for extra species generation */

#ifdef PRINTF
    printf("Dissolved- %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n", count [ DIFFCSH ],
           count [ DIFFCH ], count [ DIFFGYP ], count [ DIFFC3A ], count [ DIFFFH3 ],
           count [ DIFFETTR ], count [ DIFFAS ], count [ DIFFANH ], count [ DIFFHEM ],
           count [ DIFFCAS2 ], count [ DIFFCACL2 ], count [ DIFFCACO3 ]);
#endif

    sulf_cur = count [ DIFFGYP ] + count [ DIFFANH ] + count [ DIFFHEM ];
    /*      difffile=fopen("diffuse.out","a");
     *   fprintf(difffile,"%d %ld %f %f %f %f %f %f %f\n",cycle, sulf_cur, cs_acc, ca_acc, disprob[C3S], disprob[C3A], disprob[C4AF], dfact, dfact1);
     *
     *   fclose(difffile); */

    /* if too many diffusing gypsums already in solution */
    if ( sulf_cur > DGYPMAX ) {
        disprob [ GYPSUM ] = disprob [ GYPSUMS ] = 0.0;
    } else {
        disprob [ GYPSUM ] = disbase [ GYPSUM ];
        disprob [ ANHYDRITE ] = disbase [ ANHYDRITE ];
        disprob [ HEMIHYD ] = disbase [ HEMIHYD ];
        disprob [ GYPSUMS ] = disbase [ GYPSUMS ];
    }

#ifdef PRINTF
    printf("C3AH6 dissolved- %ld with prob. of %f \n", nhgd, disprob [ C3AH6 ]);
#endif

    fflush(stdout);
}
/* routine to add nneed one pixel elements of phase randid at random */
/* locations in microstructure */
/* Special features for addition of 1-pixel CACO3 and INERT particles */
/* added 5/26/2004 */
/* Called by main program */
/* Calls no other routines */
void CemhydMatStatus :: addrand(int randid, long int nneed) {
    int ix, iy, iz;
    long int ic;
    int success, cpores;

    /* Add number of requested phase pixels at random pore locations */
    for ( ic = 1; ic <= nneed; ic++ ) {
        success = 0;
        while ( success == 0 ) {
            ix = ( int ) ( ( float ) SYSIZE * ran1(seed) );
            iy = ( int ) ( ( float ) SYSIZE * ran1(seed) );
            iz = ( int ) ( ( float ) SYSIZE * ran1(seed) );
            if ( ix == SYSIZE ) {
                ix = 0;
            }

            if ( iy == SYSIZE ) {
                iy = 0;
            }

            if ( iz == SYSIZE ) {
                iz = 0;
            }

            if ( mic [ ix ] [ iy ] [ iz ] == POROSITY ) {
                if ( ( randid != CACO3 ) && ( randid != INERT ) ) {
                    mic [ ix ] [ iy ] [ iz ] = randid;
                    micorig [ ix ] [ iy ] [ iz ] = randid;
                    success = 1;
                } else {
                    cpores = countboxc(3, ix, iy, iz);
                    if ( cpores >= 26 ) {
                        mic [ ix ] [ iy ] [ iz ] = randid;
                        micorig [ ix ] [ iy ] [ iz ] = randid;
                        success = 1;
                    }
                }
            }
        }
    }
}
/* Routine measuresurf to measure initial surface counts for cement */
/* and for all phases (cement= C3S, C2S, C3A, C4AF, and calcium sulfates */
void CemhydMatStatus :: measuresurf(void)
{
    int sx, sy, sz, jx, jy, jz, faceid;

    for ( sx = 0; sx < SYSIZE; sx++ ) {
        for ( sy = 0; sy < SYSIZE; sy++ ) {
            for ( sz = 0; sz < SYSIZE; sz++ ) {
                if ( mic [ sx ] [ sy ] [ sz ] == POROSITY ) {
                    for ( faceid = 0; faceid < 6; faceid++ ) {
                        if ( faceid == 1 ) {
                            jx = sx - 1;
                            if ( jx < 0 ) {
                                jx = SYSIZE - 1;
                            }

                            jy = sy;
                            jz = sz;
                        } else if ( faceid == 0 ) {
                            jx = sx + 1;
                            if ( jx > ( SYSIZE - 1 ) ) {
                                jx = 0;
                            }

                            jy = sy;
                            jz = sz;
                        } else if ( faceid == 2 ) {
                            jy = sy + 1;
                            if ( jy > ( SYSIZE - 1 ) ) {
                                jy = 0;
                            }

                            jx = sx;
                            jz = sz;
                        } else if ( faceid == 3 ) {
                            jy = sy - 1;
                            if ( jy < 0 ) {
                                jy = SYSIZE - 1;
                            }

                            jx = sx;
                            jz = sz;
                        } else if ( faceid == 4 ) {
                            jz = sz + 1;
                            if ( jz > ( SYSIZE - 1 ) ) {
                                jz = 0;
                            }

                            jx = sx;
                            jy = sy;
                        } else if ( faceid == 5 ) {
                            jz = sz - 1;
                            if ( jz < 0 ) {
                                jz = SYSIZE - 1;
                            }

                            jx = sx;
                            jy = sy;
                        }

                        /* If the neighboring pixel is solid, update surface counts */
                        if ( ( mic [ jx ] [ jy ] [ jz ] == C3S ) || ( mic [ jx ] [ jy ] [ jz ] == C2S ) || ( mic [ jx ] [ jy ] [ jz ] == C3A ) || ( mic [ jx ] [ jy ] [ jz ] == C4AF ) || ( mic [ jx ] [ jy ] [ jz ] == INERT ) || ( mic [ jx ] [ jy ] [ jz ] == CACO3 ) ) {
                            scnttotal += 1;
                            if ( ( mic [ jx ] [ jy ] [ jz ] == C3S ) || ( mic [ jx ] [ jy ] [ jz ] == C2S ) || ( mic [ jx ] [ jy ] [ jz ] == C3A ) || ( mic [ jx ] [ jy ] [ jz ] == C4AF ) ) {
                                scntcement += 1;
                            }
                        }
                    }
                }
            }
        }
    }

#ifdef PRINTF
    printf("Cement surface count is %ld \n", scntcement);
    printf("Total surface count is %ld \n", scnttotal);
#endif

    surffract = ( float ) scntcement / ( float ) scnttotal;

#ifdef PRINTF
    printf("Surface fraction is %f \n", surffract);
#endif

    fflush(stdout);
}

/* Routine to resaturate all empty porosity */
/* and continue with hydration under saturated conditions */
void CemhydMatStatus :: resaturate(void)
{
    int sx, sy, sz;
    long int nresat = 0;

    for ( sx = 0; sx < SYSIZE; sx++ ) {
        for ( sy = 0; sy < SYSIZE; sy++ ) {
            for ( sz = 0; sz < SYSIZE; sz++ ) {
                if ( mic [ sx ] [ sy ] [ sz ] == EMPTYP ) {
                    mic [ sx ] [ sy ] [ sz ] = POROSITY;
                    nresat++;
                }
            }
        }
    }

    if ( nresat > 0 ) {
        porefl1 = porefl2 = porefl3 = 1;
    }

    printf("Number resaturated is %ld \n", nresat);
    fflush(stdout);
}


/*called from main loop*/
void CemhydMatStatus :: outputImageFileUnperc(char ***m) {
    FILE *imgout;
    char extension [ 10 ];
    //char outputname[80];
    char *prefix;
    prefix = ( char * ) malloc(80);

    //mkdir("unperc", 0777);//make directory (every time)
    //system("mk unperc");
    //system("mkdir unperc 2> /dev/null");
    fprintf(infoUnperc, "%d %.4f %.3f\n", icyc, alpha_cur, time_cur);
    fflush(infoUnperc);

    sprintf(extension, "%04d", icyc);
    strcpy(prefix, "unperc/out5."); //see line with mkdir in order to be the same
    strcat(prefix, extension);
    strcat(prefix, ".img");
#ifdef PRINTF
    printf("Name of output file is %s\n", prefix);
#endif

    if ( ( imgout = fopen(prefix, "w") ) == NULL ) {
        printf("File %s can not be opened\n", prefix);
        free(prefix);
        return;
    }

    for ( int dz = 0; dz < SYSIZE; dz++ ) {
        for ( int dy = 0; dy < SYSIZE; dy++ ) {
            for ( int dx = 0; dx < SYSIZE; dx++ ) {
                fprintf(imgout, "%d\n", m [ dx ] [ dy ] [ dz ]);
            }
        }
    }

    fclose(imgout);
#ifdef PRINTF
    printf("Unpercolated file %s wrote\n", prefix);
#endif
    free(prefix);
}


void CemhydMatStatus :: readhydrparam()
{
    int fidc3s, fidc2s, fidc3a, fidc4af, fidgyp, fidagg, ffac3a;
    int fidhem, fidanh, fidcaco3;
    int read_micr;
    long int valin;
    int ix, iy, iz, phtodo;
    long int nadd;
    char filei [ 80 ];
    FILE *infile;

    /* Get random number seed */
    //printf("Enter random number seed \n");
    //fscanf(in, "%d",&nseed);
    nseed = iseed;
    seed = ( & nseed );

#ifdef PRINTF
    printf("Seed %d\n", * seed);
#endif

    //check whether microstructure from file should be used
#ifdef CMLFILE
    F->get_value(13, ( long & )read_micr);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Given_microstructure", 0, read_micr);
#endif

    /* Open file and read in original cement particle microstructure if required*/
    if ( read_micr ) {
#ifdef PRINTF
        printf("Enter name of file to read initial microstructure from \n");
#endif
#ifdef CMLFILE
        F->get_value(1, filei);
#endif
#ifdef TINYXML
        QueryStringAttributeExt(xmlFile, "Input_img_file", 0, filei);
#endif
        //fscanf(in, "%s",filei);

#ifdef PRINTF
        printf("%s\n", filei);
        fflush(stdout);
#endif
    }

    /* Get phase assignments for original microstructure */
    /* to transform to needed ID values */

#ifdef PRINTF
    printf("Enter IDs in file for C3S, C2S, C3A, C4AF, Gypsum, Hemihydrate, Anhydrite, Aggregate, CaCO3\n");
#endif

    fidc3s = C3S;
    fidc2s = C2S;
    fidc3a = C3A;
    fidc4af = C4AF;
    fidgyp = GYPSUM;
    fidhem = HEMIHYD;
    fidanh = ANHYDRITE;
    fidagg = INERTAGG;
    fidcaco3 = CACO3;

    //scanf("%d %d %d %d %d %d %d %d %d",&fidc3s,&fidc2s,&fidc3a,&fidc4af,&fidgyp,&fidhem,&fidanh,&fidagg,&fidcaco3);
#ifdef PRINTF
    printf("%d %d %d %d %d %d %d %d %d\n", fidc3s, fidc2s, fidc3a, fidc4af, fidgyp, fidhem, fidanh, fidagg, fidcaco3);
    printf("Enter ID in file for C3A in fly ash (default=35)\n");
#endif
    //fscanf(in, "%d",&ffac3a);
    ffac3a = 35;
#ifdef PRINTF
    printf("%d\n", ffac3a);
#endif
    fflush(stdout);

    if ( read_micr ) {
        if ( ( infile = fopen(filei, "r") ) == NULL ) {
            printf("CEMHYD3D microstructure file %s not found\n", filei);
            exit(0);
        }

        for ( iz = 0; iz < SYSIZE; iz++ ) { //fixed cycle order
            for ( iy = 0; iy < SYSIZE; iy++ ) {
                for ( ix = 0; ix < SYSIZE; ix++ ) {
                    cshage [ ix ] [ iy ] [ iz ] = 0;
                    faces [ ix ] [ iy ] [ iz ] = 0;

                    if ( fscanf(infile, "%ld", & valin) == EOF ) {
                        printf("End of file %s reached, terminating\n", filei);
                        exit(1);
                    }

                    if ( valin < 0 ) {
                        printf("Error in the reading at x=%d y=%d z=%d, value %ld\n", ix, iy, iz, valin);
                        exit(1);
                    }

                    mic [ ix ] [ iy ] [ iz ] = valin;
                    if ( valin == fidc3s ) {
                        mic [ ix ] [ iy ] [ iz ] = C3S;
                    } else if ( valin == fidc2s ) {
                        mic [ ix ] [ iy ] [ iz ] = C2S;
                    } else if ( ( valin == fidc3a ) || ( valin == ffac3a ) ) {
                        mic [ ix ] [ iy ] [ iz ] = C3A;
                    } else if ( valin == fidc4af ) {
                        mic [ ix ] [ iy ] [ iz ] = C4AF;
                    } else if ( valin == fidgyp ) {
                        mic [ ix ] [ iy ] [ iz ] = GYPSUM;
                    } else if ( valin == fidanh ) {
                        mic [ ix ] [ iy ] [ iz ] = ANHYDRITE;
                    } else if ( valin == fidhem ) {
                        mic [ ix ] [ iy ] [ iz ] = HEMIHYD;
                    } else if ( valin == fidcaco3 ) {
                        mic [ ix ] [ iy ] [ iz ] = CACO3;
                    } else if ( valin == fidagg ) {
                        mic [ ix ] [ iy ] [ iz ] = INERTAGG;
                    }

                    micorig [ ix ] [ iy ] [ iz ] = mic [ ix ] [ iy ] [ iz ];
                }
            }
        }

        fclose(infile);
    } else {
#ifdef PRINTF
        printf("Generated microstructure from memory will be used\n");
#endif
    }

    fflush(stdout);

    /* Now read in particle IDs from file if required*/
    if ( read_micr ) {
#ifdef PRINTF
        printf("Enter name of file to read particle IDs from \n");
#endif
#ifdef CMLFILE
        F->get_value(2, filei);
#endif
#ifdef TINYXML
        QueryStringAttributeExt(xmlFile, "Input_id_file", 0, filei);
#endif
        //fscanf(in, "%s",filei);

#ifdef PRINTF
        printf("%s\n", filei);
#endif

        if ( ( infile = fopen(filei, "r") ) == NULL ) {
            printf("CEMHYD3D microstructure ID file %s not found\n", filei);
            exit(0);
        }

        for ( iz = 0; iz < SYSIZE; iz++ ) { //fixed cycle order
            for ( iy = 0; iy < SYSIZE; iy++ ) {
                for ( ix = 0; ix < SYSIZE; ix++ ) {
                    if ( fscanf(infile, "%ld", & valin) == EOF ) {
                        printf("End of file %s reached, terminating\n", filei);
                        exit(1);
                    }

                    micpart [ ix ] [ iy ] [ iz ] = valin;
                }
            }
        }

        fclose(infile);
    } else {
#ifdef PRINTF
        printf("Generated microstructure ID from memory will be used\n");
#endif
    }


    fflush(stdout);


    /* Allow user to iteratively add one pixel particles of various phases */
    /* Typical application would be for addition of silica fume */
#ifdef PRINTF
    printf("Enter number of one pixel particles to add (0 to quit) \n");
#endif
    //fscanf(in, "%ld",&nadd);
    nadd = 0;

#ifdef PRINTF
    printf("%ld\n", nadd);
#endif

    while ( nadd > 0 ) {
#ifdef PRINTF
        printf("Enter phase to add \n");
        printf(" C3S 1 \n");
        printf(" C2S 2 \n");
        printf(" C3A 3 \n");
        printf(" C4AF 4 \n");
        printf(" GYPSUM 5 \n");
        printf(" HEMIHYD 6 \n");
        printf(" ANHYDRITE 7 \n");
        printf(" POZZ 8 \n");
        printf(" INERT 9 \n");
        printf(" SLAG 10 \n");
        printf(" ASG 11 \n");
        printf(" CAS2 12 \n");
        printf(" CH 13 \n");
        printf(" CSH 14 \n");
        printf(" C3AH6 15 \n");
        printf(" Ettringite 16 \n");
        printf(" Stable Ettringite from C4AF 17 \n");
        printf(" AFM 18 \n");
        printf(" FH3 19 \n");
        printf(" POZZCSH 20 \n");
        printf(" SLAGCSH 21 \n");
        printf(" CACL2 22 \n");
        printf(" Friedels salt 23 \n");
        printf(" Stratlingite 24 \n");
        printf(" Calcium carbonate 26 \n");
#endif
        //fscanf(in, "%d",&phtodo);
        phtodo = 0;
#ifdef PRINTF
        printf("%d \n", phtodo);
#endif
        if ( ( phtodo < 0 ) || ( phtodo > CACO3 ) ) {
            printf("Error in phase input for one pixel particles \n");
            exit(1);
        }

        addrand(phtodo, nadd);

#ifdef PRINTF
        printf("Enter number of one pixel particles to add (0 to quit) \n");
#endif

        //fscanf(in, "%ld",&nadd);
        nadd = 0;
#ifdef PRINTF
        printf("%ld\n", nadd);
#endif
    }

    fflush(stdout);

    //printf("Do you wish hydration under 0) saturated or 1) sealed conditions \n");
#ifdef CMLFILE
    F->get_value(3, ( long & )sealed);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Saturated_sealed", 0, sealed);
#endif
    //fscanf(in, "%d",&sealed);
    //printf("%d \n",sealed);
    //printf("Enter max. # of diffusion steps per cycle (500) \n");
    //fscanf(in, "%d",&ntimes);
#ifdef CMLFILE
    F->get_value(23, ( long & )ntimes);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Diffusion_steps_per_cycle", 0, ntimes);
#endif
    //ntimes=500;
    //printf("%d \n",ntimes);
    //printf("Enter nuc. prob. and scale factor for CH nucleation \n");
    //fscanf(in, "%f %f",&pnucch,&pscalech);
#ifdef CMLFILE
    F->get_value(24, pnucch);
    F->get_value(25, pscalech);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "CH_nucleation_probability", 0, pnucch);
    QueryNumAttributeExt(xmlFile, "CH_scale_factor", 0, pscalech);
#endif
    //pnucch=0.0001;
    //pscalech=9000.;
    //printf("%f %f \n",pnucch,pscalech);
    //printf("Enter nuc. prob. and scale factor for gypsum nucleation \n");
    //fscanf(in, "%f %f",&pnucgyp,&pscalegyp);
#ifdef CMLFILE
    F->get_value(26, pnucgyp);
    F->get_value(27, pscalegyp);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Gypsum_nucleation_probability", 0, pnucgyp);
    QueryNumAttributeExt(xmlFile, "Gypsum_scale_factor", 0, pscalegyp);
#endif
    //pnucgyp=0.01;
    //pscalegyp=9000.;
    //printf("%f %f \n",pnucgyp,pscalegyp);
    //printf("Enter nuc. prob. and scale factor for C3AH6 nucleation \n");
    //fscanf(in, "%f %f",&pnuchg,&pscalehg);
#ifdef CMLFILE
    F->get_value(28, pnuchg);
    F->get_value(29, pscalehg);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "C3AH6_nucleation_probability", 0, pnuchg);
    QueryNumAttributeExt(xmlFile, "C3AH6_scale_factor", 0, pscalehg);
#endif
    //pnuchg=0.00002;
    //pscalehg=10000.;
    //printf("%f %f \n",pnuchg,pscalehg);
    //printf("Enter nuc. prob. and scale factor for FH3 nucleation \n");
    //fscanf(in, "%f %f",&pnucfh3,&pscalefh3);
#ifdef CMLFILE
    F->get_value(30, pnucfh3);
    F->get_value(31, pscalefh3);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "FH3_nucleation_probability", 0, pnucfh3);
    QueryNumAttributeExt(xmlFile, "FH3_scale_factor", 0, pscalefh3);
#endif
    //pnucfh3=0.002;
    //pscalefh3=2500.;
    //printf("%f %f \n",pnucfh3,pscalefh3);
    //printf("Enter cycle frequency for checking pore space percolation \n");
    //fscanf(in, "%d",&burnfreq);
#ifdef CMLFILE
    F->get_value(17, ( long & )burnfreq);
    F->get_value(33, adiabatic_curing);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Cycle_freq_perc_pore", 0, burnfreq);
    QueryNumAttributeExt(xmlFile, "Adiabatic_conditions", 0, adiabatic_curing);
#endif

#ifdef __TM_MODULE //OOFEM transport module
    adiabatic_curing = 0;
#endif

#ifdef PRINTF
    printf("Enter cycle frequency for checking percolation of solids (set) \n");
    printf("%d\n", burnfreq);
    printf("Adiabatic curing conditions %d\n", adiabatic_curing);
#endif
    //fscanf(in, "%d",&setfreq);
#ifdef CMLFILE
    F->get_value(18, ( long & )setfreq);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Cycle_freq_perc_sol", 0, setfreq);
#endif
    //setfreq=50000;
#ifdef PRINTF
    printf("Enter cycle frequency for checking hydration of particles \n");
    printf("%d\n", setfreq);
#endif
    //fscanf(in, "%d",&phydfreq);
    phydfreq = 50000; //disabled
    //printf("%d\n",phydfreq);
    //printf("Enter cycle frequency for outputting hydrating microstructure \n");
    //fscanf(in, "%d",&outfreq);
    outfreq = 50000;
    //printf("%d\n",outfreq);
    //printf("Enter the induction time in hours \n");
#ifdef CMLFILE
    F->get_value(4, ind_time);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Induction_time", 0, ind_time);
#endif
    //fscanf(in, "%lf",&ind_time);
    //printf("%f \n",ind_time);
    time_cur += ind_time;
    /* heat transfer coefficient in J/g/C/s - not used*/
    //fscanf(in, "%f",&U_coeff);
    U_coeff = 0;
#ifdef PRINTF
    printf("Enter apparent activation energy for hydration in kJ/mole \n");
#endif
#ifdef CMLFILE
    F->get_value(5, E_act);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Ea_cement", 0, E_act);
#endif
    //fscanf(in, "%lf",&E_act);
#ifdef PRINTF
    printf("%f \n", E_act);
    printf("Enter apparent activation energy for pozzolanic reactions in kJ/mole \n");
#endif
#ifdef CMLFILE
    F->get_value(6, E_act_pozz);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Ea_pozz", 0, E_act_pozz);
#endif
    //fscanf(in, "%f",&E_act_pozz);
#ifdef PRINTF
    printf("%f \n", E_act_pozz);
    printf("Enter apparent activation energy for slag reactions in kJ/mole \n");
#endif
#ifdef CMLFILE
    F->get_value(7, E_act_slag);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Ea_slag", 0, E_act_slag);
#endif
    //fscanf(in, "%f",&E_act_slag);
#ifdef PRINTF
    printf("%f \n", E_act_slag);
    printf("Enter kinetic factor to convert cycles to time for 25 C \n");
#endif
#ifdef CMLFILE
    F->get_value(8, beta);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Beta", 0, beta);
#endif
    //fscanf(in, "%lf",&beta);
#ifdef PRINTF
    printf("%f \n", beta);
    printf("Enter mass fraction of aggregate in concrete \n");
#endif
#ifdef CMLFILE
    F->get_value(9, mass_agg);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Mass_SCM_FA_CA_inert_frac", 0, mass_agg);
#endif
    //fscanf(in, "%f",&mass_agg);
#ifdef PRINTF
    printf("%f \n", mass_agg);
    printf("Enter kg of cement (+admixtures) in 1 m3\n");
#endif
#ifdef CMLFILE
    F->get_value(10, Mass_cement_concrete);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Mass_cement_concrete", 0, Mass_cement_concrete);
#endif
    //fscanf(in, "%f",&Mass_cement_concrete);
#ifdef PRINTF
    printf("%f \n", Mass_cement_concrete);
    printf("Enter heat capacity of aggregate in concrete J/g/C\n");
#endif
#ifdef CMLFILE
    F->get_value(11, Cp_agg);
    Cp_agg /= 1000.; //scale to J/g/K
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Cp_SCM_FA_CA_inert", 0, Cp_agg);
    Cp_agg /= 1000.; //scale to J/g/K
#endif
    //fscanf(in, "%f",&Cp_agg);
#ifdef PRINTF
    printf("%f \n", Cp_agg);
    printf("Enter heat capacity of cement J/g/C\n");
#endif
#ifdef CMLFILE
    F->get_value(12, Cp_cement);
    Cp_cement /= 1000.; //scale to J/g/K
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Cp_cem", 0, Cp_cement);
    Cp_cement /= 1000.; //scale to J/g/K
#endif
    //fscanf(in, "%f",&Cp_cement);
#ifdef PRINTF
    printf("%f \n", Cp_cement);
    printf("CSH to pozzolanic CSH 0) prohibited or 1) allowed \n");
#endif
    //fscanf(in, "%d",&csh2flag);
    csh2flag = 0;
#ifdef PRINTF
    printf("%d \n", csh2flag);
    printf("CH precipitation on aggregate surfaces 0) prohibited or 1) allowed \n");
#endif
    //fscanf(in, "%d",&chflag);
    chflag = 0;
#ifdef PRINTF
    printf("%d \n", chflag);
    printf("Enter number of cycles before executing total resaturation \n");
#endif
    //fscanf(in, "%d\n",&resatcyc);
    resatcyc = 0;
#ifdef PRINTF
    printf("%d\n", resatcyc);
    printf("Enter choice for C-S-H geometry 0) random or 1) plates \n");
#endif
    //fscanf(in, "%d\n",&cshgeom);
    cshgeom = 0;
#ifdef PRINTF
    printf("%d \n", cshgeom);
    printf("Does pH influence hydration kinetics 0) no or 1) yes \n");
#endif
    //fscanf(in, "%d\n",&pHactive);
    pHactive = 0;
#ifdef PRINTF
    printf("%d\n", pHactive);
#endif
#ifdef CMLFILE
    F->get_value(34, Vol_cement_clinker_gypsum);
    F->get_value(35, Vol_cement_SCM);
    F->get_value(36, Vol_water);
    F->get_value(37, Vol_FA);
    F->get_value(38, Vol_CA);
    F->get_value(39, Vol_inert_filler);
    F->get_value(40, Vol_entrained_entrapped_air);
    F->get_value(41, Grain_average_FA);
    F->get_value(42, Grain_average_CA);
    F->get_value(43, ITZ_thickness);
    F->get_value(44, ITZ_Young_red);
    F->get_value(45, Young_SCM);
    F->get_value(46, Poisson_SCM);
    F->get_value(47, Young_FA);
    F->get_value(48, Poisson_FA);
    F->get_value(49, Young_CA);
    F->get_value(50, Poisson_CA);
    F->get_value(51, Young_inert);
    F->get_value(52, Poisson_inert);
    F->get_value(53, Calculate_elastic_homogenization);
#endif
#ifdef TINYXML
    QueryNumAttributeExt(xmlFile, "Vol_cement_clinker_gypsum", 0, Vol_cement_clinker_gypsum);
    QueryNumAttributeExt(xmlFile, "Vol_cement_SCM", 0, Vol_cement_SCM);
    QueryNumAttributeExt(xmlFile, "Vol_water", 0, Vol_water);
    QueryNumAttributeExt(xmlFile, "Vol_FA", 0, Vol_FA);
    QueryNumAttributeExt(xmlFile, "Vol_CA", 0, Vol_CA);
    QueryNumAttributeExt(xmlFile, "Vol_inert_filler", 0, Vol_inert_filler);
    QueryNumAttributeExt(xmlFile, "Vol_entrained_entrapped_air", 0, Vol_entrained_entrapped_air);
    QueryNumAttributeExt(xmlFile, "Grain_average_FA", 0, Grain_average_FA);
    QueryNumAttributeExt(xmlFile, "Grain_average_CA", 0, Grain_average_CA);
    QueryNumAttributeExt(xmlFile, "ITZ_thickness", 0, ITZ_thickness);
    QueryNumAttributeExt(xmlFile, "ITZ_Young_red", 0, ITZ_Young_red);
    QueryNumAttributeExt(xmlFile, "Young_SCM", 0, Young_SCM);
    QueryNumAttributeExt(xmlFile, "Poisson_SCM", 0, Poisson_SCM);
    QueryNumAttributeExt(xmlFile, "Young_FA", 0, Young_FA);
    QueryNumAttributeExt(xmlFile, "Poisson_FA", 0, Poisson_FA);
    QueryNumAttributeExt(xmlFile, "Young_CA", 0, Young_CA);
    QueryNumAttributeExt(xmlFile, "Poisson_CA", 0, Poisson_CA);
    QueryNumAttributeExt(xmlFile, "Young_inert", 0, Young_inert);
    QueryNumAttributeExt(xmlFile, "Poisson_inert", 0, Poisson_inert);
    QueryNumAttributeExt(xmlFile, "Calculate_elastic_homogenization", 0, Calculate_elastic_homogenization);
    //needed in OOFEM for conductivity and capacity calculations
    QueryNumAttributeExt(xmlFile, "Mass_tot_concrete", 0, Mass_tot_concrete);
    QueryNumAttributeExt(xmlFile, "Cp_SCM", 0, Cp_SCM);
    QueryNumAttributeExt(xmlFile, "Cp_FA", 0, Cp_FA);
    QueryNumAttributeExt(xmlFile, "Cp_CA", 0, Cp_CA);
    QueryNumAttributeExt(xmlFile, "Cp_inert", 0, Cp_inert);
    QueryNumAttributeExt(xmlFile, "Mass_SCM_frac", 0, Mass_SCM_frac);
    QueryNumAttributeExt(xmlFile, "Mass_FA_frac", 0, Mass_FA_frac);
    QueryNumAttributeExt(xmlFile, "Mass_CA_frac", 0, Mass_CA_frac);
    QueryNumAttributeExt(xmlFile, "Mass_inert_frac", 0, Mass_inert_frac);
    QueryNumAttributeExt(xmlFile, "Concrete_thermal_conductivity", 0, Concrete_thermal_conductivity); //[W/m/K]
    QueryNumAttributeExt(xmlFile, "Concrete_bulk_density", 0, Concrete_bulk_density); //[kg/m3]
#endif
}


/* Calls init, dissolve and addrand */
void CemhydMatStatus :: disrealnew_init(void)
{
#ifdef OUTFILES
    char fileo [ 80 ];
#endif

    nseed = iseed; //set random number generator the same number
    seed = & nseed;

    ncshplategrow = 0;
    ncshplateinit = 0;
    slagemptyp = 0;
    countpore = 0;
    time_step = 0.0;
    w_to_c = 0.0;
    totfract = 1.0;
    tfractw04 = 0.438596;
    fractwithfill = 1.0;
    tfractw05 = 0.384615;
    surffract = 0.0;
    pfractw05 = 0.615385;
    scntcement = 0;
    scnttotal = 0;
    saturation = 1.0;
    cs_acc = 1.0;
    ca_acc = 1.0;
    dismin_c3a = DISMIN_C3A_0;
    dismin_c4af = DISMIN_C4AF_0;
    gsratio2 = 0.0;
    onepixelbias = 1.0;
    DIFFCHdeficit = 0;
    slaginit = 0;
    slagcum = 0;
    chgone = 0;
    nch_slag = 0;
    sulf_cur = 0;
    pHfactor = 0.0;
    conccaplus = 0.0;
    moles_syn_precip = 0.0;
    concsulfate = 0.0;

    primevalues [ 0 ] = 2;
    primevalues [ 1 ] = 3;
    primevalues [ 2 ] = 5;
    primevalues [ 3 ] = 7;
    primevalues [ 4 ] = 11;
    primevalues [ 5 ] = 13;

    ngoing = 0;
    porefl1 = porefl2 = porefl3 = 1;
    pore_off = water_off = 0;
    cycflag = 0;
    heat_old = heat_new = 0.0;
    chold = chnew = 0;     /* Current and previous cycle CH counts */
    cubesize = CUBEMAX;
    ppozz = PPOZZ;
    poregone = poretodo = 0;

    /* Initialize counters, etc. */
    npr = nasr = nslagr = 0;
    nfill = 0;
    ncsbar = 0;
    netbar = 0;
    porinit = 0;
    cyccnt = 0;
    setflag = 0;
    c3sinit = c2sinit = c3ainit = c4afinit = anhinit = heminit = slaginit = 0;

    /* Initialize structure for ants */
    headant = ( struct ants * ) malloc( sizeof( struct ants ) );
    headant->prevant = NULL;
    headant->nextant = NULL;
    headant->x = 0;
    headant->y = 0;
    headant->z = 0;
    headant->id = 100;   /* special ID indicating first ant in list */
    headant->cycbirth = 0;
    tailant = headant;

    LastHydrTime = 0.;
    LastCycHeat = 0.;
    LastTotHeat = 0.;
    LastCycCnt = 0;
    PrevHydrTime = 0.;
    PrevCycHeat = 0.;


    init();
#ifdef PRINTF
    printf("After init routine \n");
#endif
    fflush(stdout);



#ifdef OUTFILES
    fileperc = fopen("percpore.out", "w"); //in burn3d.cpp
    //fprintf(fileperc,"Cycle time(h) alpha_mass conn_por total_por frac_conn\n");
    fprintf(fileperc, "#cycle alpha #through #tot_phase\n");
    disprobfile = fopen("disprob.out", "w");
    //percfile = fopen("percset.out", "w"); //in burnset.cpp
    //fprintf(percfile, "#cycle time alpha #through #C3S+C2S+C3A+C4AF+CAS2+SLAG+ASG+POZZ+ETTR+ETTRC4AF+CSH conn_frac\n");
    strcpy(heatname, "heat.out");
    heatfile = fopen(heatname, "w");
    strcpy(chshrname, "chemshr.out");
    //chsfile=fopen(chshrname,"w");
    strcpy(adianame, "adiabatic.out");
    strcpy(pHname, "pHfile.out");
    pHfile = fopen(pHname, "w");
    fprintf(pHfile, "Cycle time(h) alpha_mass pH sigma [Na+] [K+] [Ca++] [SO4--] activityCa activityOH activitySO4 activityK molesSyngenite\n");
    strcpy(fileo, "out5.img");
    strcpy(phasname, "phases.out");
    phasfile = fopen(phasname, "w");
    strcpy(ppsname, "percpore.out");
    strcpy(phrname, "parthydr.out");
    perc_phases = fopen("perc_phases.out", "w");
    adiafile = fopen(adianame, "w");
    fprintf(adiafile, "Cyc Time(h) Temperature  Alpha  Krate   Cp_now  Mass_cem kpozz/khyd kslag/khyd DoR_Blend\n");
    elasfile = fopen("elas.out", "w");
    fprintf(elasfile, "Cyc\tTime[h] Alpha E_paste nu_paste E_paste_fil nu_paste_fil E_mortar nu_mortar E_concrete nu_concrete\n");
    CSHfile = fopen("CSH.out", "w");
    fprintf(CSHfile, "cyc alpha [h] HD CSH-total HD/CSH\n");
    system("mkdir perc 2> /dev/null");
    system("mkdir unperc 2> /dev/null");
#endif

#ifdef IMAGEFILES
    infoperc = fopen("perc/info.dat", "w");
    fprintf(infoperc, "#Cyc DoH Time[h]\n");
    infoUnperc = fopen("unperc/info.dat", "w");
    fprintf(infoUnperc, "#Cyc DoH Time[h]\n");
#endif
    /* Set initial properties of CSH */
    molarvcsh [ 0 ] = molarv [ CSH ];
    watercsh [ 0 ] = waterc [ CSH ];

    krate = exp( -( 1000. * E_act / 8.314 ) * ( ( 1. / ( temp_cur + 273.15 ) ) - ( 1. / 298.15 ) ) );
    /* Determine pozzolanic and slag reaction rate constants */
    kpozz = exp( -( 1000. * E_act_pozz / 8.314 ) * ( ( 1. / ( temp_cur + 273.15 ) ) - ( 1. / 298.15 ) ) );
    kslag = exp( -( 1000. * E_act_slag / 8.314 ) * ( ( 1. / ( temp_cur + 273.15 ) ) - ( 1. / 298.15 ) ) );
    ppozz = PPOZZ * kpozz / krate;
    /* Assume same holds for dissolution of fly ash phases */
    disprob [ ASG ] = disbase [ ASG ] * kpozz / krate;
    disprob [ CAS2 ] = disbase [ CAS2 ] * kpozz / krate;
    /* Update probability of slag dissolution */
    disprob [ SLAG ] = slagreact * disbase [ SLAG ] * kslag / krate;


    /* Determine surface counts */
    measuresurf();
}


/*perform as many hydration cycles as necessary
 * flag = 0 - controlled with hydrationTime
 * flag != 0 - perform exactly flag cycles
 * GiveTemp is effective in the first cycle of adiabatic regime only
 */
void CemhydMatStatus :: disrealnew(double GiveTemp, double hydrationTime, int flag) {
    int counter = 0;

    if ( adiabatic_curing != 1 ) {
        temp_cur = GiveTemp;
        temp_0 = GiveTemp;
    }

    InitTime = time_cur;

    if ( icyc == 1 ) { //first cycle
        temp_cur = GiveTemp;
        temp_0 = GiveTemp;
        disrealnew_init();
    }

    /****************************************************************/
    /**************************MAIN LOOP*****************************/
    /****************************************************************/

    /*determine the minimum number of cycles*/
    while ( flag == 0 ? ( time_cur < hydrationTime ) : counter < flag ) {
        alpha_last = alpha_cur;

        //printf("CNT seed=%d\n", *seed);
        //update counters of "previous", i.e. before any last cycle
        PrevHydrTime = time_cur * 3600.;
        PrevCycHeat = ( double ) ( heat_new ) * heat_cf;

        if ( ( sealed == 1 ) && ( icyc == ( resatcyc + 1 ) ) && ( resatcyc != 0 ) ) {
            resaturate();
            sealed = 0;
        }

        if ( temp_cur <= 80.0 ) {
            molarvcsh [ icyc ] = molarv [ CSH ] - 8.0 * ( ( temp_cur - 20. ) / ( 80. - 20. ) );
            watercsh [ icyc ] = waterc [ CSH ] - 1.3 * ( ( temp_cur - 20. ) / ( 80. - 20. ) );
        } else {
            molarvcsh [ icyc ] = molarv [ CSH ] - 8.0;
            watercsh [ icyc ] = waterc [ CSH ] - 1.3;
        }

        dissolve(icyc); //cyccnt is increased
#ifdef OUTFILES
        //printf("Number dissolved this pass- %ld total diffusing- %ld \n",nmade,ngoing);

        fflush(stdout);
        if ( icyc == 1 ) {
            // printf("ncsbar is %ld   netbar is %ld \n",ncsbar,netbar);
        }

#endif

        hydrate(cycflag, ntimes, pnucch, pscalech, pnuchg, pscalehg, pnucfh3, pscalefh3, pnucgyp, pscalegyp);
        temp_0 = temp_cur;

        /* Handle adiabatic case first */
        /* Cement + aggregate +water + filler=1;  that's all there is */
        mass_cement = 1. - mass_agg - mass_fill - mass_water - mass_CH;
        mass_cem_now = mass_cement;

        /* determine heat capacity of current mixture, */
        /* accounting for imbided water if necessary */
        if ( sealed == 1 ) {
            Cp_now = mass_agg * Cp_agg;
            Cp_now += Cp_pozz * mass_fill;
            Cp_now += Cp_cement * mass_cement;
            Cp_now += Cp_CH * mass_CH;
            Cp_now += ( Cp_h2o * mass_water - alpha_cur * WN * mass_cement * ( Cp_h2o - Cp_bh2o ) );
            mass_cem_now = mass_cement;
        }
        /* Else need to account for extra capillary water drawn in */
        /* Basis is WCHSH(0.06) g H2O per gram cement for chemical shrinkage */
        /* Need to adjust mass basis to account for extra imbibed H2O */
        else {
            mass_cur = 1. + WCHSH * mass_cement * alpha_cur;
            Cp_now = mass_agg * Cp_agg / mass_cur;
            Cp_now += Cp_pozz * mass_fill / mass_cur;
            Cp_now += Cp_cement * mass_cement / mass_cur;
            Cp_now += Cp_CH * mass_CH / mass_cur;
            Cp_now += ( Cp_h2o * mass_water - alpha_cur * WN * mass_cement * ( Cp_h2o - Cp_bh2o ) );
            Cp_now += ( WCHSH * Cp_h2o * alpha_cur * mass_cement );
            mass_cem_now = mass_cement / mass_cur;
        }


        /* Recall that temp_cur is in degrees Celsius */
        krate = exp( -( 1000. * E_act / 8.314 ) * ( ( 1. / ( temp_cur + 273.15 ) ) - ( 1. / 298.15 ) ) );
        /* Determine pozzolanic and slag reaction rate constant */
        kpozz = exp( -( 1000. * E_act_pozz / 8.314 ) * ( ( 1. / ( temp_cur + 273.15 ) ) - ( 1. / 298.15 ) ) );
        kslag = exp( -( 1000. * E_act_slag / 8.314 ) * ( ( 1. / ( temp_cur + 273.15 ) ) - ( 1. / 298.15 ) ) );
        /* Update probability of pozzolanic and slag reactions */
        /* based on ratio of pozzolanic (slag) reaction rate to hydration rate */
        ppozz = PPOZZ * kpozz / krate;
        disprob [ ASG ] = disbase [ ASG ] * kpozz / krate;
        disprob [ CAS2 ] = disbase [ CAS2 ] * kpozz / krate;
        disprob [ SLAG ] = slagreact * disbase [ SLAG ] * kslag / krate;

        /* Update temperature based on heat generated and current Cp only in adiabatic case*/
        if ( adiabatic_curing == 1 ) {
            if ( mass_cem_now > 0.01 ) {
                temp_cur = temp_0 + mass_cem_now * heat_cf * ( heat_new - heat_old ) / Cp_now;
            } else {
                temp_cur = temp_0 + mass_fill_pozz * heat_cf * ( heat_new - heat_old ) / Cp_now;
            }
        }

        /* Update system temperature due to heat loss/gain to/from */
        /* surroundings (semi-adiabatic case) */
        //temp_cur-=(temp_cur-T_ambient)*time_step*U_coeff/Cp_now;

        /* Update time based on simple numerical integration */
        /* simulating maturity approach */
        /* with parabolic kinetics (Knudsen model) */
        if ( cyccnt > 1 ) {
            time_cur += ( 2. * ( double ) ( cyccnt - 1 ) - 1.0 ) * beta / krate;
            time_step = ( 2. * ( double ) ( cyccnt - 1 ) - 1.0 ) * beta / krate;
        }

        gsratio2 = 0.0;
        gsratio2 += ( float ) ( count [ CH ] + count [ CSH ] + count [ C3AH6 ] + count [ ETTR ] );
        gsratio2 += ( float ) ( count [ POZZCSH ] + count [ SLAGCSH ] + count [ FH3 ] + count [ AFM ] + count [ ETTRC4AF ] );
        gsratio2 += ( float ) ( count [ FREIDEL ] + count [ STRAT ] + count [ ABSGYP ] + count [ AFMC ] );
        gsratio2 = ( gsratio2 ) / ( gsratio2 + ( float ) ( count [ POROSITY ] + count [ EMPTYP ] ) );


#ifdef OUTFILES
        fprintf(adiafile, "%d \t %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", cyccnt, time_cur, temp_cur, alpha_cur, krate, Cp_now, mass_cem_now, kpozz / krate, kslag / krate, alpha_fa_cur);
        fprintf( heatfile, "%d %f %f %f %f %f  %f \n",
                cyccnt - 1, time_cur, alpha, alpha_cur, heat_new * heat_cf, gsratio2, ( ( 0.68 * alpha_cur ) / ( 0.32 * alpha_cur + w_to_c ) ) );
        fflush(adiafile);
        fflush(heatfile);
#endif

        /*chemical shrinkage*/
        //fprintf(chsfile,"%d %f %f %f\n", cyccnt-1,time_cur,alpha_cur,chs_new);
        //fflush(chsfile);
        //pHpred();-dangerous function - can lead to zero division etc.

        /* Check percolation of pore space */
        /* Note that first variable passed corresponds to phase to check */
        /* Could easily add calls to check for percolation of CH, CSH, etc. */

        if ( ( ( icyc % burnfreq ) == 0 ) && ( ( porefl1 + porefl2 + porefl3 ) != 0 ) ) {
            porefl1 = burn3d(0, 1, 0, 0);
            porefl2 = burn3d(0, 0, 1, 0);
            porefl3 = burn3d(0, 0, 0, 1);
            /* Switch to self-desiccating conditions when porosity */
            /* disconnects */
            if ( ( ( porefl1 + porefl2 + porefl3 ) == 0 ) && ( sealed == 0 ) ) {
                water_off = water_left;
                pore_off = countkeep;
                sealed = 1;
#ifdef PRINTF
                printf("Switching to self-desiccating at cycle %d \n", cyccnt);
#endif
                fflush(stdout);
            }
        }


        /*update microstructure in mic_CSH[][][] to HDCSH and use it in*/
        /*burn_phases function and outputImageFile function*/
        //CSHbox(CSH_vicinity); - called independently, does not influence percolation


        /* Check percolation of solids (set point) */
        if ( ( ( icyc % setfreq ) == 0 ) && ( setflag == 0 ) ) {
            //             sf1 = burnset(1, 0, 0);
            //             sf2 = burnset(0, 1, 0);
            //             sf3 = burnset(0, 0, 1);
            //             setflag = sf1 * sf2 * sf3;
        }

#ifdef IMAGEFILES
        //outputImageFilePerc();//is executed from CemhydMatStatus :: burn_phases
        outputImageFileUnperc(mic);
#endif
        /* Check hydration of particles */
        if ( ( icyc % phydfreq ) == 0 ) {
            //  parthyd();
        }

        /* Calculate homogenized elastic values if turned on */
        if ( alpha_cur >= TargDoHelas && Calculate_elastic_homogenization ) {
            //increase DoH by 0.02 till 0.2 and then by 0.05 step
            if ( TargDoHelas < 0.1 ) {
                TargDoHelas += 0.01; //0.02
            } else if ( TargDoHelas < 0.2 ) {
                TargDoHelas += 0.01; //0.05
            } else {
                TargDoHelas += 0.01; //0.05
            }

            double E_paste_SCM_filler_air, nu_paste_SCM_filler_air, E_mortar, nu_mortar;

            CreateHDCSH(); //alters microstructure so it contains HDCSH - need for all homogenizations
            PercolateForOutput();
            AnalyticHomogenizationPaste(last_values [ 2 ], last_values [ 3 ], 0); //0-perc, 1-unperc
            AnalyticHomogenizationConcrete(last_values [ 2 ], last_values [ 3 ], & E_paste_SCM_filler_air, & nu_paste_SCM_filler_air, & E_mortar, & nu_mortar, last_values [ 4 ], last_values [ 5 ]);
#ifdef OUTFILES
            fprintf(elasfile, "%d \t%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", cyccnt, time_cur, alpha_cur, last_values [ 2 ], last_values [ 3 ], E_paste_SCM_filler_air, nu_paste_SCM_filler_air, E_mortar, nu_mortar, last_values [ 4 ], last_values [ 5 ]);
            fflush(elasfile);
#endif
        }

#ifdef OUTFILES
        fprintf(CSHfile, "%d %.4f %.4f %ld %ld %f\n", icyc, alpha_cur, time_cur, count [ HDCSH ], count [ CSH ], ( double ) count [ HDCSH ] / count [ CSH ]);
        fflush(CSHfile);
#endif

        icyc++;
        counter++;
    }

    //end of main while loop
}


/*Function returning the concrete released heat during the difference in given times
 * of hydration. Units are in kJ/kg_all_cement_solids.
 * All solids include cement clinker minerals, gypsum forms, [INERT],
 * [SLAG],[POZZ],[CACL2],[ASG],[CAS2], NOT aggregates!!!
 * The heat is independent on the hydrating cube size but may slightly
 * decrease as some reversible reactions exist!!! (checked)
 * ^
 **|Released heat                                            (stops now)
 |                                                            CYC_m
 |             (last stop)                      t_n+1           |
 |                CYC   ..CYC..  CYC_m-1          |             |
 |  t_n            !              |               |             |
 |  !LastTotHeat   !LastHydrTime  |PrevHydrTime   |TargTime     |time_cur
 |  !LastCallTime  !LastCycHeat   |PrevCycHeat    |LastTargTime |heat_new
 |  !              !LastCycCnt    |               |             |
 |||------------------------------------------------------------------->time
 *
 * let all time counters in seconds! (numerical precision)
 * TargTime, LastCallTime, LastTargTime are in absolute time taken from OOFEM
 * time_cur, PrevHydrTime, LastHydrTime is hydration time, ofset with the castingTime
 */

///GiveTemp in [C], TargTime in [s]
double CemhydMatStatus :: GivePower(double GiveTemp, double TargTime) {
#ifdef __TM_MODULE //OOFEM transport module
    double castingTime = this->gp->giveMaterial()->giveCastingTime();
#else
    double castingTime = 0.;
#endif

    //do not calculate anything before casting time
    if ( TargTime - castingTime <= 0 ) {
        PartHeat = 0.;
        LastCallTime = castingTime;
        return 0.;
    }

    LastTargTime = TargTime;

    if ( TargTime < LastCallTime ) {
        printf("Cannot go backwards in hydration, TargTime = %f s < LastCallTime = %f s\n", TargTime, LastCallTime);
        exit(0);
    }

#ifdef __TM_MODULE //OOFEM transport module
    CemhydMat *cemhydmat = ( CemhydMat * ) this->gp->giveMaterial();
    if ( !cemhydmat->nowarnings.at(4) && ( GiveTemp > 200. ) ) {
        printf("Temperature exceeds 200 C (file %s, line %d),\n", __FILE__, __LINE__);
    }

#else
    if ( GiveTemp > 200. ) {
        printf("Temperature exceeds 200 C (file %s, line %d),\n", __FILE__, __LINE__);
    }

#endif


    /*perform as many cycles as necessary, always more than required
     * first, calculate everything on the cement paste*/

    disrealnew(GiveTemp, ( TargTime - castingTime ) / 3600., 0); //perform hydration cycles controled with TargTime [h]

#ifdef PRINTF
    printf("Hydration heat %f [J/g of cement]\n", heat_new * heat_cf);
    printf("time_cur %f [h]\n", time_cur);
    printf("PrevHydrTime %f [h]\n", PrevHydrTime / 3600.);
#endif

    PartHeat = 0.;

    if ( TargTime != LastCallTime && icyc >  1 ) {
        if ( LastCycCnt - icyc ) { //if at least one hydration cycle has elapsed
            //interpolate between CYC_m-1 and CYC_m
            PartHeat = ( heat_new * heat_cf - PrevCycHeat ) * ( TargTime - PrevHydrTime - castingTime ) / ( 3600. * time_cur - PrevHydrTime );
            //add between CYC and CYC_m-1
            PartHeat += PrevCycHeat - LastCycHeat;
            //add between t_n and CYC
            PartHeat += LastCycHeat - LastTotHeat;
        } else { //not even one cycle elapsed
            PartHeat = ( LastCycHeat - LastTotHeat ) * ( TargTime - LastCallTime ) / ( LastHydrTime - LastCallTime + castingTime );
        }

        LastTotHeat += PartHeat;

        if ( TargTime != LastCallTime ) {
            //Mass_cement_concrete tells the amount of cement [kg] in a m3 of concrete
            //add aggregates and rescale to kJ/m3 of concrete
            PartHeat *= Mass_cement_concrete;
            //rescale from kJ/m3/time_interval to J/m3/s or W/m3
            PartHeat /= ( TargTime - LastCallTime );
            PartHeat *= 1000;
        }
    }

    LastCallTime = TargTime;
    LastHydrTime = time_cur * 3600.;
    LastCycHeat = ( double ) ( heat_new ) * heat_cf;
    LastCycCnt = icyc;

#ifdef PRINTF
    printf("PartHeat %f [W/m3 of concrete]\n", PartHeat);
#endif

    //PartHeat = 1500;
    return PartHeat;
}


//move hydration model by several steps and update all internal variables
double CemhydMatStatus :: MoveCycles(double GiveTemp, int cycles) {
    double PartHeat;
    double castingTime = 0.;
#ifdef __TM_MODULE
    this->gp->giveMaterial()->giveCastingTime();
#endif

    disrealnew(GiveTemp, -1., cycles); //perform amount of hydration cycles controlled

    //at least one cycle elapsed
    //difference in heat before and now at CYC_m
    PartHeat = heat_new * heat_cf - PrevCycHeat;
    //add between CYC and CYC_m-1
    PartHeat += PrevCycHeat - LastCycHeat;
    //add between t_n and CYC
    PartHeat += LastCycHeat - LastTotHeat;

    LastTotHeat += PartHeat;

    //Mass_cement_concrete tells the amount of cement [kg] in a m3 of concrete
    //add aggregates and rescale to kJ/m3 of concrete
    PartHeat *= Mass_cement_concrete;
    //rescale from kJ/m3/time_interval to J/m3/s
    PartHeat /= ( 3600. * time_cur - LastCallTime + castingTime );
    PartHeat *= 1000;

    LastCallTime = 3600. * time_cur + castingTime;
    LastHydrTime = 3600. * time_cur;
    LastCycHeat = ( double ) ( heat_new ) * heat_cf;
    LastCycCnt = icyc;

    LastTargTime = time_cur * 3600.;
    return time_cur; //return in hours
}

//move to desired DoH (the next DoH in CEMHYD3D cycle) limited by maxcyc
//return 1 if the number of cycles was exceeded
int CemhydMatStatus :: MoveToDoH(double GiveTemp, double DesiredDoH, int maxcyc) {
    int cycle = 1;

    while ( alpha_cur < DesiredDoH ) {
        if ( cycle > maxcyc ) {
#ifdef PRINTF
            printf("Target DoH %f (now is %f) was not reached after %d cycles (%s, line %d)\n", DesiredDoH, alpha_cur, maxcyc, __FILE__, __LINE__);
#endif
            return 1;
        }

        MoveCycles(GiveTemp, 1); //move one cycle
        cycle++;
    }

    return 0;
}

//TargTime in [h]
//move to desired Time[h] (the next time in CEMHYD3D cycle, no interpolation)
//no cycle limit is specified (any time can be reached)
int CemhydMatStatus :: MoveToTime(double GiveTemp, double TargTime) {
    GivePower(GiveTemp, 3600. * TargTime);
    return 0;
}


//return total released heat from cement in kJ/kg_cement_solids at TargTime
double CemhydMatStatus :: GiveTotCemHeat(void)
{
    return LastTotHeat;
}

//return total heat in J/m3_of_concrete at TargTime
double CemhydMatStatus :: GiveTotHeat(void)
{
    return 1000 * LastTotHeat * Mass_cement_concrete;
}

//return heat capacity in kJ/K/kg_of_everything in the last CEMHYD cycle
double CemhydMatStatus :: GiveCp(void)
{
    return Cp_now;
}


//see D. P. Bentz: Transient Plane Source Measurements of the Thermal Properties of Hydrating Cement Pastes, Materials and Structures, 1073-1080, 2007
double CemhydMatStatus :: computeConcreteCapacityBentz(void)
{
    double capacityPaste, capacityConcrete;
    capacityPaste = ( 4180. * Vol_water * 1.0 + Cp_cement * 1000. * Mass_cement_concrete ) / Mass_tot_concrete; //capacity of cement paste [J/kg/K]
    capacityPaste *= ( 1. - 0.26 * ( 1 - pow(2.71828, -2.9 * alpha_cur) ) );
    capacityConcrete = capacityPaste + Cp_SCM * Mass_SCM_frac + Cp_FA * Mass_FA_frac + Cp_CA * Mass_CA_frac + Cp_inert * Mass_inert_frac;

    return capacityConcrete;
}

double CemhydMatStatus :: GiveDensity(void)
{
    return Concrete_bulk_density;
}

//return degree of hydration in the last CEMHYD cycle, corresponding to time_cur
double CemhydMatStatus :: GiveDoHLastCyc(void)
{
    return alpha_cur;
}

//return actual degree of hydration (interpolated)
double CemhydMatStatus :: GiveDoHActual(void)
{
    double castingTime = 0;

    if ( icyc <= 1 ) { //not even one cycle has elapsed
        return 0;
    }

#ifdef __TM_MODULE
    castingTime = this->gp->giveMaterial()->giveCastingTime();
#endif

    return alpha_last + ( alpha_cur - alpha_last ) * ( LastTargTime - castingTime - PrevHydrTime ) / ( 3600. * time_cur - PrevHydrTime ); //interpolate
}


//return cycle in the last CEMHYD cycle (higher or the same as actual time)
int CemhydMatStatus :: GiveCycNum(void)
{
    return icyc;
}

//return actual time of the last cycle (higher or the same as actual time)
double CemhydMatStatus :: GiveCycTime(void)
{
    return time_cur;
}



/* routine to assess the connectivity (percolation) of a single phase */
/* Two matrices are used here: one to store the recently burnt locations */
/*         burn3d.cpp:97: implicit declaration of function `int printf(...)'
 *         the other to store the newly found burnt locations */

int CemhydMatStatus :: burn3d(int npix, int d1, int d2, int d3)
/* npix is ID of phase to perform burning on */
/* directional flags */
{
    long int ntop, nthrough, ncur, nnew, ntot, nphc;
    int i, inew, j, k;
    int xl, xh, j1, k1, px, py, pz, qx, qy, qz, xcn, ycn, zcn;
    int x1, y1, z1, igood;
    int *nmatx, *nmaty, *nmatz, *nnewx, *nnewy, *nnewz;
    int jnew, icur;
    int bflag;
    float mass_burn = 0.0, alpha_burn = 0.0;

    nmatx = new int [ SIZE2D ];
    nmaty = new int [ SIZE2D ];
    nmatz = new int [ SIZE2D ];
    nnewx = new int [ SIZE2D ];
    nnewy = new int [ SIZE2D ];
    nnewz = new int [ SIZE2D ];


    /* counters for number of pixels of phase accessible from surface #1 */
    /* and number which are part of a percolated pathway to surface #2 */
    ntop = 0;
    bflag = 0;
    nthrough = 0;
    nphc = 0;

    /* percolation is assessed from top to bottom only */
    /* and burning algorithm is periodic in other two directions */
    /* use of directional flags allow transformation of coordinates */
    /* to burn in direction of choosing (x, y, or z) */
    i = 0;

    for ( k = 0; k < SYSIZE; k++ ) {
        for ( j = 0; j < SYSIZE; j++ ) {
            igood = 0;
            ncur = 0;
            ntot = 0;
            /* Transform coordinates */
            px = cx(i, j, k, d1, d2, d3);
            py = cy(i, j, k, d1, d2, d3);
            pz = cz(i, j, k, d1, d2, d3);
            if ( mic [ px ] [ py ] [ pz ] == npix ) {
                /* Start a burn front */
                mic [ px ] [ py ] [ pz ] = BURNT;
                ntot += 1;
                ncur += 1;
                /* burn front is stored in matrices nmat* */
                /* and nnew* */
                nmatx [ ncur ] = i;
                nmaty [ ncur ] = j;
                nmatz [ ncur ] = k;
                /* Burn as long as new (fuel) pixels are found */
                do {
                    nnew = 0;
                    for ( inew = 1; inew <= ncur; inew++ ) {
                        xcn = nmatx [ inew ];
                        ycn = nmaty [ inew ];
                        zcn = nmatz [ inew ];

                        /* Check all six neighbors */
                        for ( jnew = 1; jnew <= 6; jnew++ ) {
                            x1 = xcn;
                            y1 = ycn;
                            z1 = zcn;
                            if ( jnew == 1 ) {
                                x1 -= 1;
                            }

                            if ( jnew == 2 ) {
                                x1 += 1;
                            }

                            if ( jnew == 3 ) {
                                y1 -= 1;
                            }

                            if ( jnew == 4 ) {
                                y1 += 1;
                            }

                            if ( jnew == 5 ) {
                                z1 -= 1;
                            }

                            if ( jnew == 6 ) {
                                z1 += 1;
                            }

                            /* Periodic in y and */
                            if ( y1 >= SYSIZE ) {
                                y1 -= SYSIZE;
                            } else if ( y1 < 0 ) {
                                y1 += SYSIZE;
                            }

                            /* Periodic in z direction */
                            if ( z1 >= SYSIZE ) {
                                z1 -= SYSIZE;
                            } else if ( z1 < 0 ) {
                                z1 += SYSIZE;
                            }

                            /* Nonperiodic so be sure to remain in the 3-D box */
                            if ( ( x1 >= 0 ) && ( x1 < SYSIZE ) ) {
                                /* Transform coordinates */
                                px = cx(x1, y1, z1, d1, d2, d3);
                                py = cy(x1, y1, z1, d1, d2, d3);
                                pz = cz(x1, y1, z1, d1, d2, d3);
                                if ( mic [ px ] [ py ] [ pz ] == npix ) {
                                    ntot += 1;
                                    mic [ px ] [ py ] [ pz ] = BURNT;
                                    nnew += 1;
                                    if ( nnew >= SIZE2D ) {
                                        printf("error in size of nnew \n");
                                    }

                                    nnewx [ nnew ] = x1;
                                    nnewy [ nnew ] = y1;
                                    nnewz [ nnew ] = z1;
                                }
                            }
                        }
                    }

                    if ( nnew > 0 ) {
                        ncur = nnew;
                        /* update the burn front matrices */
                        for ( icur = 1; icur <= ncur; icur++ ) {
                            nmatx [ icur ] = nnewx [ icur ];
                            nmaty [ icur ] = nnewy [ icur ];
                            nmatz [ icur ] = nnewz [ icur ];
                        }
                    }
                } while ( nnew > 0 );

                ntop += ntot;
                xl = 0;
                xh = SYSIZE - 1;
                /* See if current path extends through the microstructure */
                for ( j1 = 0; j1 < SYSIZE; j1++ ) {
                    for ( k1 = 0; k1 < SYSIZE; k1++ ) {
                        px = cx(xl, j1, k1, d1, d2, d3);
                        py = cy(xl, j1, k1, d1, d2, d3);
                        pz = cz(xl, j1, k1, d1, d2, d3);
                        qx = cx(xh, j1, k1, d1, d2, d3);
                        qy = cy(xh, j1, k1, d1, d2, d3);
                        qz = cz(xh, j1, k1, d1, d2, d3);
                        if ( ( mic [ px ] [ py ] [ pz ] == BURNT ) && ( mic [ qx ] [ qy ] [ qz ] == BURNT ) ) {
                            igood = 2;
                        }

                        if ( mic [ px ] [ py ] [ pz ] == BURNT ) {
                            mic [ px ] [ py ] [ pz ] = BURNT + 1;
                        }

                        if ( mic [ qx ] [ qy ] [ qz ] == BURNT ) {
                            mic [ qx ] [ qy ] [ qz ] = BURNT + 1;
                        }
                    }
                }

                if ( igood == 2 ) {
                    nthrough += ntot;
                }
            }
        }
    }

    /* return the burnt sites to their original phase values */
    for ( i = 0; i < SYSIZE; i++ ) {
        for ( j = 0; j < SYSIZE; j++ ) {
            for ( k = 0; k < SYSIZE; k++ ) {
                if ( mic [ i ] [ j ] [ k ] >= BURNT ) {
                    nphc += 1;
                    mic [ i ] [ j ] [ k ] = npix;
                } else if ( mic [ i ] [ j ] [ k ] == npix ) {
                    nphc += 1;
                }
            }
        }
    }

#ifdef PRINTF
    printf("Phase ID= %d \n", npix);
    printf("Number accessible from first surface = %ld \n", ntop);
    printf("Number contained in through pathways= %ld \n", nthrough);
#endif


    mass_burn += specgrav [ C3S ] * count [ C3S ];
    mass_burn += specgrav [ C2S ] * count [ C2S ];
    mass_burn += specgrav [ C3A ] * count [ C3A ];
    mass_burn += specgrav [ C4AF ] * count [ C4AF ];
    alpha_burn = 1. - ( mass_burn / cemmass );

#ifdef OUTFILES
    fprintf(fileperc, "%d %f %ld %ld \n", cyccnt, alpha_cur, nthrough, nphc);
    fflush(fileperc);
#endif

    if ( nthrough > 0 ) {
        bflag = 1;
    }

    delete [] nmatx;
    delete [] nmaty;
    delete [] nmatz;
    delete [] nnewx;
    delete [] nnewy;
    delete [] nnewz;
    return ( bflag );
}


/* routine to assess connectivity (percolation) of solids for set estimation*/
/* Definition of set is a through pathway of cement particles connected */
/* together by CSH or ettringite */
/* Two matrices are used here: one to store the recently burnt locations */
/* the other to store the newly found burnt locations */
//use perc_phases instead. It contains all solid phases, :: burnset function is outdated

int CemhydMatStatus :: burnset(int d1, int d2, int d3)
{
    long int ntop, nthrough, icur, inew, ncur, nnew, ntot, count_solid;
    int i, j, k, setyet;
    //static int nmatx[SIZESET],nmaty[SIZESET],nmatz[SIZESET];
    int *nmatx, *nmaty, *nmatz;
    int xl, xh, j1, k1, px, py, pz, qx, qy, qz;
    int xcn, ycn, zcn, x1, y1, z1, igood;
    //static int nnewx[SIZESET],nnewy[SIZESET],nnewz[SIZESET];
    int *nnewx, *nnewy, *nnewz;
    int jnew;
    float mass_burn = 0.0, alpha_burn = 0.0, con_frac;
    //static char newmat [SYSIZE] [SYSIZE] [SYSIZE];
    char ***newmat;

    alloc_char_3D(newmat, SYSIZE);
    nmatx = new int [ SIZESET ];
    nmaty = new int [ SIZESET ];
    nmatz = new int [ SIZESET ];
    nnewx = new int [ SIZESET ];
    nnewy = new int [ SIZESET ];
    nnewz = new int [ SIZESET ];

    /* counters for number of pixels of phase accessible from surface #1 */
    /* and number which are part of a percolated pathway to surface #2 */
    ntop = 0;
    nthrough = 0;
    setyet = 0;
    for ( k = 0; k < SYSIZE; k++ ) {
        for ( j = 0; j < SYSIZE; j++ ) {
            for ( i = 0; i < SYSIZE; i++ ) {
                newmat [ i ] [ j ] [ k ] = mic [ i ] [ j ] [ k ];
            }
        }
    }

    /* percolation is assessed from top to bottom only */
    /* in transformed coordinates */
    /* and burning algorithm is periodic in other two directions */
    i = 0;

    for ( k = 0; k < SYSIZE; k++ ) {
        for ( j = 0; j < SYSIZE; j++ ) {
            igood = 0;
            ncur = 0;
            ntot = 0;
            /* Transform coordinates */
            px = cx(i, j, k, d1, d2, d3);
            py = cy(i, j, k, d1, d2, d3);
            pz = cz(i, j, k, d1, d2, d3);
            /* start from a cement clinker, slag, fly ash ettringite, C3AH6, or
             * CSH pixel */
            if ( ( mic [ px ] [ py ] [ pz ] == C3S ) ||
                ( mic [ px ] [ py ] [ pz ] == C2S ) ||
                ( mic [ px ] [ py ] [ pz ] == SLAG ) ||
                ( mic [ px ] [ py ] [ pz ] == ASG ) ||
                ( mic [ px ] [ py ] [ pz ] == CAS2 ) ||
                ( mic [ px ] [ py ] [ pz ] == POZZ ) ||
                ( mic [ px ] [ py ] [ pz ] == CSH ) ||
                ( mic [ px ] [ py ] [ pz ] == C3AH6 ) ||
                ( mic [ px ] [ py ] [ pz ] == ETTR ) ||
                ( mic [ px ] [ py ] [ pz ] == ETTRC4AF ) ||
                ( mic [ px ] [ py ] [ pz ] == C3A ) ||
                ( mic [ px ] [ py ] [ pz ] == C4AF ) ) {
                /* Start a burn front */
                mic [ px ] [ py ] [ pz ] = BURNT;
                ntot += 1;
                ncur += 1;
                /* burn front is stored in matrices nmat* */
                /* and nnew* */
                nmatx [ ncur ] = i;
                nmaty [ ncur ] = j;
                nmatz [ ncur ] = k;
                /* Burn as long as new (fuel) pixels are found */
                do {
                    nnew = 0;
                    for ( inew = 1; inew <= ncur; inew++ ) {
                        xcn = nmatx [ inew ];
                        ycn = nmaty [ inew ];
                        zcn = nmatz [ inew ];
                        /* Convert to directional coordinates */
                        qx = cx(xcn, ycn, zcn, d1, d2, d3);
                        qy = cy(xcn, ycn, zcn, d1, d2, d3);
                        qz = cz(xcn, ycn, zcn, d1, d2, d3);

                        /* Check all six neighbors */
                        for ( jnew = 1; jnew <= 6; jnew++ ) {
                            x1 = xcn;
                            y1 = ycn;
                            z1 = zcn;
                            if ( jnew == 1 ) {
                                x1 -= 1;
                            }

                            if ( jnew == 2 ) {
                                x1 += 1;
                            }

                            if ( jnew == 3 ) {
                                y1 -= 1;
                            }

                            if ( jnew == 4 ) {
                                y1 += 1;
                            }

                            if ( jnew == 5 ) {
                                z1 -= 1;
                            }

                            if ( jnew == 6 ) {
                                z1 += 1;
                            }

                            /* Periodic in y and */
                            if ( y1 >= SYSIZE ) {
                                y1 -= SYSIZE;
                            } else if ( y1 < 0 ) {
                                y1 += SYSIZE;
                            }

                            /* Periodic in z direction */
                            if ( z1 >= SYSIZE ) {
                                z1 -= SYSIZE;
                            } else if ( z1 < 0 ) {
                                z1 += SYSIZE;
                            }

                            /* Nonperiodic so be sure to remain in the 3-D box */
                            if ( ( x1 >= 0 ) && ( x1 < SYSIZE ) ) {
                                px = cx(x1, y1, z1, d1, d2, d3);
                                py = cy(x1, y1, z1, d1, d2, d3);
                                pz = cz(x1, y1, z1, d1, d2, d3);
                                /* Conditions for propagation of burning */
                                /* 1) new pixel is CH, CSH or ETTR or C3AH6*/
                                if ( ( mic [ px ] [ py ] [ pz ] == CSH ) || ( mic [ px ] [ py ] [ pz ] == ETTRC4AF ) || ( mic [ px ] [ py ] [ pz ] == C3AH6 ) || ( mic [ px ] [ py ] [ pz ] == ETTR ) ) {
                                    ntot += 1;
                                    mic [ px ] [ py ] [ pz ] = BURNT;
                                    nnew += 1;
                                    if ( nnew >= SIZESET ) {
                                        printf("error in size of nnew %ld\n", nnew);
                                    }

                                    nnewx [ nnew ] = x1;
                                    nnewy [ nnew ] = y1;
                                    nnewz [ nnew ] = z1;
                                }
                                /* 2) old pixel is CSH or ETTR or C3AH6 and new pixel is one of cement clinker, slag, or fly ash phases */
                                else if ( ( ( newmat [ qx ] [ qy ] [ qz ] == CSH ) || ( newmat [ qx ] [ qy ] [ qz ] == ETTRC4AF ) || ( newmat [ qx ] [ qy ] [ qz ] == C3AH6 ) || ( newmat [ qx ] [ qy ] [ qz ] == ETTR ) ) &&
                                         ( ( mic [ px ] [ py ] [ pz ] == C3S ) ||
                                          ( mic [ px ] [ py ] [ pz ] == C2S ) ||
                                          ( mic [ px ] [ py ] [ pz ] == CAS2 ) ||
                                          ( mic [ px ] [ py ] [ pz ] == SLAG ) ||
                                          ( mic [ px ] [ py ] [ pz ] == POZZ ) ||
                                          ( mic [ px ] [ py ] [ pz ] == ASG ) ||
                                          ( mic [ px ] [ py ] [ pz ] == C3A ) ||
                                          ( mic [ px ] [ py ] [ pz ] == C4AF ) ) ) {
                                    ntot += 1;
                                    mic [ px ] [ py ] [ pz ] = BURNT;
                                    nnew += 1;
                                    if ( nnew >= SIZESET ) {
                                        printf("error in size of nnew %ld\n", nnew);
                                    }

                                    nnewx [ nnew ] = x1;
                                    nnewy [ nnew ] = y1;
                                    nnewz [ nnew ] = z1;
                                }
                                /* 3) old and new pixels belong to one of cement clinker, slag, or fly ash phases and */
                                /* are contained in the same initial cement particle */
                                /* and it is not a one-pixel particle */
                                else if ( ( micpart [ qx ] [ qy ] [ qz ] == micpart [ px ] [ py ] [ pz ] ) &&
                                         ( ( mic [ px ] [ py ] [ pz ] == C3S ) ||
                                          ( mic [ px ] [ py ] [ pz ] == C2S ) ||
                                          ( mic [ px ] [ py ] [ pz ] == POZZ ) ||
                                          ( mic [ px ] [ py ] [ pz ] == SLAG ) ||
                                          ( mic [ px ] [ py ] [ pz ] == ASG ) ||
                                          ( mic [ px ] [ py ] [ pz ] == CAS2 ) ||
                                          ( mic [ px ] [ py ] [ pz ] == C3A ) ||
                                          ( mic [ px ] [ py ] [ pz ] == C4AF ) ) && ( ( newmat [ qx ] [ qy ] [ qz ] == C3S ) ||
                                                                                     ( newmat [ qx ] [ qy ] [ qz ] == C2S ) ||
                                                                                     ( newmat [ qx ] [ qy ] [ qz ] == SLAG ) ||
                                                                                     ( newmat [ qx ] [ qy ] [ qz ] == ASG ) ||
                                                                                     ( newmat [ qx ] [ qy ] [ qz ] == POZZ ) ||
                                                                                     ( newmat [ qx ] [ qy ] [ qz ] == CAS2 ) ||
                                                                                     ( newmat [ qx ] [ qy ] [ qz ] == C3A ) ||
                                                                                     ( newmat [ qx ] [ qy ] [ qz ] == C4AF ) ) ) {
                                    ntot += 1;
                                    mic [ px ] [ py ] [ pz ] = BURNT;
                                    nnew += 1;
                                    if ( nnew >= SIZESET ) {
                                        printf("error in size of nnew %ld\n", nnew);
                                    }

                                    nnewx [ nnew ] = x1;
                                    nnewy [ nnew ] = y1;
                                    nnewz [ nnew ] = z1;
                                }
                            }

                            /* nonperiodic if delimiter */
                        }

                        /* neighbors loop */
                    }

                    /* propagators loop */
                    if ( nnew > 0 ) {
                        ncur = nnew;
                        /* update the burn front matrices */
                        for ( icur = 1; icur <= ncur; icur++ ) {
                            nmatx [ icur ] = nnewx [ icur ];
                            nmaty [ icur ] = nnewy [ icur ];
                            nmatz [ icur ] = nnewz [ icur ];
                        }
                    }
                } while ( nnew > 0 );

                ntop += ntot;
                xl = 0;
                xh = SYSIZE - 1;
                /* Check for percolated path through system */
                for ( j1 = 0; j1 < SYSIZE; j1++ ) {
                    for ( k1 = 0; k1 < SYSIZE; k1++ ) {
                        px = cx(xl, j1, k1, d1, d2, d3);
                        py = cy(xl, j1, k1, d1, d2, d3);
                        pz = cz(xl, j1, k1, d1, d2, d3);
                        qx = cx(xh, j1, k1, d1, d2, d3);
                        qy = cy(xh, j1, k1, d1, d2, d3);
                        qz = cz(xh, j1, k1, d1, d2, d3);
                        if ( ( mic [ px ] [ py ] [ pz ] == BURNT ) && ( mic [ qx ] [ qy ] [ qz ] == BURNT ) ) {
                            igood = 2;
                        }

                        if ( mic [ px ] [ py ] [ pz ] == BURNT ) {
                            mic [ px ] [ py ] [ pz ] = BURNT + 1;
                        }

                        if ( mic [ qx ] [ qy ] [ qz ] == BURNT ) {
                            mic [ qx ] [ qy ] [ qz ] = BURNT + 1;
                        }
                    }
                }

                if ( igood == 2 ) {
                    nthrough += ntot;
                }
            }
        }
    }

#ifdef PRINTF
    printf("Phase ID= Solid Phases \n");
    printf("Number accessible from first surface = %ld \n", ntop);
    printf("Number contained in through pathways= %ld \n", nthrough);
#endif

    mass_burn += specgrav [ C3S ] * count [ C3S ];
    mass_burn += specgrav [ C2S ] * count [ C2S ];
    mass_burn += specgrav [ C3A ] * count [ C3A ];
    mass_burn += specgrav [ C4AF ] * count [ C4AF ];
    alpha_burn = 1. - ( mass_burn / cemmass );

    con_frac = 0.0;
    count_solid = count [ C3S ] + count [ C2S ] + count [ C3A ] + count [ C4AF ] + count [ ETTR ] + count [ CSH ] + count [ C3AH6 ] + count [ ETTRC4AF ] + count [ POZZ ] + count [ ASG ] + count [ SLAG ] + count [ CAS2 ];

    if ( count_solid > 0 ) {
        con_frac = ( float ) nthrough / ( float ) count_solid;
    }

#ifdef OUTFILES
    //fprintf(percfile, "%d  %f %f  %ld %ld %f\n", cyccnt, time_cur + ( 2. * ( float ) ( cyccnt ) - 1.0 ) * beta / krate, alpha_cur, nthrough, count [ C3S ] + count [ C2S ] + count [ C3A ] + count [ C4AF ] + count [ CAS2 ] + count [ SLAG ] + count [ ASG ] + count [ POZZ ] + count [ ETTR ] + count [ C3AH6 ] + count [ ETTRC4AF ] + count [ CSH ], con_frac);
    //fflush(percfile);
#endif

    if ( con_frac > 0.985 ) {
        setyet = 1;
    }

    outputImageFileUnperc(mic);

    /* return the burnt sites to their original phase values */
    for ( i = 0; i < SYSIZE; i++ ) {
        for ( j = 0; j < SYSIZE; j++ ) {
            for ( k = 0; k < SYSIZE; k++ ) {
                if ( mic [ i ] [ j ] [ k ] >= BURNT ) {
                    mic [ i ] [ j ] [ k ] = newmat [ i ] [ j ] [ k ];
                }
            }
        }
    }


    dealloc_char_3D(newmat, SYSIZE);
    delete []  nmatx;
    delete []  nmaty;
    delete []  nmatz;
    delete []  nnewx;
    delete []  nnewy;
    delete []  nnewz;

    /* Return flag indicating if set has indeed occurred */
    return ( setyet );
}



/* Routine to assess relative particle hydration */
void CemhydMatStatus :: parthyd(void) {
    int norig [ 100000 ], nleft [ 100000 ];
    int ix, iy, iz;
    char valmic, valmicorig;
    int valpart, partmax;
    float alpart;
    FILE *phydfile;

    /* Initialize the particle count arrays */
    for ( ix = 0; ix < 100000; ix++ ) {
        nleft [ ix ] = norig [ ix ] = 0;
    }

    phydfile = fopen(phrname, "a");
    fprintf(phydfile, "%d %f\n", cyccnt, alpha_cur);

    partmax = 0;
    /* Scan the microstructure pixel by pixel and update counts */
    for ( ix = 0; ix < SYSIZE; ix++ ) {
        for ( iy = 0; iy < SYSIZE; iy++ ) {
            for ( iz = 0; iz < SYSIZE; iz++ ) {
                if ( micpart [ ix ] [ iy ] [ iz ] != 0 ) {
                    valpart = micpart [ ix ] [ iy ] [ iz ];
                    if ( valpart > partmax ) {
                        partmax = valpart;
                    }

                    valmic = mic [ ix ] [ iy ] [ iz ];
                    if ( ( valmic == C3S ) || ( valmic == C2S ) || ( valmic == C3A ) || ( valmic == C4AF ) ) {
                        nleft [ valpart ] += 1;
                    }

                    valmicorig = micorig [ ix ] [ iy ] [ iz ];
                    if ( ( valmicorig == C3S ) || ( valmicorig == C2S ) || ( valmicorig == C3A ) || ( valmicorig == C4AF ) ) {
                        norig [ valpart ] += 1;
                    }
                }
            }
        }
    }

    /* Output results to end of particle hydration file */
    for ( ix = 100; ix <= partmax; ix++ ) {
        alpart = 0.0;
        if ( norig [ ix ] != 0 ) {
            alpart = 1. - ( float ) nleft [ ix ] / ( float ) norig [ ix ];
        }

        fprintf(phydfile, "%d %d %d %.3f\n", ix, norig [ ix ], nleft [ ix ], alpart);
    }

    fclose(phydfile);
}


/* extpozz, movefh3, movech, extc3ah6, movec3a */
/* extfreidel, movecacl2, extstrat, moveas */
int CemhydMatStatus :: moveone(int *xloc, int *yloc, int *zloc, int *act, int sumold)
{
    int plok, sumnew, xl1, yl1, zl1, act1;

    sumnew = 1;
    /* store the input values for location */
    xl1 = ( * xloc );
    yl1 = ( * yloc );
    zl1 = ( * zloc );
    act1 = ( * act );

    /* Choose one of six directions (at random) for the new */
    /* location */
    plok = ( int ) ( 6. * ran1(seed) );
    if ( ( plok > 5 ) || ( plok < 0 ) ) {
        plok = 5;
    }

    switch ( plok ) {
    case 0:
        xl1 -= 1;
        act1 = 1;
        if ( xl1 < 0 ) {
            xl1 = ( SYSIZEM1 );
        }

        if ( sumold % 2 != 0 ) {
            sumnew = 2;
        }

        break;
    case 1:
        xl1 += 1;
        act1 = 2;
        if ( xl1 >= SYSIZE ) {
            xl1 = 0;
        }

        if ( sumold % 3 != 0 ) {
            sumnew = 3;
        }

        break;
    case 2:
        yl1 -= 1;
        act1 = 3;
        if ( yl1 < 0 ) {
            yl1 = ( SYSIZEM1 );
        }

        if ( sumold % 5 != 0 ) {
            sumnew = 5;
        }

        break;
    case 3:
        yl1 += 1;
        act1 = 4;
        if ( yl1 >= SYSIZE ) {
            yl1 = 0;
        }

        if ( sumold % 7 != 0 ) {
            sumnew = 7;
        }

        break;
    case 4:
        zl1 -= 1;
        act1 = 5;
        if ( zl1 < 0 ) {
            zl1 = ( SYSIZEM1 );
        }

        if ( sumold % 11 != 0 ) {
            sumnew = 11;
        }

        break;
    case 5:
        zl1 += 1;
        act1 = 6;
        if ( zl1 >= SYSIZE ) {
            zl1 = 0;
        }

        if ( sumold % 13 != 0 ) {
            sumnew = 13;
        }

        break;
    default:
        break;
    }

    /* Return the new location */
    * xloc = xl1;
    * yloc = yl1;
    * zloc = zl1;
    * act = act1;
    /* sumnew returns a prime number indicating that a specific direction */
    /* has been chosen */
    return ( sumnew );
}

/* routine to return count of number of neighboring pixels for pixel */
/* (xck,yck,zck) which are not phase ph1, ph2, or ph3 which are input as */
/* parameters */
/* Calls no other routines */
/* Called by extettr, extfh3, extch, extafm, extpozz, extc3ah6 */
/* extfreidel, extcsh, and extstrat */
int CemhydMatStatus :: edgecnt(int xck, int yck, int zck, int ph1, int ph2, int ph3)
{
    int ixe, iye, ize, edgeback, x2, y2, z2, check;

    /* counter for number of neighboring pixels which are not ph1, ph2, or ph3 */
    edgeback = 0;

    /* Examine all pixels in a 3*3*3 box centered at (xck,yck,zck) */
    /* except for the central pixel */
    for ( ixe = ( -1 ); ixe <= 1; ixe++ ) {
        x2 = xck + ixe;
        for ( iye = ( -1 ); iye <= 1; iye++ ) {
            y2 = yck + iye;
            for ( ize = ( -1 ); ize <= 1; ize++ ) {
                if ( ( ixe != 0 ) || ( iye != 0 ) || ( ize != 0 ) ) {
                    z2 = zck + ize;
                    /* adjust to maintain periodic boundaries */
                    if ( x2 < 0 ) {
                        x2 = ( SYSIZEM1 );
                    } else if ( x2 >= SYSIZE ) {
                        x2 = 0;
                    }

                    if ( y2 < 0 ) {
                        y2 = ( SYSIZEM1 );
                    } else if ( y2 >= SYSIZE ) {
                        y2 = 0;
                    }

                    if ( z2 < 0 ) {
                        z2 = ( SYSIZEM1 );
                    } else if ( z2 >= SYSIZE ) {
                        z2 = 0;
                    }

                    check = mic [ x2 ] [ y2 ] [ z2 ];
                    if ( ( check != ph1 ) && ( check != ph2 ) && ( check != ph3 ) ) {
                        edgeback += 1;
                    }
                }
            }
        }
    }

    /* return number of neighboring pixels which are not ph1, ph2, or ph3 */
    return ( edgeback );
}

/* routine to add extra CSH when diffusing CSH reacts */
/* Called by movecsh */
/* Calls edgecnt */
void CemhydMatStatus :: extcsh(void)
{
    int numnear, xchr, ychr, zchr, fchr, check, msface;
    long int tries;

    fchr = 0;
    tries = 0;
    /* locate CSH at random location */
    /* in pore space in contact with at least another CSH or C3S or C2S */
    while ( fchr == 0 ) {
        tries += 1;
        /* generate a random location in the 3-D system */
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if location is porosity, locate the CSH there */
        if ( check == POROSITY ) {
            numnear = edgecnt(xchr, ychr, zchr, CSH, C3S, C2S);
            /* be sure that at least one neighboring pixel */
            /* is C2S, C3S, or diffusing CSH */
            if ( ( numnear < 26 ) || ( tries > 5000 ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = CSH;
                count [ CSH ] += 1;
                count [ POROSITY ] -= 1;
                cshage [ xchr ] [ ychr ] [ zchr ] = cyccnt;
                if ( cshgeom == 1 ) {
                    msface = ( int ) ( 3. * ran1(seed) + 1. );
                    if ( msface > 3 ) {
                        msface = 1;
                    }

                    faces [ xchr ] [ ychr ] [ zchr ] = msface;
                    ncshplateinit += 1;
                }

                fchr = 1;
            }
        }
    }
}

/* routine to move a diffusing CSH species */
/* Inputs: current location (xcur,ycur,zcur) and flag indicating if final */
/* step in diffusion process */
/* Returns flag indicating action taken (reaction or diffusion/no movement) */
/* Calls moveone,extcsh */
/* Called by hydrate */
int CemhydMatStatus :: movecsh(int xcur, int ycur, int zcur, int finalstep, int cycorig)
{
    int xnew, ynew, znew, action, sumback, sumin, check;
    int msface, mstest = 0, mstest2 = 0;
    float prcsh, prcsh1, prcsh2, prtest;

    action = 0;
    /* Store current location of species */
    xnew = xcur;
    ynew = ycur;
    znew = zcur;
    sumin = 1;
    sumback = moveone(& xnew, & ynew, & znew, & action, sumin);
    if ( cshgeom == 1 ) {
        /* Determine eligible faces based on direction of move */
        if ( xnew != xcur ) {
            mstest = 1;
            mstest2 = 2;
        }

        if ( ynew != ycur ) {
            mstest = 2;
            mstest2 = 3;
        }

        if ( znew != zcur ) {
            mstest = 3;
            mstest2 = 1;
        }
    }

    if ( action == 0 ) {
        printf("Error in value of action \n");
    }

    check = mic [ xnew ] [ ynew ] [ znew ];


    /* if new location is solid CSH and plate growth is favorable, */
    /* then convert diffusing CSH species to solid CSH */
    prcsh = ran1(seed);
    if ( ( check == CSH ) && ( ( cshgeom == 0 ) || ( faces [ xnew ] [ ynew ] [ znew ] == 0 ) || ( faces [ xnew ] [ ynew ] [ znew ] == mstest ) || ( faces [ xnew ] [ ynew ] [ znew ] == mstest2 ) ) ) {
        /* decrement count of diffusing CSH species */
        count [ DIFFCSH ] -= 1;
        /* and increment count of solid CSH if needed */
        prtest = molarvcsh [ cyccnt ] / molarvcsh [ cycorig ];
        prcsh1 = ran1(seed);
        if ( prcsh1 <= prtest ) {
            mic [ xcur ] [ ycur ] [ zcur ] = CSH;
            if ( cshgeom == 1 ) {
                faces [ xcur ] [ ycur ] [ zcur ] = faces [ xnew ] [ ynew ] [ znew ];
                ncshplategrow += 1;
            }

            cshage [ xcur ] [ ycur ] [ zcur ] = cyccnt;
            count [ CSH ] += 1;
        } else {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            count [ POROSITY ] += 1;
        }

        /* May need extra solid CSH if temperature goes down with time */
        if ( prtest > 1.0 ) {
            prcsh2 = ran1(seed);
            if ( prcsh2 < ( prtest - 1.0 ) ) {
                extcsh();
            }
        }

        action = 0;
    }
    /* Changed prcsh limit from 0.1 to 0.01 for CH test 1/27/05 */
    else if ( ( check == SLAGCSH ) || ( check == POZZCSH ) || ( finalstep == 1 ) ||
             ( ( ( check == C3S ) || ( check == C2S ) ) && ( prcsh < 0.001 ) ) ||
             ( ( ( check == C3A ) || ( check == C4AF ) ) && ( prcsh < 0.2 ) ) ||
             ( ( check == CH ) && ( prcsh < 0.01 ) ) ||
             ( check == CACO3 ) || ( check == INERT ) ) {
        /* decrement count of diffusing CSH species */
        count [ DIFFCSH ] -= 1;
        /* and increment count of solid CSH if needed */
        prtest = molarvcsh [ cyccnt ] / molarvcsh [ cycorig ];
        prcsh1 = ran1(seed);
        if ( prcsh1 <= prtest ) {
            mic [ xcur ] [ ycur ] [ zcur ] = CSH;
            cshage [ xcur ] [ ycur ] [ zcur ] = cyccnt;
            if ( cshgeom == 1 ) {
                msface = ( int ) ( 2. * ran1(seed) + 1. );
                if ( msface > 2 ) {
                    msface = 1;
                }

                if ( msface == 1 ) {
                    faces [ xcur ] [ ycur ] [ zcur ] = mstest;
                } else {
                    faces [ xcur ] [ ycur ] [ zcur ] = mstest2;
                }

                ncshplateinit += 1;
            }

            count [ CSH ] += 1;
        } else {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            count [ POROSITY ] += 1;
        }

        /* May need extra solid CSH if temperature goes down with time */
        if ( prtest > 1.0 ) {
            prcsh2 = ran1(seed);
            if ( prcsh2 < ( prtest - 1.0 ) ) {
                extcsh();
            }
        }

        action = 0;
    }

    if ( action != 0 ) {
        /* if diffusion step is possible, perform it */
        if ( check == POROSITY ) {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            mic [ xnew ] [ ynew ] [ znew ] = DIFFCSH;
        } else {
            /* indicate that diffusing CSH species remained */
            /* at original location */
            action = 7;
        }
    }

    return ( action );
}

/* routine to add extra FH3 when gypsum, hemihydrate, anhydrite, CAS2, or */
/* CaCl2 reacts with C4AF at location (xpres,ypres,zpres) */
/* Called by movegyp, moveettr, movecas2, movehem, moveanh, and movecacl2 */
/* Calls moveone and edgecnt */
void CemhydMatStatus :: extfh3(int xpres, int ypres, int zpres)
{
    int multf, numnear, sump, xchr, ychr, zchr, check, fchr, i1, newact;
    long int tries;

    /* first try 6 neighboring locations until      */
    /*    a) successful                */
    /*    b) all 6 sites are tried and full or    */
    /*    c) 500 tries are made           */
    fchr = 0;
    sump = 1;
    for ( i1 = 1; ( ( i1 <= 500 ) && ( fchr == 0 ) && ( sump != 30030 ) ); i1++ ) {
        /* choose a neighbor at random */
        xchr = xpres;
        ychr = ypres;
        zchr = zpres;
        newact = 0;
        multf = moveone(& xchr, & ychr, & zchr, & newact, sump);
        if ( newact == 0 ) {
            printf("Error in value of newact in extfh3 \n");
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if neighbor is porosity   */
        /* then locate the FH3 there */
        if ( check == POROSITY ) {
            mic [ xchr ] [ ychr ] [ zchr ] = FH3;
            count [ FH3 ] += 1;
            count [ POROSITY ] -= 1;
            fchr = 1;
        } else {
            sump *= multf;
        }
    }

    /* if no neighbor available, locate FH3 at random location */
    /* in pore space in contact with at least one FH3 */
    tries = 0;
    while ( fchr == 0 ) {
        tries += 1;
        /* generate a random location in the 3-D system */
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];
        /* if location is porosity, locate the FH3 there */
        if ( check == POROSITY ) {
            numnear = edgecnt(xchr, ychr, zchr, FH3, FH3, DIFFFH3);
            /* be sure that at least one neighboring pixel */
            /* is FH3 or diffusing FH3 */
            if ( ( numnear < 26 ) || ( tries > 5000 ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = FH3;
                count [ FH3 ] += 1;
                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }
}

/* routine to add extra ettringite when gypsum, anhydrite, or hemihydrate */
/* reacts with aluminates addition adjacent to location (xpres,ypres,zpres) */
/* in a fashion to preserve needle growth */
/* etype=0 indicates primary ettringite */
/* etype=1 indicates iron-rich stable ettringite */
/* Returns flag indicating action taken */
/* Calls moveone and edgecnt */
/* Called by movegyp, movehem, moveanh, and movec3a */
int CemhydMatStatus :: extettr(int xpres, int ypres, int zpres, int etype)
{
    int check, newact, multf, numnear, sump, xchr, ychr, zchr, fchr, i1;
    int numalum, numsil;
    float pneigh, ptest;
    long int tries;

    /* first try neighboring locations until        */
    /*    a) successful                */
    /*    c) 1000 tries are made          */
    fchr = 0;
    sump = 1;
    /* Note that 30030 = 2*3*5*7*11*13 */
    /* indicating that all six sites have been tried */
    for ( i1 = 1; ( ( i1 <= 1000 ) && ( fchr == 0 ) ); i1++ ) {
        /* determine location of neighbor (using periodic boundaries) */
        xchr = xpres;
        ychr = ypres;
        zchr = zpres;
        newact = 0;
        multf = moveone(& xchr, & ychr, & zchr, & newact, sump);
        if ( newact == 0 ) {
            printf("Error in value of action \n");
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if neighbor is porosity, and conditions are favorable */
        /* based on number of neighboring ettringite, C3A, or C4AF */
        /* pixels then locate the ettringite there */
        if ( check == POROSITY ) {
            /* be sure ettringite doesn't touch C3S */
            numsil = edgecnt(xchr, ychr, zchr, C3S, C2S, C3S);
            numsil = 26 - numsil;
            if ( etype == 0 ) {
                numnear = edgecnt(xchr, ychr, zchr, ETTR, ETTR, ETTR);
                numalum = edgecnt(xchr, ychr, zchr, C3A, C3A, C3A);
                numalum = 26 - numalum;
            } else {
                numnear = edgecnt(xchr, ychr, zchr, ETTRC4AF, ETTRC4AF, ETTRC4AF);
                numalum = edgecnt(xchr, ychr, zchr, C4AF, C4AF, C4AF);
                numalum = 26 - numalum;
            }

            pneigh = ( float ) ( numnear + 1 ) / 26.0;
            pneigh *= pneigh;
            if ( numalum <= 1 ) {
                pneigh = 0.0;
            }

            if ( numalum >= 2 ) {
                pneigh += 0.5;
            }

            if ( numalum >= 3 ) {
                pneigh += 0.25;
            }

            if ( numalum >= 5 ) {
                pneigh += 0.25;
            }

            ptest = ran1(seed);
            if ( numsil < 1 ) {
                if ( pneigh >= ptest ) {
                    if ( etype == 0 ) {
                        mic [ xchr ] [ ychr ] [ zchr ] = ETTR;
                        count [ ETTR ] += 1;
                    } else {
                        mic [ xchr ] [ ychr ] [ zchr ] = ETTRC4AF;
                        count [ ETTRC4AF ] += 1;
                    }

                    fchr = 1;
                    count [ POROSITY ] -= 1;
                }
            }
        }
    }

    /* if no neighbor available, locate ettringite at random location */
    /* in pore space in contact with at least another ettringite */
    /* or aluminate surface  */
    tries = 0;
    while ( fchr == 0 ) {
        tries += 1;
        newact = 7;
        /* generate a random location in the 3-D system */
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];
        /* if location is porosity, locate the ettringite there */
        if ( check == POROSITY ) {
            numsil = edgecnt(xchr, ychr, zchr, C3S, C2S, C3S);
            numsil = 26 - numsil;
            if ( etype == 0 ) {
                numnear = edgecnt(xchr, ychr, zchr, ETTR, C3A, C4AF);
            } else {
                numnear = edgecnt(xchr, ychr, zchr, ETTRC4AF, C3A, C4AF);
            }

            /* be sure that at least one neighboring pixel */
            /* is ettringite, or aluminate clinker */
            if ( ( tries > 5000 ) || ( ( numnear < 26 ) && ( numsil < 1 ) ) ) {
                if ( etype == 0 ) {
                    mic [ xchr ] [ ychr ] [ zchr ] = ETTR;
                    count [ ETTR ] += 1;
                } else {
                    mic [ xchr ] [ ychr ] [ zchr ] = ETTRC4AF;
                    count [ ETTRC4AF ] += 1;
                }

                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }

    return ( newact );
}

/* routine to add extra CH when gypsum, hemihydrate, anhydrite, CaCl2, or */
/* diffusing CAS2  reacts with C4AF */
/* Called by movegyp, movehem, moveanh, moveettr, movecas2, and movecacl2 */
/* Calls edgecnt */
void CemhydMatStatus :: extch(void)
{
    int numnear, xchr, ychr, zchr, fchr, check;
    long int tries;

    fchr = 0;
    tries = 0;
    /* locate CH at random location */
    /* in pore space in contact with at least another CH */
    while ( fchr == 0 ) {
        tries += 1;
        /* generate a random location in the 3-D system */
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if location is porosity, locate the CH there */
        if ( check == POROSITY ) {
            numnear = edgecnt(xchr, ychr, zchr, CH, DIFFCH, CH);
            /* be sure that at least one neighboring pixel */
            /* is CH or diffusing CH */
            if ( ( numnear < 26 ) || ( tries > 5000 ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = CH;
                count [ CH ] += 1;
                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }
}

/* routine to add extra gypsum when hemihydrate or anhydrite hydrates */
/* Called by movehem and moveanh */
/* Calls moveone and edgecnt */
void CemhydMatStatus :: extgyps(int xpres, int ypres, int zpres)
{
    int multf, numnear, sump, xchr, ychr, zchr, check, fchr, i1, newact;
    long int tries;

    /* first try 6 neighboring locations until      */
    /*    a) successful                */
    /*    b) all 6 sites are tried and full or    */
    /*    c) 500 tries are made           */
    fchr = 0;
    sump = 1;
    for ( i1 = 1; ( ( i1 <= 500 ) && ( fchr == 0 ) && ( sump != 30030 ) ); i1++ ) {
        /* choose a neighbor at random */
        xchr = xpres;
        ychr = ypres;
        zchr = zpres;
        newact = 0;
        multf = moveone(& xchr, & ychr, & zchr, & newact, sump);
        if ( newact == 0 ) {
            printf("Error in value of newact in extfh3 \n");
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if neighbor is porosity   */
        /* then locate the GYPSUMS there */
        if ( check == POROSITY ) {
            mic [ xchr ] [ ychr ] [ zchr ] = GYPSUMS;
            count [ GYPSUMS ] += 1;
            count [ POROSITY ] -= 1;
            fchr = 1;
        } else {
            sump *= multf;
        }
    }

    /* if no neighbor available, locate GYPSUMS at random location */
    /* in pore space in contact with at least one GYPSUMS */
    tries = 0;
    while ( fchr == 0 ) {
        tries += 1;
        /* generate a random location in the 3-D system */
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];
        /* if location is porosity, locate the GYPSUMS there */
        if ( check == POROSITY ) {
            numnear = edgecnt(xchr, ychr, zchr, HEMIHYD, GYPSUMS, ANHYDRITE);
            /* be sure that at least one neighboring pixel */
            /* is Gypsum in some form */
            if ( ( numnear < 26 ) || ( tries > 5000 ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = GYPSUMS;
                count [ GYPSUMS ] += 1;
                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }
}

/* routine to move a diffusing ANHYDRITE species */
/* Inputs: current location (xcur,ycur,zcur) and flag indicating if final */
/* step in diffusion process */
/* Returns flag indicating action taken (reaction or diffusion/no movement) */
/* Calls moveone */
/* Called by hydrate */
int CemhydMatStatus :: moveanh(int xcur, int ycur, int zcur, int finalstep, float nucprgyp)
{
    int xnew, ynew, znew, action, sumback, sumin, check = -1;
    int nexp, iexp, xexp, yexp, zexp, newact, ettrtype;
    float pgen, pexp, pext, p2diff;

    action = 0;
    /* first check for nucleation */
    pgen = ran1(seed);
    p2diff = ran1(seed);
    if ( ( nucprgyp >= pgen ) || ( finalstep == 1 ) ) {
        action = 0;
        mic [ xcur ] [ ycur ] [ zcur ] = GYPSUMS;
        count [ DIFFANH ] -= 1;
        count [ GYPSUMS ] += 1;
        pexp = ran1(seed);
        if ( pexp < 0.4 ) {
            extgyps(xcur, ycur, zcur);
        }
    } else {
        /* Store current location of species */
        xnew = xcur;
        ynew = ycur;
        znew = zcur;
        sumin = 1;
        sumback = moveone(& xnew, & ynew, & znew, & action, sumin);

        if ( action == 0 ) {
            printf("Error in value of action \n");
        }

        check = mic [ xnew ] [ ynew ] [ znew ];

        /* if new location is solid GYPSUM(S) or diffusing GYPSUM, then convert */
        /* diffusing ANHYDRITE species to solid GYPSUM */
        if ( ( check == GYPSUM ) || ( check == GYPSUMS ) || ( check == DIFFGYP ) ) {
            mic [ xcur ] [ ycur ] [ zcur ] = GYPSUMS;
            /* decrement count of diffusing ANHYDRITE species */
            /* and increment count of solid GYPSUMS */
            count [ DIFFANH ] -= 1;
            count [ GYPSUMS ] += 1;
            action = 0;
            /* Add extra gypsum as necessary */
            pexp = ran1(seed);
            if ( pexp < 0.4 ) {
                extgyps(xnew, ynew, znew);
            }
        }
        /* if new location is C3A or diffusing C3A, execute conversion */
        /* to ettringite (including necessary volumetric expansion) */
        else if ( ( ( check == C3A ) && ( p2diff < SOLIDC3AGYP ) ) || ( ( check == DIFFC3A ) && ( p2diff < C3AGYP ) ) || ( ( check == DIFFC4A ) && ( p2diff < C3AGYP ) ) ) {
            /* Convert diffusing gypsum to an ettringite pixel */
            ettrtype = 0;
            mic [ xcur ] [ ycur ] [ zcur ] = ETTR;
            if ( check == DIFFC4A ) {
                ettrtype = 1;
                mic [ xcur ] [ ycur ] [ zcur ] = ETTRC4AF;
            }

            action = 0;
            count [ DIFFANH ] -= 1;
            count [ check ] -= 1;

            /* determine if C3A should be converted to ettringite */
            /* 1 unit of hemihydrate requires 0.569 units of C3A */
            /* and should form 4.6935 units of ettringite */
            pexp = ran1(seed);
            nexp = 3;
            if ( pexp <= 0.569 ) {
                if ( ettrtype == 0 ) {
                    mic [ xnew ] [ ynew ] [ znew ] = ETTR;
                    count [ ETTR ] += 1;
                } else {
                    mic [ xnew ] [ ynew ] [ znew ] = ETTRC4AF;
                    count [ ETTRC4AF ] += 1;
                }

                nexp = 2;
            } else {
                /* maybe someday, use a new FIXEDC3A here */
                /* so it won't dissolve later */
                if ( check == C3A ) {
                    mic [ xnew ] [ ynew ] [ znew ] = C3A;
                    count [ C3A ] += 1;
                } else {
                    if ( ettrtype == 0 ) {
                        count [ DIFFC3A ] += 1;
                        mic [ xnew ] [ ynew ] [ znew ] = DIFFC3A;
                    } else {
                        count [ DIFFC4A ] += 1;
                        mic [ xnew ] [ ynew ] [ znew ] = DIFFC4A;
                    }
                }

                nexp = 3;
            }

            /* create extra ettringite pixels to maintain volume stoichiometry */
            /* xexp, yexp, and zexp hold coordinates of most recently added ettringite */
            /* species as we attempt to grow a needle like structure */
            xexp = xcur;
            yexp = ycur;
            zexp = zcur;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extettr(xexp, yexp, zexp, ettrtype);
                /* update xexp, yexp and zexp as needed */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last ettringite pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.6935 ) {
                newact = extettr(xexp, yexp, zexp, ettrtype);
            }
        }

        /* if new location is C4AF execute conversion */
        /* to ettringite (including necessary volumetric expansion) */
        if ( ( check == C4AF ) && ( p2diff < SOLIDC4AFGYP ) ) {
            mic [ xcur ] [ ycur ] [ zcur ] = ETTRC4AF;
            count [ ETTRC4AF ] += 1;
            count [ DIFFANH ] -= 1;

            /* determine if C4AF should be converted to ettringite */
            /* 1 unit of gypsum requires 0.8174 units of C4AF */
            /* and should form 4.6935 units of ettringite */
            pexp = ran1(seed);
            nexp = 3;
            if ( pexp <= 0.8174 ) {
                mic [ xnew ] [ ynew ] [ znew ] = ETTRC4AF;
                count [ ETTRC4AF ] += 1;
                count [ C4AF ] -= 1;
                nexp = 2;
                pext = ran1(seed);
                /* Addition of extra CH */
                if ( pext < 0.2584 ) {
                    extch();
                }

                pext = ran1(seed);
                /* Addition of extra FH3 */
                if ( pext < 0.5453 ) {
                    extfh3(xnew, ynew, znew);
                }
            } else {
                /* maybe someday, use a new FIXEDC4AF here */
                /* so it won't dissolve later */
                mic [ xnew ] [ ynew ] [ znew ] = C4AF;
                nexp = 3;
            }

            /* create extra ettringite pixels to maintain volume stoichiometry */
            /* xexp, yexp and zexp hold coordinates of most recently added ettringite */
            /* species as we attempt to grow a needle like structure */
            xexp = xcur;
            yexp = ycur;
            zexp = zcur;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extettr(xexp, yexp, zexp, 1);
                /* update xexp, yexp and zexp as needed */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last ettringite pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.6935 ) {
                newact = extettr(xexp, yexp, zexp, 1);
            }

            action = 0;
        }
    }

    if ( action != 0 ) {
        /* if diffusion step is possible, perform it */
        if ( check == POROSITY ) {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            mic [ xnew ] [ ynew ] [ znew ] = DIFFANH;
        } else {
            /* indicate that diffusing ANHYDRITE species remained */
            /* at original location */
            action = 7;
        }
    }

    return ( action );
}

/* routine to move a diffusing HEMIHYDRATE species */
/* Inputs: current location (xcur,ycur,zcur) and flag indicating if final */
/* step in diffusion process */
/* Returns flag indicating action taken (reaction or diffusion/no movement) */
/* Calls moveone, extettr, extch, and extfh3 */
/* Called by hydrate */
int CemhydMatStatus :: movehem(int xcur, int ycur, int zcur, int finalstep, float nucprgyp)
{
    int xnew, ynew, znew, action, sumback, sumin, check = -1;
    int nexp, iexp, xexp, yexp, zexp, newact, ettrtype;
    float pgen, pexp, pext, p2diff;

    action = 0;
    /* first check for nucleation */
    pgen = ran1(seed);
    p2diff = ran1(seed);
    if ( ( nucprgyp >= pgen ) || ( finalstep == 1 ) ) {
        action = 0;
        mic [ xcur ] [ ycur ] [ zcur ] = GYPSUMS;
        count [ DIFFHEM ] -= 1;
        count [ GYPSUMS ] += 1;
        /* Add extra gypsum as necessary */
        pexp = ran1(seed);
        if ( pexp < 0.4 ) {
            extgyps(xcur, ycur, zcur);
        }
    } else {
        /* Store current location of species */
        xnew = xcur;
        ynew = ycur;
        znew = zcur;
        sumin = 1;
        sumback = moveone(& xnew, & ynew, & znew, & action, sumin);

        if ( action == 0 ) {
            printf("Error in value of action \n");
        }

        check = mic [ xnew ] [ ynew ] [ znew ];

        /* if new location is solid GYPSUM(S) or diffusing GYPSUM, then convert */
        /* diffusing HEMIHYDRATE species to solid GYPSUM */
        if ( ( check == GYPSUM ) || ( check == GYPSUMS ) || ( check == DIFFGYP ) ) {
            mic [ xcur ] [ ycur ] [ zcur ] = GYPSUMS;
            /* decrement count of diffusing HEMIHYDRATE species */
            /* and increment count of solid GYPSUMS */
            count [ DIFFHEM ] -= 1;
            count [ GYPSUMS ] += 1;
            action = 0;
            /* Add extra gypsum as necessary */
            pexp = ran1(seed);
            if ( pexp < 0.4 ) {
                extgyps(xnew, ynew, znew);
            }
        }
        /* if new location is C3A or diffusing C3A, execute conversion */
        /* to ettringite (including necessary volumetric expansion) */
        else if ( ( ( check == C3A ) && ( p2diff < SOLIDC3AGYP ) ) || ( ( check == DIFFC3A ) && ( p2diff < C3AGYP ) ) || ( ( check == DIFFC4A ) && ( p2diff < C3AGYP ) ) ) {
            /* Convert diffusing gypsum to an ettringite pixel */
            ettrtype = 0;
            mic [ xcur ] [ ycur ] [ zcur ] = ETTR;
            if ( check == DIFFC4A ) {
                ettrtype = 1;
                mic [ xcur ] [ ycur ] [ zcur ] = ETTRC4AF;
            }

            action = 0;
            count [ DIFFHEM ] -= 1;
            count [ check ] -= 1;

            /* determine if C3A should be converted to ettringite */
            /* 1 unit of hemihydrate requires 0.5583 units of C3A */
            /* and should form 4.6053 units of ettringite */
            pexp = ran1(seed);
            nexp = 3;
            if ( pexp <= 0.5583 ) {
                if ( ettrtype == 0 ) {
                    mic [ xnew ] [ ynew ] [ znew ] = ETTR;
                    count [ ETTR ] += 1;
                } else {
                    mic [ xnew ] [ ynew ] [ znew ] = ETTRC4AF;
                    count [ ETTRC4AF ] += 1;
                }

                nexp = 2;
            } else {
                /* maybe someday, use a new FIXEDC3A here */
                /* so it won't dissolve later */
                if ( check == C3A ) {
                    mic [ xnew ] [ ynew ] [ znew ] = C3A;
                    count [ C3A ] += 1;
                } else {
                    if ( ettrtype == 0 ) {
                        count [ DIFFC3A ] += 1;
                        mic [ xnew ] [ ynew ] [ znew ] = DIFFC3A;
                    } else {
                        count [ DIFFC4A ] += 1;
                        mic [ xnew ] [ ynew ] [ znew ] = DIFFC4A;
                    }
                }

                nexp = 3;
            }

            /* create extra ettringite pixels to maintain volume stoichiometry */
            /* xexp, yexp, and zexp hold coordinates of most recently added ettringite */
            /* species as we attempt to grow a needle like structure */
            xexp = xcur;
            yexp = ycur;
            zexp = zcur;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extettr(xexp, yexp, zexp, ettrtype);
                /* update xexp, yexp and zexp as needed */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last ettringite pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.6053 ) {
                newact = extettr(xexp, yexp, zexp, ettrtype);
            }
        }

        /* if new location is C4AF execute conversion */
        /* to ettringite (including necessary volumetric expansion) */
        if ( ( check == C4AF ) && ( p2diff < SOLIDC4AFGYP ) ) {
            mic [ xcur ] [ ycur ] [ zcur ] = ETTRC4AF;
            count [ ETTRC4AF ] += 1;
            count [ DIFFHEM ] -= 1;

            /* determine if C4AF should be converted to ettringite */
            /* 1 unit of gypsum requires 0.802 units of C4AF */
            /* and should form 4.6053 units of ettringite */
            pexp = ran1(seed);
            nexp = 3;
            if ( pexp <= 0.802 ) {
                mic [ xnew ] [ ynew ] [ znew ] = ETTRC4AF;
                count [ ETTRC4AF ] += 1;
                count [ C4AF ] -= 1;
                nexp = 2;
                pext = ran1(seed);
                /* Addition of extra CH */
                if ( pext < 0.2584 ) {
                    extch();
                }

                pext = ran1(seed);
                /* Addition of extra FH3 */
                if ( pext < 0.5453 ) {
                    extfh3(xnew, ynew, znew);
                }
            } else {
                /* maybe someday, use a new FIXEDC4AF here */
                /* so it won't dissolve later */
                mic [ xnew ] [ ynew ] [ znew ] = C4AF;
                nexp = 3;
            }

            /* create extra ettringite pixels to maintain volume stoichiometry */
            /* xexp, yexp and zexp hold coordinates of most recently added ettringite */
            /* species as we attempt to grow a needle like structure */
            xexp = xcur;
            yexp = ycur;
            zexp = zcur;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extettr(xexp, yexp, zexp, 1);
                /* update xexp, yexp and zexp as needed */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last ettringite pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.6053 ) {
                newact = extettr(xexp, yexp, zexp, 1);
            }

            action = 0;
        }
    }

    if ( action != 0 ) {
        /* if diffusion step is possible, perform it */
        if ( check == POROSITY ) {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            mic [ xnew ] [ ynew ] [ znew ] = DIFFHEM;
        } else {
            /* indicate that diffusing HEMIHYDRATE species */
            /* remained at original location */
            action = 7;
        }
    }

    return ( action );
}

/* routine to add extra Freidel's salt when CaCl2 reacts with */
/* C3A or C4AF at location (xpres,ypres,zpres) */
/* Called by movecacl2 and movec3a */
/* Calls moveone and edgecnt */
int CemhydMatStatus :: extfreidel(int xpres, int ypres, int zpres)
{
    int multf, numnear, sump, xchr, ychr, zchr, check, fchr, i1, newact;
    long int tries;

    /* first try 6 neighboring locations until      */
    /*    a) successful                */
    /*    b) all 6 sites are tried and full or    */
    /*    c) 500 tries are made           */
    fchr = 0;
    sump = 1;
    for ( i1 = 1; ( ( i1 <= 500 ) && ( fchr == 0 ) && ( sump != 30030 ) ); i1++ ) {
        /* choose a neighbor at random */
        xchr = xpres;
        ychr = ypres;
        zchr = zpres;
        newact = 0;
        multf = moveone(& xchr, & ychr, & zchr, & newact, sump);
        if ( newact == 0 ) {
            printf("Error in value of newact in extfreidel \n");
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if neighbor is porosity   */
        /* then locate the freidel's salt there */
        if ( check == POROSITY ) {
            mic [ xchr ] [ ychr ] [ zchr ] = FREIDEL;
            count [ FREIDEL ] += 1;
            count [ POROSITY ] -= 1;
            fchr = 1;
        } else {
            sump *= multf;
        }
    }

    /* if no neighbor available, locate FREIDEL at random location */
    /* in pore space in contact with at least one FREIDEL */
    tries = 0;
    while ( fchr == 0 ) {
        tries += 1;
        newact = 7;
        /* generate a random location in the 3-D system */
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];
        /* if location is porosity, locate the FREIDEL there */
        if ( check == POROSITY ) {
            numnear = edgecnt(xchr, ychr, zchr, FREIDEL, FREIDEL, DIFFCACL2);
            /* be sure that at least one neighboring pixel */
            /* is FREIDEL or diffusing CACL2 */
            if ( ( numnear < 26 ) || ( tries > 5000 ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = FREIDEL;
                count [ FREIDEL ] += 1;
                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }

    return ( newact );
}

/* routine to add extra stratlingite when AS reacts with */
/* CH at location (xpres,ypres,zpres) */
/* or when diffusing CAS2 reacts with aluminates */
/* Called by moveas, movech, and movecas2 */
/* Calls moveone and edgecnt */
int CemhydMatStatus :: extstrat(int xpres, int ypres, int zpres)
{
    int multf, numnear, sump, xchr, ychr, zchr, check, fchr, i1, newact;
    long int tries;

    /* first try 6 neighboring locations until      */
    /*    a) successful                */
    /*    b) all 6 sites are tried and full or    */
    /*    c) 500 tries are made           */
    fchr = 0;
    sump = 1;
    for ( i1 = 1; ( ( i1 <= 500 ) && ( fchr == 0 ) && ( sump != 30030 ) ); i1++ ) {
        /* choose a neighbor at random */
        xchr = xpres;
        ychr = ypres;
        zchr = zpres;
        newact = 0;
        multf = moveone(& xchr, & ychr, & zchr, & newact, sump);
        if ( newact == 0 ) {
            printf("Error in value of newact in extstrat \n");
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if neighbor is porosity   */
        /* then locate the stratlingite there */
        if ( check == POROSITY ) {
            mic [ xchr ] [ ychr ] [ zchr ] = STRAT;
            count [ STRAT ] += 1;
            count [ POROSITY ] -= 1;
            fchr = 1;
        } else {
            sump *= multf;
        }
    }

    /* if no neighbor available, locate STRAT at random location */
    /* in pore space in contact with at least one STRAT */
    tries = 0;
    while ( fchr == 0 ) {
        tries += 1;
        newact = 7;
        /* generate a random location in the 3-D system */
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];
        /* if location is porosity, locate the STRAT there */
        if ( check == POROSITY ) {
            numnear = edgecnt(xchr, ychr, zchr, STRAT, DIFFCAS2, DIFFAS);
            /* be sure that at least one neighboring pixel */
            /* is STRAT, diffusing CAS2, or diffusing AS */
            if ( ( numnear < 26 ) || ( tries > 5000 ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = STRAT;
                count [ STRAT ] += 1;
                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }

    return ( newact );
}

/* routine to move a diffusing gypsum species */
/* from current location (xcur,ycur,zcur) */
/* Returns action flag indicating response taken */
/* Called by hydrate */
/* Calls moveone, extettr, extch, and extfh3 */
int CemhydMatStatus :: movegyp(int xcur, int ycur, int zcur, int finalstep)
{
    int check, xnew, ynew, znew, action, nexp, iexp;
    int xexp, yexp, zexp, newact, sumold, sumgarb, ettrtype;
    float pexp, pext, p2diff;

    sumold = 1;

    /* First be sure that a diffusing gypsum species is located at xcur,ycur,zcur */
    /* if not, return to calling routine */
    if ( mic [ xcur ] [ ycur ] [ zcur ] != DIFFGYP ) {
        action = 0;
        return ( action );
    }

    /* Determine new coordinates (periodic boundaries are used) */
    xnew = xcur;
    ynew = ycur;
    znew = zcur;
    action = 0;
    sumgarb = moveone(& xnew, & ynew, & znew, & action, sumold);
    if ( action == 0 ) {
        printf("Error in value of action in movegyp \n");
    }

    check = mic [ xnew ] [ ynew ] [ znew ];
    p2diff = ran1(seed);
    /* if new location is CSH, check for absorption of gypsum */
    if ( ( check == CSH ) && ( ( float ) count [ ABSGYP ] < ( gypabsprob * ( float ) count [ CSH ] ) ) ) {
        pexp = ran1(seed);
        if ( pexp < AGRATE ) {
            /* update counts for absorbed and diffusing gypsum */
            count [ ABSGYP ] += 1;
            count [ DIFFGYP ] -= 1;
            mic [ xcur ] [ ycur ] [ zcur ] = ABSGYP;
            action = 0;
        }
    }
    /* if new location is C3A or diffusing C3A, execute conversion */
    /* to ettringite (including necessary volumetric expansion) */
    /* Use p2diff to try to favor formation of ettringite on */
    /* aluminate surfaces as opposed to in solution */
    else if ( ( ( check == C3A ) && ( p2diff < SOLIDC3AGYP ) ) || ( ( check == DIFFC3A ) && ( p2diff < C3AGYP ) ) || ( ( check == DIFFC4A ) && ( p2diff < C3AGYP ) ) ) {
        /* Convert diffusing gypsum to an ettringite pixel */
        ettrtype = 0;
        mic [ xcur ] [ ycur ] [ zcur ] = ETTR;
        if ( check == DIFFC4A ) {
            ettrtype = 1;
            mic [ xcur ] [ ycur ] [ zcur ] = ETTRC4AF;
        }

        action = 0;
        count [ DIFFGYP ] -= 1;
        count [ check ] -= 1;

        /* determine if C3A should be converted to ettringite */
        /* 1 unit of gypsum requires 0.40 units of C3A */
        /* and should form 3.30 units of ettringite */
        pexp = ran1(seed);
        nexp = 2;
        if ( pexp <= 0.40 ) {
            if ( ettrtype == 0 ) {
                mic [ xnew ] [ ynew ] [ znew ] = ETTR;
                count [ ETTR ] += 1;
            } else {
                mic [ xnew ] [ ynew ] [ znew ] = ETTRC4AF;
                count [ ETTRC4AF ] += 1;
            }

            nexp = 1;
        } else {
            /* maybe someday, use a new FIXEDC3A here */
            /* so it won't dissolve later */
            if ( check == C3A ) {
                mic [ xnew ] [ ynew ] [ znew ] = C3A;
                count [ C3A ] += 1;
            } else {
                if ( ettrtype == 0 ) {
                    count [ DIFFC3A ] += 1;
                    mic [ xnew ] [ ynew ] [ znew ] = DIFFC3A;
                } else {
                    count [ DIFFC4A ] += 1;
                    mic [ xnew ] [ ynew ] [ znew ] = DIFFC4A;
                }
            }

            nexp = 2;
        }

        /* create extra ettringite pixels to maintain volume stoichiometry */
        /* xexp, yexp, and zexp hold coordinates of most recently added ettringite */
        /* species as we attempt to grow a needle like structure */
        xexp = xcur;
        yexp = ycur;
        zexp = zcur;
        for ( iexp = 1; iexp <= nexp; iexp++ ) {
            newact = extettr(xexp, yexp, zexp, ettrtype);
            /* update xexp, yexp and zexp as needed */
            switch ( newact ) {
            case 1:
                xexp -= 1;
                if ( xexp < 0 ) {
                    xexp = ( SYSIZEM1 );
                }

                break;
            case 2:
                xexp += 1;
                if ( xexp >= SYSIZE ) {
                    xexp = 0;
                }

                break;
            case 3:
                yexp -= 1;
                if ( yexp < 0 ) {
                    yexp = ( SYSIZEM1 );
                }

                break;
            case 4:
                yexp += 1;
                if ( yexp >= SYSIZE ) {
                    yexp = 0;
                }

                break;
            case 5:
                zexp -= 1;
                if ( zexp < 0 ) {
                    zexp = ( SYSIZEM1 );
                }

                break;
            case 6:
                zexp += 1;
                if ( zexp >= SYSIZE ) {
                    zexp = 0;
                }

                break;
            default:
                break;
            }
        }

        /* probabilistic-based expansion for last ettringite pixel */
        pexp = ran1(seed);
        if ( pexp <= 0.30 ) {
            newact = extettr(xexp, yexp, zexp, ettrtype);
        }
    }

    /* if new location is C4AF execute conversion */
    /* to ettringite (including necessary volumetric expansion) */
    if ( ( check == C4AF ) && ( p2diff < SOLIDC4AFGYP ) ) {
        mic [ xcur ] [ ycur ] [ zcur ] = ETTRC4AF;
        count [ ETTRC4AF ] += 1;
        count [ DIFFGYP ] -= 1;

        /* determine if C4AF should be converted to ettringite */
        /* 1 unit of gypsum requires 0.575 units of C4AF */
        /* and should form 3.30 units of ettringite */
        pexp = ran1(seed);
        nexp = 2;
        if ( pexp <= 0.575 ) {
            mic [ xnew ] [ ynew ] [ znew ] = ETTRC4AF;
            count [ ETTRC4AF ] += 1;
            count [ C4AF ] -= 1;
            nexp = 1;
            pext = ran1(seed);
            /* Addition of extra CH */
            if ( pext < 0.2584 ) {
                extch();
            }

            pext = ran1(seed);
            /* Addition of extra FH3 */
            if ( pext < 0.5453 ) {
                extfh3(xnew, ynew, znew);
            }
        } else {
            /* maybe someday, use a new FIXEDC4AF here */
            /* so it won't dissolve later */
            mic [ xnew ] [ ynew ] [ znew ] = C4AF;
            nexp = 2;
        }

        /* create extra ettringite pixels to maintain volume stoichiometry */
        /* xexp, yexp and zexp hold coordinates of most recently added ettringite */
        /* species as we attempt to grow a needle like structure */
        xexp = xcur;
        yexp = ycur;
        zexp = zcur;
        for ( iexp = 1; iexp <= nexp; iexp++ ) {
            newact = extettr(xexp, yexp, zexp, 1);
            /* update xexp, yexp and zexp as needed */
            switch ( newact ) {
            case 1:
                xexp -= 1;
                if ( xexp < 0 ) {
                    xexp = ( SYSIZEM1 );
                }

                break;
            case 2:
                xexp += 1;
                if ( xexp >= SYSIZE ) {
                    xexp = 0;
                }

                break;
            case 3:
                yexp -= 1;
                if ( yexp < 0 ) {
                    yexp = ( SYSIZEM1 );
                }

                break;
            case 4:
                yexp += 1;
                if ( yexp >= SYSIZE ) {
                    yexp = 0;
                }

                break;
            case 5:
                zexp -= 1;
                if ( zexp < 0 ) {
                    zexp = ( SYSIZEM1 );
                }

                break;
            case 6:
                zexp += 1;
                if ( zexp >= SYSIZE ) {
                    zexp = 0;
                }

                break;
            default:
                break;
            }
        }

        /* probabilistic-based expansion for last ettringite pixel */
        pexp = ran1(seed);
        if ( pexp <= 0.30 ) {
            newact = extettr(xexp, yexp, zexp, 1);
        }

        action = 0;
    }

    /* if last diffusion step and no reaction, convert back to */
    /* primary solid gypsum */
    if ( ( action != 0 ) && ( finalstep == 1 ) ) {
        action = 0;
        count [ DIFFGYP ] -= 1;
        count [ GYPSUM ] += 1;
        mic [ xcur ] [ ycur ] [ zcur ] = GYPSUM;
    }

    if ( action != 0 ) {
        /* if diffusion is possible, execute it */
        if ( check == POROSITY ) {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            mic [ xnew ] [ ynew ] [ znew ] = DIFFGYP;
        } else {
            /* indicate that diffusing gypsum remained at */
            /* original location */
            action = 7;
        }
    }

    return ( action );
}

/* routine to move a diffusing CaCl2 species */
/* from current location (xcur,ycur,zcur) */
/* Returns action flag indicating response taken */
/* Called by hydrate */
/* Calls moveone, extfreidel, extch, and extfh3 */
int CemhydMatStatus :: movecacl2(int xcur, int ycur, int zcur, int finalstep)
{
    int check, xnew, ynew, znew, action, nexp, iexp;
    int xexp, yexp, zexp, newact, sumold, sumgarb, keep;
    float pexp, pext;

    sumold = 1;
    keep = 0;

    /* First be sure that a diffusing CaCl2 species is located at xcur,ycur,zcur */
    /* if not, return to calling routine */
    if ( mic [ xcur ] [ ycur ] [ zcur ] != DIFFCACL2 ) {
        action = 0;
        return ( action );
    }

    /* Determine new coordinates (periodic boundaries are used) */
    xnew = xcur;
    ynew = ycur;
    znew = zcur;
    action = 0;
    sumgarb = moveone(& xnew, & ynew, & znew, & action, sumold);
    if ( action == 0 ) {
        printf("Error in value of action in movecacl2 \n");
    }

    check = mic [ xnew ] [ ynew ] [ znew ];

    /* if new location is C3A or diffusing C3A, execute conversion */
    /* to freidel's salt (including necessary volumetric expansion) */
    if ( ( check == C3A ) || ( check == DIFFC3A ) || ( check == DIFFC4A ) ) {
        /* Convert diffusing C3A or C3A to a freidel's salt pixel */
        action = 0;
        mic [ xnew ] [ ynew ] [ znew ] = FREIDEL;
        count [ FREIDEL ] += 1;
        count [ check ] -= 1;

        /* determine if diffusing CaCl2 should be converted to FREIDEL */
        /* 0.5793 unit of CaCl2 requires 1 unit of C3A */
        /* and should form 3.3295 units of FREIDEL */
        pexp = ran1(seed);
        nexp = 2;
        if ( pexp <= 0.5793 ) {
            mic [ xcur ] [ ycur ] [ zcur ] = FREIDEL;
            count [ FREIDEL ] += 1;
            count [ DIFFCACL2 ] -= 1;
            nexp = 1;
        } else {
            keep = 1;
            nexp = 2;
        }

        /* create extra Freidel's salt pixels to maintain volume stoichiometry */
        /* xexp, yexp, and zexp hold coordinates of most recently added FREIDEL */
        xexp = xcur;
        yexp = ycur;
        zexp = zcur;
        for ( iexp = 1; iexp <= nexp; iexp++ ) {
            newact = extfreidel(xexp, yexp, zexp);
            /* update xexp, yexp and zexp as needed */
            switch ( newact ) {
            case 1:
                xexp -= 1;
                if ( xexp < 0 ) {
                    xexp = ( SYSIZEM1 );
                }

                break;
            case 2:
                xexp += 1;
                if ( xexp >= SYSIZE ) {
                    xexp = 0;
                }

                break;
            case 3:
                yexp -= 1;
                if ( yexp < 0 ) {
                    yexp = ( SYSIZEM1 );
                }

                break;
            case 4:
                yexp += 1;
                if ( yexp >= SYSIZE ) {
                    yexp = 0;
                }

                break;
            case 5:
                zexp -= 1;
                if ( zexp < 0 ) {
                    zexp = ( SYSIZEM1 );
                }

                break;
            case 6:
                zexp += 1;
                if ( zexp >= SYSIZE ) {
                    zexp = 0;
                }

                break;
            default:
                break;
            }
        }

        /* probabilistic-based expansion for last FREIDEL pixel */
        pexp = ran1(seed);
        if ( pexp <= 0.3295 ) {
            newact = extfreidel(xexp, yexp, zexp);
        }
    }
    /* if new location is C4AF execute conversion */
    /* to freidel's salt (including necessary volumetric expansion) */
    else if ( check == C4AF ) {
        mic [ xnew ] [ ynew ] [ znew ] = FREIDEL;
        count [ FREIDEL ] += 1;
        count [ C4AF ] -= 1;

        /* determine if CACL2 should be converted to FREIDEL */
        /* 0.4033 unit of CaCl2 requires 1 unit of C4AF */
        /* and should form 2.3176 units of FREIDEL */
        /* Also 0.6412 units of CH and 1.3522 units of FH3 */
        /* per unit of CACL2 */
        pexp = ran1(seed);
        nexp = 1;
        if ( pexp <= 0.4033 ) {
            mic [ xcur ] [ ycur ] [ zcur ] = FREIDEL;
            count [ FREIDEL ] += 1;
            count [ DIFFCACL2 ] -= 1;
            nexp = 0;
            pext = ran1(seed);
            /* Addition of extra CH */
            if ( pext < 0.6412 ) {
                extch();
            }

            pext = ran1(seed);
            /* Addition of extra FH3 */
            if ( pext < 0.3522 ) {
                extfh3(xnew, ynew, znew);
            }

            extfh3(xnew, ynew, znew);
        } else {
            nexp = 1;
            keep = 1;
        }

        /* create extra freidel's salt pixels to maintain volume stoichiometry */
        /* xexp, yexp and zexp hold coordinates of most recently added FREIDEL */
        xexp = xcur;
        yexp = ycur;
        zexp = zcur;
        for ( iexp = 1; iexp <= nexp; iexp++ ) {
            newact = extfreidel(xexp, yexp, zexp);
            /* update xexp, yexp and zexp as needed */
            switch ( newact ) {
            case 1:
                xexp -= 1;
                if ( xexp < 0 ) {
                    xexp = ( SYSIZEM1 );
                }

                break;
            case 2:
                xexp += 1;
                if ( xexp >= SYSIZE ) {
                    xexp = 0;
                }

                break;
            case 3:
                yexp -= 1;
                if ( yexp < 0 ) {
                    yexp = ( SYSIZEM1 );
                }

                break;
            case 4:
                yexp += 1;
                if ( yexp >= SYSIZE ) {
                    yexp = 0;
                }

                break;
            case 5:
                zexp -= 1;
                if ( zexp < 0 ) {
                    zexp = ( SYSIZEM1 );
                }

                break;
            case 6:
                zexp += 1;
                if ( zexp >= SYSIZE ) {
                    zexp = 0;
                }

                break;
            default:
                break;
            }
        }

        /* probabilistic-based expansion for last FREIDEL pixel */
        pexp = ran1(seed);
        if ( pexp <= 0.3176 ) {
            newact = extfreidel(xexp, yexp, zexp);
        }

        action = 0;
    }

    /* if last diffusion step and no reaction, convert back to */
    /* solid CaCl2 */
    if ( ( action != 0 ) && ( finalstep == 1 ) ) {
        action = 0;
        count [ DIFFCACL2 ] -= 1;
        count [ CACL2 ] += 1;
        mic [ xcur ] [ ycur ] [ zcur ] = CACL2;
    }

    if ( action != 0 ) {
        /* if diffusion is possible, execute it */
        if ( check == POROSITY ) {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            mic [ xnew ] [ ynew ] [ znew ] = DIFFCACL2;
        } else {
            /* indicate that diffusing CACL2 remained at */
            /* original location */
            action = 7;
        }
    }

    if ( keep == 1 ) {
        action = 7;
    }

    return ( action );
}

/* routine to move a diffusing CAS2 species */
/* from current location (xcur,ycur,zcur) */
/* Returns action flag indicating response taken */
/* Called by hydrate */
/* Calls moveone, extstrat, extch, and extfh3 */
int CemhydMatStatus :: movecas2(int xcur, int ycur, int zcur, int finalstep)
{
    int check, xnew, ynew, znew, action, nexp, iexp;
    int xexp, yexp, zexp, newact, sumold, sumgarb, keep;
    float pexp, pext;

    sumold = 1;
    keep = 0;

    /* First be sure that a diffusing CAS2 species is located at xcur,ycur,zcur */
    /* if not, return to calling routine */
    if ( mic [ xcur ] [ ycur ] [ zcur ] != DIFFCAS2 ) {
        action = 0;
        return ( action );
    }

    /* Determine new coordinates (periodic boundaries are used) */
    xnew = xcur;
    ynew = ycur;
    znew = zcur;
    action = 0;
    sumgarb = moveone(& xnew, & ynew, & znew, & action, sumold);
    if ( action == 0 ) {
        printf("Error in value of action in movecas2 \n");
    }

    check = mic [ xnew ] [ ynew ] [ znew ];

    /* if new location is C3A or diffusing C3A, execute conversion */
    /* to stratlingite (including necessary volumetric expansion) */
    if ( ( check == C3A ) || ( check == DIFFC3A ) || ( check == DIFFC4A ) ) {
        /* Convert diffusing CAS2 to a stratlingite pixel */
        action = 0;
        mic [ xcur ] [ ycur ] [ zcur ] = STRAT;
        count [ STRAT ] += 1;
        count [ DIFFCAS2 ] -= 1;

        /* determine if diffusing or solid C3A should be converted to STRAT*/
        /* 1 unit of CAS2 requires 0.886 units of C3A */
        /* and should form 4.286 units of STRAT */
        pexp = ran1(seed);
        nexp = 3;
        if ( pexp <= 0.886 ) {
            mic [ xnew ] [ ynew ] [ znew ] = STRAT;
            count [ STRAT ] += 1;
            count [ check ] -= 1;
            nexp = 2;
        }

        /* create extra stratlingite pixels to maintain volume stoichiometry */
        /* xexp, yexp, and zexp hold coordinates of most recently added STRAT */
        xexp = xcur;
        yexp = ycur;
        zexp = zcur;
        for ( iexp = 1; iexp <= nexp; iexp++ ) {
            newact = extstrat(xexp, yexp, zexp);
            /* update xexp, yexp and zexp as needed */
            switch ( newact ) {
            case 1:
                xexp -= 1;
                if ( xexp < 0 ) {
                    xexp = ( SYSIZEM1 );
                }

                break;
            case 2:
                xexp += 1;
                if ( xexp >= SYSIZE ) {
                    xexp = 0;
                }

                break;
            case 3:
                yexp -= 1;
                if ( yexp < 0 ) {
                    yexp = ( SYSIZEM1 );
                }

                break;
            case 4:
                yexp += 1;
                if ( yexp >= SYSIZE ) {
                    yexp = 0;
                }

                break;
            case 5:
                zexp -= 1;
                if ( zexp < 0 ) {
                    zexp = ( SYSIZEM1 );
                }

                break;
            case 6:
                zexp += 1;
                if ( zexp >= SYSIZE ) {
                    zexp = 0;
                }

                break;
            default:
                break;
            }
        }

        /* probabilistic-based expansion for last STRAT pixel */
        pexp = ran1(seed);
        if ( pexp <= 0.286 ) {
            newact = extstrat(xexp, yexp, zexp);
        }
    }
    /* if new location is C4AF execute conversion */
    /* to stratlingite (including necessary volumetric expansion) */
    else if ( check == C4AF ) {
        mic [ xnew ] [ ynew ] [ znew ] = STRAT;
        count [ STRAT ] += 1;
        count [ C4AF ] -= 1;

        /* determine if CAS2 should be converted to STRAT */
        /* 0.786 units of CAS2 requires 1 unit of C4AF */
        /* and should form 3.37 units of STRAT */
        /* Also 0.2586 units of CH and 0.5453 units of FH3 */
        /* per unit of C4AF */
        pexp = ran1(seed);
        nexp = 2;
        if ( pexp <= 0.786 ) {
            mic [ xcur ] [ ycur ] [ zcur ] = STRAT;
            count [ STRAT ] += 1;
            count [ DIFFCAS2 ] -= 1;
            nexp = 1;
            pext = ran1(seed);
            /* Addition of extra CH */
            /* 0.329= 0.2586/0.786 */
            if ( pext < 0.329 ) {
                extch();
            }

            pext = ran1(seed);
            /* Addition of extra FH3 */
            /* 0.6938= 0.5453/0.786 */
            if ( pext < 0.6938 ) {
                extfh3(xnew, ynew, znew);
            }
        } else {
            nexp = 2;
            keep = 1;
        }

        /* create extra stratlingite pixels to maintain volume stoichiometry */
        /* xexp, yexp and zexp hold coordinates of most recently added STRAT */
        xexp = xcur;
        yexp = ycur;
        zexp = zcur;
        for ( iexp = 1; iexp <= nexp; iexp++ ) {
            newact = extstrat(xexp, yexp, zexp);
            /* update xexp, yexp and zexp as needed */
            switch ( newact ) {
            case 1:
                xexp -= 1;
                if ( xexp < 0 ) {
                    xexp = ( SYSIZEM1 );
                }

                break;
            case 2:
                xexp += 1;
                if ( xexp >= SYSIZE ) {
                    xexp = 0;
                }

                break;
            case 3:
                yexp -= 1;
                if ( yexp < 0 ) {
                    yexp = ( SYSIZEM1 );
                }

                break;
            case 4:
                yexp += 1;
                if ( yexp >= SYSIZE ) {
                    yexp = 0;
                }

                break;
            case 5:
                zexp -= 1;
                if ( zexp < 0 ) {
                    zexp = ( SYSIZEM1 );
                }

                break;
            case 6:
                zexp += 1;
                if ( zexp >= SYSIZE ) {
                    zexp = 0;
                }

                break;
            default:
                break;
            }
        }

        /* probabilistic-based expansion for last STRAT pixel */
        pexp = ran1(seed);
        if ( pexp <= 0.37 ) {
            newact = extstrat(xexp, yexp, zexp);
        }

        action = 0;
    }

    /* if last diffusion step and no reaction, convert back to */
    /* solid CAS2 */
    if ( ( action != 0 ) && ( finalstep == 1 ) ) {
        action = 0;
        count [ DIFFCAS2 ] -= 1;
        count [ CAS2 ] += 1;
        mic [ xcur ] [ ycur ] [ zcur ] = CAS2;
    }

    if ( action != 0 ) {
        /* if diffusion is possible, execute it */
        if ( check == POROSITY ) {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            mic [ xnew ] [ ynew ] [ znew ] = DIFFCAS2;
        } else {
            /* indicate that diffusing CAS2 remained at */
            /* original location */
            action = 7;
        }
    }

    if ( keep == 1 ) {
        action = 7;
    }

    return ( action );
}

/* routine to move a diffusing AS species */
/* from current location (xcur,ycur,zcur) */
/* Returns action flag indicating response taken */
/* Called by hydrate */
/* Calls moveone, extstrat */
int CemhydMatStatus :: moveas(int xcur, int ycur, int zcur, int finalstep)
{
    int check, xnew, ynew, znew, action, nexp, iexp;
    int xexp, yexp, zexp, newact, sumold, sumgarb, keep;
    float pexp;

    sumold = 1;
    keep = 0;

    /* First be sure that a diffusing AS species is located at xcur,ycur,zcur */
    /* if not, return to calling routine */
    if ( mic [ xcur ] [ ycur ] [ zcur ] != DIFFAS ) {
        action = 0;
        return ( action );
    }

    /* Determine new coordinates (periodic boundaries are used) */
    xnew = xcur;
    ynew = ycur;
    znew = zcur;
    action = 0;
    sumgarb = moveone(& xnew, & ynew, & znew, & action, sumold);
    if ( action == 0 ) {
        printf("Error in value of action in moveas \n");
    }

    check = mic [ xnew ] [ ynew ] [ znew ];

    /* if new location is CH or diffusing CH, execute conversion */
    /* to stratlingite (including necessary volumetric expansion) */
    if ( ( check == CH ) || ( check == DIFFCH ) ) {
        /* Convert diffusing CH or CH to a stratlingite pixel */
        action = 0;
        mic [ xnew ] [ ynew ] [ znew ] = STRAT;
        count [ STRAT ] += 1;
        count [ check ] -= 1;

        /* determine if diffusing AS should be converted to STRAT */
        /* 0.7538 unit of AS requires 1 unit of CH */
        /* and should form 3.26 units of STRAT */
        pexp = ran1(seed);
        nexp = 2;
        if ( pexp <= 0.7538 ) {
            mic [ xcur ] [ ycur ] [ zcur ] = STRAT;
            count [ STRAT ] += 1;
            count [ DIFFAS ] -= 1;
            nexp = 1;
        } else {
            keep = 1;
            nexp = 2;
        }

        /* create extra stratlingite pixels to maintain volume stoichiometry */
        /* xexp, yexp, and zexp hold coordinates of most recently added STRAT */
        xexp = xcur;
        yexp = ycur;
        zexp = zcur;
        for ( iexp = 1; iexp <= nexp; iexp++ ) {
            newact = extstrat(xexp, yexp, zexp);
            /* update xexp, yexp and zexp as needed */
            switch ( newact ) {
            case 1:
                xexp -= 1;
                if ( xexp < 0 ) {
                    xexp = ( SYSIZEM1 );
                }

                break;
            case 2:
                xexp += 1;
                if ( xexp >= SYSIZE ) {
                    xexp = 0;
                }

                break;
            case 3:
                yexp -= 1;
                if ( yexp < 0 ) {
                    yexp = ( SYSIZEM1 );
                }

                break;
            case 4:
                yexp += 1;
                if ( yexp >= SYSIZE ) {
                    yexp = 0;
                }

                break;
            case 5:
                zexp -= 1;
                if ( zexp < 0 ) {
                    zexp = ( SYSIZEM1 );
                }

                break;
            case 6:
                zexp += 1;
                if ( zexp >= SYSIZE ) {
                    zexp = 0;
                }

                break;
            default:
                break;
            }
        }

        /* probabilistic-based expansion for last stratlingite pixel */
        pexp = ran1(seed);
        if ( pexp <= 0.326 ) {
            newact = extstrat(xexp, yexp, zexp);
        }
    }

    /* if last diffusion step and no reaction, convert back to */
    /* solid ASG */
    if ( ( action != 0 ) && ( finalstep == 1 ) ) {
        action = 0;
        count [ DIFFAS ] -= 1;
        count [ ASG ] += 1;
        mic [ xcur ] [ ycur ] [ zcur ] = ASG;
    }

    if ( action != 0 ) {
        /* if diffusion is possible, execute it */
        if ( check == POROSITY ) {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            mic [ xnew ] [ ynew ] [ znew ] = DIFFAS;
        } else {
            /* indicate that diffusing AS remained at */
            /* original location */
            action = 7;
        }
    }

    if ( keep == 1 ) {
        action = 7;
    }

    return ( action );
}

/* routine to move a diffusing CACO3 species */
/* from current location (xcur,ycur,zcur) */
/* Returns action flag indicating response taken */
/* Called by hydrate */
/* Calls moveone, extettr */
int CemhydMatStatus :: movecaco3(int xcur, int ycur, int zcur, int finalstep)
{
    int check, xnew, ynew, znew, action;
    int xexp, yexp, zexp, newact, sumold, sumgarb, keep;
    float pexp;

    sumold = 1;
    keep = 0;

    /* First be sure that a diffusing CACO3 species is located at xcur,ycur,zcur */
    /* if not, return to calling routine */
    if ( mic [ xcur ] [ ycur ] [ zcur ] != DIFFCACO3 ) {
        action = 0;
        return ( action );
    }

    /* Determine new coordinates (periodic boundaries are used) */
    xnew = xcur;
    ynew = ycur;
    znew = zcur;
    action = 0;
    sumgarb = moveone(& xnew, & ynew, & znew, & action, sumold);
    if ( action == 0 ) {
        printf("Error in value of action in moveas \n");
    }

    check = mic [ xnew ] [ ynew ] [ znew ];

    /* if new location is AFM execute conversion */
    /* to carboaluminate and ettringite (including necessary */
    /* volumetric expansion) */
    if ( check == AFM ) {
        /* Convert AFM to a carboaluminate or ettringite pixel */
        action = 0;
        pexp = ran1(seed);
        if ( pexp <= 0.479192 ) {
            mic [ xnew ] [ ynew ] [ znew ] = AFMC;
            count [ AFMC ] += 1;
        } else {
            mic [ xnew ] [ ynew ] [ znew ] = ETTR;
            count [ ETTR ] += 1;
        }

        count [ check ] -= 1;

        /* determine if diffusing CACO3 should be converted to AFMC */
        /* 0.078658 unit of AS requires 1 unit of AFM */
        /* and should form 0.55785 units of AFMC */
        pexp = ran1(seed);
        if ( pexp <= 0.078658 ) {
            mic [ xcur ] [ ycur ] [ zcur ] = AFMC;
            count [ AFMC ] += 1;
            count [ DIFFCACO3 ] -= 1;
        } else {
            keep = 1;
        }

        /* create extra ettringite pixels to maintain volume stoichiometry */
        /* xexp, yexp, and zexp hold coordinates of most recently added ETTR */
        xexp = xnew;
        yexp = ynew;
        zexp = znew;

        /* probabilistic-based expansion for new ettringite pixel */
        pexp = ran1(seed);
        if ( pexp <= 0.26194 ) {
            newact = extettr(xexp, yexp, zexp, 0);
        }
    }

    /* if last diffusion step and no reaction, convert back to */
    /* solid CACO3 */
    if ( ( action != 0 ) && ( finalstep == 1 ) ) {
        action = 0;
        count [ DIFFCACO3 ] -= 1;
        count [ CACO3 ] += 1;
        mic [ xcur ] [ ycur ] [ zcur ] = CACO3;
    }

    if ( action != 0 ) {
        /* if diffusion is possible, execute it */
        if ( check == POROSITY ) {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            mic [ xnew ] [ ynew ] [ znew ] = DIFFCACO3;
        } else {
            /* indicate that diffusing CACO3 remained at */
            /* original location */
            action = 7;
        }
    }

    if ( keep == 1 ) {
        action = 7;
    }

    return ( action );
}

/* routine to add extra AFm phase when diffusing ettringite reacts */
/* with C3A (diffusing or solid) at location (xpres,ypres,zpres) */
/* Called by moveettr and movec3a */
/* Calls moveone and edgecnt */
void CemhydMatStatus :: extafm(int xpres, int ypres, int zpres)
{
    int check, sump, xchr, ychr, zchr, fchr, i1, newact, numnear;
    long int tries;

    /* first try 6 neighboring locations until      */
    /*    a) successful                */
    /*    b) all 6 sites are tried or             */
    /*    c) 100 tries are made           */
    fchr = 0;
    sump = 1;
    for ( i1 = 1; ( ( i1 <= 100 ) && ( fchr == 0 ) && ( sump != 30030 ) ); i1++ ) {
        /* determine location of neighbor (using periodic boundaries) */
        xchr = xpres;
        ychr = ypres;
        zchr = zpres;
        newact = 0;
        sump *= moveone(& xchr, & ychr, & zchr, & newact, sump);
        if ( newact == 0 ) {
            printf("Error in value of newact in extafm \n");
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if neighbor is porosity, locate the AFm phase there */
        if ( check == POROSITY ) {
            mic [ xchr ] [ ychr ] [ zchr ] = AFM;
            count [ AFM ] += 1;
            count [ POROSITY ] -= 1;
            fchr = 1;
        }
    }

    /* if no neighbor available, locate AFm phase at random location */
    /* in pore space */
    tries = 0;
    while ( fchr == 0 ) {
        tries += 1;
        /* generate a random location in the 3-D system */
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if location is porosity, locate the extra AFm there */
        if ( check == POROSITY ) {
            numnear = edgecnt(xchr, ychr, zchr, AFM, C3A, C4AF);
            /* Be sure that at least one neighboring pixel is */
            /* Afm phase, C3A, or C4AF */
            if ( ( tries > 5000 ) || ( numnear < 26 ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = AFM;
                count [ AFM ] += 1;
                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }
}

/* routine to move a diffusing ettringite species */
/* currently located at (xcur,ycur,zcur) */
/* Called by hydrate */
/* Calls moveone, extch, extfh3, and extafm */
int CemhydMatStatus :: moveettr(int xcur, int ycur, int zcur, int finalstep)
{
    int check, xnew, ynew, znew, action;
    int sumold, sumgarb;
    float pexp, pafm, pgrow;

    /* First be sure a diffusing ettringite species is located at xcur,ycur,zcur */
    /* if not, return to calling routine */
    if ( mic [ xcur ] [ ycur ] [ zcur ] != DIFFETTR ) {
        action = 0;
        return ( action );
    }

    /* Determine new coordinates (periodic boundaries are used) */
    xnew = xcur;
    ynew = ycur;
    znew = zcur;
    action = 0;
    sumold = 1;
    sumgarb = moveone(& xnew, & ynew, & znew, & action, sumold);
    if ( action == 0 ) {
        printf("Error in value of action in moveettr \n");
    }

    check = mic [ xnew ] [ ynew ] [ znew ];

    /* if new location is C4AF, execute conversion */
    /* to AFM phase (including necessary volumetric expansion) */
    if ( check == C4AF ) {
        /* Convert diffusing ettringite to AFM phase */
        mic [ xcur ] [ ycur ] [ zcur ] = AFM;
        count [ AFM ] += 1;
        count [ DIFFETTR ] -= 1;

        /* determine if C4AF should be converted to Afm */
        /* or FH3- 1 unit of ettringite requires 0.348 units */
        /* of C4AF to form 1.278 units of Afm, */
        /* 0.0901 units of CH and 0.1899 units of FH3 */
        pexp = ran1(seed);

        if ( pexp <= 0.278 ) {
            mic [ xnew ] [ ynew ] [ znew ] = AFM;
            count [ AFM ] += 1;
            count [ C4AF ] -= 1;
            pafm = ran1(seed);
            /* 0.3241= 0.0901/0.278 */
            if ( pafm <= 0.3241 ) {
                extch();
            }

            pafm = ran1(seed);
            /* 0.4313= ((.1899-(.348-.278))/.278)   */
            if ( pafm <= 0.4313 ) {
                extfh3(xnew, ynew, znew);
            }
        } else if ( pexp <= 0.348 ) {
            mic [ xnew ] [ ynew ] [ znew ] = FH3;
            count [ FH3 ] += 1;
            count [ C4AF ] -= 1;
        }

        action = 0;
    }
    /* if new location is C3A or diffusing C3A, execute conversion */
    /* to AFM phase (including necessary volumetric expansion) */
    else if ( ( check == C3A ) || ( check == DIFFC3A ) ) {
        /* Convert diffusing ettringite to AFM phase */
        action = 0;
        mic [ xcur ] [ ycur ] [ zcur ] = AFM;
        count [ DIFFETTR ] -= 1;
        count [ AFM ] += 1;
        count [ check ] -= 1;

        /* determine if C3A should be converted to AFm */
        /* 1 unit of ettringite requires 0.2424 units of C3A */
        /* and should form 1.278 units of AFm phase */
        pexp = ran1(seed);
        if ( pexp <= 0.2424 ) {
            mic [ xnew ] [ ynew ] [ znew ] = AFM;
            count [ AFM ] += 1;
            pafm = ( -0.1 );
        } else {
            /* maybe someday, use a new FIXEDC3A here */
            /* so it won't dissolve later */
            if ( check == C3A ) {
                mic [ xnew ] [ ynew ] [ znew ] = C3A;
                count [ C3A ] += 1;
            } else {
                count [ DIFFC3A ] += 1;
                mic [ xnew ] [ ynew ] [ znew ] = DIFFC3A;
            }

            /*                      pafm=(0.278-0.2424)/(1.0-0.2424);  */
            pafm = 0.04699;
        }

        /* probabilistic-based expansion for new AFm phase pixel */
        pexp = ran1(seed);
        if ( pexp <= pafm ) {
            extafm(xcur, ycur, zcur);
        }
    }
    /* Check for conversion back to solid ettringite */
    else if ( check == ETTR ) {
        pgrow = ran1(seed);
        if ( pgrow <= ETTRGROW ) {
            mic [ xcur ] [ ycur ] [ zcur ] = ETTR;
            count [ ETTR ] += 1;
            action = 0;
            count [ DIFFETTR ] -= 1;
        }
    }

    /* if last diffusion step and no reaction, convert back to */
    /* solid ettringite */
    if ( ( action != 0 ) && ( finalstep == 1 ) ) {
        action = 0;
        count [ DIFFETTR ] -= 1;
        count [ ETTR ] += 1;
        mic [ xcur ] [ ycur ] [ zcur ] = ETTR;
    }

    if ( action != 0 ) {
        /* if diffusion is possible, execute it */
        if ( check == POROSITY ) {
            mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
            mic [ xnew ] [ ynew ] [ znew ] = DIFFETTR;
        } else {
            /* indicate that diffusing ettringite remained at */
            /* original location */
            action = 7;
        }
    }

    return ( action );
}

/* routine to add extra pozzolanic CSH when CH reacts at */
/* pozzolanic surface (e.g. silica fume) located at (xpres,ypres,zpres) */
/* Called by movech */
/* Calls moveone and edgecnt */
void CemhydMatStatus :: extpozz(int xpres, int ypres, int zpres)
{
    int check, sump, xchr, ychr, zchr, fchr, i1, action, numnear;
    long int tries;

    /* first try 6 neighboring locations until      */
    /*    a) successful                */
    /*    b) all 6 sites are tried or             */
    /*    c) 100 tries are made           */
    fchr = 0;
    sump = 1;
    for ( i1 = 1; ( ( i1 <= 100 ) && ( fchr == 0 ) && ( sump != 30030 ) ); i1++ ) {
        /* determine location of neighbor (using periodic boundaries) */
        xchr = xpres;
        ychr = ypres;
        zchr = zpres;
        action = 0;
        sump *= moveone(& xchr, & ychr, & zchr, & action, sump);
        if ( action == 0 ) {
            printf("Error in value of action in extpozz \n");
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if neighbor is porosity, locate the pozzolanic CSH there */
        if ( check == POROSITY ) {
            mic [ xchr ] [ ychr ] [ zchr ] = POZZCSH;
            count [ POZZCSH ] += 1;
            count [ POROSITY ] -= 1;
            fchr = 1;
        }
    }

    /* if no neighbor available, locate pozzolanic CSH at random location */
    /* in pore space */
    tries = 0;
    while ( fchr == 0 ) {
        tries += 1;
        /* generate a random location in the 3-D system */
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];
        /* if location is porosity, locate the extra pozzolanic CSH there */
        if ( check == POROSITY ) {
            numnear = edgecnt(xchr, ychr, zchr, POZZ, CSH, POZZCSH);
            /* Be sure that one neighboring species is CSH or */
            /* pozzolanic material */
            if ( ( tries > 5000 ) || ( numnear < 26 ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = POZZCSH;
                count [ POZZCSH ] += 1;
                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }
}

/* routine to move a diffusing FH3 species */
/* from location (xcur,ycur,zcur) with nucleation probability nucprob */
/* Called by hydrate */
/* Calls moveone */
int CemhydMatStatus :: movefh3(int xcur, int ycur, int zcur, int finalstep, float nucprob)
{
    int check, xnew, ynew, znew, action, sumold, sumgarb;
    float pgen;

    /* first check for nucleation */
    pgen = ran1(seed);

    if ( ( nucprob >= pgen ) || ( finalstep == 1 ) ) {
        action = 0;
        mic [ xcur ] [ ycur ] [ zcur ] = FH3;
        count [ FH3 ] += 1;
        count [ DIFFFH3 ] -= 1;
    } else {
        /* determine new location (using periodic boundaries) */
        xnew = xcur;
        ynew = ycur;
        znew = zcur;
        action = 0;
        sumold = 1;
        sumgarb = moveone(& xnew, & ynew, & znew, & action, sumold);
        if ( action == 0 ) {
            printf("Error in value of action in movefh3 \n");
        }

        check = mic [ xnew ] [ ynew ] [ znew ];

        /* check for growth of FH3 crystal */
        if ( check == FH3 ) {
            mic [ xcur ] [ ycur ] [ zcur ] = FH3;
            count [ FH3 ] += 1;
            count [ DIFFFH3 ] -= 1;
            action = 0;
        }

        if ( action != 0 ) {
            /* if diffusion is possible, execute it */
            if ( check == POROSITY ) {
                mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
                mic [ xnew ] [ ynew ] [ znew ] = DIFFFH3;
            } else {
                /* indicate that diffusing FH3 species */
                /* remained at original location */
                action = 7;
            }
        }
    }

    return ( action );
}

/* routine to move a diffusing CH species */
/* from location (xcur,ycur,zcur) with nucleation probability nucprob */
/* Called by hydrate */
/* Calls moveone and extpozz */
int CemhydMatStatus :: movech(int xcur, int ycur, int zcur, int finalstep, float nucprob)
{
    int check, xnew, ynew, znew, action, sumgarb, sumold;
    float pexp, pgen, pfix;

    /* first check for nucleation */
    pgen = ran1(seed);
    if ( ( nucprob >= pgen ) || ( finalstep == 1 ) ) {
        action = 0;
        mic [ xcur ] [ ycur ] [ zcur ] = CH;
        count [ DIFFCH ] -= 1;
        count [ CH ] += 1;
    } else {
        /* determine new location (using periodic boundaries) */
        xnew = xcur;
        ynew = ycur;
        znew = zcur;
        action = 0;
        sumold = 1;
        sumgarb = moveone(& xnew, & ynew, & znew, & action, sumold);
        if ( action == 0 ) {
            printf("Error in value of action in movech \n");
        }

        check = mic [ xnew ] [ ynew ] [ znew ];

        /* check for growth of CH crystal */
        if ( ( check == CH ) && ( pgen <= CHGROW ) ) {
            mic [ xcur ] [ ycur ] [ zcur ] = CH;
            count [ DIFFCH ] -= 1;
            count [ CH ] += 1;
            action = 0;
        }
        /* check for growth of CH crystal on aggregate or CaCO3 surface */
        /* re suggestion of Sidney Diamond */
        else if ( ( ( check == INERTAGG ) || ( check == CACO3 ) || ( check == INERT ) ) && ( pgen <= CHGROWAGG ) && ( chflag == 1 ) ) {
            mic [ xcur ] [ ycur ] [ zcur ] = CH;
            count [ DIFFCH ] -= 1;
            count [ CH ] += 1;
            action = 0;
        }
        /* check for pozzolanic reaction */
        /* 36.41 units CH can react with 27 units of S */
        else if ( ( pgen <= ppozz ) && ( check == POZZ ) && ( npr <= ( int ) ( ( float ) nfill * 1.35 ) ) ) {
            action = 0;
            mic [ xcur ] [ ycur ] [ zcur ] = POZZCSH;
            count [ POZZCSH ] += 1;
            /* update counter of number of diffusing CH */
            /* which have reacted pozzolanically */
            npr += 1;
            count [ DIFFCH ] -= 1;
            /* Convert pozzolan to pozzolanic CSH as needed */
            pfix = ran1(seed);
            if ( pfix <= ( 1. / 1.35 ) ) {
                mic [ xnew ] [ ynew ] [ znew ] = POZZCSH;
                count [ POZZ ] -= 1;
                count [ POZZCSH ] += 1;
            }

            /* allow for extra pozzolanic CSH as needed */
            pexp = ran1(seed);
            /* should form 101.81 units of pozzolanic CSH for */
            /* each 36.41 units of CH and 27 units of S */
            /* 1.05466=(101.81-36.41-27)/36.41 */
            extpozz(xcur, ycur, zcur);
            if ( pexp <= 0.05466 ) {
                extpozz(xcur, ycur, zcur);
            }
        } else if ( check == DIFFAS ) {
            action = 0;
            mic [ xcur ] [ ycur ] [ zcur ] = STRAT;
            count [ STRAT ] += 1;
            /* update counter of number of diffusing CH */
            /* which have reacted to form stratlingite */
            nasr += 1;
            count [ DIFFCH ] -= 1;
            /* Convert DIFFAS to STRAT as needed */
            pfix = ran1(seed);
            if ( pfix <= 0.7538 ) {
                mic [ xnew ] [ ynew ] [ znew ] = STRAT;
                count [ STRAT ] += 1;
                count [ DIFFAS ] -= 1;
            }

            /* allow for extra stratlingite as needed */
            /* 1.5035=(215.63-66.2-49.9)/66.2 */
            extstrat(xcur, ycur, zcur);
            pexp = ran1(seed);
            if ( pexp <= 0.5035 ) {
                extstrat(xcur, ycur, zcur);
            }
        }

        if ( action != 0 ) {
            /* if diffusion is possible, execute it */
            if ( check == POROSITY ) {
                mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
                mic [ xnew ] [ ynew ] [ znew ] = DIFFCH;
            } else {
                /* indicate that diffusing CH species */
                /* remained at original location */
                action = 7;
            }
        }
    }

    return ( action );
}

/* routine to add extra C3AH6 when diffusing C3A nucleates or reacts at */
/* C3AH6 surface at location (xpres,ypres,zpres) */
/* Called by movec3a */
/* Calls moveone and edgecnt */
void CemhydMatStatus :: extc3ah6(int xpres, int ypres, int zpres)
{
    int check, sump, xchr, ychr, zchr, fchr, i1, action, numnear;
    long int tries;

    /* First try 6 neighboring locations until      */
    /*  a) successful                */
    /*    b) all 6 sites are tried or             */
    /*    c) 100 random attempts are made     */
    fchr = 0;
    sump = 1;
    for ( i1 = 1; ( ( i1 <= 100 ) && ( fchr == 0 ) && ( sump != 30030 ) ); i1++ ) {
        /* determine new coordinates (using periodic boundaries) */
        xchr = xpres;
        ychr = ypres;
        zchr = zpres;
        action = 0;
        sump *= moveone(& xchr, & ychr, & zchr, & action, sump);
        if ( action == 0 ) {
            printf("Error in action value in extc3ah6 \n");
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        /* if neighbor is pore space, convert it to C3AH6 */
        if ( check == POROSITY ) {
            mic [ xchr ] [ ychr ] [ zchr ] = C3AH6;
            count [ C3AH6 ] += 1;
            count [ POROSITY ] -= 1;
            fchr = 1;
        }
    }

    /* if unsuccessful, add C3AH6 at random location in pore space */
    tries = 0;
    while ( fchr == 0 ) {
        tries += 1;
        xchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        ychr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        zchr = ( int ) ( ( float ) SYSIZE * ran1(seed) );
        if ( xchr >= SYSIZE ) {
            xchr = 0;
        }

        if ( ychr >= SYSIZE ) {
            ychr = 0;
        }

        if ( zchr >= SYSIZE ) {
            zchr = 0;
        }

        check = mic [ xchr ] [ ychr ] [ zchr ];

        if ( check == POROSITY ) {
            numnear = edgecnt(xchr, ychr, zchr, C3AH6, C3A, C3AH6);
            /* Be sure that new C3AH6 is in contact with */
            /* at least one C3AH6 or C3A */
            if ( ( tries > 5000 ) || ( numnear < 26 ) ) {
                mic [ xchr ] [ ychr ] [ zchr ] = C3AH6;
                count [ C3AH6 ] += 1;
                count [ POROSITY ] -= 1;
                fchr = 1;
            }
        }
    }
}

/* routine to move a diffusing C3A species */
/* from location (xcur,ycur,zcur) with nucleation probability of nucprob */
/* Called by hydrate */
/* Calls extc3ah6, moveone, extettr, and extafm */
int CemhydMatStatus :: movec3a(int xcur, int ycur, int zcur, int finalstep, float nucprob)
{
    int check, xnew, ynew, znew, action, sumgarb, sumold;
    int xexp, yexp, zexp, nexp, iexp, newact;
    float pgen, pexp, pafm, pgrow, p2diff;

    /* First be sure that a diffusing C3A species is at (xcur,ycur,zcur) */
    if ( mic [ xcur ] [ ycur ] [ zcur ] != DIFFC3A ) {
        action = 0;
        return ( action );
    }

    /* Check for nucleation into solid C3AH6 */
    pgen = ran1(seed);
    p2diff = ran1(seed);

    if ( ( nucprob >= pgen ) || ( finalstep == 1 ) ) {
        action = 0;
        mic [ xcur ] [ ycur ] [ zcur ] = C3AH6;
        count [ C3AH6 ] += 1;
        /* decrement count of diffusing C3A species */
        count [ DIFFC3A ] -= 1;
        /* allow for probabilistic-based expansion of C3AH6 */
        /* crystal to account for volume stoichiometry */
        pexp = ran1(seed);
        if ( pexp <= 0.69 ) {
            extc3ah6(xcur, ycur, zcur);
        }
    } else {
        /* determine new coordinates (using periodic boundaries) */
        xnew = xcur;
        ynew = ycur;
        znew = zcur;
        action = 0;
        sumold = 1;
        sumgarb = moveone(& xnew, & ynew, & znew, & action, sumold);
        if ( action == 0 ) {
            printf("Error in value of action in movec3a \n");
        }

        check = mic [ xnew ] [ ynew ] [ znew ];

        /* check for growth of C3AH6 crystal */
        if ( check == C3AH6 ) {
            pgrow = ran1(seed);
            /* Try to slow down growth of C3AH6 crystals to */
            /* promote ettringite and Afm formation */
            if ( pgrow <= C3AH6GROW ) {
                mic [ xcur ] [ ycur ] [ zcur ] = C3AH6;
                count [ C3AH6 ] += 1;
                count [ DIFFC3A ] -= 1;
                action = 0;
                /* allow for probabilistic-based expansion of C3AH6 */
                /* crystal to account for volume stoichiometry */
                pexp = ran1(seed);
                if ( pexp <= 0.69 ) {
                    extc3ah6(xcur, ycur, zcur);
                }
            }
        }
        /* examine reaction with diffusing gypsum to form ettringite */
        /* Only allow reaction with diffusing gypsum */
        else if ( ( check == DIFFGYP ) && ( p2diff < C3AGYP ) ) {
            /* convert diffusing gypsum to ettringite */
            mic [ xnew ] [ ynew ] [ znew ] = ETTR;
            count [ ETTR ] += 1;
            /* decrement counts of diffusing gypsum */
            count [ DIFFGYP ] -= 1;
            action = 0;

            /* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
            pexp = ran1(seed);
            nexp = 2;
            if ( pexp <= 0.40 ) {
                mic [ xcur ] [ ycur ] [ zcur ] = ETTR;
                count [ ETTR ] += 1;
                count [ DIFFC3A ] -= 1;
                nexp = 1;
            } else {
                /* indicate that diffusing species remains in current location */
                action = 7;
                nexp = 2;
            }

            /* Perform expansion that occurs when ettringite is formed */
            /* xexp, yexp and zexp are the coordinates of the last ettringite */
            /* pixel to be added */
            xexp = xnew;
            yexp = ynew;
            zexp = znew;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extettr(xexp, yexp, zexp, 0);
                /* update xexp, yexp and zexp */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last ettringite pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.30 ) {
                newact = extettr(xexp, yexp, zexp, 0);
            }
        }
        /* examine reaction with diffusing hemihydrate to form ettringite */
        /* Only allow reaction with diffusing hemihydrate */
        else if ( ( check == DIFFHEM ) && ( p2diff < C3AGYP ) ) {
            /* convert diffusing hemihydrate to ettringite */
            mic [ xnew ] [ ynew ] [ znew ] = ETTR;
            count [ ETTR ] += 1;
            /* decrement counts of diffusing hemihydrate */
            count [ DIFFHEM ] -= 1;
            action = 0;

            /* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
            pexp = ran1(seed);
            nexp = 3;
            if ( pexp <= 0.5583 ) {
                mic [ xcur ] [ ycur ] [ zcur ] = ETTR;
                count [ ETTR ] += 1;
                count [ DIFFC3A ] -= 1;
                nexp = 2;
            } else {
                /* indicate that diffusing species remains in current location */
                action = 7;
                nexp = 3;
            }

            /* Perform expansion that occurs when ettringite is formed */
            /* xexp, yexp and zexp are the coordinates of the last ettringite */
            /* pixel to be added */
            xexp = xnew;
            yexp = ynew;
            zexp = znew;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extettr(xexp, yexp, zexp, 0);
                /* update xexp, yexp and zexp */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last ettringite pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.6053 ) {
                newact = extettr(xexp, yexp, zexp, 0);
            }
        }
        /* examine reaction with diffusing anhydrite to form ettringite */
        /* Only allow reaction with diffusing anhydrite */
        else if ( ( check == DIFFANH ) && ( p2diff < C3AGYP ) ) {
            /* convert diffusing anhydrite to ettringite */
            mic [ xnew ] [ ynew ] [ znew ] = ETTR;
            count [ ETTR ] += 1;
            /* decrement counts of diffusing anhydrite */
            count [ DIFFANH ] -= 1;
            action = 0;

            /* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
            pexp = ran1(seed);
            nexp = 3;
            if ( pexp <= 0.569 ) {
                mic [ xcur ] [ ycur ] [ zcur ] = ETTR;
                count [ ETTR ] += 1;
                count [ DIFFC3A ] -= 1;
                nexp = 2;
            } else {
                /* indicate that diffusing species remains in current location */
                action = 7;
                nexp = 3;
            }

            /* Perform expansion that occurs when ettringite is formed */
            /* xexp, yexp and zexp are the coordinates of the last ettringite */
            /* pixel to be added */
            xexp = xnew;
            yexp = ynew;
            zexp = znew;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extettr(xexp, yexp, zexp, 0);
                /* update xexp, yexp and zexp */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last ettringite pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.6935 ) {
                newact = extettr(xexp, yexp, zexp, 0);
            }
        }
        /* examine reaction with diffusing CaCl2 to form FREIDEL */
        /* Only allow reaction with diffusing CaCl2 */
        else if ( check == DIFFCACL2 ) {
            /* convert diffusing C3A to Freidel's salt */
            mic [ xcur ] [ ycur ] [ zcur ] = FREIDEL;
            count [ FREIDEL ] += 1;
            /* decrement counts of diffusing C3A and CaCl2 */
            count [ DIFFC3A ] -= 1;
            action = 0;

            /* convert diffusing CACL2 to solid FREIDEL or else leave as a diffusing CACL2 */
            pexp = ran1(seed);
            nexp = 2;
            if ( pexp <= 0.5793 ) {
                mic [ xnew ] [ ynew ] [ znew ] = FREIDEL;
                count [ FREIDEL ] += 1;
                count [ DIFFCACL2 ] -= 1;
                nexp = 1;
            } else {
                nexp = 2;
            }

            /* Perform expansion that occurs when Freidel's salt is formed */
            /* xexp, yexp and zexp are the coordinates of the last FREIDEL */
            /* pixel to be added */
            xexp = xnew;
            yexp = ynew;
            zexp = znew;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extfreidel(xexp, yexp, zexp);
                /* update xexp, yexp and zexp */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last FREIDEL pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.3295 ) {
                newact = extfreidel(xexp, yexp, zexp);
            }
        }
        /* examine reaction with diffusing CAS2 to form STRAT */
        /* Only allow reaction with diffusing (not solid) CAS2 */
        else if ( check == DIFFCAS2 ) {
            /* convert diffusing CAS2 to stratlingite */
            mic [ xnew ] [ ynew ] [ znew ] = STRAT;
            count [ STRAT ] += 1;
            /* decrement counts of diffusing C3A and CAS2 */
            count [ DIFFCAS2 ] -= 1;
            action = 0;

            /* convert diffusing C3A to solid STRAT or else leave as a diffusing C3A */
            pexp = ran1(seed);
            nexp = 3;
            if ( pexp <= 0.886 ) {
                mic [ xcur ] [ ycur ] [ zcur ] = STRAT;
                count [ STRAT ] += 1;
                count [ DIFFC3A ] -= 1;
                nexp = 2;
            } else {
                action = 7;
                nexp = 3;
            }

            /* Perform expansion that occurs when stratlingite is formed */
            /* xexp, yexp and zexp are the coordinates of the last STRAT */
            /* pixel to be added */
            xexp = xnew;
            yexp = ynew;
            zexp = znew;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extstrat(xexp, yexp, zexp);
                /* update xexp, yexp and zexp */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last STRAT pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.286 ) {
                newact = extstrat(xexp, yexp, zexp);
            }
        }

        /* check for reaction with diffusing or solid ettringite to form AFm */
        /* reaction at solid ettringite only possible if ettringite is soluble */
        /* and even then on a limited bases to avoid a great formation of AFm */
        /* when ettringite first becomes soluble */
        pgrow = ran1(seed);
        if ( ( check == DIFFETTR ) || ( ( check == ETTR ) && ( soluble [ ETTR ] == 1 ) && ( pgrow <= C3AETTR ) ) ) {
            /* convert diffusing or solid ettringite to AFm */
            mic [ xnew ] [ ynew ] [ znew ] = AFM;
            count [ AFM ] += 1;
            /* decrement count of ettringite */
            count [ check ] -= 1;
            action = 0;

            /* convert diffusing C3A to AFm or leave as diffusing C3A */
            pexp = ran1(seed);
            if ( pexp <= 0.2424 ) {
                mic [ xcur ] [ ycur ] [ zcur ] = AFM;
                count [ AFM ] += 1;
                count [ DIFFC3A ] -= 1;
                pafm = ( -0.1 );
            } else {
                action = 7;
                pafm = 0.04699;
            }

            /* probabilistic-based expansion for new AFm pixel */
            pexp = ran1(seed);
            if ( pexp <= pafm ) {
                extafm(xnew, ynew, znew);
            }
        }

        if ( ( action != 0 ) && ( action != 7 ) ) {
            /* if diffusion is possible, execute it */
            if ( check == POROSITY ) {
                mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
                mic [ xnew ] [ ynew ] [ znew ] = DIFFC3A;
            } else {
                /* indicate that diffusing C3A remained */
                /* at original location */
                action = 7;
            }
        }
    }

    return ( action );
}

/* routine to move a diffusing C4A species */
/* from location (xcur,ycur,zcur) with nucleation probability of nucprob */
/* Called by hydrate */
/* Calls extc3ah6, moveone, extettr, and extafm */
int CemhydMatStatus :: movec4a(int xcur, int ycur, int zcur, int finalstep, float nucprob)
{
    int check, xnew, ynew, znew, action, sumgarb, sumold;
    int xexp, yexp, zexp, nexp, iexp, newact;
    float pgen, pexp, pafm, pgrow, p2diff;

    /* First be sure that a diffusing C4A species is at (xcur,ycur,zcur) */
    if ( mic [ xcur ] [ ycur ] [ zcur ] != DIFFC4A ) {
        action = 0;
        return ( action );
    }

    /* Check for nucleation into solid C3AH6 */
    pgen = ran1(seed);
    p2diff = ran1(seed);

    if ( ( nucprob >= pgen ) || ( finalstep == 1 ) ) {
        action = 0;
        mic [ xcur ] [ ycur ] [ zcur ] = C3AH6;
        count [ C3AH6 ] += 1;
        /* decrement count of diffusing C3A species */
        count [ DIFFC4A ] -= 1;
        /* allow for probabilistic-based expansion of C3AH6 */
        /* crystal to account for volume stoichiometry */
        pexp = ran1(seed);
        if ( pexp <= 0.69 ) {
            extc3ah6(xcur, ycur, zcur);
        }
    } else {
        /* determine new coordinates (using periodic boundaries) */
        xnew = xcur;
        ynew = ycur;
        znew = zcur;
        action = 0;
        sumold = 1;
        sumgarb = moveone(& xnew, & ynew, & znew, & action, sumold);
        if ( action == 0 ) {
            printf("Error in value of action in movec4a \n");
        }

        check = mic [ xnew ] [ ynew ] [ znew ];

        /* check for growth of C3AH6 crystal */
        if ( check == C3AH6 ) {
            pgrow = ran1(seed);
            /* Try to slow down growth of C3AH6 crystals to */
            /* promote ettringite and Afm formation */
            if ( pgrow <= C3AH6GROW ) {
                mic [ xcur ] [ ycur ] [ zcur ] = C3AH6;
                count [ C3AH6 ] += 1;
                count [ DIFFC4A ] -= 1;
                action = 0;
                /* allow for probabilistic-based expansion of C3AH6 */
                /* crystal to account for volume stoichiometry */
                pexp = ran1(seed);
                if ( pexp <= 0.69 ) {
                    extc3ah6(xcur, ycur, zcur);
                }
            }
        }
        /* examine reaction with diffusing gypsum to form ettringite */
        /* Only allow reaction with diffusing gypsum */
        else if ( ( check == DIFFGYP ) && ( p2diff < C3AGYP ) ) {
            /* convert diffusing gypsum to ettringite */
            mic [ xnew ] [ ynew ] [ znew ] = ETTRC4AF;
            count [ ETTRC4AF ] += 1;
            /* decrement counts of diffusing gypsum */
            count [ DIFFGYP ] -= 1;
            action = 0;

            /* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
            pexp = ran1(seed);
            nexp = 2;
            if ( pexp <= 0.40 ) {
                mic [ xcur ] [ ycur ] [ zcur ] = ETTRC4AF;
                count [ ETTRC4AF ] += 1;
                count [ DIFFC4A ] -= 1;
                nexp = 1;
            } else {
                /* indicate that diffusing species remains in current location */
                action = 7;
                nexp = 2;
            }

            /* Perform expansion that occurs when ettringite is formed */
            /* xexp, yexp and zexp are the coordinates of the last ettringite */
            /* pixel to be added */
            xexp = xnew;
            yexp = ynew;
            zexp = znew;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extettr(xexp, yexp, zexp, 1);
                /* update xexp, yexp and zexp */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last ettringite pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.30 ) {
                newact = extettr(xexp, yexp, zexp, 1);
            }
        }
        /* examine reaction with diffusing hemi to form ettringite */
        /* Only allow reaction with diffusing hemihydrate */
        else if ( ( check == DIFFHEM ) && ( p2diff < C3AGYP ) ) {
            /* convert diffusing hemihydrate to ettringite */
            mic [ xnew ] [ ynew ] [ znew ] = ETTRC4AF;
            count [ ETTRC4AF ] += 1;
            /* decrement counts of diffusing hemihydrate */
            count [ DIFFHEM ] -= 1;
            action = 0;

            /* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
            pexp = ran1(seed);
            nexp = 3;
            if ( pexp <= 0.5583 ) {
                mic [ xcur ] [ ycur ] [ zcur ] = ETTRC4AF;
                count [ ETTRC4AF ] += 1;
                count [ DIFFC4A ] -= 1;
                nexp = 2;
            } else {
                /* indicate that diffusing species remains in current location */
                action = 7;
                nexp = 3;
            }

            /* Perform expansion that occurs when ettringite is formed */
            /* xexp, yexp and zexp are the coordinates of the last ettringite */
            /* pixel to be added */
            xexp = xnew;
            yexp = ynew;
            zexp = znew;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extettr(xexp, yexp, zexp, 1);
                /* update xexp, yexp and zexp */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last ettringite pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.6053 ) {
                newact = extettr(xexp, yexp, zexp, 1);
            }
        }
        /* examine reaction with diffusing anhydrite to form ettringite */
        /* Only allow reaction with diffusing anhydrite */
        else if ( ( check == DIFFANH ) && ( p2diff < C3AGYP ) ) {
            /* convert diffusing anhydrite to ettringite */
            mic [ xnew ] [ ynew ] [ znew ] = ETTRC4AF;
            count [ ETTRC4AF ] += 1;
            /* decrement counts of diffusing anhydrite */
            count [ DIFFANH ] -= 1;
            action = 0;

            /* convert diffusing C3A to solid ettringite or else leave as a diffusing C3A */
            pexp = ran1(seed);
            nexp = 3;
            if ( pexp <= 0.569 ) {
                mic [ xcur ] [ ycur ] [ zcur ] = ETTRC4AF;
                count [ ETTRC4AF ] += 1;
                count [ DIFFC4A ] -= 1;
                nexp = 2;
            } else {
                /* indicate that diffusing species remains in current location */
                action = 7;
                nexp = 3;
            }

            /* Perform expansion that occurs when ettringite is formed */
            /* xexp, yexp and zexp are the coordinates of the last ettringite */
            /* pixel to be added */
            xexp = xnew;
            yexp = ynew;
            zexp = znew;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extettr(xexp, yexp, zexp, 1);
                /* update xexp, yexp and zexp */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last ettringite pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.6935 ) {
                newact = extettr(xexp, yexp, zexp, 1);
            }
        }
        /* examine reaction with diffusing CaCl2 to form FREIDEL */
        /* Only allow reaction with diffusing CaCl2 */
        else if ( check == DIFFCACL2 ) {
            /* convert diffusing C3A to Freidel's salt */
            mic [ xcur ] [ ycur ] [ zcur ] = FREIDEL;
            count [ FREIDEL ] += 1;
            /* decrement counts of diffusing C3A and CaCl2 */
            count [ DIFFC4A ] -= 1;
            action = 0;

            /* convert diffusing CACL2 to solid FREIDEL or else leave as a diffusing CACL2 */
            pexp = ran1(seed);
            nexp = 2;
            if ( pexp <= 0.5793 ) {
                mic [ xnew ] [ ynew ] [ znew ] = FREIDEL;
                count [ FREIDEL ] += 1;
                count [ DIFFCACL2 ] -= 1;
                nexp = 1;
            } else {
                nexp = 2;
            }

            /* Perform expansion that occurs when Freidel's salt is formed */
            /* xexp, yexp and zexp are the coordinates of the last FREIDEL */
            /* pixel to be added */
            xexp = xnew;
            yexp = ynew;
            zexp = znew;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extfreidel(xexp, yexp, zexp);
                /* update xexp, yexp and zexp */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last FREIDEL pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.3295 ) {
                newact = extfreidel(xexp, yexp, zexp);
            }
        }
        /* examine reaction with diffusing CAS2 to form STRAT */
        /* Only allow reaction with diffusing (not solid) CAS2 */
        else if ( check == DIFFCAS2 ) {
            /* convert diffusing CAS2 to stratlingite */
            mic [ xnew ] [ ynew ] [ znew ] = STRAT;
            count [ STRAT ] += 1;
            /* decrement counts of diffusing CAS2 */
            count [ DIFFCAS2 ] -= 1;
            action = 0;

            /* convert diffusing C3A to solid STRAT or else leave as a diffusing C3A */
            pexp = ran1(seed);
            nexp = 3;
            if ( pexp <= 0.886 ) {
                mic [ xcur ] [ ycur ] [ zcur ] = STRAT;
                count [ STRAT ] += 1;
                count [ DIFFC4A ] -= 1;
                nexp = 2;
            } else {
                action = 7;
                nexp = 3;
            }

            /* Perform expansion that occurs when stratlingite is formed */
            /* xexp, yexp and zexp are the coordinates of the last STRAT */
            /* pixel to be added */
            xexp = xnew;
            yexp = ynew;
            zexp = znew;
            for ( iexp = 1; iexp <= nexp; iexp++ ) {
                newact = extstrat(xexp, yexp, zexp);
                /* update xexp, yexp and zexp */
                switch ( newact ) {
                case 1:
                    xexp -= 1;
                    if ( xexp < 0 ) {
                        xexp = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xexp += 1;
                    if ( xexp >= SYSIZE ) {
                        xexp = 0;
                    }

                    break;
                case 3:
                    yexp -= 1;
                    if ( yexp < 0 ) {
                        yexp = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    yexp += 1;
                    if ( yexp >= SYSIZE ) {
                        yexp = 0;
                    }

                    break;
                case 5:
                    zexp -= 1;
                    if ( zexp < 0 ) {
                        zexp = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zexp += 1;
                    if ( zexp >= SYSIZE ) {
                        zexp = 0;
                    }

                    break;
                default:
                    break;
                }
            }

            /* probabilistic-based expansion for last STRAT pixel */
            pexp = ran1(seed);
            if ( pexp <= 0.286 ) {
                newact = extstrat(xexp, yexp, zexp);
            }
        }

        /* check for reaction with diffusing or solid ettringite to form AFm */
        /* reaction at solid ettringite only possible if ettringite is soluble */
        /* and even then on a limited bases to avoid a great formation of AFm */
        /* when ettringite first becomes soluble */
        pgrow = ran1(seed);
        if ( ( check == DIFFETTR ) || ( ( check == ETTR ) && ( soluble [ ETTR ] == 1 ) && ( pgrow <= C3AETTR ) ) ) {
            /* convert diffusing or solid ettringite to AFm */
            mic [ xnew ] [ ynew ] [ znew ] = AFM;
            count [ AFM ] += 1;
            /* decrement count of ettringite */
            count [ check ] -= 1;
            action = 0;

            /* convert diffusing C4A to AFm or leave as diffusing C4A */
            pexp = ran1(seed);
            if ( pexp <= 0.2424 ) {
                mic [ xcur ] [ ycur ] [ zcur ] = AFM;
                count [ AFM ] += 1;
                count [ DIFFC4A ] -= 1;
                pafm = ( -0.1 );
            } else {
                action = 7;
                pafm = 0.04699;
            }

            /* probabilistic-based expansion for new AFm pixel */
            pexp = ran1(seed);
            if ( pexp <= pafm ) {
                extafm(xnew, ynew, znew);
            }
        }

        if ( ( action != 0 ) && ( action != 7 ) ) {
            /* if diffusion is possible, execute it */
            if ( check == POROSITY ) {
                mic [ xcur ] [ ycur ] [ zcur ] = POROSITY;
                mic [ xnew ] [ ynew ] [ znew ] = DIFFC4A;
            } else {
                /* indicate that diffusing C4A remained */
                /* at original location */
                action = 7;
            }
        }
    }

    return ( action );
}

/* routine to oversee hydration by updating position of all */
/* remaining diffusing species */
/* Calls movech, movec3a, movefh3, moveettr, movecsh, and movegyp */

void CemhydMatStatus :: hydrate(int fincyc, int stepmax, float chpar1, float chpar2, float hgpar1, float hgpar2, float fhpar1, float fhpar2, float gypar1, float gypar2)
{
    int xpl, ypl, zpl, phpl, agepl, xpnew, ypnew, zpnew;
    float chprob, c3ah6prob, fh3prob, gypprob;
    long int nleft, ntodo, ndale;
    int istep, termflag, reactf = -1;
    float beterm;
    struct ants *curant, *antgone;

    ntodo = nmade;
    nleft = nmade;
    termflag = 0;

    /* Perform diffusion until all reacted or max. # of diffusion steps reached */
    for ( istep = 1; ( ( istep <= stepmax ) && ( nleft > 0 ) ); istep++ ) {
        if ( ( fincyc == 1 ) && ( istep == stepmax ) ) {
            termflag = 1;
        }

        nleft = 0;
        ndale = 0;

        /* determine probabilities for CH and C3AH6 nucleation */
        beterm = exp(-( double ) ( count [ DIFFCH ] ) * 1000000. / SYSIZE_POW3 / chpar2);       //fixed
        chprob = chpar1 * ( 1. - beterm );
        beterm = exp(-( double ) ( count [ DIFFC3A ] ) * 1000000. / SYSIZE_POW3 / hgpar2);       //fixed
        c3ah6prob = hgpar1 * ( 1. - beterm );
        beterm = exp(-( double ) ( count [ DIFFFH3 ] ) * 1000000. / SYSIZE_POW3 / fhpar2);       //fixed
        fh3prob = fhpar1 * ( 1. - beterm );
        beterm = exp(-( double ) ( count [ DIFFANH ] + count [ DIFFHEM ] ) * 1000000. / SYSIZE_POW3 / gypar2);       //fixed
        gypprob = gypar1 * ( 1. - beterm );

        /* Process each diffusing species in turn */
        curant = headant->nextant;
        while ( curant != NULL ) {
            ndale += 1;
            xpl = curant->x;
            ypl = curant->y;
            zpl = curant->z;
            phpl = curant->id;
            agepl = curant->cycbirth;

            /* based on ID, call appropriate routine to process diffusing species */
            if ( phpl == DIFFCSH ) {
                /* printf("Calling movecsh \n");
                 * fflush(stdout); */
                reactf = movecsh(xpl, ypl, zpl, termflag, agepl);
            } else if ( phpl == DIFFANH ) {
                /* printf("Calling moveanh \n");
                 * fflush(stdout); */
                reactf = moveanh(xpl, ypl, zpl, termflag, gypprob);
            } else if ( phpl == DIFFHEM ) {
                /* printf("Calling movehem \n");
                 * fflush(stdout); */
                reactf = movehem(xpl, ypl, zpl, termflag, gypprob);
            } else if ( phpl == DIFFCH ) {
                /* printf("Calling movech \n");
                 * fflush(stdout); */
                reactf = movech(xpl, ypl, zpl, termflag, chprob);
            } else if ( phpl == DIFFFH3 ) {
                /* printf("Calling movefh3 \n");
                 * fflush(stdout); */
                reactf = movefh3(xpl, ypl, zpl, termflag, fh3prob);
            } else if ( phpl == DIFFGYP ) {
                /* printf("Calling movegyp \n");
                 * fflush(stdout); */
                reactf = movegyp(xpl, ypl, zpl, termflag);
            } else if ( phpl == DIFFC3A ) {
                /* printf("Calling movec3a \n");
                 * fflush(stdout); */
                reactf = movec3a(xpl, ypl, zpl, termflag, c3ah6prob);
            } else if ( phpl == DIFFC4A ) {
                /* printf("Calling movec4a \n");
                 * fflush(stdout); */
                reactf = movec4a(xpl, ypl, zpl, termflag, c3ah6prob);
            } else if ( phpl == DIFFETTR ) {
                /* printf("Calling moveettr \n");
                 * fflush(stdout); */
                reactf = moveettr(xpl, ypl, zpl, termflag);
            } else if ( phpl == DIFFCACL2 ) {
                /* printf("Calling movecacl2 \n");
                 * fflush(stdout); */
                reactf = movecacl2(xpl, ypl, zpl, termflag);
            } else if ( phpl == DIFFCAS2 ) {
                /* printf("Calling movecas2 \n");
                 * fflush(stdout); */
                reactf = movecas2(xpl, ypl, zpl, termflag);
            } else if ( phpl == DIFFAS ) {
                /* printf("Calling moveas \n");
                 * fflush(stdout); */
                reactf = moveas(xpl, ypl, zpl, termflag);
            } else if ( phpl == DIFFCACO3 ) {
                /* printf("Calling movecaco3 \n");
                 * fflush(stdout); */
                reactf = movecaco3(xpl, ypl, zpl, termflag);
            } else {
                printf("Error in ID of phase \n");
            }

            /* if no reaction */
            if ( reactf != 0 ) {
                nleft += 1;
                xpnew = xpl;
                ypnew = ypl;
                zpnew = zpl;

                /* update location of diffusing species */
                switch ( reactf ) {
                case 1:
                    xpnew -= 1;
                    if ( xpnew < 0 ) {
                        xpnew = ( SYSIZEM1 );
                    }

                    break;
                case 2:
                    xpnew += 1;
                    if ( xpnew >= SYSIZE ) {
                        xpnew = 0;
                    }

                    break;
                case 3:
                    ypnew -= 1;
                    if ( ypnew < 0 ) {
                        ypnew = ( SYSIZEM1 );
                    }

                    break;
                case 4:
                    ypnew += 1;
                    if ( ypnew >= SYSIZE ) {
                        ypnew = 0;
                    }

                    break;
                case 5:
                    zpnew -= 1;
                    if ( zpnew < 0 ) {
                        zpnew = ( SYSIZEM1 );
                    }

                    break;
                case 6:
                    zpnew += 1;
                    if ( zpnew >= SYSIZE ) {
                        zpnew = 0;
                    }

                    break;
                default:
                    break;
                }

                /* store new location of diffusing species */
                curant->x = xpnew;
                curant->y = ypnew;
                curant->z = zpnew;
                curant->id = phpl;
                curant = curant->nextant;
            } /* end of reactf!=0 loop */
              /* else remove ant from list */
            else {
                if ( ndale == 1 ) {
                    headant->nextant = curant->nextant;
                } else {
                    ( curant->prevant )->nextant = curant->nextant;
                }

                if ( curant->nextant != NULL ) {
                    ( curant->nextant )->prevant = curant->prevant;
                } else {
                    tailant = curant->prevant;
                }

                antgone = curant;
                curant = curant->nextant;
                //printf("DEMEM antgone addr %p\n", (const void *)antgone);
                free(antgone);
                ngoing -= 1;
            }
        }

        /* end of curant loop */
        ntodo = nleft;
    }

    /* end of istep loop */
}


void CemhydMatStatus :: laguer(fcomplex_cem a[], int m, fcomplex_cem *x, float eps, int polish) {
    int j, iter;
    float err, dxold, cdx, abx;
    fcomplex_cem sq, h, gp, gm, g2, g, b, d, dx, f, x1;

    dxold = Cabs(* x);
    for ( iter = 1; iter <= MAXIT; iter++ ) {
        b = a [ m ];
        err = Cabs(b);
        d = f = ComplexCemhyd(0.0, 0.0);
        abx = Cabs(* x);
        for ( j = m - 1; j >= 0; j-- ) {
            f = Cadd(Cmul(* x, f), d);
            d = Cadd(Cmul(* x, d), b);
            b = Cadd(Cmul(* x, b), a [ j ]);
            err = Cabs(b) + abx * err;
        }

        err *= EPSS;
        if ( Cabs(b) <= err ) {
            return;
        }

        g = Cdiv(d, b);
        g2 = Cmul(g, g);
        h = Csub( g2, RCmul( 2.0, Cdiv(f, b) ) );
        sq = Csqrt( RCmul( ( float ) ( m - 1 ), Csub(RCmul( ( float ) m, h ), g2) ) );
        gp = Cadd(g, sq);
        gm = Csub(g, sq);
        if ( Cabs(gp) < Cabs(gm) ) {
            gp = gm;
        }

        dx = Cdiv(ComplexCemhyd( ( float ) m, 0.0 ), gp);
        x1 = Csub(* x, dx);
        if ( x->r == x1.r && x->i == x1.i ) {
            return;
        }

        * x = x1;
        cdx = Cabs(dx);
        if ( iter > 6 && cdx >= dxold ) {
            return;
        }

        dxold = cdx;
        if ( !polish ) {
            if ( cdx <= eps * Cabs(* x) ) {
                return;
            }
        }
    }

    nrerror("Too many iterations in routine CemhydMat - LAGUER\n");
}



void CemhydMatStatus :: zroots(fcomplex_cem a[], int m, fcomplex_cem roots[], int polish) {
    int jj, j, i;
    fcomplex_cem x, b, c;
    fcomplex_cem *ad;
    //void laguer();fixed-no influence
    ad = new fcomplex_cem [ MAXM ];

    for ( j = 0; j <= m; j++ ) {
        ad [ j ] = a [ j ];
    }

    for ( j = m; j >= 1; j-- ) {
        x = ComplexCemhyd(0.0, 0.0);
        laguer(ad, j, & x, EPSP, 0);
        if ( fabs(x.i) <= ( 2.0 * EPSP * fabs(x.r) ) ) {
            x.i = 0.0;
        }

        roots [ j ] = x;
        b = ad [ j ];
        for ( jj = j - 1; jj >= 0; jj-- ) {
            c = ad [ jj ];
            ad [ jj ] = b;
            b = Cadd(Cmul(x, b), c);
        }
    }

    if ( polish ) {
        for ( j = 1; j <= m; j++ ) {
            laguer(a, m, & roots [ j ], EPSP, 1);
        }
    }

    for ( j = 2; j <= m; j++ ) {
        x = roots [ j ];
        for ( i = j - 1; i >= 1; i-- ) {
            if ( roots [ i ].r <= x.r ) {
                break;
            }

            roots [ i + 1 ] = roots [ i ];
        }

        roots [ i + 1 ] = x;
    }

    delete [] ad;
}



void CemhydMatStatus :: pHpred(void) { //dangerous function - can lead to zero division etc.
    int j, syngen_change = 0, syn_old = 0, counter = 0;
    double concnaplus, conckplus;
    double concohminus = 0., A, B, C, conctest, concsulfate1;
    double volpore, grams_cement;
    double releasedna, releasedk, activitySO4, activityK, test_precip;
    double activityCa, activityOH, Istrength, Anow, Bnow, Inew;
    double conductivity = 0.0;
    fcomplex_cem coef [ 5 ], roots [ 5 ];
    float sumbest, sumtest, pozzreact, KspCH;

    /* Update CH activity product based on current system temperature */
    /* Factors derived from fitting CH solubility vs. temperature */
    /* data in Taylor book (p. 117) */
    KspCH = KspCH25C * ( 1.534385 - 0.02057 * temp_cur );

    if ( conccaplus > 1.0 ) {
        conccaplus = 0.0;
    }

    /* Calculate volume of pore solution in the concrete in Liters */
    volpore = ( double ) count [ POROSITY ] + 0.38 * ( double ) count [ CSH ] + 0.20 * ( double ) count [ POZZCSH ] + 0.20 * count [ SLAGCSH ]; //relative amount, no need to rescale
    /* Convert from pixels (cubic micrometers) to liters (cubic decimeters) */
    volpore *= VOLFACTOR;
    volpore *= VOLFACTOR;
    volpore *= VOLFACTOR;
    /* Compute pore volume per gram of cement */
    grams_cement = cemmasswgyp * MASSFACTOR * MASSFACTOR * MASSFACTOR;
    /* Compute pore volume per gram of cement */
    volpore /= grams_cement;
    /* Compute grams of pozzolan which have reacted */
    pozzreact = ( ( float ) npr / 1.35 ) * MASSFACTOR * MASSFACTOR * MASSFACTOR * specgrav [ POZZ ];

    /* Compute moles of released potassium and sodium per gram of cement*/
    if ( time_cur > 1.0 ) {
        releasedk = ( 2. * ( rspotassium + ( totpotassium - rspotassium ) * alpha_cur ) );
        releasedk /= MMK2O;
        releasedna = ( 2. * ( rssodium + ( totsodium - rssodium ) * alpha_cur ) );
        releasedna /= MMNa2O;
    } else {
        /* Proportion the early time release over the first hour */
        /* 90% immediately and the remaining 10% over the first hour */
        /* based on limited data from Davide Zampini (MBT) */
        releasedk = ( 2. * ( ( 0.9 + 0.1 * time_cur ) * ( rspotassium ) + ( totpotassium - rspotassium ) * alpha_cur ) );
        releasedk /= MMK2O;
        releasedna = ( 2. * ( ( 0.9 + 0.1 * time_cur ) * ( rssodium ) + ( totsodium - rssodium ) * alpha_cur ) );
        releasedna /= MMNa2O;
    }

    /* Compute concentrations of K+ and Na+ in pore solution currently */
    /* Remember to decrease K+ by KperSyn*moles of syngenite precipitated */
    /* Units must be in moles/gram for both */
    conckplus = ( ( releasedk - moles_syn_precip * KperSyn ) / ( volpore + BK * alpha_cur + BprimeK * pozzreact ) );
    concnaplus = ( releasedna ) / ( volpore + BNa * alpha_cur + BprimeNa * pozzreact );

    do { /* while Syngenite precipitating loop */
         /* Now compute the activities (estimated) of Ca++ and OH- */
        activityCa = activityOH = activitySO4 = activityK = 1.0;
        Inew = 0.0;
        if ( ( ( concnaplus + conckplus ) > 0.0 ) && ( soluble [ ETTR ] == 0 ) ) {
            /* Factor of 1000 to convert from M to mmol/L */
            Istrength = 1000. * ( zK * zK * conckplus + zNa * zNa * concnaplus + zCa * zCa * conccaplus );
            if ( Istrength < 1.0 ) {
                Istrength = 1.0;
            }

            while ( ( fabs(Istrength - Inew) / Istrength ) > 0.10 ) {
                Istrength = 1000. * ( zK * zK * conckplus + zNa * zNa * concnaplus + zCa * zCa * conccaplus );
                if ( Istrength < 1.0 ) {
                    Istrength = 1.0;
                }

                Anow = activeA0 * 295. * sqrt(295.) / ( ( temp_cur + 273.15 ) * sqrt(temp_cur + 273.15) );
                Bnow = activeB0 * sqrt(295.) / ( sqrt(temp_cur + 273.15) );
                /* Equations from papers of Marchand et al. */
                activityCa = ( -Anow * zCa * zCa * sqrt(Istrength) ) / ( 1. + aCa * Bnow * sqrt(Istrength) );
                activityCa += ( 0.2 - 0.0000417 * Istrength ) * Anow * zCa * zCa * Istrength / sqrt(1000.);
                activityCa = exp(activityCa);
                activityOH = ( -Anow * zOH * zOH * sqrt(Istrength) ) / ( 1. + aOH * Bnow * sqrt(Istrength) );
                activityOH += ( 0.2 - 0.0000417 * Istrength ) * Anow * zOH * zOH * Istrength / sqrt(1000.);
                activityOH = exp(activityOH);
                activityK = ( -Anow * zK * zK * sqrt(Istrength) ) / ( 1. + aK * Bnow * sqrt(Istrength) );
                activityK += ( 0.2 - 0.0000417 * Istrength ) * Anow * zK * zK * Istrength / sqrt(1000.);
                activityK = exp(activityK);
                activitySO4 = ( -Anow * zSO4 * zSO4 * sqrt(Istrength) ) / ( 1. + aSO4 * Bnow * sqrt(Istrength) );
                activitySO4 += ( 0.2 - 0.0000417 * Istrength ) * Anow * zSO4 * zSO4 * Istrength / sqrt(1000.);
                activitySO4 = exp(activitySO4);
                /* Now try to find roots of fourth degree polynomial */
                /* to determine sulfate, OH-, and calcium ion concentrations */
                /*    A=(-KspCH); */
                /* Now with activities */
                A = ( -KspCH / ( activityCa * activityOH * activityOH ) );
                B = conckplus + concnaplus;
                C = ( -2. * KspGypsum / ( activityCa * activitySO4 ) );
                concohminus = conckplus + concnaplus;
                coef [ 0 ] = ComplexCemhyd(C, 0.0);
                coef [ 1 ] = ComplexCemhyd( ( A + 2. * B * C ) / C, 0.0 );
                coef [ 2 ] = ComplexCemhyd(B * B / C + 4., 0.0);
                coef [ 3 ] = ComplexCemhyd(4. * B / C, 0.0);
                coef [ 4 ] = ComplexCemhyd(4. / C, 0.0);
                /*         printf("coef 0 is (%f,%f)\n",coef[0].r,coef[0].i);
                 *     printf("coef 1 is (%f,%f)\n",coef[1].r,coef[1].i);
                 *     printf("coef 2 is (%f,%f)\n",coef[2].r,coef[2].i);
                 *     printf("coef 3 is (%f,%f)\n",coef[3].r,coef[3].i);
                 *     printf("coef 4 is (%f,%f)\n",coef[4].r,coef[4].i); */
                roots [ 1 ] = ComplexCemhyd(0.0, 0.0);
                roots [ 2 ] = ComplexCemhyd(0.0, 0.0);
                roots [ 3 ] = ComplexCemhyd(0.0, 0.0);
                roots [ 4 ] = ComplexCemhyd(0.0, 0.0);
                zroots(coef, 4, roots, 1);
                sumbest = 100;
                /* Find the best real root for electoneutrality */
                for ( j = 1; j <= 4; j++ ) {
                    //printf("pH root %d is (%f,%f)\n",j,roots[j].r,roots[j].i);
                    fflush(stdout);
                    if ( ( ( roots [ j ].i ) == 0.0 ) && ( ( roots [ j ].r ) > 0.0 ) ) {
                        conctest = sqrt( KspCH / ( roots [ j ].r * activityCa * activityOH * activityOH ) );
                        concsulfate1 = KspGypsum / ( roots [ j ].r * activityCa * activitySO4 );
                        sumtest = concnaplus + conckplus + 2. * roots [ j ].r - conctest - 2. * concsulfate1;
                        if ( fabs(sumtest) < sumbest ) {
                            sumbest = fabs(sumtest);
                            concohminus = conctest;
                            conccaplus = roots [ j ].r;
                            concsulfate = concsulfate1;
                        }
                    }
                }

                /* Update ionic strength */
                Inew = 1000. * ( zK * zK * conckplus + zNa * zNa * concnaplus + zCa * zCa * conccaplus );
            }

            /* end of while loop for Istrength-Inew */
        } else {
            /* Factor of 1000 to convert from M to mmol/L */
            Istrength = 1000. * ( zK * zK * conckplus + zNa * zNa * concnaplus + zCa * zCa * conccaplus );
            if ( Istrength < 1.0 ) {
                Istrength = 1.0;
            }

            while ( ( fabs(Istrength - Inew) / Istrength ) > 0.10 && counter < 4 ) {
                //add counter to limit to 3 iterations, normally it suffices but when high wcr it hangs - smilauer@cml.fsv.cvut.cz
                Istrength = 1000. * ( zK * zK * conckplus + zNa * zNa * concnaplus + zCa * zCa * conccaplus );
                Anow = activeA0 * 295. * sqrt(295.) / ( ( temp_cur + 273.15 ) * sqrt(temp_cur + 273.15) );
                Bnow = activeB0 * sqrt(295.) / ( sqrt(temp_cur + 273.15) );
                /* Equations from papers of Marchand et al. */
                activityCa = ( -Anow * zCa * zCa * sqrt(Istrength) ) / ( 1. + aCa * Bnow * sqrt(Istrength) );
                activityCa += ( 0.2 - 0.0000417 * Istrength ) * Anow * zCa * zCa * Istrength / sqrt(1000.);
                activityCa = exp(activityCa);
                activityOH = ( -Anow * zOH * zOH * sqrt(Istrength) ) / ( 1. + aOH * Bnow * sqrt(Istrength) );
                activityOH += ( 0.2 - 0.0000417 * Istrength ) * Anow * zOH * zOH * Istrength / sqrt(1000.);
                activityOH = exp(activityOH);
                activityK = ( -Anow * zK * zK * sqrt(Istrength) ) / ( 1. + aK * Bnow * sqrt(Istrength) );
                activityK += ( 0.2 - 0.0000417 * Istrength ) * Anow * zK * zK * Istrength / sqrt(1000.);
                activityK = exp(activityK);
                /* Calculate pH assuming simply that OH- balances sum of Na+ and K+ */
                concohminus = conckplus + concnaplus;
                if ( ( conccaplus ) > ( 0.1 * ( concohminus ) ) ) {
                    concohminus += ( 2. * conccaplus );
                }

                conccaplus = ( KspCH / ( activityCa * activityOH * activityOH * concohminus * concohminus ) );
                concsulfate = 0.0;
                /* Update ionic strength */
                Inew = 1000. * ( zK * zK * conckplus + zNa * zNa * concnaplus + zCa * zCa * conccaplus );
                counter++;
            }

            /* end of while loop for Istrength-Inew */
        }

        /* Check for syngenite precipitation */
        syngen_change = 0;
        if ( syn_old != 2 ) {
            test_precip = conckplus * conckplus * activityK * activityK;
            test_precip *= conccaplus * activityCa;
            test_precip *= concsulfate * concsulfate * activitySO4 * activitySO4;
            if ( test_precip > KspSyngenite ) {
#ifdef PRINTF
                printf("Syngenite precipitating at cycle %d\n", icyc);
#endif
                syngen_change = syn_old = 1;
                /* Units of moles_syn_precip are moles per gram of cement */
                if ( conckplus > 0.002 ) {
                    conckplus -= 0.001;
                    moles_syn_precip += 0.001 * volpore / KperSyn;
                } else if ( conckplus > 0.0002 ) {
                    conckplus -= 0.0001;
                    moles_syn_precip += 0.0001 * volpore / KperSyn;
                } else {
                    moles_syn_precip += conckplus * volpore / KperSyn;
                    conckplus = 0.0;
                }
            }

            /* Check for syngenite dissolution */
            /* How to control dissolution rates??? */
            /* Have 0.001*KperSyn increase in conckplus each cycle */
            /* Only one dissolution per cycle --- purpose of syn_old */
            /* and no dissolution if some precipitation in that cycle */
            if ( ( syn_old == 0 ) && ( moles_syn_precip > 0.0 ) ) {
                syngen_change = syn_old = 2;
                /*               conckplus+=(moles_syn_precip/10.0)/volpore; */
                if ( ( moles_syn_precip / volpore ) > 0.001 ) {
                    conckplus += 0.001 * KperSyn;
                    moles_syn_precip -= ( 0.001 * volpore );
                } else {
                    conckplus += ( moles_syn_precip * KperSyn / volpore );
                    moles_syn_precip = 0.0;
                }
            }
        }
    } while ( syngen_change != 0 );

    if ( concohminus < ( 0.0000001 ) ) {
        concohminus = 0.0000001;
        conccaplus = ( KspCH / ( activityCa * activityOH * activityOH * concohminus * concohminus ) );
    }

    pH_cur = 14.0 + log10(concohminus * activityOH);
    /* Calculation of solution conductivity (Snyder and Feng basis) */
    /* First convert ionic strength back to M units */
    Istrength /= 1000.;
    conductivity += zCa * conccaplus * ( lambdaCa_0 / ( 1. + GCa * sqrt(Istrength) ) );
    conductivity += zOH * concohminus * ( lambdaOH_0 / ( 1. + GOH * sqrt(Istrength) ) );
    conductivity += zNa * concnaplus * ( lambdaNa_0 / ( 1. + GNa * sqrt(Istrength) ) );
    conductivity += zK * conckplus * ( lambdaK_0 / ( 1. + GK * sqrt(Istrength) ) );
    conductivity += zSO4 * concsulfate * ( lambdaSO4_0 / ( 1. + GSO4 * sqrt(Istrength) ) );
    conductivity *= cm2perL2m;

    /* Output results to logging file */
#ifdef OUTFILES
    fprintf(pHfile, "%d %.4f %f %.4f %f %f %f %f %f %.4f %.4f %.4f %.4f %f\n", cyccnt - 1, time_cur, alpha_cur, pH_cur, conductivity, concnaplus, conckplus, conccaplus, concsulfate, activityCa, activityOH, activitySO4, activityK, moles_syn_precip);
    fflush(pHfile);
#endif
}


int CemhydMatStatus :: IsSolidPhase(int phase) {
    if ( ( phase >= C3S && phase <= ABSGYP ) || ( phase == HDCSH ) ) {
        return 1;
    } else {
        return 0;
    }
}

void CemhydMatStatus :: burn_phases(int d1, int d2, int d3) {
    long int icur, inew, ncur, nnew;
    int i, j, k, cnt_perc, cnt_tot, kh;
    int *nmatx, *nmaty, *nmatz;
    int xl, xh, j1, k1, px, py, pz, qx, qy, qz;
    int xcn, ycn, zcn, x1, y1, z1, igood;
    int *nnewx, *nnewy, *nnewz;
    int jnew;
    int phase_temp [ 51 ];
    char ***newmat;


    alloc_char_3D(newmat, SYSIZE);
    nmatx = new int [ SIZESET ];
    nmaty = new int [ SIZESET ];
    nmatz = new int [ SIZESET ];
    nnewx = new int [ SIZESET ];
    nnewy = new int [ SIZESET ];
    nnewz = new int [ SIZESET ];

    last = NULL;

    /* counters for number of pixels of phase accessible from surface #1 */
    /* and number which are part of a percolated pathway to surface #2 */
    /* create copy of original microstructure to newmat[][][]*/
    /* array mic_CSH[][][] will not be changed through this routine */
    /* phase_temp[] and phase[] are for phases storage in percolated pathway */
    //  ntop=0;
    //  nthrough=0;
    for ( k = 0; k < SYSIZE; k++ ) {
        for ( j = 0; j < SYSIZE; j++ ) {
            for ( i = 0; i < SYSIZE; i++ ) {
                newmat [ i ] [ j ] [ k ] = mic_CSH [ i ] [ j ] [ k ];
                //assign 0 or EMPTYP to ArrPerc[][][]
                if ( mic_CSH [ i ] [ j ] [ k ] == EMPTYP ) {
                    ArrPerc [ i ] [ j ] [ k ] = EMPTYP;
                } else {
                    ArrPerc [ i ] [ j ] [ k ] = 0;
                }
            }
        }
    }

    for ( k = 0; k < 51; k++ ) {
        phase [ k ] = 0;
    }

    /* percolation is assessed from top to bottom only */
    /* in transformed coordinates */
    /* and burning algorithm is periodic in other two directions */
    i = 0;

    for ( k = 0; k < SYSIZE; k++ ) {
        for ( j = 0; j < SYSIZE; j++ ) {
            igood = 0;
            ncur = 0;
            for ( kh = 0; kh < 51; kh++ ) {
                phase_temp [ kh ] = 0;
            }


            //delete all array in case of next cycle
            while ( last != NULL ) { //go through all voxel coordinates
                current = last;
                last = last->prev;
                /*deallocate memory*/
                delete current;
            }

            /* Transform coordinates */
            px = cx(i, j, k, d1, d2, d3);
            py = cy(i, j, k, d1, d2, d3);
            pz = cz(i, j, k, d1, d2, d3);
            /* start from any solid phase */
            if ( IsSolidPhase(newmat [ px ] [ py ] [ pz ]) == 1 ) {
                /* Start a burn front */
                phase_temp [ ( int ) newmat [ px ] [ py ] [ pz ] ]++;

                //store in unnamed list
                WriteUnsortedList(px, py, pz);
                newmat [ px ] [ py ] [ pz ] = BURNT;
                ncur += 1;



                /* burn front is stored in matrices nmat* */
                /* and nnew* */
                nmatx [ ncur ] = i;
                nmaty [ ncur ] = j;
                nmatz [ ncur ] = k;
                /* Burn as long as new (fuel) pixels are found */
                do {
                    nnew = 0;
                    for ( inew = 1; inew <= ncur; inew++ ) {
                        xcn = nmatx [ inew ];
                        ycn = nmaty [ inew ];
                        zcn = nmatz [ inew ];
                        /* Convert to directional coordinates */
                        qx = cx(xcn, ycn, zcn, d1, d2, d3);
                        qy = cy(xcn, ycn, zcn, d1, d2, d3);
                        qz = cz(xcn, ycn, zcn, d1, d2, d3);

                        /* Check all six neighbors */
                        for ( jnew = 1; jnew <= 6; jnew++ ) {
                            x1 = xcn;
                            y1 = ycn;
                            z1 = zcn;
                            if ( jnew == 1 ) {
                                x1 -= 1;
                            }

                            if ( jnew == 2 ) {
                                x1 += 1;
                            }

                            if ( jnew == 3 ) {
                                y1 -= 1;
                            }

                            if ( jnew == 4 ) {
                                y1 += 1;
                            }

                            if ( jnew == 5 ) {
                                z1 -= 1;
                            }

                            if ( jnew == 6 ) {
                                z1 += 1;
                            }

                            /* Periodic in y and */
                            if ( y1 >= SYSIZE ) {
                                y1 -= SYSIZE;
                            } else if ( y1 < 0 ) {
                                y1 += SYSIZE;
                            }

                            /* Periodic in z direction */
                            if ( z1 >= SYSIZE ) {
                                z1 -= SYSIZE;
                            } else if ( z1 < 0 ) {
                                z1 += SYSIZE;
                            }

                            /* Nonperiodic so be sure to remain in the 3-D box */
                            if ( ( x1 >= 0 ) && ( x1 < SYSIZE ) ) {
                                px = cx(x1, y1, z1, d1, d2, d3);
                                py = cy(x1, y1, z1, d1, d2, d3);
                                pz = cz(x1, y1, z1, d1, d2, d3);
                                /* Conditions for propagation of burning */
                                /* 1) new pixel is any solid except clinker phases */

                                if ( ( IsSolidPhase(newmat [ px ] [ py ] [ pz ]) == 1 ) &&
                                    ( newmat [ px ] [ py ] [ pz ] != C3S ) &&
                                    ( newmat [ px ] [ py ] [ pz ] != C2S ) &&
                                    ( newmat [ px ] [ py ] [ pz ] != C3A ) &&
                                    ( newmat [ px ] [ py ] [ pz ] != C4AF ) ) {
                                    /*
                                     * if((newmat[px][py][pz]==CSH)||(newmat[px][py][pz]==ETTRC4AF)||(newmat[px][py][pz]==ETTR)){
                                     */
                                    phase_temp [ ( int ) newmat [ px ] [ py ] [ pz ] ]++;
                                    WriteUnsortedList(px, py, pz);
                                    newmat [ px ] [ py ] [ pz ] = BURNT;
                                    nnew += 1;
                                    if ( nnew >= SIZESET ) {
                                        printf("error in size of nnew %ld\n", nnew);
                                    }

                                    nnewx [ nnew ] = x1;
                                    nnewy [ nnew ] = y1;
                                    nnewz [ nnew ] = z1;
                                }
                                /* 2) old pixel is solid except clinker and new pixel is one of cement clinker phases */
                                else if ( ( ( IsSolidPhase(mic_CSH [ qx ] [ qy ] [ qz ]) == 1 ) &&
                                           ( mic_CSH [ qx ] [ qy ] [ qz ] != C3S ) &&
                                           ( mic_CSH [ qx ] [ qy ] [ qz ] != C2S ) &&
                                           ( mic_CSH [ qx ] [ qy ] [ qz ] != C3A ) &&
                                           ( mic_CSH [ qx ] [ qy ] [ qz ] != C4AF ) ) &&
                                         ( ( newmat [ px ] [ py ] [ pz ] == C3S ) ||
                                          ( newmat [ px ] [ py ] [ pz ] == C2S ) ||
                                          ( newmat [ px ] [ py ] [ pz ] == C3A ) ||
                                          ( newmat [ px ] [ py ] [ pz ] == C4AF ) ) ) {
                                    phase_temp [ ( int ) newmat [ px ] [ py ] [ pz ] ]++;
                                    WriteUnsortedList(px, py, pz);
                                    newmat [ px ] [ py ] [ pz ] = BURNT;
                                    nnew += 1;
                                    if ( nnew >= SIZESET ) {
                                        printf("error in size of nnew %ld\n", nnew);
                                    }

                                    nnewx [ nnew ] = x1;
                                    nnewy [ nnew ] = y1;
                                    nnewz [ nnew ] = z1;
                                }
                                /* 3) old and new pixels belong to one of cement clinker phases and */
                                /* are contained in the same initial cement particle */
                                else if ( ( micpart [ qx ] [ qy ] [ qz ] == micpart [ px ] [ py ] [ pz ] ) &&
                                         ( ( newmat [ px ] [ py ] [ pz ] == C3S ) ||
                                          ( newmat [ px ] [ py ] [ pz ] == C2S ) ||
                                          ( newmat [ px ] [ py ] [ pz ] == C3A ) ||
                                          ( newmat [ px ] [ py ] [ pz ] == C4AF ) ) &&
                                         ( ( mic_CSH [ qx ] [ qy ] [ qz ] == C3S ) ||
                                          ( mic_CSH [ qx ] [ qy ] [ qz ] == C2S ) ||
                                          ( mic_CSH [ qx ] [ qy ] [ qz ] == C3A ) ||
                                          ( mic_CSH [ qx ] [ qy ] [ qz ] == C4AF ) ) ) {
                                    //          ntot+=1;
                                    phase_temp [ ( int ) newmat [ px ] [ py ] [ pz ] ]++;
                                    WriteUnsortedList(px, py, pz);
                                    newmat [ px ] [ py ] [ pz ] = BURNT;
                                    nnew += 1;
                                    if ( nnew >= SIZESET ) {
                                        printf("error in size of nnew %ld\n", nnew);
                                    }

                                    nnewx [ nnew ] = x1;
                                    nnewy [ nnew ] = y1;
                                    nnewz [ nnew ] = z1;
                                }
                            }

                            /* nonperiodic if delimiter */
                        }

                        /* neighbors loop */
                    }

                    /* propagators loop */
                    if ( nnew > 0 ) {
                        ncur = nnew;
                        /* update the burn front matrices */
                        for ( icur = 1; icur <= ncur; icur++ ) {
                            nmatx [ icur ] = nnewx [ icur ];
                            nmaty [ icur ] = nnewy [ icur ];
                            nmatz [ icur ] = nnewz [ icur ];
                        }
                    }
                } while ( nnew > 0 );

                //    ntop+=ntot;
                xl = 0;
                xh = SYSIZE - 1;
                /* Check for percolated path through system */
                for ( j1 = 0; j1 < SYSIZE; j1++ ) {
                    for ( k1 = 0; k1 < SYSIZE; k1++ ) {
                        px = cx(xl, j1, k1, d1, d2, d3);
                        py = cy(xl, j1, k1, d1, d2, d3);
                        pz = cz(xl, j1, k1, d1, d2, d3);
                        qx = cx(xh, j1, k1, d1, d2, d3);
                        qy = cy(xh, j1, k1, d1, d2, d3);
                        qz = cz(xh, j1, k1, d1, d2, d3);
                        if ( ( newmat [ px ] [ py ] [ pz ] == BURNT ) && ( newmat [ qx ] [ qy ] [ qz ] == BURNT ) ) {
                            igood = 2;
                        }

                        if ( newmat [ px ] [ py ] [ pz ] == BURNT ) {
                            newmat [ px ] [ py ] [ pz ] = BURNT + 1;
                        }

                        if ( newmat [ qx ] [ qy ] [ qz ] == BURNT ) {
                            newmat [ qx ] [ qy ] [ qz ] = BURNT + 1;
                        }
                    }
                }

                if ( igood == 2 ) {
                    //      nthrough+=ntot;

                    while ( last != NULL ) { //go through all voxel coordinates
                        ArrPerc [ last->x ] [ last->y ] [ last->z ] = mic_CSH [ last->x ] [ last->y ] [ last->z ];
                        //    printf("AAA%d %d %d\n", last->x, last->y, last->z);
                        current = last;
                        last = last->prev;
                        /*deallocate memory*/
                        delete current;
                    }


                    for ( kh = 0; kh < 51; kh++ ) {
                        phase [ kh ] += phase_temp [ kh ];
                    }
                }
            }

            //if(IsSolidPhase(newmat [px] [py] [pz])==1)
        }

        //loop j
    }

    //loop k

    //GenerateConnNumbers();//only for FEM analysis

    //output only solid phases in percolation path
#ifdef IMAGEFILES
    outputImageFilePerc();
#endif

    cnt_perc = 0;
    for ( kh = 0; kh < 51; kh++ ) {
        cnt_perc += phase [ kh ];
    }

#ifdef OUTFILES
    fprintf(perc_phases, "%d %f %f %d ", cyccnt, time_cur, alpha_cur, cnt_perc);
#endif
    cnt_tot = 0;
    phase [ EMPTYP ] = count [ EMPTYP ];

    /* calculate sum of all solid phases*/
    for ( kh = 1; kh <= HDCSH; kh++ ) {
        if ( ( kh >= DIFFCSH && kh <= EMPTYP ) || kh == HDCSH ) {
            continue;
        }

        cnt_tot += count [ kh ];
    }

#ifdef OUTFILES
    fprintf(perc_phases, "%d %f | ", cnt_tot, ( double ) cnt_perc / cnt_tot);
    for ( kh = 0; kh <= HDCSH; kh++ ) {
        if ( kh > ABSGYP && kh < EMPTYP ) {
            continue;
        }

        //assign EMPTYP to percolated path (used in homogenization)
        phase [ EMPTYP ] = count [ EMPTYP ];
        fprintf(perc_phases, "%ld %ld  ",
                phase [ kh ], kh == CSH ? ( count [ CSH ] - count [ HDCSH ] ) : count [ kh ]);
    }

    fprintf(perc_phases, "\n");
    fflush(perc_phases);
#endif

    dealloc_char_3D(newmat, SYSIZE);
    delete []  nmatx;
    delete []  nmaty;
    delete []  nmatz;
    delete []  nnewx;
    delete []  nnewy;
    delete []  nnewz;
}


/*Check whether two adjacent voxels by side are/aren't connected
 * i.e. if they are in the percolated pathway.
 * Return function will always terminate loops.
 * Return 0 if one of phases is not solid.
 * Return 1 if phases are solid and are disconnected.
 * Return 2 if phases are solid and connected.
 */

int CemhydMatStatus :: IsConnected(int cx, int cy, int cz, int dx, int dy, int dz) {
    int CentPhase, NeighPhase;

    CentPhase = ArrPerc [ AdjCoord(cx) ] [ AdjCoord(cy) ] [ AdjCoord(cz) ];
    if ( IsSolidPhase(CentPhase) ) { //is non-zero
        NeighPhase = ArrPerc [ AdjCoord(cx + dx) ] [ AdjCoord(cy + dy) ] [ AdjCoord(cz + dz) ];

        //skip non-solid phases
        if ( IsSolidPhase(NeighPhase) ) { //is non-zero
            /* 1) old voxel is any solid, new voxel is any
             * solid except clinker phases */
            if ( ( NeighPhase != C3S ) &&
                NeighPhase != C2S &&
                NeighPhase != C3A &&
                NeighPhase != C4AF ) {
                return 2;
            }
            /*2) old voxel is solid except clinker and new voxel is one
             * of cement clinker phases */
            else if ( ( CentPhase != C3S &&
                        CentPhase != C2S &&
                        CentPhase != C3A &&
                        CentPhase != C4AF ) &&
                     ( NeighPhase == C3S ||
                       NeighPhase == C2S ||
                       NeighPhase == C3A ||
                       NeighPhase == C4AF ) ) {
                return 2;
            }
            /* 3) old and new voxels belong to one of cement clinker phases and
             * are contained in the same initial cement particle
             * one-voxel particles have also non-zero ID number*/
            else if ( ( micpart [ AdjCoord(cx) ] [ AdjCoord(cy) ] [ AdjCoord(cz) ] ==
                        micpart [ AdjCoord(cx + dx) ] [ AdjCoord(cy + dy) ] [ AdjCoord(cz + dz) ] ) &&
                     ( CentPhase == C3S ||
                       CentPhase == C2S ||
                       CentPhase == C3A ||
                       CentPhase == C4AF ) &&
                     ( NeighPhase == C3S ||
                       NeighPhase == C2S ||
                       NeighPhase == C3A ||
                       NeighPhase == C4AF ) ) {
                return 2;
            } else {
                return 1;
            }
        }

        //IsSolidPhase(NeighPhase)
        return 0;
    }

    //IsSolidPhase(CentPhase)
    return 0;
}



/*Routine to perform check of connected and disconnected phases
 * on percolated microstructure stored in ArrPerc[dx][dy][dz].
 * Can not be executed during percolation algorithm since whole microstructure
 * is unknown.
 * Each vertex has possible 8x multiple node.
 * Voxel coordinates are the closest to the vertex coordinates
 * in the positive octant.
 * Disconnected node is added to the array as the sum of vertexes
 * by the following scheme (Big Voxel):
 *
 * ---------4096
 * /.       /|
 * / .      32|
 * /  . 4   /  256
 * /   .    / 2 |
 **|--16---512  |
 |   .....|...2048
 |  .     8   /
 | .  1   | 128
 |||.       | /
 * ----64---1024
 *
 * e.g. if number 1 is present, that means that solid voxel below is disconnected
 * number 1024 means that voxel in front, right and down is solid and disconnected
 *
 * z^
 |
 |
 * ----->y
 * /
 * /x
 * coordinates transformation is irrelevant, used x,y,z loop
 */
void CemhydMatStatus :: GenerateConnNumbers(void) {
    int cx, cy, cz;
    int CentPhase;


    for ( cz = 0; cz < SYSIZE; cz++ ) {
        for ( cy = 0; cy < SYSIZE; cy++ ) {
            for ( cx = 0; cx < SYSIZE; cx++ ) {
                //set zero values
                ConnNumbers [ cx ] [ cy ] [ cz ] = 0;

                //if voxel is any solid phase
                CentPhase = ArrPerc [ cx ] [ cy ] [ cz ];
                if ( IsSolidPhase(CentPhase) ) { //is non-zero value
                    /*Each vertex has 7 surrouning neighbors, go through them
                     * Return 0 if one of phases is not solid.
                     * Return 1 if phases are solid and are disconnected.
                     * Return 2 if phases are solid and connected.*/

                    //(1) x+1
                    if ( IsConnected(cx, cy, cz, 1, 0, 0) == 1 ) {
                        ConnNumbers [ cx ] [ cy ] [ cz ] += 1;
                    }

                    //(2) y+1
                    if ( IsConnected(cx, cy, cz, 0, 1, 0) == 1 ) {
                        ConnNumbers [ cx ] [ cy ] [ cz ] += 2;
                    }

                    //(4) z+1
                    if ( IsConnected(cx, cy, cz, 0, 0, 1) == 1 ) {
                        ConnNumbers [ cx ] [ cy ] [ cz ] += 4;
                    }

                    /*(8) x+1, y+1, 2 possible ways of connection
                     * voxel may not be connected to master (CentPhase), i.e.
                     * not connected in first way && not connected in second way
                     */
                    if ( IsSolidPhase(ArrPerc [ AdjCoord(cx + 1) ] [ AdjCoord(cy + 1) ] [ AdjCoord(cz) ]) ) {
                        if ( ( IsConnected(cx, cy, cz, 1, 0, 0) != 2 ||
                               IsConnected(cx + 1, cy, cz, 0, 1, 0) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz, 1, 0, 0) != 2 ) ) {
                            ConnNumbers [ cx ] [ cy ] [ cz ] += 8;
                        }
                    }

                    /*(16) x+1, z+1, 2 possible ways of connection
                     * voxel may not be connected to master (CentPhase), i.e.
                     * not connected in first way && not connected in second way
                     */
                    if ( IsSolidPhase(ArrPerc [ AdjCoord(cx + 1) ] [ AdjCoord(cy) ] [ AdjCoord(cz + 1) ]) ) {
                        if ( ( IsConnected(cx, cy, cz, 1, 0, 0) != 2 ||
                               IsConnected(cx + 1, cy, cz, 0, 0, 1) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 0, 1) != 2 ||
                              IsConnected(cx, cy, cz + 1, 1, 0, 0) != 2 ) ) {
                            ConnNumbers [ cx ] [ cy ] [ cz ] += 16;
                        }
                    }

                    /*(32) y+1, z+1, 2 possible ways of connection
                     * voxel may not be connected to master (CentPhase), i.e.
                     * not connected in first way && not connected in second way
                     */
                    if ( IsSolidPhase(ArrPerc [ AdjCoord(cx) ] [ AdjCoord(cy + 1) ] [ AdjCoord(cz + 1) ]) ) {
                        if ( ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                               IsConnected(cx, cy + 1, cz, 0, 0, 1) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 0, 1) != 2 ||
                              IsConnected(cx, cy, cz + 1, 0, 1, 0) != 2 ) ) {
                            ConnNumbers [ cx ] [ cy ] [ cz ] += 32;
                        }
                    }

                    /*(64) x+1, z-1, 2 possible ways of connection
                     * AdjCoord MUST be present at z direction
                     * voxel may not be connected to master (CentPhase), i.e.
                     * not connected in first way && not connected in second way
                     */
                    if ( IsSolidPhase(ArrPerc [ AdjCoord(cx + 1) ] [ AdjCoord(cy) ] [ AdjCoord(cz - 1) ]) ) {
                        if ( ( IsConnected(cx, cy, cz, 1, 0, 0) != 2 ||
                               IsConnected(cx + 1, cy, cz, 0, 0, -1) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 0, -1) != 2 ||
                              IsConnected(cx, cy, cz - 1, 1, 0, 0) != 2 ) ) {
                            ConnNumbers [ cx ] [ cy ] [ cz ] += 64;
                        }
                    }

                    /*(128) y+1, z-1, 2 possible ways of connection
                     * AdjCoord MUST be present at z direction
                     * voxel may not be connected to master (CentPhase), i.e.
                     * not connected in first way && not connected in second way
                     */
                    if ( IsSolidPhase(ArrPerc [ AdjCoord(cx) ] [ AdjCoord(cy + 1) ] [ AdjCoord(cz - 1) ]) ) {
                        if ( ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                               IsConnected(cx, cy + 1, cz, 0, 0, -1) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 0, -1) != 2 ||
                              IsConnected(cx, cy, cz - 1, 0, 1, 0) != 2 ) ) {
                            ConnNumbers [ cx ] [ cy ] [ cz ] += 128;
                        }
                    }


                    /*(256) x-1, y+1,  2 possible ways of connection
                     * AdjCoord MUST be present at x direction
                     * voxel may not be connected to master (CentPhase), i.e.
                     * not connected in first way && not connected in second way
                     */
                    if ( IsSolidPhase(ArrPerc [ AdjCoord(cx - 1) ] [ AdjCoord(cy + 1) ] [ AdjCoord(cz) ]) ) {
                        if ( ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                               IsConnected(cx, cy + 1, cz, -1, 0, 0) != 2 ) &&
                            ( IsConnected(cx, cy, cz, -1, 0, 0) != 2 ||
                              IsConnected(cx - 1, cy, cz, 0, 1, 0) != 2 ) ) {
                            ConnNumbers [ cx ] [ cy ] [ cz ] += 256;
                        }
                    }

                    /*(512) x+1, y+1, z+1, 6 possible ways of connection
                     * voxel may not be connected to master (CentPhase), i.e.
                     * not connected in first way && not connected in second way && not in
                     * third way
                     */
                    if ( IsSolidPhase(ArrPerc [ AdjCoord(cx + 1) ] [ AdjCoord(cy + 1) ] [ AdjCoord(cz + 1) ]) ) {
                        //        printf("%d %d %d\n", IsConnected(cx, cy, cz, 1, 0, 0), IsConnected(cx+1, cy, cz, 1, 0, 1), IsConnected(cx+1, cy, cz+1, 1, 1, 1));

                        if ( ( IsConnected(cx, cy, cz, 1, 0, 0) != 2 ||
                               IsConnected(cx + 1, cy, cz, 0, 0, 1) != 2 ||
                               IsConnected(cx + 1, cy, cz + 1, 0, 1, 0) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 1, 0, 0) != 2 ||
                              IsConnected(cx + 1, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx + 1, cy + 1, cz, 0, 0, 1) != 2 ) &&

                            ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz, 1, 0, 0) != 2 ||
                              IsConnected(cx + 1, cy + 1, cz, 0, 0, 1) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz, 0, 0, 1) != 2 ||
                              IsConnected(cx, cy + 1, cz + 1, 1, 0, 0) != 2 ) &&

                            ( IsConnected(cx, cy, cz, 0, 0, 1) != 2 ||
                              IsConnected(cx, cy, cz + 1, 1, 0, 0) != 2 ||
                              IsConnected(cx + 1, cy, cz + 1, 0, 1, 0) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 0, 1) != 2 ||
                              IsConnected(cx, cy, cz + 1, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz + 1, 1, 0, 0) != 2 ) ) {
                            ConnNumbers [ cx ] [ cy ] [ cz ] += 512;
                        }
                    }


                    /*(1024) x+1, y+1, z-1, 6 possible ways of connection
                     * voxel may not be connected to master (CentPhase), i.e.
                     * not connected in first way && not connected in second way && not in
                     * third way
                     */
                    if ( IsSolidPhase(ArrPerc [ AdjCoord(cx + 1) ] [ AdjCoord(cy + 1) ] [ AdjCoord(cz - 1) ]) ) {
                        if ( ( IsConnected(cx, cy, cz, 1, 0, 0) != 2 ||
                               IsConnected(cx + 1, cy, cz, 0, 0, -1) != 2 ||
                               IsConnected(cx + 1, cy, cz - 1, 0, 1, 0) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 1, 0, 0) != 2 ||
                              IsConnected(cx + 1, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx + 1, cy + 1, cz, 0, 0, -1) != 2 ) &&

                            ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz, 1, 0, 0) != 2 ||
                              IsConnected(cx + 1, cy + 1, cz, 0, 0, -1) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz, 0, 0, -1) != 2 ||
                              IsConnected(cx, cy + 1, cz - 1, 1, 0, 0) != 2 ) &&

                            ( IsConnected(cx, cy, cz, 0, 0, -1) != 2 ||
                              IsConnected(cx, cy, cz - 1, 1, 0, 0) != 2 ||
                              IsConnected(cx + 1, cy, cz - 1, 0, 1, 0) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 0, -1) != 2 ||
                              IsConnected(cx, cy, cz - 1, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz - 1, 1, 0, 0) != 2 ) ) {
                            ConnNumbers [ cx ] [ cy ] [ cz ] += 1024;
                        }
                    }


                    /*(2048) x-1, y+1, z-1, 6 possible ways of connection
                     * voxel may not be connected to master (CentPhase), i.e.
                     * not connected in first way && not connected in second way && not in
                     * third way
                     */
                    if ( IsSolidPhase(ArrPerc [ AdjCoord(cx - 1) ] [ AdjCoord(cy + 1) ] [ AdjCoord(cz - 1) ]) ) {
                        if ( ( IsConnected(cx, cy, cz, -1, 0, 0) != 2 ||
                               IsConnected(cx - 1, cy, cz, 0, 0, -1) != 2 ||
                               IsConnected(cx - 1, cy, cz - 1, 0, 1, 0) != 2 ) &&
                            ( IsConnected(cx, cy, cz, -1, 0, 0) != 2 ||
                              IsConnected(cx - 1, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx - 1, cy + 1, cz, 0, 0, -1) != 2 ) &&

                            ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz, -1, 0, 0) != 2 ||
                              IsConnected(cx - 1, cy + 1, cz, 0, 0, -1) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz, 0, 0, -1) != 2 ||
                              IsConnected(cx, cy + 1, cz - 1, -1, 0, 0) != 2 ) &&

                            ( IsConnected(cx, cy, cz, 0, 0, -1) != 2 ||
                              IsConnected(cx, cy, cz - 1, -1, 0, 0) != 2 ||
                              IsConnected(cx - 1, cy, cz - 1, 0, 1, 0) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 0, -1) != 2 ||
                              IsConnected(cx, cy, cz - 1, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz - 1, -1, 0, 0) != 2 ) ) {
                            ConnNumbers [ cx ] [ cy ] [ cz ] += 2048;
                        }
                    }

                    /*(4096) x-1, y+1, z+1, 6 possible ways of connection
                     * voxel may not be connected to master (CentPhase), i.e.
                     * not connected in first way && not connected in second way && not in
                     * third way
                     */
                    if ( IsSolidPhase(ArrPerc [ AdjCoord(cx - 1) ] [ AdjCoord(cy + 1) ] [ AdjCoord(cz + 1) ]) ) {
                        if ( ( IsConnected(cx, cy, cz, -1, 0, 0) != 2 ||
                               IsConnected(cx - 1, cy, cz, 0, 0, 1) != 2 ||
                               IsConnected(cx - 1, cy, cz + 1, 0, 1, 0) != 2 ) &&
                            ( IsConnected(cx, cy, cz, -1, 0, 0) != 2 ||
                              IsConnected(cx - 1, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx - 1, cy + 1, cz, 0, 0, 1) != 2 ) &&

                            ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz, -1, 0, 0) != 2 ||
                              IsConnected(cx - 1, cy + 1, cz, 0, 0, 1) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz, 0, 0, 1) != 2 ||
                              IsConnected(cx, cy + 1, cz + 1, -1, 0, 0) != 2 ) &&

                            ( IsConnected(cx, cy, cz, 0, 0, 1) != 2 ||
                              IsConnected(cx, cy, cz + 1, -1, 0, 0) != 2 ||
                              IsConnected(cx - 1, cy, cz + 1, 0, 1, 0) != 2 ) &&
                            ( IsConnected(cx, cy, cz, 0, 0, 1) != 2 ||
                              IsConnected(cx, cy, cz + 1, 0, 1, 0) != 2 ||
                              IsConnected(cx, cy + 1, cz + 1, -1, 0, 0) != 2 ) ) {
                            ConnNumbers [ cx ] [ cy ] [ cz ] += 4096;
                        }
                    }
                }

                //IsSolidPhase(CentPhase)
            }

            //loop cx
        }

        //loop cy
    }

    //loop cz
}

void CemhydMatStatus :: outputImageFilePerc(void) {
    FILE *perc_img;
    char extension [ 10 ];
    //char outputname[80];
    char *prefix;
    prefix = ( char * ) malloc(80);

    //mkdir("perc", 0777);//make directory (every time)
    //system("mk perc");
    //system("mkdir perc 2> /dev/null");
    fprintf(infoperc, "%d %.4f %.3f\n", icyc, alpha_cur, time_cur);
    fflush(infoperc);
    sprintf(extension, "%04d", icyc);
    strcpy(prefix, "perc/out5."); //see line with mkdir in order to be the same
    strcat(prefix, extension);
    strcat(prefix, ".p.img");
#ifdef PRINTF
    printf("Name of percolated output file is %s\n", prefix);
#endif


    if ( ( perc_img = fopen(prefix, "w") ) == NULL ) {
        printf("\nFile %s can not be opened\n", prefix);
        free(prefix);
        return;
    }

    for ( int dz = 0; dz < SYSIZE; dz++ ) {
        for ( int dy = 0; dy < SYSIZE; dy++ ) {
            for ( int dx = 0; dx < SYSIZE; dx++ ) {
                //fprintf(perc_img, "%d %d\n", ArrPerc[dx][dy][dz], ConnNumbers[dx][dy][dz]);
                fprintf(perc_img, "%d\n", ArrPerc [ dx ] [ dy ] [ dz ]);
            }
        }
    }

    fclose(perc_img);
#ifdef PRINTF
    printf("Percolated file %s wrote\n", prefix);
#endif
    free(prefix);
}



void CemhydMatStatus :: WriteUnsortedList(int px, int py, int pz) {
    current = new percolatedpath;
    if ( last != NULL ) {
        last->next = current;
    } else {
        current->prev = NULL;
    }

    current->x = px;
    current->y = py;
    current->z = pz;
    current->prev = last;
    last = current;
}



inline int CemhydMatStatus :: AdjCoord(int coord) {
    if ( coord < 0 ) {
        coord += SYSIZE;
    }

    if ( coord >= SYSIZE ) {
        coord -= SYSIZE;
    }

    return coord;
}


int CemhydMatStatus :: NumSol(int cx, int cy, int cz) {
    int cnt = 0;
    int *p_arr;

    /*check if box is eligible for CSH transformation*/
    for ( int dx = -BoxSize; dx <= BoxSize; dx++ ) {
        for ( int dy = -BoxSize; dy <= BoxSize; dy++ ) {
            for ( int dz = -BoxSize; dz <= BoxSize; dz++ ) {
                if ( mic_CSH [ AdjCoord(cx + dx) ] [ AdjCoord(cy + dy) ] [ AdjCoord(cz + dz) ] != POROSITY  && mic_CSH [ AdjCoord(cx + dx) ] [ AdjCoord(cy + dy) ] [ AdjCoord(cz + dz) ] != EMPTYP ) {
                    cnt++;
                }
            }
        }
    }

    /*update count*/
    CSH_vicinity [ cnt ]++;


    /*if enough phases found, change all CSH in the box to HDCSH*/
    if ( cnt >= SolidLimit ) {
        cnt = 0;

        for ( int dx = -BoxSize; dx <= BoxSize; dx++ ) {
            for ( int dy = -BoxSize; dy <= BoxSize; dy++ ) {
                for ( int dz = -BoxSize; dz <= BoxSize; dz++ ) {
                    p_arr = & mic_CSH [ AdjCoord(cx + dx) ] [ AdjCoord(cy + dy) ] [ AdjCoord(cz + dz) ];
                    //printf("%d at %d %d %d \n", *p_arr, dx, dy, dz);
                    if ( * p_arr == CSH ) {
                        * p_arr = HDCSH;
                        cnt++;
                        //   printf("ph %d\n", mic_CSH[AdjCoord(cx+dx)][AdjCoord(cy+dy)][AdjCoord(cz+dz)]);
                    }
                }
            }
        }
    } else {
        cnt = 0;
    }

    return cnt;
}

/*use pointer to array as argument*/
/*CSH_vicinity[(2*BoxSize+1)*(2*BoxSize+1)*(2*BoxSize+1)]*/

void CemhydMatStatus :: CSHbox(unsigned int *CSH_vicinity) {
    int TotPhase = 0;
    int cx, cy, cz;

    for ( int i = 0; i <= ( ( 2 * BoxSize + 1 ) * ( 2 * BoxSize + 1 ) * ( 2 * BoxSize + 1 ) ); i++ ) {
        CSH_vicinity [ i ] = 0;
    }

    /*make copy of microstructure*/
    for ( cx = 0; cx < SYSIZE; cx++ ) {
        for ( cy = 0; cy < SYSIZE; cy++ ) {
            for ( cz = 0; cz < SYSIZE; cz++ ) {
                mic_CSH [ cx ] [ cy ] [ cz ] = mic [ cx ] [ cy ] [ cz ];
            }
        }
    }


    /*routine to scan microstructure and find eligible voxels*/
    for ( cx = 0; cx < SYSIZE; cx++ ) {
        for ( cy = 0; cy < SYSIZE; cy++ ) {
            for ( cz = 0; cz < SYSIZE; cz++ ) {
                /*CSH may be already in the CSHbox*/
                if ( mic_CSH [ cx ] [ cy ] [ cz ] == CSH || mic_CSH [ cx ] [ cy ] [ cz ] == HDCSH ) {
                    TotPhase += NumSol(cx, cy, cz);
                }
            }
        }
    }

    //  printf("TotPhase:%d", TotPhase);
    count [ HDCSH ] = TotPhase;
}



void CemhydMatStatus :: nrerror(const char *error_text) {
    printf("\nNumerical Recipes run-time error...\n");
    printf("%s\n", error_text);
    printf("...now exiting to system...\n");
    exit(1);
}


float *CemhydMatStatus :: vector(int nl, int nh) {
    float *v;

    v = ( float * ) malloc( ( unsigned ) ( nh - nl + 1 ) * sizeof( float ) );
    if ( !v ) {
        nrerror("allocation failure in vector()");
    }

    return v - nl;
}

int *CemhydMatStatus :: ivector(int nl, int nh) {
    int *v;

    v = ( int * ) malloc( ( unsigned ) ( nh - nl + 1 ) * sizeof( int ) );
    if ( !v ) {
        nrerror("allocation failure in ivector()");
    }

    return v - nl;
}

double *CemhydMatStatus :: dvector(int nl, int nh) {
    double *v;

    v = ( double * ) malloc( ( unsigned ) ( nh - nl + 1 ) * sizeof( double ) );
    if ( !v ) {
        nrerror("allocation failure in dvector()");
    }

    return v - nl;
}



//float** CemhydMat::matrix(int nrl,int nrh,int ncl,int nch){
float **CemhydMatStatus :: matrix_cem(int nrl, int nrh, int ncl, int nch) {
    int i;
    float **m;

    m = ( float ** ) malloc( ( unsigned ) ( nrh - nrl + 1 ) * sizeof( float * ) );
    if ( !m ) {
        nrerror("allocation failure 1 in matrix_cem()");
    }

    m -= nrl;

    for ( i = nrl; i <= nrh; i++ ) {
        m [ i ] = ( float * ) malloc( ( unsigned ) ( nch - ncl + 1 ) * sizeof( float ) );
        if ( !m [ i ] ) {
            nrerror("allocation failure 2 in matrix_cem()");
        }

        m [ i ] -= ncl;
    }

    return m;
}

double **CemhydMatStatus :: dmatrix(int nrl, int nrh, int ncl, int nch) {
    int i;
    double **m;

    m = ( double ** ) malloc( ( unsigned ) ( nrh - nrl + 1 ) * sizeof( double * ) );
    if ( !m ) {
        nrerror("allocation failure 1 in dmatrix()");
    }

    m -= nrl;

    for ( i = nrl; i <= nrh; i++ ) {
        m [ i ] = ( double * ) malloc( ( unsigned ) ( nch - ncl + 1 ) * sizeof( double ) );
        if ( !m [ i ] ) {
            nrerror("allocation failure 2 in dmatrix()");
        }

        m [ i ] -= ncl;
    }

    return m;
}

int **CemhydMatStatus :: imatrix(int nrl, int nrh, int ncl, int nch) {
    int i, **m;

    m = ( int ** ) malloc( ( unsigned ) ( nrh - nrl + 1 ) * sizeof( int * ) );
    if ( !m ) {
        nrerror("allocation failure 1 in imatrix()");
    }

    m -= nrl;

    for ( i = nrl; i <= nrh; i++ ) {
        m [ i ] = ( int * ) malloc( ( unsigned ) ( nch - ncl + 1 ) * sizeof( int ) );
        if ( !m [ i ] ) {
            nrerror("allocation failure 2 in imatrix()");
        }

        m [ i ] -= ncl;
    }

    return m;
}



float **CemhydMatStatus :: submatrix(float **a, int oldrl, int oldrh, int oldcl, int, int newrl, int newcl) {
    int i, j;
    float **m;

    m = ( float ** ) malloc( ( unsigned ) ( oldrh - oldrl + 1 ) * sizeof( float * ) );
    if ( !m ) {
        nrerror("allocation failure in submatrix()");
    }

    m -= newrl;

    for ( i = oldrl, j = newrl; i <= oldrh; i++, j++ ) {
        m [ j ] = a [ i ] + oldcl - newcl;
    }

    return m;
}



void free_vector(float *v, int nl) {
    free( ( char * ) ( v + nl ) );
}

void free_ivector(int *v, int nl) {
    free( ( char * ) ( v + nl ) );
}

void free_dvector(double *v, int nl) {
    free( ( char * ) ( v + nl ) );
}



void free_matrix(float **m, int nrl, int nrh, int ncl) {
    int i;

    for ( i = nrh; i >= nrl; i-- ) {
        free( ( char * ) ( m [ i ] + ncl ) );
    }

    free( ( char * ) ( m + nrl ) );
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl) {
    int i;

    for ( i = nrh; i >= nrl; i-- ) {
        free( ( char * ) ( m [ i ] + ncl ) );
    }

    free( ( char * ) ( m + nrl ) );
}

void free_imatrix(int **m, int nrl, int nrh, int ncl) {
    int i;

    for ( i = nrh; i >= nrl; i-- ) {
        free( ( char * ) ( m [ i ] + ncl ) );
    }

    free( ( char * ) ( m + nrl ) );
}

void free_submatrix(float *b, int nrl) {
    free( ( char * ) ( b + nrl ) );
}

float **CemhydMatStatus :: convert_matrix(float *a, int nrl, int nrh, int ncl, int nch) {
    int i, j, nrow, ncol;
    float **m;

    nrow = nrh - nrl + 1;
    ncol = nch - ncl + 1;
    m = ( float ** ) malloc( ( unsigned ) ( nrow ) * sizeof( float * ) );
    if ( !m ) {
        nrerror("allocation failure in convert_matrix()");
    }

    m -= nrl;
    for ( i = 0, j = nrl; i <= nrow - 1; i++, j++ ) {
        m [ j ] = a + ncol * i - ncl;
    }

    return m;
}

void free_convert_matrix(float **b, int nrl) {
    free( ( char * ) ( b + nrl ) );
}

fcomplex_cem CemhydMatStatus :: Cadd(fcomplex_cem a, fcomplex_cem b) {
    fcomplex_cem c;
    c.r = a.r + b.r;
    c.i = a.i + b.i;
    return c;
}

fcomplex_cem CemhydMatStatus :: Csub(fcomplex_cem a, fcomplex_cem b) {
    fcomplex_cem c;
    c.r = a.r - b.r;
    c.i = a.i - b.i;
    return c;
}

fcomplex_cem CemhydMatStatus :: Cmul(fcomplex_cem a, fcomplex_cem b) {
    fcomplex_cem c;
    c.r = a.r * b.r - a.i * b.i;
    c.i = a.i * b.r + a.r * b.i;
    return c;
}

fcomplex_cem CemhydMatStatus :: ComplexCemhyd(float re, float im) {
    fcomplex_cem c;
    c.r = re;
    c.i = im;
    return c;
}

fcomplex_cem CemhydMatStatus :: Conjg(fcomplex_cem z) {
    fcomplex_cem c;
    c.r = z.r;
    c.i = -z.i;
    return c;
}

fcomplex_cem CemhydMatStatus :: Cdiv(fcomplex_cem a, fcomplex_cem b) {
    fcomplex_cem c;
    float r, den;
    if ( fabs(b.r) >= fabs(b.i) ) {
        r = b.i / b.r;
        den = b.r + r * b.i;
        c.r = ( a.r + r * a.i ) / den;
        c.i = ( a.i - r * a.r ) / den;
    } else {
        r = b.r / b.i;
        den = b.i + r * b.r;
        c.r = ( a.r * r + a.i ) / den;
        c.i = ( a.i * r - a.r ) / den;
    }

    return c;
}

float CemhydMatStatus :: Cabs(fcomplex_cem z) {
    float x, y, ans, temp;
    x = fabs(z.r);
    y = fabs(z.i);
    if ( x == 0.0 ) {
        ans = y;
    } else if ( y == 0.0 ) {
        ans = x;
    } else if ( x > y ) {
        temp = y / x;
        ans = x * sqrt(1.0 + temp * temp);
    } else {
        temp = x / y;
        ans = y * sqrt(1.0 + temp * temp);
    }

    return ans;
}

fcomplex_cem CemhydMatStatus :: Csqrt(fcomplex_cem z) {
    fcomplex_cem c;
    float x, y, w, r;
    if ( ( z.r == 0.0 ) && ( z.i == 0.0 ) ) {
        c.r = 0.0;
        c.i = 0.0;
        return c;
    } else {
        x = fabs(z.r);
        y = fabs(z.i);
        if ( x >= y ) {
            r = y / x;
            w = sqrt(x) * sqrt( 0.5 * ( 1.0 + sqrt(1.0 + r * r) ) );
        } else {
            r = x / y;
            w = sqrt(y) * sqrt( 0.5 * ( r + sqrt(1.0 + r * r) ) );
        }

        if ( z.r >= 0.0 ) {
            c.r = w;
            c.i = z.i / ( 2.0 * w );
        } else {
            c.i = ( z.i >= 0 ) ? w : -w;
            c.r = z.i / ( 2.0 * c.i );
        }

        return c;
    }
}

fcomplex_cem CemhydMatStatus :: RCmul(float x, fcomplex_cem a) {
    fcomplex_cem c;
    c.r = x * a.r;
    c.i = x * a.i;
    return c;
}


//public function to call the routine for HDCSH creation
void CemhydMatStatus :: CreateHDCSH(void) {
    CSHbox(CSH_vicinity);
}

//public function to call the routine for solid percolation
void CemhydMatStatus :: PercolateForOutput(void) {
    burn_phases(1, 0, 0);
}

//public function to access w/c ratio
double CemhydMatStatus :: GiveWcr(void) {
    return w_to_c;
}

/*function read into the string
 * int rand_seed_num, double fineness, double dihydrate, double hemihydrate, double anhydrite, double C3S_mass, double C2S_mass, double C3A_mass, double C4AF_mass, double wcr, int satur_sealed
 */

void CemhydMatStatus :: GetInputParams(char *my_string) {
    sprintf(my_string, "%d", iseed);
}

long CemhydMatStatus :: cx(int x, int y, int z, int a, int b, int c) {
    return ( 1 - b - c ) * x + ( 1 - a - c ) * y + ( 1 - a - b ) * z;
}


long CemhydMatStatus :: cy(int x, int y, int z, int a, int b, int c) {
    return ( 1 - a - b ) * x + ( 1 - b - c ) * y + ( 1 - a - c ) * z;
}

long CemhydMatStatus :: cz(int x, int y, int z, int a, int b, int c) {
    return ( 1 - a - c ) * x + ( 1 - a - b ) * y + ( 1 - b - c ) * z;
}

//homogenization of plain cement paste without filler
void CemhydMatStatus :: AnalyticHomogenizationPaste(double &E, double &nu, int perc_unperc_flag) {
    //flag 0 - percolated phases are is in phase[]
    //flag 1 - unpercolated phases are in count[]
    int index, x, i;
    double sum = 0.;
    FloatMatrix PhaseMatrix(32, 3), LevelI(2, 3);
    double E_CSH_hmg, nu_CSH_hmg;
    Homogenize CSH_level, Paste_level;


    //determine phase fractions
    //0-23, PhaseFrac[0] = none, [1]=POROSITY, [2]=C3S ..  [31]=EMPTYP [32]=HDCSH
    for ( x = 0; x < 34; x++ ) {
        PhaseFrac [ x ] = 0.;
    }

    for ( x = 0; x < 51; x++ ) {
        index = x;
        if ( x == HDCSH ) {
            index = 31;
        } else if ( x == EMPTYP ) { //EMPTYP
            index = 30;
        } else if ( x >= DIFFCSH ) { //all dissolved phases
            continue;
        }


        if ( perc_unperc_flag == 0 ) { //percolated
            PhaseFrac [ index + 1 ] += phase [ x ];
        } else { //unpercolated
            if ( x == CSH ) {
                PhaseFrac [ index + 1 ] += ( count [ CSH ] - count [ HDCSH ] );
            } else {
                PhaseFrac [ index + 1 ] += count [ x ];
            }
        }
    }

    for ( x = 2; x < 34; x++ ) {
        PhaseFrac [ x ] /= SYSIZE_POW3;
        sum += PhaseFrac [ x ];
    }


    //water-filled porosity is what is missing
    PhaseFrac [ 1 ] = 1 - sum;

    for ( x = 0; x < 34; x++ ) { //34
        //printf("Phase %d PhaseFrac %f\n", x, PhaseFrac[x]);
        //printf("Phase %d perc / unper %ld %ld\n", x, phase[x], count[x]);
    }

    PhaseMatrix(0, 0) = PhaseFrac [ 1  ]; //POROSITY    0
    PhaseMatrix(0, 1) = 0.001;
    PhaseMatrix(0, 2) = 0.499924;
    PhaseMatrix(1, 0) = PhaseFrac [ 2  ]; //C3S         1
    PhaseMatrix(1, 1) = 135.0;
    PhaseMatrix(1, 2) = 0.300000;
    PhaseMatrix(2, 0) = PhaseFrac [ 3  ]; //C2S         2
    PhaseMatrix(2, 1) = 130.0;
    PhaseMatrix(2, 2) = 0.300000;
    PhaseMatrix(3, 0) = PhaseFrac [ 4  ]; //C3A         3
    PhaseMatrix(3, 1) = 145.0;
    PhaseMatrix(3, 2) = 0.300000;
    PhaseMatrix(4, 0) = PhaseFrac [ 5  ]; //C4AF        4
    PhaseMatrix(4, 1) = 125.0;
    PhaseMatrix(4, 2) = 0.300000;
    PhaseMatrix(5, 0) = PhaseFrac [ 6  ]; //GYPSUM      5
    PhaseMatrix(5, 1) = 30.00;
    PhaseMatrix(5, 2) = 0.300000;
    PhaseMatrix(6, 0) = PhaseFrac [ 7  ]; //HEMIHYD     6
    PhaseMatrix(6, 1) = 62.90;
    PhaseMatrix(6, 2) = 0.300000;
    PhaseMatrix(7, 0) = PhaseFrac [ 8  ]; //ANHYDRITE   7
    PhaseMatrix(7, 1) = 73.60;
    PhaseMatrix(7, 2) = 0.275000;
    PhaseMatrix(8, 0) = PhaseFrac [ 9  ]; //POZZ,SiO2   8
    PhaseMatrix(8, 1) = 72.80;
    PhaseMatrix(8, 2) = 0.167000;
    PhaseMatrix(9, 0) = PhaseFrac [ 10 ]; //INERT       9
    PhaseMatrix(9, 1) = 79.60;
    PhaseMatrix(9, 2) = 0.309900;
    PhaseMatrix(10, 0) = PhaseFrac [ 11 ]; //SLAG       10
    PhaseMatrix(10, 1) = 72.80;
    PhaseMatrix(10, 2) = 0.16700;
    PhaseMatrix(11, 0) = PhaseFrac [ 12 ]; //ASG        11
    PhaseMatrix(11, 1) = 72.80;
    PhaseMatrix(11, 2) = 0.16700;
    PhaseMatrix(12, 0) = PhaseFrac [ 13 ]; //CAS2 slag  12
    PhaseMatrix(12, 1) = 72.80;
    PhaseMatrix(12, 2) = 0.16700;
    PhaseMatrix(13, 0) = PhaseFrac [ 14 ]; //CH         13
    PhaseMatrix(13, 1) = 38.00;
    PhaseMatrix(13, 2) = 0.30500;
    PhaseMatrix(14, 0) = -10.; //CSH upsc.  14
    PhaseMatrix(14, 1) = -10.0;
    PhaseMatrix(14, 2) = -10.000;
    PhaseMatrix(15, 0) = PhaseFrac [ 16 ]; //C3AH6      15
    PhaseMatrix(15, 1) = 22.40;
    PhaseMatrix(15, 2) = 0.25000;
    PhaseMatrix(16, 0) = PhaseFrac [ 17 ]; //ETTR       16
    PhaseMatrix(16, 1) = 22.40;
    PhaseMatrix(16, 2) = 0.25000;
    PhaseMatrix(17, 0) = PhaseFrac [ 18 ]; //ETTRC4AF   17
    PhaseMatrix(17, 1) = 22.40;
    PhaseMatrix(17, 2) = 0.25000;
    PhaseMatrix(18, 0) = PhaseFrac [ 19 ]; //AFM        18
    PhaseMatrix(18, 1) = 42.30;
    PhaseMatrix(18, 2) = 0.32380;
    PhaseMatrix(19, 0) = PhaseFrac [ 20 ]; //FH3        19
    PhaseMatrix(19, 1) = 22.40;
    PhaseMatrix(19, 2) = 0.24590;
    PhaseMatrix(20, 0) = PhaseFrac [ 21 ]; //POZZCSH    20
    PhaseMatrix(20, 1) = 22.40;
    PhaseMatrix(20, 2) = 0.24590;
    PhaseMatrix(21, 0) = PhaseFrac [ 22 ]; //SLAGCSH    21
    PhaseMatrix(21, 1) = 22.40;
    PhaseMatrix(21, 2) = 0.24590;
    PhaseMatrix(22, 0) = PhaseFrac [ 23 ]; //CACL2      22
    PhaseMatrix(22, 1) = 42.30;
    PhaseMatrix(22, 2) = 0.32380;
    PhaseMatrix(23, 0) = PhaseFrac [ 24 ]; //FREIDEL    23
    PhaseMatrix(23, 1) = 22.40;
    PhaseMatrix(23, 2) = 0.24590;
    PhaseMatrix(24, 0) = PhaseFrac [ 25 ]; //STRAT      24
    PhaseMatrix(24, 1) = 22.40;
    PhaseMatrix(24, 2) = 0.25000;
    PhaseMatrix(25, 0) = PhaseFrac [ 26 ]; //GYPSUMS    25
    PhaseMatrix(25, 1) = 30.00;
    PhaseMatrix(25, 2) = 0.30000;
    PhaseMatrix(26, 0) = PhaseFrac [ 27 ]; //CaCO3      26
    PhaseMatrix(26, 1) = 60.00;
    PhaseMatrix(26, 2) = 0.30000;
    PhaseMatrix(27, 0) = PhaseFrac [ 28 ]; //AFMC       27
    PhaseMatrix(27, 1) = 42.30;
    PhaseMatrix(27, 2) = 0.32380;
    PhaseMatrix(28, 0) = PhaseFrac [ 29 ]; //INERTAGG   28
    PhaseMatrix(28, 1) = 79.60;
    PhaseMatrix(28, 2) = 0.30990;
    PhaseMatrix(29, 0) = PhaseFrac [ 30 ]; //ABSGYP     29
    PhaseMatrix(29, 1) = 30.00;
    PhaseMatrix(29, 2) = 0.30000;
    PhaseMatrix(30, 0) = PhaseFrac [ 31 ]; //EMPTYP     30
    PhaseMatrix(30, 1) = 0.001;
    PhaseMatrix(30, 2) = 0.00100;

    //CSH level
    if ( PhaseFrac [ 15 ] == 0. && PhaseFrac [ 32 ] == 0. ) { //beginning of hydration
        LevelI(0, 0) = 1.;
        LevelI(1, 0) = 0.;
    } else {
        LevelI(0, 0) = PhaseFrac [ 15 ] / ( PhaseFrac [ 15 ] + PhaseFrac [ 32 ] ); //LDCSH
        LevelI(1, 0) = PhaseFrac [ 32 ] / ( PhaseFrac [ 15 ] + PhaseFrac [ 32 ] ); //HDCSH
    }

    LevelI(0, 1) = 21.7; //21.7 LDCSH
    LevelI(0, 2) = 0.24;

    LevelI(1, 1) = 29.4; //29.4 HDCSH
    LevelI(1, 2) = 0.24;

    //constructor without parameters
    CSH_level.moriTanaka(LevelI, 0);

    //fill actual values into PhaseMatrix
    for ( i = 0; i < 31; i++ ) {
        PhaseMatrix(i, 0) =  PhaseFrac [ i + 1 ];
    }

    //cement paste level

    PhaseMatrix(14, 0) = PhaseFrac [ 15 ] + PhaseFrac [ 32 ];
    PhaseMatrix(14, 1) = CSH_level.E_hmg;
    PhaseMatrix(14, 2) = CSH_level.nu_hmg;

    E_CSH_hmg = CSH_level.E_hmg; //E from CSH level
    nu_CSH_hmg = CSH_level.nu_hmg; //nu from CSH level

    Paste_level.selfConsistent(PhaseMatrix);

    E = Paste_level.E_hmg;
    nu =  Paste_level.nu_hmg;
}

//use results from cement paste level and add SCM, inert filler, entrained air, ITZ, FA, CA
void CemhydMatStatus :: AnalyticHomogenizationConcrete(double E_paste_inp, double nu_paste_inp, double *E_paste, double *nu_paste, double *E_mortar, double *nu_mortar, double &E_concrete, double &nu_concrete) {
    Homogenize Paste_level, Mortar_level, Concrete_level;

    //CEMENT PASTE LEVEL - homogenize cement paste with SCM, inert filler, entrained and entrapped air using Mori-Tanaka with cement paste as reference via MT scheme
    double vol_tot_paste = Vol_cement_clinker_gypsum + Vol_cement_SCM + Vol_water + Vol_inert_filler + Vol_entrained_entrapped_air;
    FloatMatrix PhaseMatrix(4, 3), Mortar(4, 3), Concrete(4, 3);

    PhaseMatrix(0, 0) = ( Vol_cement_clinker_gypsum + Vol_water ) / vol_tot_paste; //cement paste (reference medium)
    PhaseMatrix(0, 1) = E_paste_inp;
    PhaseMatrix(0, 2) = nu_paste_inp;
    PhaseMatrix(1, 0) = Vol_cement_SCM / vol_tot_paste; //SCM
    PhaseMatrix(1, 1) = Young_SCM;
    PhaseMatrix(1, 2) = Poisson_SCM;
    PhaseMatrix(2, 0) = Vol_inert_filler / vol_tot_paste; //inert filler
    PhaseMatrix(2, 1) = Young_inert;
    PhaseMatrix(2, 2) = Poisson_inert;
    PhaseMatrix(3, 0) = Vol_entrained_entrapped_air / vol_tot_paste; //entrained + entrapped air
    PhaseMatrix(3, 1) = 0.001;
    PhaseMatrix(3, 2) = 0.001;

    Paste_level.moriTanaka(PhaseMatrix, 0);
    * E_paste = Paste_level.E_hmg;
    * nu_paste = Paste_level.nu_hmg;

    //MORTAR LEVEL - homogenize cement paste with fine aggregates including ITZ via Herve-Zaoui scheme
    //ITZ covers all aggregates, therefore its fraction associated with mortar and concrete has to be determined
    //ITZ occupies additional space
    double n_FA, n_CA, vol_ITZ_FA, vol_ITZ_CA;
    n_FA = Vol_FA / ( 4. / 3. * M_PI * pow(0.01 * Grain_average_FA / 2., 3.) ); //amount of FA particles
    vol_ITZ_FA = n_FA * 4. / 3. * M_PI * ( pow(0.01 * Grain_average_FA / 2. + 0.00001 * ITZ_thickness, 3) - pow(0.01 * Grain_average_FA / 2., 3) ); //volume occupied by ITZ in FA [l]
    n_CA = Vol_CA / ( 4. / 3. * M_PI * pow(0.01 * Grain_average_CA / 2., 3) ); //amount of CA particles
    vol_ITZ_CA = n_CA * 4. / 3. * M_PI * ( pow(0.01 * Grain_average_CA / 2. + 0.00001 * ITZ_thickness, 3) - pow(0.01 * Grain_average_CA / 2., 3) ); //volume occupied by ITZ in CA [l]
    double vol_tot_mortar = vol_tot_paste + Vol_FA;

    Mortar(0, 0) = Vol_FA / vol_tot_mortar; //Fine aggregates
    Mortar(0, 1) = Young_FA;
    Mortar(0, 2) = Poisson_FA;
    Mortar(1, 0) = vol_ITZ_FA / vol_tot_mortar; //ITZ
    Mortar(1, 1) = * E_paste * ITZ_Young_red;
    Mortar(1, 2) = * nu_paste;
    Mortar(2, 0) = ( vol_tot_paste - vol_ITZ_FA ) / vol_tot_mortar; //cement paste
    Mortar(2, 1) = * E_paste;
    Mortar(2, 2) = * nu_paste;
    Mortar(3, 0) = 0.; //vol. fraction must be zero - reference medium
    Mortar(3, 1) = 10.; //arbitrary value
    Mortar(3, 2) = 0.3; //arbitrary value

    Mortar_level.herveZaoui(Mortar);
    * E_mortar = Mortar_level.E_hmg;
    * nu_mortar = Mortar_level.nu_hmg;

    //CONCRETE LEVEL
    double vol_tot_concrete = vol_tot_mortar + Vol_CA;

    Concrete(0, 0) = Vol_CA / vol_tot_concrete; //Coarse aggregates
    Concrete(0, 1) = Young_CA;
    Concrete(0, 2) = Poisson_CA;
    Concrete(1, 0) = vol_ITZ_CA / vol_tot_concrete; //ITZ
    Concrete(1, 1) = * E_paste * ITZ_Young_red;
    Concrete(1, 2) = * nu_paste;
    Concrete(2, 0) = ( vol_tot_mortar - vol_ITZ_CA ) / vol_tot_concrete; //mortar
    Concrete(2, 1) = * E_mortar;
    Concrete(2, 2) = * nu_mortar;
    Concrete(3, 0) = 0.; //vol. fraction must be zero - reference medium
    Concrete(3, 1) = 10.; //arbitrary value
    Concrete(3, 2) = 0.3; //arbitrary value

    Concrete_level.herveZaoui(Concrete);
    E_concrete = Concrete_level.E_hmg;
    nu_concrete = Concrete_level.nu_hmg;
}

//at least one cycle has to pass, volume fraction
void CemhydMatStatus :: GetInitClinkerPhases(double &c3s, double &c2s, double &c3a, double &c4af, double &gypsum, double &hemi, double &anh) {
    double sum = c3sinit + c2sinit + c3ainit + c4afinit + ncsbar + heminit + anhinit;

    if ( icyc > 1 ) {
        c3s = ( double ) c3sinit / sum;
        c2s = ( double ) c2sinit / sum;
        c3a = ( double ) c3ainit / sum;
        c4af = ( double ) c4afinit / sum;
        gypsum = ( double ) ncsbar / sum;
        hemi = ( double ) heminit / sum;
        anh = ( double ) anhinit / sum;
    }
}


#ifdef __TM_MODULE  //transport module for OOFEM
void
CemhydMatStatus :: updateYourself(TimeStep *atTime) {
    TransportMaterialStatus :: updateYourself(atTime);
};

double CemhydMatStatus :: giveAverageTemperature(void) {
    CemhydMat *cemhydmat = ( CemhydMat * ) this->gp->giveMaterial();
    CemhydMatStatus *ms = ( CemhydMatStatus * ) cemhydmat->giveStatus(gp);
    double ret;

    if ( !cemhydmat->eachGP ) {
        ret = cemhydmat->MasterCemhydMatStatus->averageTemperature;
    } else {
        ret = ms->giveStateVector().at(1);
    }

    return ret;
}


void
CemhydMatStatus :: printOutputAt(FILE *file, TimeStep *atTime)
{
    CemhydMat *cemhydmat = ( CemhydMat * ) this->gp->giveMaterial();
    CemhydMatStatus *ms = this;
    if ( cemhydmat->MasterCemhydMatStatus ) {
        ms = cemhydmat->MasterCemhydMatStatus;
    }

    TransportMaterialStatus :: printOutputAt(file, atTime);
    fprintf(file, "   status {");
    fprintf( file, "cyc %d  time %e  DoH %f  conductivity %f  capacity %f  density %f", cemhydmat->giveCycleNumber(this->gp), 3600 * cemhydmat->giveTimeOfCycle(this->gp), cemhydmat->giveDoHActual(this->gp), cemhydmat->giveConcreteConductivity(this->gp), cemhydmat->giveConcreteCapacity(this->gp), cemhydmat->giveConcreteDensity(this->gp) );
    if ( ms->Calculate_elastic_homogenization ) {
        fprintf(file, " EVirginPaste %f NuVirginPaste %f EConcrete %f NuConcrete %f", ms->last_values [ 2 ], ms->last_values [ 3 ], ms->last_values [ 4 ], ms->last_values [ 5 ]);
    }


    if ( !cemhydmat->eachGP ) {
        if ( cemhydmat->giveStatus(gp) == cemhydmat->MasterCemhydMatStatus ) {
            fprintf( file, " master material %d", cemhydmat->giveNumber() );
        } else {
            fprintf( file, " slave of material %d", cemhydmat->giveNumber() );
        }
    } else {
        fprintf( file, " independent microstructure %p from material %d", this, cemhydmat->giveNumber() );
    }

    fprintf(file, "}\n");
}
#endif


#ifdef CEMPY
///constructors
CemhydMatStatus :: CemhydMatStatus()
{
    PartHeat = 0.;
    temp_cur = 0.;
 #ifdef PRINTF
    printf("Constructor of CemhydMatStatus called\n");
    fflush(stdout);
 #endif
}


void
CemhydMatStatus :: InitializePy(const char *inp) {
    initializeMicrostructure();
    readInputFileAndInitialize(inp, ( bool ) 1);
}

#endif //CEMPY
} //end of namespace oofem



#ifdef CEMPY
 #include <boost/python.hpp>
BOOST_PYTHON_MODULE(cemhydmodule)
{
    //boost::python::class_< oofem::CemhydMatStatus >("cemhydmatstatus", init<std::string>() ) - does not work, why ?
    boost :: python :: class_< oofem :: CemhydMatStatus >("cemhydmatstatus") //use default class constructor instead
    .def("InitializePy", & oofem :: CemhydMatStatus :: InitializePy)
    .def("MoveCycles", & oofem :: CemhydMatStatus :: MoveCycles)
    ;
}
#endif
