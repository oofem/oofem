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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef nlisomoisturemat_h
#define nlisomoisturemat_h

#include "tm/Materials/isomoisturemat.h"
#include "floatarray.h"
#include "scalarfunction.h"
#include "function.h"

///@name Input fields for NlIsoMoistureMaterial
//@{
#define _IFT_NlIsoMoistureMaterial_Name "nlisomoisturemat"
#define _IFT_NlIsoMoistureMaterial_isothermtype "isothermtype"
#define _IFT_NlIsoMoistureMaterial_permeabilitytype "permeabilitytype"
#define _IFT_NlIsoMoistureMaterial_capillarytransporttype "capillarytransporttype"
#define _IFT_NlIsoMoistureMaterial_rhodry "rhodry"
#define _IFT_NlIsoMoistureMaterial_capa "capa"
#define _IFT_NlIsoMoistureMaterial_hx "hx"
#define _IFT_NlIsoMoistureMaterial_dx "dx"
#define _IFT_NlIsoMoistureMaterial_iso_offset "isooffset"
#define _IFT_NlIsoMoistureMaterial_iso_h "iso_h"
#define _IFT_NlIsoMoistureMaterial_iso_wh "iso_w(h)"
#define _IFT_NlIsoMoistureMaterial_dd "dd"
#define _IFT_NlIsoMoistureMaterial_wf "wf"
#define _IFT_NlIsoMoistureMaterial_b "b"
#define _IFT_NlIsoMoistureMaterial_uh "uh"
#define _IFT_NlIsoMoistureMaterial_a "a"
#define _IFT_NlIsoMoistureMaterial_nn "nn"
#define _IFT_NlIsoMoistureMaterial_c "c"
#define _IFT_NlIsoMoistureMaterial_k "k"
#define _IFT_NlIsoMoistureMaterial_vm "vm"
#define _IFT_NlIsoMoistureMaterial_perm_h "perm_h"
#define _IFT_NlIsoMoistureMaterial_perm_ch "perm_c(h)"
#define _IFT_NlIsoMoistureMaterial_hc "hc"
#define _IFT_NlIsoMoistureMaterial_alpha0 "alpha0"
#define _IFT_NlIsoMoistureMaterial_c1 "c1"
#define _IFT_NlIsoMoistureMaterial_n "n"
#define _IFT_NlIsoMoistureMaterial_alphah "alphah"
#define _IFT_NlIsoMoistureMaterial_betah "betah"
#define _IFT_NlIsoMoistureMaterial_gammah "gammah"
#define _IFT_NlIsoMoistureMaterial_rhoh2o "rhoh2o"
#define _IFT_NlIsoMoistureMaterial_capperm_h "capperm_h"
#define _IFT_NlIsoMoistureMaterial_capperm_dwh "capperm_dw(h)"
#define _IFT_NlIsoMoistureMaterial_capperm_wv "capperm_wv"
#define _IFT_NlIsoMoistureMaterial_capperm_dwwv "capperm_dw(wv)"
#define _IFT_NlIsoMoistureMaterial_abs "abs"
#define _IFT_NlIsoMoistureMaterial_pl "pl"
#define _IFT_NlIsoMoistureMaterial_mu "mu"
#define _IFT_NlIsoMoistureMaterial_timescale "timescale"
#define _IFT_NlIsoMoistureMaterial_wn "wn"
#define _IFT_NlIsoMoistureMaterial_alpha "alpha"
#define _IFT_NlIsoMoistureMaterial_capil_coef "capil_coef"
#define _IFT_NlIsoMoistureMaterial_t "t"
#define _IFT_NlIsoMoistureMaterial_ttf "ttf"
#define _IFT_NlIsoMoistureMaterial_vg_b "vg_b"
#define _IFT_NlIsoMoistureMaterial_vg_m "vg_m"
//@}

namespace oofem {
/**
 * This class implements various functions for concrete moisture permeability and moisture capacity
 */
class NlIsoMoistureMaterial : public IsotropicMoistureTransferMaterial
{
protected:
    enum isothermType { linear, multilinear, Ricken, Kuenzel, Hansen, BSB, bilinear, vanGenuchten } Isotherm;

    /// density of the dry solid phase
    double rhodry = 0.;

    /// values of the linear isotherm
    double moistureCapacity = 0.;

    /// values of the multilinear isotherm
    FloatArray iso_h;
    FloatArray iso_wh;

    /// parameters of the Ricken isotherm
    double dd = 0.;

    /// parameters of the Kuenzel isotherm
    double wf = 0., b = 0.;

    /// parameters of the isotherm proposed by P. Freiesleben Hansen (Coupled moisture/heat transport in cross sections of structures, Beton og Konstruktionsinstituttet, 1985)
    double uh = 0., A = 0., nn = 0.;

    /// parameters of the BSB isotherm
    double c = 0., k = 0., Vm = 0.;

    /// values of the bilinear isotherm
    double hx = 0., dx = 0.;
    double iso_offset = 0.;
    double c1 = 0., c2 = 0., capa2 = 0.;

    /// parameters of vanGenuchten isotherm
    double vG_b = 0., vG_m = 0.;

    /// Nonevaporable water content per m3 of concrete at complete hydration
    double wn = 0.;

    /// Function of degree of hydration
    ScalarFunction alpha;

    enum permeabilityType { multilin, Bazant, Xi, KunzelPerm } Permeability;

    /// values of the multilinear permeability
    FloatArray perm_h;
    FloatArray perm_ch;

    /// "permeability" according to Bazant
    double C1 = 0., n = 0., alpha0 = 0., hC = 0.;

    /// permeability parameters according to Xi, Bazant & Jennings
    double alphah = 0., betah = 0., gammah = 0.;

    enum capillaryTransportType { Multilin_h, Multilin_wV, KunzelCT } CapillaryTransport;
    /// water absorption coefficient [kg m^-2 s^-0.5]
    double Abs = 0.;

    /// water vapor diffusion resistance [-]
    double mu = 0.;

    /// ambient atmospheric pressure [Pa]
    double PL = 101325.;

    /// = 1 for analysis in seconds, = 86400 for analysis in days, etc.
    double timeScale = 1.;

    /// parameter in liquid conduction
    double capillary_transport_coef = 1000.;

    /// constant temperature [K]
    double T = 0.;
    /// explicitly prescribed evolution of temperature by a time function (e.g. piecewise-linear dfined externally) in [K]
    int T_TF = 0;


    /// values of the multilinear capillary transport function
    FloatArray capPerm_h;
    FloatArray capPerm_Dwh;
    FloatArray capPerm_wV;
    FloatArray capPerm_DwwV;
    double rhoH2O = 0.;

public:
    NlIsoMoistureMaterial(int n, Domain *d) : IsotropicMoistureTransferMaterial(n, d) { }

    void initializeFrom(InputRecord &ir) override;

    /// evaluates slope of the sorption isotherm
    double giveMoistureCapacity(GaussPoint *gp, TimeStep *tStep) const override;
    double giveMoistureContent(double humidity) const override;
    double givePermeability(GaussPoint *gp, TimeStep *tStep) const override;
    double computeCapTranspCoeff(double humidity) const;

    /// compute vapor diffusion coefficient in air [kg m^-1 s^-1 Pa^-1]
    double computeVaporDiffusionCoeff(GaussPoint *gp, TimeStep *tStep) const;
    /// compute saturation water vapor pressure
    double computeSaturationWaterVaporPressure(GaussPoint *gp, TimeStep *tStep) const;
    /// evaluate temperature effect on water viscosity - liquid water capillary conduction
    double computeTemperatureEffectOnViscosity(GaussPoint *gp, TimeStep *tStep) const;
    /// returns temperature in [K]
    double giveTemperature(GaussPoint *gp, TimeStep *tStep) const;

    const char *giveInputRecordName() const override { return _IFT_NlIsoMoistureMaterial_Name; }
    const char *giveClassName() const override { return "NlIsoMoistureMaterial"; }

    double giveHumidity(GaussPoint *gp, ValueModeType mode) const override;


    bool hasInternalSource() const override;
    void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;
};
} // end namespace oofem
#endif // nlisomoisturemat_h
