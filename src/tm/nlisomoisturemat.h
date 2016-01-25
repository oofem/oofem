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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#include "isomoisturemat.h"
#include "floatarray.h"

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
#define _IFT_NlIsoMoistureMaterial_iso_offset "isoOffset"
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
#define _IFT_NlIsoMoistureMaterial_t "t"
#define _IFT_NlIsoMoistureMaterial_timescale "timescale"


//@}

namespace oofem {
/**
 * This class implements various functions for concrete moisture permeability and moisture capacity
 */
class NlIsoMoistureMaterial : public IsotropicMoistureTransferMaterial
{
protected:
  enum isothermType { linear, multilinear, Ricken, Kuenzel, Hansen, BSB, bilinear } Isotherm;

    /// density of the dry solid phase
    double rhodry;

    /// values of the linear isotherm
    double moistureCapacity;

    /// values of the multilinear isotherm
    FloatArray iso_h;
    FloatArray iso_wh;

    /// parameters of the Ricken isotherm
    double dd;

    /// parameters of the Kuenzel isotherm
    double wf, b;

    /// parameters of the isotherm proposed by P. Freiesleben Hansen (Coupled moisture/heat transport in cross sections of structures, Beton og Konstruktionsinstituttet, 1985)
    double uh, A, nn;

    /// parameters of the BSB isotherm
    double c, k, Vm;

    /// values of the bilinear isotherm
    double hx, dx;
    double iso_offset;
    double c1, c2, capa2;


    enum permeabilityType { multilin, Bazant, Xi, KunzelPerm } Permeability;

    /// values of the multilinear permeability
    FloatArray perm_h;
    FloatArray perm_ch;

    /// "permeability" according to Bazant
    double C1, n, alpha0, hC;

    /// permeability parameters according to Xi, Bazant & Jennings
    double alphah, betah, gammah;

    /// permeability parameters according to Kunzel
    double deltap, p_sat; // permation
    enum capillaryTransportType {Multilin_h, Multilin_wV, KunzelCT} CapillaryTransport;
    double Abs;            ///< water absorption coefficient [kg m^-2 s^-0.5]

    /// values of the multilinear capillary transport function
    FloatArray capPerm_h;
    FloatArray capPerm_Dwh;
    FloatArray capPerm_wV;
    FloatArray capPerm_DwwV;
    double rhoH2O;

public:
    NlIsoMoistureMaterial(int n, Domain * d) : IsotropicMoistureTransferMaterial(n, d) { }
    virtual ~NlIsoMoistureMaterial() { }

    /// evaluates slope of the sorption isotherm
    virtual double giveMoistureCapacity(GaussPoint *gp, TimeStep *tStep);
    virtual double sorptionIsotherm(double humidity);
    virtual double givePermeability(GaussPoint *gp, TimeStep *tStep);
    virtual double computeCapTranspCoeff(double humidity);

    virtual const char *giveInputRecordName() const { return _IFT_NlIsoMoistureMaterial_Name; }
    virtual const char *giveClassName() const { return "NlIsoMoistureMaterial"; }

    virtual double giveHumidity(GaussPoint *gp, ValueModeType mode);

    virtual IRResultType initializeFrom(InputRecord *ir);
};
} // end namespace oofem
#endif // nlisomoisturemat_h
