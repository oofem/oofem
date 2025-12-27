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

#include "tm/Materials/nlisomoisturemat.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(NlIsoMoistureMaterial);

void
NlIsoMoistureMaterial::initializeFrom(InputRecord &ir)
{
    IsotropicMoistureTransferMaterial::initializeFrom(ir);

    int type = 0;
    IR_GIVE_FIELD(ir, type, _IFT_NlIsoMoistureMaterial_isothermtype);

    if ( type > 7 ) {
        throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_isothermtype, "must be equal to 0, 1, 2 ... 7");
    }

    this->Isotherm = ( isothermType ) type;

    IR_GIVE_FIELD(ir, type, _IFT_NlIsoMoistureMaterial_permeabilitytype);

    if ( type >= 4 ) {
        throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_permeabilitytype, "must be equal to 0, 1, 2, 3");
    }

    this->Permeability = ( permeabilityType ) type;

    if ( Permeability == KunzelPerm ) {
        IR_GIVE_FIELD(ir, type, _IFT_NlIsoMoistureMaterial_capillarytransporttype);
        this->CapillaryTransport = ( capillaryTransportType ) type;
        if ( type >= 3 ) {
            throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_capillarytransporttype, "must be equal to 0, 1, 2");
        }
    }

    if ( this->Isotherm == linear ) { // linear isotherm = type 0
        IR_GIVE_FIELD(ir, moistureCapacity, _IFT_NlIsoMoistureMaterial_capa);
        this->iso_offset = 0.;
        IR_GIVE_OPTIONAL_FIELD(ir, iso_offset, _IFT_NlIsoMoistureMaterial_iso_offset);
    } else if ( this->Isotherm == multilinear ) { // multilinear isotherm = type 1
        IR_GIVE_FIELD(ir, iso_h, _IFT_NlIsoMoistureMaterial_iso_h);
        IR_GIVE_FIELD(ir, iso_wh, _IFT_NlIsoMoistureMaterial_iso_wh);

        if ( iso_h.giveSize() != iso_wh.giveSize() ) {
            throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_iso_h, "the size of 'iso_h' and 'iso_w(h)' must be the same");
        }

        for ( int i = 1; i < iso_h.giveSize(); i++ ) {
            if ( ( iso_h.at(i) < 0. ) || ( iso_h.at(i) > 1. ) ) {
                throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_iso_h, "iso_h must be in the range <0; 1>");
            }
        }
    } else if ( this->Isotherm == Ricken ) { // reference mentioned in Kuenzel isotherm = type 2
        IR_GIVE_FIELD(ir, dd, _IFT_NlIsoMoistureMaterial_dd);
    } else if ( this->Isotherm == Kuenzel ) { // isotherm = type 3
        IR_GIVE_FIELD(ir, wf, _IFT_NlIsoMoistureMaterial_wf);
        IR_GIVE_FIELD(ir, b, _IFT_NlIsoMoistureMaterial_b);
    } else if ( this->Isotherm == Hansen ) { // isotherm = type 4
        IR_GIVE_FIELD(ir, uh, _IFT_NlIsoMoistureMaterial_uh);
        IR_GIVE_FIELD(ir, A, _IFT_NlIsoMoistureMaterial_a);
        IR_GIVE_FIELD(ir, nn, _IFT_NlIsoMoistureMaterial_nn);
        IR_GIVE_FIELD(ir, rhodry, _IFT_NlIsoMoistureMaterial_rhodry);
    } else if ( this->Isotherm == BSB ) { // isotherm = type 5
        IR_GIVE_FIELD(ir, c, _IFT_NlIsoMoistureMaterial_c);
        IR_GIVE_FIELD(ir, k, _IFT_NlIsoMoistureMaterial_k);
        IR_GIVE_FIELD(ir, Vm, _IFT_NlIsoMoistureMaterial_vm);
        IR_GIVE_FIELD(ir, rhodry, _IFT_NlIsoMoistureMaterial_rhodry);
    } else if ( this->Isotherm == bilinear ) { // bilinear isotherm = type 6
        IR_GIVE_FIELD(ir, moistureCapacity, _IFT_NlIsoMoistureMaterial_capa);
        IR_GIVE_FIELD(ir, wf, _IFT_NlIsoMoistureMaterial_wf);
        IR_GIVE_FIELD(ir, hx, _IFT_NlIsoMoistureMaterial_hx);
        IR_GIVE_FIELD(ir, dx, _IFT_NlIsoMoistureMaterial_dx);

        this->iso_offset = 0.;
        IR_GIVE_OPTIONAL_FIELD(ir, iso_offset, _IFT_NlIsoMoistureMaterial_iso_offset);

        this->capa2 = ( wf - hx * moistureCapacity - iso_offset ) / ( 1 - hx );
        double wa = ( hx - dx ) * moistureCapacity + iso_offset;
        double wb = hx * moistureCapacity + dx * capa2 + iso_offset;

        this->c1 = ( moistureCapacity * ( 2 * dx ) + capa2 * ( 2 * dx ) + 2 * wa - 2 * wb ) / ( 8 * dx * dx * dx );
        this->c2 = ( -3 * c1 * ( 2 * dx ) * ( 2 * dx ) - moistureCapacity + capa2 ) / ( 2 * ( 2 * dx ) );
    } else if ( this->Isotherm == vanGenuchten ) { // isotherm = type 7
        IR_GIVE_FIELD(ir, wf, _IFT_NlIsoMoistureMaterial_wf);
        IR_GIVE_FIELD(ir, vG_b, _IFT_NlIsoMoistureMaterial_vg_b);
        IR_GIVE_FIELD(ir, vG_m, _IFT_NlIsoMoistureMaterial_vg_m);
    } else {
        throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_isothermtype, "unknown isotherm type");
    }

    if ( this->Permeability == multilin ) {
        IR_GIVE_FIELD(ir, perm_h, _IFT_NlIsoMoistureMaterial_perm_h);
        IR_GIVE_FIELD(ir, perm_ch, _IFT_NlIsoMoistureMaterial_perm_ch);

        if ( perm_h.giveSize() != perm_ch.giveSize() ) {
            throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_perm_ch, "the size of 'perm_h' and 'perm_c(h)' must be the same");
        }

        for ( int i = 1; i < perm_h.giveSize(); i++ ) {
            if ( ( perm_h.at(i) < 0. ) || ( perm_h.at(i) > 1. ) ) {
                throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_perm_h, "must be in the range <0; 1>");
            }
        }
    } else if ( this->Permeability == Bazant ) {
        IR_GIVE_FIELD(ir, C1, _IFT_NlIsoMoistureMaterial_c1);
        IR_GIVE_FIELD(ir, n, _IFT_NlIsoMoistureMaterial_n);
        IR_GIVE_FIELD(ir, alpha0, _IFT_NlIsoMoistureMaterial_alpha0);
        IR_GIVE_FIELD(ir, hC, _IFT_NlIsoMoistureMaterial_hc);
    } else if ( this->Permeability == Xi ) {
        IR_GIVE_FIELD(ir, alphah, _IFT_NlIsoMoistureMaterial_alphah);
        IR_GIVE_FIELD(ir, betah, _IFT_NlIsoMoistureMaterial_betah);
        IR_GIVE_FIELD(ir, gammah, _IFT_NlIsoMoistureMaterial_gammah);
    } else if ( this->Permeability == KunzelPerm ) {
        IR_GIVE_FIELD(ir, mu, _IFT_NlIsoMoistureMaterial_mu);

        // read temperature - either constant or time-dependent
        if ( ir.hasField(_IFT_NlIsoMoistureMaterial_ttf) ) {
            IR_GIVE_FIELD(ir, T_TF, _IFT_NlIsoMoistureMaterial_ttf);
        } else {
            IR_GIVE_FIELD(ir, T, _IFT_NlIsoMoistureMaterial_t);
        }

        IR_GIVE_OPTIONAL_FIELD(ir, capillary_transport_coef, _IFT_NlIsoMoistureMaterial_capil_coef);

        IR_GIVE_OPTIONAL_FIELD(ir, PL, _IFT_NlIsoMoistureMaterial_pl);

        IR_GIVE_OPTIONAL_FIELD(ir, timeScale, _IFT_NlIsoMoistureMaterial_timescale);


        if ( CapillaryTransport == Multilin_h ) {
            IR_GIVE_FIELD(ir, capPerm_h, _IFT_NlIsoMoistureMaterial_capperm_h);
            IR_GIVE_FIELD(ir, capPerm_Dwh, _IFT_NlIsoMoistureMaterial_capperm_dwh);

            if ( !( capPerm_h.giveSize() == capPerm_Dwh.giveSize() ) ) {
                throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_capperm_dwh, "size of 'capPerm_h' and 'capPerm_dw(h)' must be the same");
            }

            for ( int i = 1; i < capPerm_h.giveSize(); i++ ) {
                if ( ( capPerm_h.at(i) < 0. ) || ( capPerm_h.at(i) > 1. ) ) {
                    throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_capperm_h, "must be in the range <0; 1>");
                }
            }
        } else if ( this->CapillaryTransport == Multilin_wV ) {
            rhoH2O = 1000.;
            IR_GIVE_OPTIONAL_FIELD(ir, rhoH2O, _IFT_NlIsoMoistureMaterial_rhoh2o);

            IR_GIVE_FIELD(ir, capPerm_wV, _IFT_NlIsoMoistureMaterial_capperm_wv);
            IR_GIVE_FIELD(ir, capPerm_DwwV, _IFT_NlIsoMoistureMaterial_capperm_dwwv);

            if ( !( capPerm_wV.giveSize() == capPerm_DwwV.giveSize() ) ) {
                throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_capperm_dwwv, "size of 'capPerm_wV' and 'capPerm_Dw(wV)' must be the same");
            }

            for ( int i = 1; i < capPerm_wV.giveSize(); i++ ) {
                if ( ( capPerm_wV.at(i) < 0. ) || ( capPerm_wV.at(i) > 1. ) ) {
                    throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_capperm_wv, "must be in the range <0; 1>");
                }
            }
        } else { // according to Kunzel
            if ( this->Isotherm == Hansen  ) {
                this->wf = rhodry * uh;
            } else {
                IR_GIVE_FIELD(ir, wf, _IFT_NlIsoMoistureMaterial_wf);
            }
            IR_GIVE_FIELD(ir, Abs, _IFT_NlIsoMoistureMaterial_abs);
        }
    } else {
        throw ValueInputException(ir, _IFT_NlIsoMoistureMaterial_permeabilitytype, "unknown permeability type");
    }

    wn = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, wn, _IFT_NlIsoMoistureMaterial_wn);
    IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_NlIsoMoistureMaterial_alpha);
}

double
NlIsoMoistureMaterial::giveMoistureCapacity(GaussPoint *gp, TimeStep *tStep) const
{
    double humidity = this->giveHumidity(gp, VM_Total);

    if ( this->Isotherm == linear ) {
        return moistureCapacity;
    } else if ( this->Isotherm == multilinear ) {
        double tol = 1.e-10;
        for ( int i = 2; i <= iso_h.giveSize(); i++ ) {
            if ( ( humidity - iso_h.at(i) ) < tol ) {
                return ( iso_wh.at(i) - iso_wh.at(i - 1) ) / ( iso_h.at(i) - iso_h.at(i - 1) );
            }
        }
    } else if ( this->Isotherm == Ricken ) {
        return 1. / ( dd * ( 1. - humidity ) );
    } else if ( this->Isotherm == Kuenzel ) {
        return wf * ( b - 1. ) * b / ( ( b - humidity ) * ( b - humidity ) );
    } else if ( this->Isotherm == Hansen ) {
        return rhodry * uh / ( A * nn * humidity * pow( ( 1 - log(humidity) / A ), ( 1 / nn + 1 ) ) );
    } else if ( this->Isotherm == BSB ) {
        double nominator, denominator;
        nominator = c * k * Vm * rhodry * ( 1. + k * k * humidity * humidity * c - k * k * humidity * humidity );
        denominator = ( 1. - k * humidity ) * ( 1. - k * humidity ) * ( 1. + ( c - 1. ) * k * humidity ) * ( 1. + ( c - 1. ) * k * humidity );
        return nominator / denominator;
    } else if ( this->Isotherm == bilinear ) {
        if ( humidity <= ( hx - dx ) ) {
            return moistureCapacity;
        } else if ( humidity >= ( hx + dx ) ) {
            return capa2;
        } else {
            return moistureCapacity + 2 * c2 * ( humidity - hx + dx ) + 3 * c1 * pow(humidity - hx + dx, 2);
        }
    } else if ( this->Isotherm == vanGenuchten ) {
        return wf * vG_m * vG_b / humidity * pow( ( pow( ( -vG_b * log(humidity) ), ( 1. / ( 1. - vG_m ) ) ) + 1. ), ( -vG_m - 1 ) ) * ( 1. / ( 1. - vG_m ) * pow( ( -vG_b * log(humidity) ), ( vG_m / ( 1. - vG_m ) ) ) );
    } else {
        OOFEM_ERROR("unknown isotherm type");
    }

    return 0.;
}

double
NlIsoMoistureMaterial::giveMoistureContent(double humidity) const
{
    if ( this->Isotherm == linear ) {
        return max(moistureCapacity * humidity + iso_offset, 0.);
    } else if ( this->Isotherm == multilinear ) {
        double tol = 1.e-10;
        for ( int i = 1; i <= iso_h.giveSize(); i++ ) {
            if ( ( humidity - iso_h.at(i) ) < tol ) {
                return iso_wh.at(i - 1) +  ( iso_wh.at(i) - iso_wh.at(i - 1) ) / ( iso_h.at(i) - iso_h.at(i - 1) ) * ( humidity - iso_h.at(i - 1) );
            }
        }
    } else if ( this->Isotherm == Ricken ) {
        return wf - log(1. - humidity) / dd;
    } else if ( this->Isotherm == Kuenzel ) {
        return wf * ( b - 1. ) * humidity / ( b - humidity );
    } else if ( this->Isotherm == Hansen ) {
        return rhodry * uh * pow( ( 1. - log(humidity) / A ), ( -1. / nn ) );
    } else if ( this->Isotherm == BSB ) {
        return rhodry * c * k * Vm * humidity / ( ( 1. - k * humidity ) * ( c - 1. ) * k * humidity );
    } else if ( this->Isotherm == bilinear ) {
        if ( humidity <= ( hx - dx ) ) {
            return max(moistureCapacity * humidity + iso_offset, 0.);
        } else if ( humidity >= ( hx + dx ) ) {
            return hx * moistureCapacity + ( humidity - hx ) * capa2  + iso_offset;
        } else {
            return ( hx - dx ) * moistureCapacity + moistureCapacity * ( humidity - hx + dx ) + c2 * pow(humidity - hx + dx, 2) + c1 * pow(humidity - hx + dx, 3)  + iso_offset;
        }
    } else if ( this->Isotherm == vanGenuchten ) {
        return wf * pow( ( pow( ( -vG_b * log(humidity) ), ( 1. / ( 1. - vG_m ) ) ) + 1. ), -vG_m);
    } else {
        OOFEM_ERROR("unknown isotherm type");
    }

    return 0.;
}

double
NlIsoMoistureMaterial::givePermeability(GaussPoint *gp, TimeStep *tStep) const
{
    double humidity = this->giveHumidity(gp, VM_Total);

    if ( this->Permeability == multilin ) {
        double tol = 1.e-10;
        for ( int i = 1; i <= perm_h.giveSize(); i++ ) {
            if ( ( humidity - perm_h.at(i) ) < tol ) {
                return perm_ch.at(i - 1) + ( perm_ch.at(i) - perm_ch.at(i - 1) ) * ( humidity - perm_h.at(i - 1) ) / ( perm_h.at(i) - perm_h.at(i - 1) );
            }
        }
    } else if ( this->Permeability == Bazant ) {
        return C1 * ( alpha0 + ( 1. - alpha0 ) / ( 1. + pow( ( 1. - humidity ) / ( 1. - hC ), n ) ) );
    } else if ( this->Permeability == Xi ) {
        double power = pow( 10., gammah * ( humidity - 1. ) );
        return alphah + betah * ( 1. - pow(2., -power) );
    } else if ( this->Permeability == KunzelPerm ) {
        // [kg/m^3]
        double dw_dh = this->giveMoistureCapacity(gp, tStep);
        // capillary transport coefficient [m^2/s]
        double Dw = this->computeCapTranspCoeff(humidity);

        double temper_factor = this->computeTemperatureEffectOnViscosity(gp, tStep);
        double p_sat = this->computeSaturationWaterVaporPressure(gp, tStep);
        double deltap = this->computeVaporDiffusionCoeff(gp, tStep);

        return temper_factor * Dw * dw_dh + deltap * p_sat;
    } else {
        OOFEM_ERROR("unknown permeability type");
    }

    return 0.;
}


double
NlIsoMoistureMaterial::computeCapTranspCoeff(double humidity) const
{
    double Dw = 0.;

    if ( this->CapillaryTransport == Multilin_h ) {
        double tol = 1.e-10;
        for ( int i = 1; i <= capPerm_h.giveSize(); i++ ) {
            if ( ( humidity - capPerm_h.at(i) ) < tol ) {
                Dw = capPerm_Dwh.at(i - 1) + ( capPerm_Dwh.at(i) - capPerm_Dwh.at(i - 1) ) * ( humidity - capPerm_h.at(i - 1) ) / ( capPerm_h.at(i) - capPerm_h.at(i - 1) );
                break;
            }
        }
    } else if ( this->CapillaryTransport == Multilin_wV ) {
        double wV = this->giveMoistureContent(humidity) / rhoH2O;
        double tol = 1.e-10;
        for ( int i = 1; i <= capPerm_wV.giveSize(); i++ ) {
            if ( ( wV - capPerm_wV.at(i) ) < tol ) {
                Dw = capPerm_DwwV.at(i - 1) + ( capPerm_DwwV.at(i) - capPerm_DwwV.at(i - 1) ) * ( wV - capPerm_wV.at(i - 1) ) / ( capPerm_wV.at(i) - capPerm_wV.at(i - 1) );
                break;
            }
        }
    } else if ( this->CapillaryTransport == KunzelCT ) {
        double w = this->giveMoistureContent(humidity);
        Dw = 3.8 * ( Abs / wf ) * ( Abs / wf ) * pow(this->capillary_transport_coef, w / wf - 1.);
    } else {
        OOFEM_ERROR("unknown capillary transport type");
    }




    return Dw;
}


double
NlIsoMoistureMaterial::computeVaporDiffusionCoeff(GaussPoint *gp, TimeStep *tStep) const
{
    double temperature = this->giveTemperature(gp, tStep);

    double delta = this->timeScale * 2.0 * 1.e-7 * pow(temperature, 0.81) / this->PL;
    return delta / this->mu;
}


double
NlIsoMoistureMaterial::computeSaturationWaterVaporPressure(GaussPoint *gp, TimeStep *tStep) const
{
    double temperature  = this->giveTemperature(gp, tStep);

    temperature -= 273.15;
    double T0, a;

    if ( temperature < 0. ) {
        a = 22.44;
        T0 = 272.44;
    } else {
        a = 17.08;
        T0 = 234.18;
    }

    return 611. * exp( a * temperature / ( T0 + temperature ) );
}


double
NlIsoMoistureMaterial::computeTemperatureEffectOnViscosity(GaussPoint *gp, TimeStep *tStep) const
{
    double temperature = this->giveTemperature(gp, tStep);

    // ratio of water viscosity at given temperature and reference temperature, taken as 20 C
    // viscosity is comptued according to Gawin
    // eta(T) = 0.6612 * (T-229)**-1.532
    return pow(64. / ( temperature - 229. ), -1.532);
}


double
NlIsoMoistureMaterial::giveTemperature(GaussPoint *gp, TimeStep *tStep) const
{
    double temperature;

    if ( this->T_TF !=  0 ) {
        temperature = domain->giveFunction(this->T_TF)->evaluateAtTime(tStep->giveTargetTime() );
    } else {
        temperature = this->T;
    }

    return temperature;
}





double
NlIsoMoistureMaterial::giveHumidity(GaussPoint *gp, ValueModeType mode) const
{
    auto tempState = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) )->giveTempField();
    if ( tempState > 1.0 || tempState < 0.0 ) {
        OOFEM_ERROR("Relative humidity %.3f is out of range", tempState);
    } else {
        return tempState;
    }
}

bool
NlIsoMoistureMaterial::hasInternalSource() const
{
    return this->wn != 0.;
}

void
NlIsoMoistureMaterial::computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
{
    val.resize(1);
    if ( mode == VM_Total || mode == VM_TotalIntrinsic ) {
        val.at(1) = -wn * ( this->alpha.eval( { { "t", tStep->giveTargetTime() } }, this->giveDomain() ) - this->alpha.eval( { { "t", tStep->giveTargetTime() - tStep->giveTimeIncrement() } }, this->giveDomain() ) ) / tStep->giveTimeIncrement();
    } else {
        OOFEM_ERROR("Undefined mode %s\n", __ValueModeTypeToString(mode) );
    }
}
} // end namespace oofem
