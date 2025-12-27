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

#ifndef linkslip_h
#define linkslip_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for LatticeSlip
//@{
#define _IFT_LinkSlip_Name "linkslip"
#define _IFT_LinkSlip_talpha "talpha"
#define _IFT_LinkSlip_kn "kn"
#define _IFT_LinkSlip_kl "kl"
#define _IFT_LinkSlip_type "type"
#define _IFT_LinkSlip_t0 "t0"
#define _IFT_LinkSlip_s1 "s1"
#define _IFT_LinkSlip_s2 "s2"
#define _IFT_LinkSlip_s3 "s3"
#define _IFT_LinkSlip_tf "tf"
#define _IFT_LinkSlip_alpha "alpha"
//@}

namespace oofem {

 /**
 * This class implements a constitutive model for a bond link for connecting beam (frame) and continuum elements in unstructured meshes.
 * The main idea is to use the rotation of the beam element and the rigid arm from the beam node to the continuum element node
 * to compute the displacement jump along the rebar element (and two components, which are perpendicular to each other and lie 
 * in a plane for which the direction along the rebar is normal to.
 * This constitutive model differs from the standard bond model, because only the first component of the jump is used to determine the bond stress.  
 *
 * @author: Peter Grassl
*/


class LinkSlipStatus : public StructuralInterfaceMaterialStatus
{
protected:

  double kappa = 0;

  double tempKappa = 0;

public:

    /// Constructor
    LinkSlipStatus(GaussPoint *g);
    
    double giveKappa(){ return this->kappa; }

    void letTempKappaBe(const double &v)
    { this->tempKappa = v; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "LinkSlipStatus"; }
    
    void initTempStatus() override;

    void updateYourself(TimeStep *) override;

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;
};

/**
 * This class implements a slip model for concrete for lattice elements.
 */
class LinkSlip : public StructuralInterfaceMaterial
    //
{
protected:

    ///Normal modulus
    double kNormal = 0.;

    ///Lateral modulus
    double kLateral = 0.;

    int type = 0;
    
    ///Strength for slip component
    double tauMax = 0., tauFinal = 0.;
    
    double s1 = 0., s2 = 0., s3 = 0.;

    double alpha = 0.;

public:

    LinkSlip(int n, Domain *d);

    bool hasAnalyticalTangentStiffness() const override { return true; }
    
    const char *giveInputRecordName() const override { return _IFT_LinkSlip_Name; }
    const char *giveClassName() const override { return "LinkSlip"; }

    void initializeFrom(InputRecord &ir) override;

    //  virtual void computeStressIndependentStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, ValueModeType mode);

    // virtual void  giveThermalDilatationVector(FloatArray &answer,  GaussPoint *gp,  TimeStep *tStep); // unused for interface materials


    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }
    
    double evaluateBondStress(const double kappa) const;

    FloatArrayF<3> giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<3,3> give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;
    
    Interface *giveInterface(InterfaceType) override;
                                  

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime) override;
};
} // end namespace oofem

#endif
