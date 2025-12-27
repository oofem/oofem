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

#ifndef concretefcmviscoelastic_h
#define concretefcmviscoelastic_h

#include "concretefcm.h"
// for the unique_ptr
#include <memory>

///@name Input fields for ConcreteFCMViscoElastic
//@{
#define _IFT_ConcreteFCMViscoElastic_Name "concretefcmviscoelastic"
#define _IFT_ConcreteFCMViscoElastic_viscoMat "viscomat"
#define _IFT_ConcreteFCMViscoElastic_timedepfracturing "timedepfracturing"
#define _IFT_ConcreteFCMViscoElastic_fib_s "fib_s"
#define _IFT_ConcreteFCMViscoElastic_fcm28 "fcm28"
#define _IFT_ConcreteFCMViscoElastic_timeFactor "timefactor"
#define _IFT_ConcreteFCMViscoElastic_stiffnessFactor "stiffnessfactor"
#define _IFT_ConcreteFCMViscoElastic_gf28 "gf28"
#define _IFT_ConcreteFCMViscoElastic_ft28 "ft28"
//@}

namespace oofem {
/**
 * This class manages the status of ConcreteFCMViscoElastic
 */
class ConcreteFCMViscoElasticStatus : public ConcreteFCMStatus
{

protected:
    std :: unique_ptr< GaussPoint >slaveGpVisco;
  
    /// hydration-degree dependent tensile strength
    double var_ft = 0.;
    /// hydration-degree dependent fracture energy
    double var_gf = 0.;
  
public:
    ConcreteFCMViscoElasticStatus(GaussPoint *g);

    double giveFractureEnergy() const { return var_gf; }
    void setFractureEnergy(double new_Gf) { var_gf = new_Gf; }

    double giveTensileStrength() const { return var_ft; }
    void setTensileStrength(double new_ft) { var_ft = new_ft; }

    GaussPoint *giveSlaveGaussPointVisco() { return this->slaveGpVisco.get(); }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "ConcreteFCMViscoElasticStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;  
};


/**
 * add description
 */
class ConcreteFCMViscoElastic : public ConcreteFCM
{
public:
    ConcreteFCMViscoElastic(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;
    const char *giveClassName() const override { return "ConcreteFCMViscoElastic"; }
    const char *giveInputRecordName() const override { return _IFT_ConcreteFCMViscoElastic_Name; }

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<ConcreteFCMViscoElasticStatus>(gp); }

    double give(int aProperty, GaussPoint *gp) const override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep) const override;

    FloatArray computeStressIndependentStrainVector(GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

    
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    
    MaterialStatus* giveStatus(GaussPoint *gp) const override;

protected:
    /// number of the viscoelastic material
    int viscoMat = 0;

    bool fib = false;
    double fib_s = 0.;
    double fib_fcm28 = 0.;
    /* scaling time factor, 1 day expressed in the time units of the analysis
     *  e.g. if the time runs in days it is 1, if in seconds it is 86400 */
    double timeFactor = 0.;
    /* scaling factor transforming PREDICTED strength and fracture energy
     *  e.g. if the stiffness should be in MPa, then stiffnessFactor = 1.e6 */
    double stiffnessFactor = 0.;    


    double giveTensileStrength(GaussPoint *gp, TimeStep *tStep) const override;
    double giveFractureEnergy(GaussPoint *gp, TimeStep *tStep) const override;    

   
    double computeOverallElasticStiffness(GaussPoint *gp, TimeStep *tStep) const override;
    double computeOverallElasticShearModulus(GaussPoint *gp, TimeStep *tStep) const override;
    

    int checkConsistency(void) override; 
   
    /// returns equivalent time (used to compute time-dependent ft and gf)
    virtual double giveEquivalentTime(GaussPoint *gp, TimeStep *tStep) const;

};
} // end namespace oofem
#endif // concretefcmviscoelastic_h
