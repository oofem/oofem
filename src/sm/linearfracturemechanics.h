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

#ifndef linearfracturemechanics_h
#define linearfracturemechanics_h

#include "linearstatic.h"
#include "meshpackagetype.h"

///@name Input fields for LinearFractureMechanics
//@{
#define _IFT_LinearFractureMechanics_Name "linearfracturemechanics"
#define _IFT_LinearFractureMechanics_Meshpackage "meshpackage"
#define _IFT_LinearFractureMechanics_FractureCriterion "fracturecriterion"
#define _IFT_LinearFractureMechanics_CrackTipCoordinates "cracktipcoordinates"
#define _IFT_LinearFractureMechanics_FractureThougness "fracturethougness"
#define _IFT_LinearFractureMechanics_FractureEnergy "fractureenergy"
#define _IFT_LinearFractureMechanics_rmax "rmax"
#define PI 3.14159265358979
//@}


/** Type determining fracture criterion */
enum FractureCriterion { FC_SIF, FC_JIntegral };


namespace oofem {
/**
 * This class implements a linear fracture mechanics.
 */
class LinearFractureMechanics : public LinearStatic
{
protected:
    /// Meshing package used for refinements.
    MeshPackageType meshPackage;
	FractureCriterion fractureCriterionType;
	double rmax;
	double fractureTougness;
	double fractureEnergy;
	FloatArray crackTipCoordinates;

public:
    LinearFractureMechanics(int i, EngngModel *_master = NULL) : LinearStatic(i, _master) { }
    virtual ~LinearFractureMechanics() { }

    virtual void updateYourself(TimeStep *tStep);

    
    virtual void terminate(TimeStep *tStep);

    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void updateDomainLinks();

    virtual IRResultType initializeFrom(InputRecord *ir);

    // identification
    virtual const char *giveClassName() const { return "LinearFractureMechanics"; }
    virtual const char *giveInputRecordName() const { return _IFT_LinearFractureMechanics_Name; }

	// Computes fracture SIF or other fracture criteria depending on fracturecriteria parameter
	virtual bool evaluateFractureCriterion(FractureCriterion fractureCriterionType, TimeStep *tStep);

	void giveElementStressIntensityFactor(FloatArray &answer, Element *elem, TimeStep *tStep);
	void giveStressIntensityFactor(double &answer, TimeStep *tStep);



};
} // end namespace oofem
#endif // linearfracturemechanics_h
