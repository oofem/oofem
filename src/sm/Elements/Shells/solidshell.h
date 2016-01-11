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

#ifndef solidshell_h
#define solidshell_h

#include "Elements/3D/lspace.h"


#define _IFT_SolidShell_Name "solidshell"

// EAS parameters
#define _IFT_SolidShell_EAS_type "eas_type"


namespace oofem {
class FEI3dHexaLin;

/**
 * This class implements a Linear 8-noded shell like solid with ANS and EAS to remove and reduce certain locking aspects
 * Each node has 3 degrees of freedom.
 *
 *@author Jim Brouzoulis 
 */
class SolidShell  : public LSpace
{
protected:
    static FEI3dHexaLin interpolation;

public:
    SolidShell(int n, Domain * d);
    virtual ~SolidShell() { }
    virtual FEInterpolation *giveInterpolation() const;
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeBHmatrixAt(FloatArray &lCoords, FloatMatrix &answer);
    
    virtual void computeBEmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep);
    virtual void computeGaussPoints();
    virtual void computeEASBmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    void computeGeometricStiffness(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep);

    void computeFVector(FloatArray &answer, FloatArray &lCoords, FloatArray &ae);
    void computeEVector(FloatArray &answer, FloatArray &lCoords, FloatArray &ae);
    void computeEASfield(FloatArray &answer, TimeStep *tStep);
    
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_LSpace_Name; }
    virtual const char *giveClassName() const { return "SolidShell"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    
    // variables associated with EAS
    int EAS_type; 
    FloatArray tempAlpha;
    FloatArray alpha;
    FloatArray u_k;
    FloatMatrix invKEE;
    FloatMatrix KEC;
    FloatArray fE;
    void computeAlpha(FloatArray &answer, FloatArray &u);
    
    void computeBondTransformationMatrix(FloatMatrix &answer, FloatMatrix &base);
    virtual void postInitialize();
    virtual int checkConsistency(){ return 1; };
    
    
    private:
    void x(int arg1);
};
} // end namespace oofem
#endif
