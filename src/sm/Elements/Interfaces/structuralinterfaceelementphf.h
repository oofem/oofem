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

#ifndef structuralinterfaceelementphf_h
#define structuralinterfaceelementphf_h

#include "Elements/Interfaces/structuralinterfaceelement.h"
#include "../sm/CrossSections/structuralinterfacecrosssection.h"
#include "element.h"
#include "floatmatrix.h"
#include "function.h"
#include "matresponsemode.h"
#include "valuemodetype.h"
#include "integrationdomain.h"
#include "dofmantransftype.h"

namespace oofem {
class TimeStep;
class Node;
class StructuralInterfaceMaterial;
class GaussPoint;
class FloatArray;
class IntArray;
class FEInterpolation;


/**
 * Interface element class with phase field (PhF) modeling of damage
 */
class StructuralInterfaceElementPhF : public StructuralInterfaceElement
{
protected:
    /// Initial displacement vector, describes the initial nodal displacements when element has been casted.
    FloatArray initialDisplacements;
    FEInterpolation *interpolation;
    /// Flag indicating if geometrical nonlinearities apply.
    int nlGeometry;

    
    IntArray loc_u, loc_d;
    
public:
    /**
     * Constructor. Creates structural element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    StructuralInterfaceElementPhF(int n, Domain * d);
    /// Destructor.
    virtual ~StructuralInterfaceElementPhF();

    //virtual FEInterpolation *giveInterpolation() const { return interpolation; };

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    
    //Phf
    void computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    void computeStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    void computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    void computeStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    void giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    void giveInternalForcesVector_d(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);    

    virtual void giveDofManDofIDMask_u(IntArray &answer) = 0;
    virtual void giveDofManDofIDMask_d(IntArray &answer) = 0;
    virtual void computeLocationArrayOfDofIDs( const IntArray &dofIdArray, IntArray &answer );
    
    virtual void getLocationArray_u( IntArray &answer )=0;
    virtual void getLocationArray_d( IntArray &answer )=0;
    
    //remove
    int computeNumberOfDofs();
    
    ///@todo move these to a cross section model later
    double internalLength;
    double giveInternalLength( ) { return internalLength; };
    //double criticalEnergy;
    //double giveCriticalEnergy() { return criticalEnergy; };
    //double relaxationTime;
    //double giveRelaxationTime( ) { return relaxationTime; };
    //double penaltyParameter;
    //double givePenaltyParameter() { return penaltyParameter; };
    //double psiBar0;
    //double givePsiBar0() { return psiBar0; };
    
    
    ///@todo put in material model
    double computeFreeEnergy( GaussPoint *gp, TimeStep *tStep );
    
    virtual double computeDamageAt(GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    

    void computeDisplacementUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *stepN);
    void computeDamageUnknowns(FloatArray &answer, ValueModeType valueMode, TimeStep *stepN);
    
    
    
    void computeBStress_u(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord);
    void computeNStress_d( FloatArray &answer, GaussPoint *gp, TimeStep *tStep, int useUpdatedGpRecord );
    

    //Interpolation matrices for phase field
    virtual void computeBd_matrixAt(GaussPoint *, FloatMatrix &);
    virtual void computeNd_matrixAt(const FloatArray &lCoords, FloatMatrix &N);
    
    
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual void computeTraction(FloatArray &traction, IntegrationPoint *ip, FloatArray &jump, TimeStep *tStep);
    //virtual void computeSpatialJump(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    //virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual void computeCovarBaseVectorsAt(GaussPoint *gp, FloatMatrix &G) = 0;
    

    //@}

    // Overloaded methods.
    virtual void updateInternalState(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);
    //virtual int checkConsistency();
    //virtual IRResultType initializeFrom(InputRecord *ir);
    //virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveClassName() const { return "StructuralInterfaceElementPhF"; };

    //StructuralInterfaceCrossSection *giveInterfaceCrossSection();

    //virtual methods that should be overloaded by the elements
    virtual void giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const double damage, TimeStep *tStep)
    {
        OOFEM_ERROR("not implemented for the current element");
    }
   


protected:


};
} // end namespace oofem
#endif // structuralinterfaceelement_h
