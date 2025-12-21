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
#ifndef structuralpenaltycontactbc_h
#define structuralpenaltycontactbc_h


#include "Contact/contactbc.h"
#include "Contact/contactpoint.h"
#include "ContactSurface/structuralfecontactsurface.h"

///@name Input fields for _IFT_ContactElement
//@{
#define _IFT_StructuralPenaltyContactBoundaryCondition_Name "structuralpenaltycontactbc"
#define _IFT_StructuralPenaltyContactBoundaryCondition_penaltyNormal "pn"
#define _IFT_StructuralPenaltyContactBoundaryCondition_penaltyTangential "pt"
#define _IFT_StructuralPenaltyContactBoundaryCondition_friction "friction"
////

#define _IFT_StructuralPenaltyContactBoundaryCondition_masterSurfaceNum "mastersurface"
#define _IFT_StructuralPenaltyContactBoundaryCondition_slaveSurfaceNum "slavesurface"
///
#define _IFT_StructuralPenaltyContactBoundaryCondition_nsd "nsd"
///
#define _IFT_StructuralPenaltyContactBoundaryCondition_algo "algo"

//@}

namespace oofem {
class Domain;
class SparseMtrx;
class TimeStep;
class DofManager;
class UnknownNumberingScheme;
class FloatMatrix;

/**
 * Boundary condition class for structural penalty contact. Maintains set of corresponding pairs of nodes and segments that are checked for contact. However, the sets should be move to contact search algorithm in the near future. 
 * @author Martin Horák, nitramkaroh@seznam.cz
 * @Note References: Tod A. Laursen: Computational Contact and Impact Mechanics: Fundamentals of Modeling Interfacial Phenomena in Nonlinear Finite Element Analysis
 *
 * Tasks:
 * calculating contact contributions to the stiffness matrix and to the force vector
 * Currently, only normal contact is working
 * Frictional contact is under development
 */


enum class ContactProcess{ None, Sticking, Sliding};
  

class OOFEM_EXPORT StructuralPenaltyContactBoundaryCondition : public ContactBoundaryCondition
{
private:
    double penalty_normal;  ///< penalty in the normal direction.
    double penalty_tangential; ///< penalty in the tangent directions.
    double friction; ///< coefficient of friction
    int surface_dimension; ///< dimension of the surface, i.e., nsd - 1
    int algo;  ///< contact search algorithm
    
    //@todo: move to contact search algorithm
    int masterSurfaceNumber;
    int slaveSurfaceNumber;
    IntArray masterSurfaceElements;
    IntArray slaveSurfaceElements;
    // should be generalized to structural contact surface, i.e., a surface that is not necessarly discretized into finite elements
    StructuralFEContactSurface *masterContactSurface;
    StructuralFEContactSurface *slaveContactSurface;

public:
    /// Constructor.
    StructuralPenaltyContactBoundaryCondition(int n, Domain *d) : ContactBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~StructuralPenaltyContactBoundaryCondition() {};
    // initialization, i.e., reading input filex
    void initializeFrom(InputRecord &ir) override;
    void postInitialize() override;
    virtual const char *giveClassName() const override { return "StructuralPenaltyContactBoundaryCondition"; }
    virtual const char *giveInputRecordName() const override { return _IFT_StructuralPenaltyContactBoundaryCondition_Name; }
    //
    void setupContactSearchAlgorithm() override;
private:
    // compute tangent stiffness matrix
    void  computeTangentFromContact(FloatMatrix &answer, ContactPair *cp, TimeStep *tStep) override;
    // compute internal forces
    void computeInternalForcesFromContact(FloatArray &answer, ContactPair *cp, TimeStep *tStep) override;
    void computeTractions(double& normalTraction, FloatArray &tangentialTraction, FloatArray &tangentialTractionTrial, ContactProcess& mode, ContactPair* contactPair, TimeStep *tStep);
    // location arrays
    void giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;
 private:
    FloatMatrix computeContravariantMetric(const std::vector<FloatArray> &tangent_vectors);
    FloatMatrix computeCovariantMetric(const std::vector<FloatArray> &tangent_vectors);


    
};
} // end namespace oofem
#endif // node2segmentpenaltycontact_h
