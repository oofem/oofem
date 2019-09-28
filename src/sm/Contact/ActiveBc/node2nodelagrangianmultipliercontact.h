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

#ifndef node2nodelagrangianmultipliercontact_h
#define node2nodelagrangianmultipliercontact_h


#include "activebc.h"


 ///@name Input fields for _IFT_ContactElement
 //@{
#define _IFT_Node2NodeLagrangianMultiplierContact_Name "n2nlagrangianmultipliercontact"
#define _IFT_Node2NodeLagrangianMultiplierContact_useTangent "usetangent"

#define _IFT_Node2NodeLagrangianMultiplierContact_masterSet "masterset"
#define _IFT_Node2NodeLagrangianMultiplierContact_slaveSet "slaveset"




//@}

namespace oofem {
    class Domain;
    class SparseMtrx;
    class TimeStep;
    class DofManager;
    class GaussPoint;
    class UnknownNumberingScheme;
    class FloatMatrix;
    class IntegrationRule;
    class ContactElement;
    class Node;

    class OOFEM_EXPORT Node2NodeLagrangianMultiplierContact : public ActiveBoundaryCondition
    {
    private:
        bool useTangent; ///< Determines if tangent should be used.
        IntArray slaveSet;
        IntArray masterSet;
        std::vector<  DofManager * >lmdm;
    public:

		/**
		* Constructor. Creates Node2NodeLagrangianMultiplierContact with given number belonging to domain aDomain.
		* @param n Node2NodeLagrangianMultiplierContact's number in domain
		* @param d reference to Node2NodeLagrangianMultiplierContact's domain
		*/
        Node2NodeLagrangianMultiplierContact(int n, Domain *d);
        /// Destructor.
        virtual ~Node2NodeLagrangianMultiplierContact() {};

		/**
		 * Initializes boundary condition from a given InputRecord
		 * @param ir Pointer to the InputRecord
		 * @return An IRResultType object.
		 */
        virtual IRResultType initializeFrom(InputRecord *ir);

		/**
		 * Assembles contact contribution to the stiffness matrix. All contact pairs are
		 * processed iteratively.
		 *
		 * @param answer [out] A SparseMtrx object for assembly of the answer
		 * @param tStep The current time step
		 * @param type Requested type of matrix. Function is only implemented for tangential stiffness
		 * @param r_s
		 * @param c_s
		 * @param scale
		 */
        virtual void assemble(SparseMtrx &answer, TimeStep *tStep,
            CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0) override;

		/**
		 * Assembles contact contribution to force vector. All contact pairs are
		 * processed iteratively.
		 *
		 * @param answer [out] A SparseMtrx object for assembly of the answer
		 * @param tStep The current time step
		 * @param type Requested type of matrix. Function is only implemented for external forces vector
		 * @param mode
		 * @param s
		 * @param eNorms
		 */
        virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
            CharType type, ValueModeType mode,
            const UnknownNumberingScheme &s, FloatArray *eNorms = NULL) override;


        virtual const char *giveClassName() const { return "Node2NodeLagrangianMultiplierContact"; }
        virtual const char *giveInputRecordName() const { return _IFT_Node2NodeLagrangianMultiplierContact_Name; }

        int giveNumberOfInternalDofManagers() override { return masterSet.giveSize(); }
        DofManager *giveInternalDofManager(int i) override { return this->lmdm.at(i - 1); }


    private:

		/**
		 * Computes the tangential stiffness matrix from contact of two nodes, given as:
		 * K_c = [ 0  Nv 
				  Nv  0 ]
		 * where Nv is the normal matrix of contact
		 *
		 * @param answer [out] a column matrix
		 * @param masterNode the master node concerned
		 * @param slaveNode the slave node concerned
		 * @param tStep the current time step
		 */
		double computeTangentFromContact(FloatMatrix &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep);
        
		/**
		 * Computes the size of gap between two nodes. Gap is negative if the nodes are on the opposite
		 * side of each other than in the undeformed configuration. Precise computation involves
		 * computing the dot product of undeformed and deformed distance vector between the nodes
		 *
		 * @param answer [out] the size of the gap
		 * @param masterNode the master node concerned
		 * @param slaveNode the slave node concerned
		 * @param tStep the current time step
		 */
		void computeGap(double &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep);

		/**
		 * Computes the normal matrix of contact. Computes the projection (distance vector), normalizes it
		 * and assembles into a column matrix in the form N = {n -n}^T
		 *
		 * @param answer [out] a column matrix
		 * @param masterNode the master node concerned
		 * @param slaveNode the slave node concerned
		 * @param tStep the current time step
		 */
        void computeNvMatrixAt(FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *TimeStep);

		/**
		 * Computes the vector of external forces from contact of two nodes, given as:
		 * f_ext = g_c
		 * where g_c is the gap
		 *
		 * @param answer [out] a column matrix (of only one member)
		 * @param masterNode the master node concerned
		 * @param slaveNode the slave node concerned
		 * @param tStep the current time step
		 */
        void computeExternalForcesFromContact(FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep);

		/**
		 * Gives the location array of the created LMs
		 * 
		 * @param answer [out] a column matrix (of only one member)
		 * @param rs a numbering scheme to be used
		 */
        void giveLagrangianMultiplierLocationArray(const UnknownNumberingScheme &r_s, std::vector< IntArray > &answer);

    public:

		/**
		 * Computes which DOFs interact with which through the contact condition. Useful for
		 * creation of global sparse matrices.
		 *
		 * @param rows [out] - row indices
		 * @param cols [out] - column indices
		 * @param type
		 * @param r_s unknown numbering scheme for rows
		 * @param c_s unknown numbering scheme for columns
		 */
        void giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;
    };
} // end namespace oofem
#endif // node2nodelagrangianmultipliercontact_h
