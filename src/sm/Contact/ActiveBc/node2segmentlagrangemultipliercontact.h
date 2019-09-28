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



#pragma once

#include "activebc.h"
#include "Contact/ContactSegment/contactsegment.h"
#include "classfactory.h"
#include "masterdof.h"
#include "domain.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "unknownnumberingscheme.h"


 ///@name Input fields for _IFT_ContactElement
 //@{
#define _IFT_Node2SegmentLagrangianMultiplierContact_Name "n2slagrangianmultipliercontact"
#define _IFT_Node2SegmentLagrangianMultiplierContact_useTangent "usetangent"

#define _IFT_Node2SegmentLagrangianMultiplierContact_nodeSet "nodeset"
#define _IFT_Node2SegmentLagrangianMultiplierContact_segmentSet "segmentset"

#define _IFT_Node2SegmentLagrangianMultiplierContact_prescribedNormal "normal"




//@}

namespace oofem {
    /*class Domain;
    class SparseMtrx;
    class TimeStep;
    class DofManager;
    class GaussPoint;
    class UnknownNumberingScheme;
    class FloatMatrix;
    class IntegrationRule;
    class ContactElement;
    class Node;
*/
    class OOFEM_EXPORT Node2SegmentLagrangianMultiplierContact : public ActiveBoundaryCondition
    {
    private:
        bool useTangent; ///< Determines if tangent should be used.
        IntArray nodeSet;
        IntArray segmentSet;
        std::vector<  DofManager * >lmdm;
        int lm_num; ///< Determines the number of Lagrange multiplier DOFs (because all nodes are checked against all segments)
	FloatArray prescribedNormal;
    public:

        /// Constructor.
        Node2SegmentLagrangianMultiplierContact(int n, Domain *d) : ActiveBoundaryCondition(n, d) { }
        /// Destructor.
        virtual ~Node2SegmentLagrangianMultiplierContact() {};

        virtual IRResultType initializeFrom(InputRecord *ir);

        virtual void assemble(SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0) override;

        virtual void assembleVector(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms = NULL) override;

        virtual const char *giveClassName() const { return "Node2SegmentLagrangianMultiplierContact"; }
        virtual const char *giveInputRecordName() const { return _IFT_Node2SegmentLagrangianMultiplierContact_Name; }

        int giveNumberOfInternalDofManagers() override { return lm_num; }
        DofManager *giveInternalDofManager(int i) override { return this->lmdm.at(i - 1); }

        void giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;
        

    private:
        double computeTangentFromContact(FloatMatrix &answer, Node *node, ContactSegment *segment, TimeStep *tStep);

        void computeGap(double &answer, Node *node, ContactSegment *segment, TimeStep *tStep);

        void computeNvMatrixAt(FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *TimeStep);

        void computeExternalForcesFromContact(FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep);

        void giveLagrangianMultiplierLocationArray(const UnknownNumberingScheme &r_s, std::vector< IntArray > &answer);

    };
} // end namespace oofem
