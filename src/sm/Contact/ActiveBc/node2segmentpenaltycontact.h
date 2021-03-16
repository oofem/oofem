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
#ifndef node2nodepenaltycontact_h
#define node2nodepenatycontact_h


#include "activebc.h"
#include "sm/Contact/ContactSegment/contactsegment.h"


///@name Input fields for _IFT_ContactElement
//@{
#define _IFT_Node2SegmentPenaltyContact_Name "n2spenaltycontact"
#define _IFT_Node2SegmentPenaltyContact_penalty "penalty"
#define _IFT_Node2SegmentPenaltyContact_useTangent "usetangent"

#define _IFT_Node2SegmentPenaltyContact_segmentSet "segmentset"
#define _IFT_Node2SegmentPenaltyContact_nodeSet "nodeset"

#define _IFT_Node2SegmentPenaltyContact_prescribedNormal "normal"

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

class OOFEM_EXPORT Node2SegmentPenaltyContact : public ActiveBoundaryCondition
{
private:
  bool useTangent; ///< Determines if tangent should be used.
  double penalty;
  IntArray nodeSet;
  IntArray segmentSet;
  FloatArray prescribedNormal;
public:

    /// Constructor.
    Node2SegmentPenaltyContact(int n, Domain *d) : ActiveBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~Node2SegmentPenaltyContact() {};

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void assemble(SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0) override;

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms = NULL) override;


    virtual const char *giveClassName() const { return "Node2SegmentPenaltyContact"; }
    virtual const char *giveInputRecordName() const { return _IFT_Node2SegmentPenaltyContact_Name; }



    void giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;

private:
    void  computeTangentFromContact(FloatMatrix & answer, Node * node, ContactSegment * segment, TimeStep * tStep);
    void computeGap(double &answer, Node *node, ContactSegment *segment, TimeStep *tStep);
    void computeExternalForcesFromContact(FloatArray & answer, Node * node, ContactSegment * segment, TimeStep * tStep);
    void giveLocationArray(IntArray &loc, const UnknownNumberingScheme &ns, const Node *node, const ContactSegment *segment);
    void computeNvMatrixAt(FloatArray & answer, Node * node, ContactSegment * segment, TimeStep * tStep);
    void computeTvMatrixAt( FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep );
    void computeBvMatrixAt( FloatArray &answer, Node *node, ContactSegment *segment, TimeStep *tStep );

};
} // end namespace oofem
#endif // node2nodecontact_h
