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

#ifndef SRC_OOFEMLIB_XFEM_NUCLEATIONCRITERION_H_
#define SRC_OOFEMLIB_XFEM_NUCLEATIONCRITERION_H_


#include <memory>
#include <vector>

namespace oofem {

class EnrichmentItem;
class Domain;
class DataReader;
class DynamicDataReader;
class InputRecord;
class EnrichmentFunction;

class NucleationCriterion
{
public:
    NucleationCriterion(Domain *ipDomain);
    virtual ~NucleationCriterion();

    virtual std::vector<std::unique_ptr<EnrichmentItem>> nucleateEnrichmentItems();

    virtual void initializeFrom(InputRecord &ir);
    virtual int instanciateYourself(DataReader &dr);
    virtual void postInitialize() {}

    virtual void appendInputRecords(DynamicDataReader &oDR);

    /// @return Class name of the receiver.
    virtual const char *giveClassName() const = 0;
    /// @return Input record name of the receiver.
    virtual const char *giveInputRecordName() const = 0;

protected:
    Domain *mpDomain;
    std::unique_ptr<EnrichmentFunction> mpEnrichmentFunc;
};

} /* namespace oofem */

#endif /* SRC_OOFEMLIB_XFEM_NUCLEATIONCRITERION_H_ */
