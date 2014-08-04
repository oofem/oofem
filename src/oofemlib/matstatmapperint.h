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

#ifndef MATSTATMAPPERINT_H_
#define MATSTATMAPPERINT_H_

namespace oofem {
class MaterialStatus;
class StructuralInterfaceMaterialStatus;
class MaterialMappingAlgorithm;
class GaussPoint;
class Domain;
class TimeStep;
class Set;

/**
 * matstatmapperint.h
 *
 * MaterialStatusMapperInterface:
 * An interface class for MaterialStatus. The purpose is to allow
 * mapping of state variables (Gauss point variables) in a
 * generic way. This is useful when the mesh changes, e.g. due to
 * propagating cracks or in adaptive analysis.
 *
 * To add adaptivity capability to a new material model, only
 * copyStateVariables() and addStateVariables need to be overloaded
 * in the corresponding material status class. Hence, mapping
 * functionality can be implemented for new material models with
 * minimum effort.
 *
 * @note bp: this is however renundant to existing MaterialModelMapperInterface (MMI_map, MMI_update),
 * which is perhaps more flexible (several mapping algorithms, etc).
 *
 * @author Erik Svenning
 *  Created on: Nov 6, 2013
 */
class MaterialStatusMapperInterface
{
public:
    MaterialStatusMapperInterface();
    virtual ~MaterialStatusMapperInterface();

    virtual void copyStateVariables(const MaterialStatus &iStatus) = 0;
    virtual void addStateVariables(const MaterialStatus &iStatus) = 0;
    //    virtual void callCopyStateVariables(MaterialStatusMapperInterface &oStatus) = 0;


    /**
     * Maps all internal state variables from
     * the old domain to the given gp status.
     * @param iGP Integration point belonging to the new domain.
     * @param iOldDom Old domain.
     * @param iTStep Time step.
     * @return Nonzero if o.k.
     */
    virtual int MSMI_map(const GaussPoint &iGP, const Domain &iOldDom, Set &sourceSet, const TimeStep &iTStep, MaterialStatus &oStatus);

    virtual int MSMI_map_cz(const GaussPoint &iGP, const Domain &iOldDom, Set &sourceSet, const TimeStep &iTStep, MaterialStatus &oStatus);
    /**
     * Updates the internal state variables from previously mapped values.
     * @param iGP Integration point belonging to the new domain.
     * @param iTStep Time step.
     * @return Nonzero if o.k.
     */
    virtual int MSMI_update(const GaussPoint &iGP, const TimeStep &iTStep);
    /**
     * Finishes the mapping for given time step. Used to perform cleanup.
     */
    virtual int MSMI_finish(const TimeStep &iTStep);

protected:
    MaterialMappingAlgorithm *mpMaterialMapper;
};
} /* namespace oofem */
#endif /* MATSTATMAPPERINT_H_ */
