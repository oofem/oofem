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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#ifndef gausspoint_h
#define gausspoint_h

#include "oofemenv.h"
#include "integrationrule.h"
#include "integrationpointstatus.h"
#include "element.h"
#include "floatarray.h"
#include "materialmode.h"
#include <map>

namespace oofem {
class Material;
class LayeredCrossSection;
class MicroplaneMaterial;
class FiberedCrossSection;
class CrossSection;
class IntegrationRule;

/**
 * Class representing integration point in finite element program.
 * Integration point maintain its space position, integration
 * weight and corresponding material mode.
 * Link to related element which integration point belongs to
 * is also maintained.
 * Integration point generally can contain list of slave integration points
 * therefore is called as master point. Slaves are used for example to implement
 * layered or fibered cross sections by cross section class. Then in one
 * "macro" master gauss point, cross section creates few slaves (one per layer)
 * and puts them into master list. When cross sections completes requests for
 * particular master integration point, it performs integration over layers.
 * It therefore calls material class for each layer, sending corresponding
 * slave as parameter and integrates results.
 * Generally, every integration point must hold
 * its own copy of history variables (which are related to corresponding
 * material model used). These material type dependent history variables
 * are stored in material type related material status, which can be
 * managed by integration point.
 * Each material model class should introduce related material status class
 * (derived from material status class or from its children), where necessary
 * history variables are kept and can be accessed by material.
 * Material class then creates unique copy of related status in all necessary
 * integration points. Because integration point is compulsory parameter of
 * all member functions of material class, particular material then can
 * easily access its associated status from integration point and therefore its
 * history variables for particular integration point.
 *
 * To provide support to integrate element contribution by parts from element subVolumes
 * (possibly with different material parameters etc), the integration point can maintain not
 * only its element natural coordinates, but also its subVolume local coordinates, that are necessary
 * to compute its jacobian, for example. These coordinates are stored in localCoordinates attribute.
 *
 */
class OOFEM_EXPORT GaussPoint
{
private:
    /// Number.
    int number;
    /// Reference to parent integration rule.
    IntegrationRule *irule;
    /// Natural Element Coordinates of receiver.
    FloatArray naturalCoordinates;
    /// Optional local sub-patch (sub-patches form element volume) coordinates of the receiver.
    std::unique_ptr<FloatArray> subPatchCoordinates;
    /// Optional global (Cartesian) coordinates
    std::unique_ptr<FloatArray> globalCoordinates;
    /// Integration weight.
    double weight;
    /// Material mode of receiver.
    MaterialMode materialMode;

protected:
    // layer and fibered material support
    /// List of slave integration points.
    std::vector< GaussPoint * > gaussPoints;
    /// Status of e.g. material in point
    std::map<int, std::unique_ptr<IntegrationPointStatus>> materialStatuses;

public:
    /**
     * Creates integration point belonging to given integration rule,
     * with given number, integration weight, coordinates and material mode.
     * @param ir Integration rule to which integration point belongs to.
     * @param n Integration point number.
     * @param iNaturalCoord Natueral coordinates.
     * @param w Integration weight.
     * @param mode Material mode.
     */
    GaussPoint(IntegrationRule *ir, int n, FloatArray iNaturalCoord, double w, MaterialMode mode);

    GaussPoint(IntegrationRule *ir, int n, double w, MaterialMode mode);

    ~GaussPoint();

    /// Returns i-th natural element coordinate of receiver
    double giveNaturalCoordinate(int i) const { return naturalCoordinates.at(i); }
    /// Returns coordinate array of receiver.
    const FloatArray &giveNaturalCoordinates() const { return naturalCoordinates; }
    void setNaturalCoordinates(const FloatArray &c) { naturalCoordinates = c; }

    /// Returns local sub-patch coordinates of the receiver
    const FloatArray &giveSubPatchCoordinates() const
    {
        if ( subPatchCoordinates ) {
            return *subPatchCoordinates;
        } else {
            return naturalCoordinates;
        }
    }
    void setSubPatchCoordinates(const FloatArray &c)
    {
        if ( subPatchCoordinates ) {
            * subPatchCoordinates = c;
        } else {
            subPatchCoordinates = std::make_unique<FloatArray>(c);
        }
    }

    inline const FloatArray &giveGlobalCoordinates()
    {
        if ( globalCoordinates ) {
            return *globalCoordinates;
        } else {
            globalCoordinates = std::make_unique<FloatArray>();
            this->giveElement()->computeGlobalCoordinates(*globalCoordinates, naturalCoordinates);
            return *globalCoordinates;
        }
    }

    void setGlobalCoordinates(const FloatArray &iCoord)
    {
        if ( globalCoordinates ) {
            *globalCoordinates = iCoord;
        } else {
            globalCoordinates = std::make_unique<FloatArray>(iCoord);
        }
    }

    /// Returns  integration weight of receiver.
    double giveWeight() { return weight; }
    void setWeight(double w) { weight = w; }
    /// Returns number of receiver.
    int giveNumber() { return number; }
    /// Returns corresponding integration rule to receiver.
    IntegrationRule *giveIntegrationRule() { return irule; }
    /// Returns corresponding element to receiver.
    Element *giveElement() { return irule->giveElement(); }

    /// Returns corresponding material mode of receiver.
    MaterialMode giveMaterialMode() { return this->materialMode; }
    /// Sets material mode of receiver.
    void setMaterialMode(MaterialMode newMode) { this->materialMode = newMode; }
    ///@todo giveMaterial routine most be removed from gauss-points, it doesn't fit with different types of cross-sections.

    /// Returns reference to material associated to related element of receiver.
    Material *giveMaterial() { return giveElement()->giveMaterial(); }

    /// Returns reference to cross section associated to related element of receiver.
    CrossSection *giveCrossSection() { return giveElement()->giveCrossSection(); }

    /**
     * Returns reference to associated material status (NULL if not defined).
     */
    IntegrationPointStatus *giveMaterialStatus(IntegrationPointStatusIDType key=IPSID_Default) { return this->materialStatuses.at(key).get(); }
    const IntegrationPointStatus *giveMaterialStatus(IntegrationPointStatusIDType key=IPSID_Default) const { return this->materialStatuses.at(key).get(); }
    bool hasMaterialStatus(IntegrationPointStatusIDType key=IPSID_Default) const { return this->materialStatuses.find(key) != this->materialStatuses.end(); }
    
    /**
     * Sets Material status managed by receiver.
     * @param ptr Pointer to new status of receiver.
     * @return Pointer to new status.
     */
    IntegrationPointStatus *setMaterialStatus(std::unique_ptr<IntegrationPointStatus> ptr, IntegrationPointStatusIDType key=IPSID_Default)
    {
        if ( this->materialStatuses.find(key) != this->materialStatuses.end() ) {
            OOFEM_ERROR("status already exist");
        }
        this->materialStatuses[key]=std::move(ptr);
        return this->materialStatuses[key].get();
    }
    /**
     * Returns index-th slave gauss point of receiver.
     * @param index Index of returned slave.
     * @return Slave gp.
     */
    GaussPoint *giveSlaveGaussPoint(int index);
    std::vector< GaussPoint * > &giveSlaveGaussPoints() {return this->gaussPoints;}

    /**
     * True if gauss point has slave points. Otherwise false.
     */
    bool hasSlaveGaussPoint();
    /**
     * Finds index of slave point in an array. Returns position. If not found, returns error.
     */
    size_t findFirstIndexOfSlaveGaussPoint(GaussPoint *gp);
     /**
     * Prints output of receiver to file. Corresponding printOutputAt  function for
     * associated status is called. The same function is also invoked for all available
     * slaves of receiver.
     */
    void printOutputAt(FILE *file, TimeStep *tStep, const char* indent="");
    /**
     * Updates internal state of receiver after finishing time step.
     * Material::updateYourself (receiver, tStep) function is called to
     * update material status. Same function is also invoked for
     * all receiver's slaves.
     */
    void updateYourself(TimeStep *tStep);

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "GaussPoint"; }

    /**
     * Sets Material status managed by receiver, this method is intended only for PyBind11 interface.
     * @param ptr Pointer to new status of receiver. Object will take ownership of the pointer.
     * @return Pointer to new status.
     */
    IntegrationPointStatus *__setMaterialStatus(IntegrationPointStatus* ptr, IntegrationPointStatusIDType key=IPSID_Default)
    {
        if ( this->materialStatuses.find(key) != this->materialStatuses.end() ) {
            OOFEM_ERROR("status already exist");
        }
        this->materialStatuses[key]=std::unique_ptr<IntegrationPointStatus>(ptr);
        return this->materialStatuses[key].get();
    }

    friend class LayeredCrossSection;
    friend class FiberedCrossSection;
};

typedef GaussPoint IntegrationPoint;
} // end namespace oofem
#endif // gausspoint_h
