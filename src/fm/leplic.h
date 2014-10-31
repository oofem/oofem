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

#ifndef leplic_h
#define leplic_h

#include "materialinterface.h"
#include "floatarray.h"
#include "mathfem.h"
#include "interface.h"

///@name Input fields for LEPLIC
//@{
#define _IFT_LEPLIC_refVol "refvol"
//@}

namespace oofem {
class LEPlic;
class Polygon;
class Element;

/**
 * Element interface for LEPlic class representing Lagrangian-Eulerian (moving) material interface.
 * The elements should provide specific functionality in order to collaborate with LEPlic and this
 * required functionality is declared in this interface.
 */
class LEPlicElementInterface : public Interface
{
protected:
    bool permanentVofFlag;
    /// Volume fraction of reference fluid in element.
    double vof, temp_vof;
    /// Line constant of line segment representing interface.
    double p, temp_p;
    /// Interface segment normal.
    FloatArray normal, temp_normal;

public:
    LEPlicElementInterface() : permanentVofFlag(false), vof(0.), temp_vof(0.) { }
    /**
     * @name The element interface required by LEPlicElementInterface
     */
    //@{
    /// Computes corresponding volume fraction to given interface position.
    virtual double computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag) = 0;
    /// Assembles the true element material polygon (takes receiver vof into accout).
    virtual void formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                        const FloatArray &normal, const double p, bool updFlag) = 0;
    /// Assembles receiver material polygon based solely on given interface line.
    virtual void formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                         const FloatArray &normal, const double p, bool updFlag) = 0;
    /// Truncates given material polygon to receiver.
    virtual double truncateMatVolume(const Polygon &matvolpoly, double &volume) = 0;
    /// Computes the receiver center (in updated Lagrangian configuration).
    virtual void giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool updFlag) = 0;
    /// Assembles receiver volume.
    virtual void formMyVolumePoly(Polygon &myPoly, LEPlic *mat_interface, bool updFlag) = 0;
    /// Computes the volume of receiver.
    virtual double computeMyVolume(LEPlic *matInterface, bool updFlag) = 0;
    /// Returns true if cell is boundary.
    bool isBoundary();
    /// Return number of receiver's element.
    virtual Element *giveElement() = 0;
    /// Computes critical time step.
    virtual double computeCriticalLEPlicTimeStep(TimeStep *tStep) = 0;
    //@}

    void setTempLineConstant(double tp) { temp_p = tp; }
    void setTempInterfaceNormal(FloatArray tg) { temp_normal = std::move(tg); }
    void setTempVolumeFraction(double v) {
        if ( !permanentVofFlag ) {
            temp_vof = v;
        }
    }
    void setPermanentVolumeFraction(double v) {
        temp_vof = vof = v;
        permanentVofFlag = true;
    }
    void addTempVolumeFraction(double v) {
        if ( !permanentVofFlag ) {
            temp_vof += v;
            if ( temp_vof > 1.0 ) {
                temp_vof = 1.0;
            }
        }
    }
    double giveVolumeFraction() { return vof; }
    double giveTempVolumeFraction() { return temp_vof; }
    void giveTempInterfaceNormal(FloatArray &n) { n = temp_normal; }
    double giveTempLineConstant() { return temp_p; }
    void updateYourself(TimeStep *tStep) {
        vof = temp_vof;
        p = temp_p;
        normal = temp_normal;
    }

    /**
     * Stores context of receiver into given stream.
     * Only non-temp internal history variables are stored.
     * @param stream stream where to write data
     * @param mode determines amount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream.
     * @param stream stream where to read data
     * @param mode determines amount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};

/**
 * Abstract base class representing Lagrangian-Eulerian (moving) material interfaces.
 * Its typical use to model moving interface (such as free surface)
 * in a fixed-grid methods (as typically used in CFD).
 * The basic tasks are representation of interface and its updating.
 */
class LEPlic : public MaterialInterface
{
protected:
    /// Array used to store updated x-coordinates of nodes as moved along streamlines
    FloatArray updated_XCoords;
    /// Array used to store updated y-coordinates of nodes as moved along streamlines
    FloatArray updated_YCoords;
    double orig_reference_fluid_volume;
public:
    /**
     * Constructor. Takes two two arguments. Creates
     * MaterialInterface instance with given number and belonging to given domain.
     * @param n Component number in particular domain. For instance, can represent
     * node number in particular domain.
     * @param d Domain to which component belongs to.
     */
    LEPlic(int n, Domain * d) : MaterialInterface(n, d) {
        orig_reference_fluid_volume = 0.0;
    }

    virtual void updatePosition(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep) { }
    virtual void giveMaterialMixtureAt(FloatArray &answer, FloatArray &position);
    virtual void giveElementMaterialMixture(FloatArray &answer, int ielem);
    virtual double giveNodalScalarRepresentation(int);
    virtual double computeCriticalTimeStep(TimeStep *tStep);

    /**
     * Returns updated nodal positions.
     */
    void giveUpdatedCoordinate(FloatArray &answer, int num) {
        answer.resize(2);
        answer.at(1) = updated_XCoords.at(num);
        answer.at(2) = updated_YCoords.at(num);
    }

    double giveUpdatedXCoordinate(int num) { return updated_XCoords.at(num); }
    double giveUpdatedYCoordinate(int num) { return updated_YCoords.at(num); }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    // identification
    virtual const char *giveClassName() const { return "LEPlic"; }

protected:
    void doLagrangianPhase(TimeStep *tStep);
    void doInterfaceReconstruction(TimeStep *tStep, bool coord_upd, bool temp_vof);
    void doInterfaceRemapping(TimeStep *tStep);
    void doCellDLS(FloatArray &fvgrad, int ie, bool coord_upd, bool temp_vof_flag);
    void findCellLineConstant(double &p, FloatArray &fvgrad, int ie, bool coord_upd, bool temp_vof_flag);


    class computeLEPLICVolumeFractionWrapper
    {
protected:
        LEPlicElementInterface *iface;
        LEPlic *minterf;
        FloatArray normal;
        double target_vof;
        bool upd;
public:
        computeLEPLICVolumeFractionWrapper(LEPlicElementInterface * i, LEPlic * mi, const FloatArray & n, const double target_vof_val, bool upd_val) :
            iface(i), minterf(mi), normal(n), target_vof(target_vof_val), upd(upd_val) { }
        void setNormal(FloatArray n) { normal = std::move(n); }
        double eval(double x) { return fabs(iface->computeLEPLICVolumeFraction(normal, x, minterf, upd) - target_vof); }
    };
};
} // end namespace oofem
#endif // leplic_h
