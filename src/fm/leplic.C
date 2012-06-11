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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "leplic.h"
#include "mathfem.h"
#include "timestep.h"
#include "geotoolbox.h"
#include "node.h"
#include "conTable.h"
#include "datastream.h"
#include "spatiallocalizer.h"
#include "contextioerr.h"
#include "element.h"

namespace oofem {
#define LEPLIC_ZERO_VOF  1.e-12
#define LEPLIC_BRENT_EPS 1.e-12



bool
LEPlicElementInterface :: isBoundary()
{
    int i, nneighbr, ineighbr;
    double fvk, fvi = this->giveTempVolumeFraction();
    IntArray currCell(1), neighborList;
    LEPlicElementInterface *ineghbrInterface;
    Domain *domain = this->giveElement()->giveDomain();
    ConnectivityTable *contable = domain->giveConnectivityTable();
    if ( ( fvi > 0. ) && ( fvi <= 1.0 ) ) {
        // potentially boundary cell
        if ( ( fvi > 0. ) && ( fvi < 1.0 ) ) {
            return true;
        }

        currCell.at(1) = this->giveElement()->giveNumber();
        contable->giveElementNeighbourList(neighborList, currCell);
        // loop over neighbors to assemble normal equations
        nneighbr = neighborList.giveSize();
        for ( i = 1; i <= nneighbr; i++ ) {
            ineighbr = neighborList.at(i);
            if ( ( ineghbrInterface =
                      ( LEPlicElementInterface * ) ( domain->giveElement(ineighbr)->giveInterface(LEPlicElementInterfaceType) ) ) ) {
                fvk = ineghbrInterface->giveTempVolumeFraction();
                if ( fvk < 1.0 ) {
                    return true;
                }
            }
        }
    }

    return false;
}


contextIOResultType
LEPlicElementInterface :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // write a raw data
    if ( !stream->write(& vof, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& p, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = normal.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

contextIOResultType
LEPlicElementInterface :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read raw data
    if ( !stream->read(& vof, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    temp_vof = vof;
    if ( !stream->read(& p, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    temp_p = p;
    if ( ( iores = normal.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    temp_normal = normal;
    return CIO_OK;
}





void
LEPlic :: updatePosition(TimeStep *atTime)
{
#ifdef __OOFEG
    //deleteLayerGraphics(OOFEG_DEBUG_LAYER);
    EVFastRedraw(myview);
#endif
    this->doLagrangianPhase(atTime);
    this->doInterfaceReconstruction(atTime, true, false);
#ifdef __OOFEG
    //ESIEventLoop (YES, "doInterfaceReconstruction Finished; Press Ctrl-p to continue");
    deleteLayerGraphics(OOFEG_DEBUG_LAYER);
#endif
    this->doInterfaceRemapping(atTime);
    // here the new VOF values are determined, now we call doInterfaceReconstruction
    // to reconstruct interface (normal, constant) on original grid
    this->doInterfaceReconstruction(atTime, false, true);
#ifdef __OOFEG
    ESIEventLoop( NO, oofem_tmpstr("doInterfaceReconstruction Finished; Press Ctrl-p to continue") );
    //ESIEventLoop (YES, "doInterfaceReconstruction Finished; Press Ctrl-p to continue");
#endif
}

void
LEPlic :: doLagrangianPhase(TimeStep *atTime)
{
    //Maps element nodes along trajectories using basic Runge-Kutta method (midpoint rule)
    int i, ci, ndofman = domain->giveNumberOfDofManagers();
    int nsd = 2;
    double dt = atTime->giveTimeIncrement();
    DofManager *dman;
    Node *inode;
    IntArray velocityMask(2);
    FloatArray x, x2(nsd), v_t, v_tn1;
    FloatMatrix t;
#if 1
    EngngModel *emodel = domain->giveEngngModel();
    int err;
#endif
    velocityMask.at(1) = V_u;
    velocityMask.at(2) = V_v;

    updated_XCoords.resize(ndofman);
    updated_YCoords.resize(ndofman);


    for ( i = 1; i <= ndofman; i++ ) {
        dman = domain->giveDofManager(i);
        // skip dofmanagers with no position information
        if ( ( dman->giveClassID() != NodeClass ) && ( dman->giveClassID() != RigidArmNodeClass ) && ( dman->giveClassID() != HangingNodeClass ) ) {
            continue;
        }

        inode = ( Node * ) dman;
        // get node coordinates
        x = * ( inode->giveCoordinates() );
        // get velocity field v(tn, x(tn)) for dof manager

#if 1
        /* Original version */
        dman->giveUnknownVector( v_t, velocityMask, EID_MomentumBalance, VM_Total, atTime->givePreviousStep() );
        /* Modified version */
        //dman->giveUnknownVector(v_t, velocityMask, EID_MomentumBalance, VM_Total, atTime);

        // Original version
        // compute updated position x(tn)+0.5*dt*v(tn,x(tn))
        for ( ci = 1; ci <= nsd; ci++ ) {
            x2.at(ci) = x.at(ci) + 0.5 *dt *v_t.at(ci);
        }

        // compute interpolated velocity field at x2 [ v(tn+1, x(tn)+0.5*dt*v(tn,x(tn))) = v(tn+1, x2) ]
        Field *vfield;
        vfield = emodel->giveContext()->giveFieldManager()->giveField(FT_Velocity);
        if ( vfield == NULL ) {
            _error("doLagrangianPhase: Velocity field not available");
        }

        err = vfield->evaluateAt(v_tn1, x2, VM_Total, atTime);
        if ( err == 1 ) {
            // point outside domain -> be explicit
            v_tn1 = v_t;
        } else if ( err != 0 ) {
            _error2("doLagrangianPhase: vfield->evaluateAt failed, error code %d", err);
        }

        // compute final updated position
        for ( ci = 1; ci <= nsd; ci++ ) {
            x2.at(ci) = x.at(ci) + dt *v_tn1.at(ci);
        }

#else
        // pure explicit version
        dman->giveUnknownVector(v_t, velocityMask, EID_MomentumBalance, VM_Total, atTime);

        for ( ci = 1; ci <= nsd; ci++ ) {
            x2.at(ci) = x.at(ci) + dt *v_t.at(ci);
        }

#endif
        // store updated node position
        updated_XCoords.at(i) = x2.at(1);
        updated_YCoords.at(i) = x2.at(2);
    }
}

void
LEPlic :: doInterfaceReconstruction(TimeStep *atTime, bool coord_upd, bool temp_vof)
{
    /* Here volume materials are reconstructed on the new Lagrangian grid */

    int ie, nelem = domain->giveNumberOfElements();
    double p;
    FloatMatrix lhs(2, 2);

    FloatArray rhs(2), fvgrad(2);
    LEPlicElementInterface *interface;


    // loop over elements
    for ( ie = 1; ie <= nelem; ie++ ) {
        /* STEP 1: first do DLS (Differential least square reconstruction) */
        this->doCellDLS(fvgrad, ie, coord_upd, temp_vof);
        /* STEP 2: Finding the line constant */
        this->findCellLineConstant(p, fvgrad, ie, coord_upd, temp_vof);

        interface = ( LEPlicElementInterface * ) domain->giveElement(ie)->giveInterface(LEPlicElementInterfaceType);
        interface->setTempLineConstant(p);
        interface->setTempInterfaceNormal(fvgrad);
    } // end loop over all elements

}


void
LEPlic :: doInterfaceRemapping(TimeStep *atTime)
{
    /*
     * Final step: deposition of volume materials truncated on Lagrangian (updated)
     * grid to the target grid, which is the original one in our Eulerian case.
     */
    int in = 0, ie = 0, neighbrNum, nelem = domain->giveNumberOfElements();
    double in_vof, total_volume = 0.0, in_vol;
    IntArray neighbours, elNum(1);
    FloatArray normal;
    Polygon matvolpoly, elemPoly;
    Graph g;

    double matVol = 0.0, matVolSum = 0.0;
    double __vol;

    LEPlicElementInterface *interface, *neghbrInterface;
    // loop over elements
    for ( ie = 1; ie <= nelem; ie++ ) {
        if ( ( interface = ( LEPlicElementInterface * ) ( domain->giveElement(ie)->giveInterface(LEPlicElementInterfaceType) ) ) ) {
            interface->setTempVolumeFraction(0.0);
        }
    }

    // loop over elements
    for ( ie = 1; ie <= nelem; ie++ ) {
        //fprintf (stderr, "doInterfaceRemapping: processing elem %d\n", ie);

        // examine only neighbours -> this is the limit on time step
        elNum.at(1) = ie;
        domain->giveConnectivityTable()->giveElementNeighbourList(neighbours, elNum);
        // form polygon of material volume on Lagrangian element
        if ( ( interface = ( LEPlicElementInterface * ) ( domain->giveElement(ie)->giveInterface(LEPlicElementInterfaceType) ) ) ) {
            if ( interface->giveVolumeFraction() > LEPLIC_ZERO_VOF ) {
                interface->giveTempInterfaceNormal(normal);
                interface->formMaterialVolumePoly(matvolpoly, this, normal, interface->giveTempLineConstant(), true);
                matVol = matvolpoly.computeVolume();

                // internal test of incompressibility
                __vol = interface->computeMyVolume(this, false);

#ifdef __OOFEG
                //EASValsSetColor(::gc[OOFEG_DEBUG_LAYER].getElementColor());
                //GraphicObj *go = matvolpoly.draw(::gc[OOFEG_DEBUG_LAYER],true);
                //EVHiliteGraphics (myview, go);
                //ESIEventLoop (YES, "Press Ctrl-p to continue");
                //EVFastRedraw(myview);
#endif


                try {
                    matVolSum = 0.0;
                    // loop over neighbours to truncate material volume on target (original) grid
                    for ( in = 1; in <= neighbours.giveSize(); in++ ) {
                        neighbrNum = neighbours.at(in);
                        if ( ( neghbrInterface = ( LEPlicElementInterface * )
                                                 ( domain->giveElement(neighbrNum)->giveInterface(LEPlicElementInterfaceType) ) ) ) {
                            in_vof = neghbrInterface->truncateMatVolume(matvolpoly, in_vol);
                            neghbrInterface->addTempVolumeFraction(in_vof);
                            total_volume += in_vol;
                            matVolSum += in_vol;
                        }
                    }
                } catch(GT_Exception & c) {
                    c.print();

                    neighbrNum = neighbours.at(in);
                    if ( ( neghbrInterface = ( LEPlicElementInterface * )
                                             ( domain->giveElement(neighbrNum)->giveInterface(LEPlicElementInterfaceType) ) ) ) {
                        in_vof = neghbrInterface->truncateMatVolume(matvolpoly, in_vol);
                        neghbrInterface->addTempVolumeFraction(in_vof);
                    }
                }

#ifdef __OOFEG
                //ESIEventLoop (YES, "Step Finished; Press Ctrl-p to continue");
#endif
                double err = fabs(matVol - matVolSum) / matVol;
                if ( ( err > 1.e-12 ) && ( fabs(matVol - matVolSum) > 1.e-4 ) && ( matVol > 1.e-6 ) ) {
                    OOFEM_WARNING4("LEPlic::doInterfaceRemapping:  volume inconsistency %5.2f%%\n\ttstep %d, element %d\n", err * 100, atTime->giveNumber(), ie);
                }

#if 0
                if ( ( err > 2.e-3 ) && ( fabs(matVol - matVolSum) > 2.e-3 ) ) {
                    //debug

 #ifdef __OOFEG
                    //ESIEventLoop (YES, "Press Ctrl-p to continue");
                    deleteLayerGraphics(OOFEG_DEBUG_LAYER);
                    EASValsSetColor( :: gc [ OOFEG_DEBUG_LAYER ].getElementColor() );
                    //GraphicObj *go = matvolpoly.draw(::gc[OOFEG_DEBUG_LAYER],true);
                    matvolpoly.draw(:: gc [ OOFEG_DEBUG_LAYER ], true);
                    //EVHiliteGraphics (myview, go);
                    //ESIEventLoop (YES, "Press Ctrl-p to continue");
                    EVFastRedraw(myview);
 #endif

                    matVolSum = 0.0;
                    // loop over neighbours to truncate material volume on target (original) grid
                    for ( in = 1; in <= neighbours.giveSize(); in++ ) {
                        neighbrNum = neighbours.at(in);
                        if ( neghbrInterface = ( LEPlicElementInterface * )
                                               ( domain->giveElement(neighbrNum)->giveInterface(LEPlicElementInterfaceType) ) ) {
                            in_vof = neghbrInterface->truncateMatVolume(matvolpoly, in_vol);
                            neghbrInterface->addTempVolumeFraction(in_vof);
                            total_volume += in_vol;
                            matVolSum += in_vol;
                        }
                    }

                    matVolSum = 0.0;
                    // loop over neighbours to truncate material volume on target (original) grid
                    for ( in = 1; in <= neighbours.giveSize(); in++ ) {
                        neighbrNum = neighbours.at(in);
                        if ( neghbrInterface = ( LEPlicElementInterface * )
                                               ( domain->giveElement(neighbrNum)->giveInterface(LEPlicElementInterfaceType) ) ) {
                            in_vof = neghbrInterface->truncateMatVolume(matvolpoly, in_vol);
                            neghbrInterface->addTempVolumeFraction(in_vof);
                            total_volume += in_vol;
                            matVolSum += in_vol;
                        }
                    }

                    matVolSum = 0.0;
                    // loop over neighbours to truncate material volume on target (original) grid
                    for ( in = 1; in <= neighbours.giveSize(); in++ ) {
                        neighbrNum = neighbours.at(in);
                        if ( neghbrInterface = ( LEPlicElementInterface * )
                                               ( domain->giveElement(neighbrNum)->giveInterface(LEPlicElementInterfaceType) ) ) {
                            in_vof = neghbrInterface->truncateMatVolume(matvolpoly, in_vol);
                            neghbrInterface->addTempVolumeFraction(in_vof);
                            total_volume += in_vol;
                            matVolSum += in_vol;
                        }
                    }

                    /*
                     * ESIEventLoop (YES, "Press Ctrl-p to exit");
                     * ESIEventLoop (YES, "Press Ctrl-p to exit");
                     */
                    exit(1);
                }

#endif
            }
        } else {
            OOFEM_ERROR("LEPlic::doInterfaceRemapping: Element with no LEPlicInterface support encountered");
        }
    } // end loop over elements

    // loop over elements
    for ( ie = 1; ie <= nelem; ie++ ) {
        if ( ( interface = ( LEPlicElementInterface * ) ( domain->giveElement(ie)->giveInterface(LEPlicElementInterfaceType) ) ) ) {
            if ( interface->giveTempVolumeFraction() > 1.0 ) {
                OOFEM_LOG_INFO("LEPlic::doInterfaceRemapping - Element %d: vof out of range, vof =%e",  ie, interface->giveTempVolumeFraction() );
            }

            if ( interface->giveTempVolumeFraction() >= 0.99999999 ) {
                interface->setTempVolumeFraction(1.0);
            }
        }
    }

    OOFEM_LOG_INFO("LEPlic::doInterfaceRemapping: Total volume is %e", total_volume);
    if ( orig_reference_fluid_volume > 0.0 ) {
        OOFEM_LOG_INFO("LEPlic::doInterfaceRemapping: Volume error is %5.2f%%",
                       ( ( total_volume - orig_reference_fluid_volume ) / orig_reference_fluid_volume ) * 100.0);
    }


#ifdef __OOFEG
    //ESIEventLoop (YES, "Step Finished; Press Ctrl-p to continue");
    //deleteLayerGraphics(OOFEG_DEBUG_LAYER);
#endif
}


void
LEPlic :: doCellDLS(FloatArray &fvgrad, int ie, bool coord_upd, bool vof_temp_flag)
{
    int i, ineighbr, nneighbr;
    double fvi, fvk, wk, dx, dy;
    bool isBoundaryCell = false;
    LEPlicElementInterface *interface, *ineghbrInterface;
    FloatMatrix lhs(2, 2);
    FloatArray rhs(2), xi(2), xk(2);
    IntArray currCell(1), neighborList;
    ConnectivityTable *contable = domain->giveConnectivityTable();

    if ( ( interface = ( LEPlicElementInterface * ) ( domain->giveElement(ie)->giveInterface(LEPlicElementInterfaceType) ) ) ) {
        if ( vof_temp_flag ) {
            fvi = interface->giveTempVolumeFraction();
        } else {
            fvi = interface->giveVolumeFraction();
        }

        if ( ( fvi > 0. ) && ( fvi <= 1.0 ) ) {
            // potentially boundary cell

            if ( ( fvi > 0. ) && ( fvi < 1.0 ) ) {
                isBoundaryCell = true;
            }

            /* DLS (Differential least square reconstruction)
             *
             * In the DLS method, volume fraction Taylor series expansion of vf (volume fraction)
             * is formed from each reference cell volume fraction vf at element center x(i) to each
             * cell neighbor at point x(k). The sum (vf(i)-vf(k))^2 over all immediate neighbors
             * is then minimized inthe least square sense.
             */
            // get list of neighbours to current cell including current cell
            currCell.at(1) = ie;
            contable->giveElementNeighbourList(neighborList, currCell);
            // loop over neighbors to assemble normal equations
            nneighbr = neighborList.giveSize();
            interface->giveElementCenter(this, xi, coord_upd);
            lhs.zero();
            rhs.zero();
            for ( i = 1; i <= nneighbr; i++ ) {
                ineighbr = neighborList.at(i);
                if ( ineighbr == ie ) {
                    continue;         // skip itself
                }

                if ( ( ineghbrInterface =
                          ( LEPlicElementInterface * ) ( domain->giveElement(ineighbr)->giveInterface(LEPlicElementInterfaceType) ) ) ) {
                    if ( vof_temp_flag ) {
                        fvk = ineghbrInterface->giveTempVolumeFraction();
                    } else {
                        fvk = ineghbrInterface->giveVolumeFraction();
                    }

                    if ( fvk < 1.0 ) {
                        isBoundaryCell = true;
                    }

                    ineghbrInterface->giveElementCenter(this, xk, coord_upd);
                    wk = xk.distance(xi);
                    dx = ( xk.at(1) - xi.at(1) ) / wk;
                    dy = ( xk.at(2) - xi.at(2) ) / wk;
                    lhs.at(1, 1) += dx * dx;
                    lhs.at(1, 2) += dx * dy;
                    lhs.at(2, 2) += dy * dy;

                    rhs.at(1) += ( fvi - fvk ) * dx / wk;
                    rhs.at(2) += ( fvi - fvk ) * dy / wk;
                }
            }

            if ( isBoundaryCell ) {
                // symmetry
                lhs.at(2, 1) = lhs.at(1, 2);

                // solve normal equation for volume fraction gradient
                lhs.solveForRhs(rhs, fvgrad);

                // compute unit normal
                fvgrad.normalize();
                fvgrad.negated();
#ifdef __OOFEG
                /*
                 * EASValsSetLayer(OOFEG_DEBUG_LAYER);
                 * WCRec p[2];
                 * double tx = -fvgrad.at(2), ty=fvgrad.at(1);
                 * p[0].x=xi.at(1)-tx*0.1;
                 * p[0].y=xi.at(2)-ty*0.1;
                 * p[1].x=xi.at(1)+tx*0.1;
                 * p[1].y=xi.at(2)+ty*0.1;
                 * p[0].z = p[1].z = 0.0;
                 * GraphicObj *go = CreateLine3D(p);
                 * EGWithMaskChangeAttributes(LAYER_MASK, go);
                 * EMAddGraphicsToModel(ESIModel(), go);
                 * ESIEventLoop (YES, "Cell DLS finished; Press Ctrl-p to continue");
                 */
#endif
            } else {
                fvgrad.zero();
            }
        }
    }
}

void
LEPlic :: findCellLineConstant(double &p, FloatArray &fvgrad, int ie, bool coord_upd, bool temp_vof_flag)
{
    /* The line constatnt p is solved from the general non-linear function
     * F(p) = V(p)-V = 0,
     * where V(p) is the truncated volume resulting from the intersection between
     *            assumed line segment and the reference cell
     *       V    is the given Volume of Fluid ratio
     * To find zero of this function, Brent's method is been used.
     */
    Element *elem = domain->giveElement(ie);
    LEPlicElementInterface *interface = ( LEPlicElementInterface * ) elem->giveInterface(LEPlicElementInterfaceType);
    int ivert, nelemnodes = elem->giveNumberOfNodes();
    double ivof, fvi;
    double ivx, ivy, pp, target_vof;
    if ( temp_vof_flag ) {
        target_vof = interface->giveTempVolumeFraction();
    } else {
        target_vof = interface->giveVolumeFraction();
    }

    computeLEPLICVolumeFractionWrapper wrapper(interface, this, fvgrad, target_vof, coord_upd);
    /*
     * Initial part: find lower and uper bounds for Brent algorithm
     * Here lines with given interface normal are passed through each vertex
     * and corresponding volume fractions are computed. The two lines forming
     * truncation volumes that bound the actual material in the cell provide
     * upper and lower bounds for line constant
     */
    double upper_vof = 10.0, lower_vof = -10.0;
    double upper_p = 0.0, lower_p = 0.0;

    if ( temp_vof_flag ) {
        fvi = interface->giveTempVolumeFraction();
    } else {
        fvi = interface->giveVolumeFraction();
    }

    if ( ( fvi > LEPLIC_ZERO_VOF ) && ( fvi < 1.0 ) ) {
        // boundary cell


        for ( ivert = 1; ivert <= nelemnodes; ivert++ ) {
            if ( coord_upd ) {
                ivx = giveUpdatedXCoordinate( elem->giveNode(ivert)->giveNumber() );
                ivy = giveUpdatedYCoordinate( elem->giveNode(ivert)->giveNumber() );
            } else {
                ivx = elem->giveNode(ivert)->giveCoordinate(1);
                ivy = elem->giveNode(ivert)->giveCoordinate(2);
            }

            // determine line constant for vertex ivert
            pp = -( fvgrad.at(1) * ivx + fvgrad.at(2) * ivy );
            ivof = interface->computeLEPLICVolumeFraction(fvgrad, pp, this, coord_upd);
            if ( ( ( ivof - target_vof ) >= 0. ) && ( ivof < upper_vof ) ) {
                upper_vof = ivof;
                upper_p = pp;
            } else if ( ( ( target_vof - ivof ) >= 0. ) && ( ivof > lower_vof ) ) {
                lower_vof = ivof;
                lower_p = pp;
            }
        }

        if ( ( lower_vof >= 0. ) && ( upper_vof <= 1.00000000001 ) ) {
            // now use brent's method to find minima of V(p)-V function
            brent(lower_p, 0.5 * ( lower_p + upper_p ), upper_p,
                  mem_fun< computeLEPLICVolumeFractionWrapper >(& wrapper, & computeLEPLICVolumeFractionWrapper :: eval),
                  LEPLIC_BRENT_EPS, p);
            //interface->setTempLineConstant (p);


#ifdef __OOFEG
            /*
             * Polygon grp, cell;
             * // check here
             * //Polygon testvolpoly;
             * //interface->formMaterialVolumePoly(testvolpoly, this, fvgrad, p, true);
             * //double debug_vof = interface->truncateMatVolume(testvolpoly);
             *
             * //printf ("Cell %d-> target_vof is %e, debug val is %e\n", ie, target_vof, debug_vof);
             * interface->formMyVolumePoly (cell, this, true);
             * GraphicObj *goc = cell.draw(::gc[0], false);
             * EASValsSetLayer(OOFEG_DEBUG_LAYER);
             * EASValsSetColor(::gc[0].getBcIcColor());
             * EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK, goc);
             * EMDrawGraphics (ESIModel(), goc);
             *
             * interface->formVolumeInterfacePoly (grp, this, fvgrad, p, true);
             * GraphicObj *go = grp.draw(::gc[0], true);
             * EASValsSetColor(::gc[0].getActiveCrackColor());
             * EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK, go);
             * EMDrawGraphics (ESIModel(), go);
             * //ESIEventLoop (YES, "findCellLineConstant -> Press Ctrl-p to continue");
             */
#endif
        } else {
            OOFEM_ERROR("LEPlic::findCellLineConstant: finding lower and uper bounds of line constant value failed");
        }
    }
}

IRResultType
LEPlic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;

    orig_reference_fluid_volume = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, orig_reference_fluid_volume, IFT_LEPLIC_refVol, "refvol");
    return IRRT_OK;
}

double
LEPlic :: computeCriticalTimeStep(TimeStep *tStep)
{
    int ie, nelem = domain->giveNumberOfElements();
    double dt = 1.e6;
    LEPlicElementInterface *interface;

    for ( ie = 1; ie <= nelem; ie++ ) {
        if ( ( interface = ( LEPlicElementInterface * ) ( domain->giveElement(ie)->giveInterface(LEPlicElementInterfaceType) ) ) ) {
            dt = min( dt, interface->computeCriticalLEPlicTimeStep(tStep) );
        }
    }

    return 0.9 * dt;
}

void
LEPlic :: giveMaterialMixtureAt(FloatArray &answer, FloatArray &position)
{
    answer.resize(2);
    Element *elem = domain->giveSpatialLocalizer()->giveElementContainingPoint(position);
    LEPlicElementInterface *interface = ( LEPlicElementInterface * ) elem->giveInterface(LEPlicElementInterfaceType);
    if ( interface ) {
        Polygon pg;
        FloatArray n;
        interface->giveTempInterfaceNormal(n);
        interface->formVolumeInterfacePoly(pg, this, n, interface->giveTempLineConstant(), false);
        if ( pg.testPoint( position.at(1), position.at(2) ) ) {
            answer.at(1) = 1.0;
            answer.at(2) = 0.0;
        } else {
            answer.at(1) = 0.0;
            answer.at(2) = 1.0;
        }
    } else {
        answer.at(1) = 1.0;
        answer.at(2) = 0.0;
    }
}

void
LEPlic :: giveElementMaterialMixture(FloatArray &answer, int ie)
{
    answer.resize(2);
    Element *elem = domain->giveElement(ie);
    LEPlicElementInterface *interface = ( LEPlicElementInterface * ) elem->giveInterface(LEPlicElementInterfaceType);
    if ( interface ) {
        answer.at(1) = interface->giveTempVolumeFraction();
        answer.at(2) = 1. - answer.at(1);
    } else {
        answer.at(1) = 1.0;
        answer.at(2) = 0.0;
    }
}

double
LEPlic :: giveNodalScalarRepresentation(int inode)
{
    bool vof_1 = false, vof_0 = false;
    double vof, vofsum = 0.0;
    const IntArray *shelem = domain->giveConnectivityTable()->giveDofManConnectivityArray(inode);

    for ( int i = 1; i <= shelem->giveSize(); i++ ) {
        LEPlicElementInterface *interface = ( LEPlicElementInterface * ) domain->giveElement( shelem->at(i) )->giveInterface(LEPlicElementInterfaceType);
        if ( interface ) {
            vof = interface->giveTempVolumeFraction();
            if ( vof == 0.0 ) {
                vof_0 = true;
            } else if ( vof == 1.0 ) {
                vof_1 = true;
            }

            vofsum += vof;
        }
    }

    if ( vof_0 && vof_1 ) {
        return 0.5;
    } else if ( vof_0 ) {
        return 0.0;
    } else if ( vof_1 ) {
        return 1.0;
    } else {
        return vofsum / shelem->giveSize();
    }
}
} // end namespace oofem
