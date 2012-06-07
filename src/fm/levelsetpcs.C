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

#include "levelsetpcs.h"
#include "mathfem.h"
#include "timestep.h"
#include "node.h"
#include "element.h"
#include "datastream.h"
#include "conTable.h"
#include "spatiallocalizer.h"
#include "geotoolbox.h"
#include "fastmarchingmethod.h"
#include "error.h"
#include "contextioerr.h"

namespace oofem {
void
LevelSetPCS :: initialize()
{
    if ( 0 ) {
        if ( previousLevelSetValues.giveSize() != domain->giveNumberOfDofManagers() ) {
            OOFEM_ERROR("LevelSetPCS::initialize size of levelSetValues does not match number of dof managers");
        }
    } else {
        if ( initialRefMatFlag ) {
            int nnodes = domain->giveNumberOfDofManagers();
            previousLevelSetValues.resize(nnodes);
            for ( int i = 1; i <= nnodes; i++ ) {
                previousLevelSetValues.at(i) = ( -1.0 ) * initialRefMatVol.pointDistance( domain->giveNode(i)->giveCoordinate(ci1),
                                                                                         domain->giveNode(i)->giveCoordinate(ci2) );
            }
        }

        levelSetValues = previousLevelSetValues;
        //previousLevelSetValues.printYourself();
    }

    levelSetVersion++;
}


IRResultType
LevelSetPCS :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;

    IR_GIVE_OPTIONAL_FIELD(ir, previousLevelSetValues, IFT_LSPCS_levelSetValues, "levelset");
    if ( !previousLevelSetValues.giveSize() ) {
        FloatArray refmatpoly_x, refmatpoly_y;
        IR_GIVE_OPTIONAL_FIELD(ir, refmatpoly_x, IFT_LSPCS_refmatpoly_x, "refmatpolyx");
        IR_GIVE_OPTIONAL_FIELD(ir, refmatpoly_y, IFT_LSPCS_refmatpoly_y, "refmatpolyy");
        ci1 = 1, ci2 = 2;
        IR_GIVE_OPTIONAL_FIELD(ir, ci1, IFT_LSPCS_ci1, "ci1");
        IR_GIVE_OPTIONAL_FIELD(ir, ci2, IFT_LSPCS_ci2, "ci2");

        int nvert = refmatpoly_x.giveSize();
        if ( nvert ) {
            Vertex v;
            for ( int i = 1; i <= nvert; i++ ) {
                // create polygonal representation
                v.setCoords( refmatpoly_x.at(i), refmatpoly_y.at(i) );
                initialRefMatVol.addVertex(v);
            }

            // close polygon (add first vertex at the end
            v.setCoords( refmatpoly_x.at(1), refmatpoly_y.at(1) );
            initialRefMatVol.addVertex(v);
        }

        initialRefMatFlag = true;
    }

    reinit_alg = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, reinit_alg, IFT_LSPCS_reinit_alg, "lsra");

    reinit_dt = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, reinit_dt, IFT_LSPCS_reinit_dt, "rdt");
    if ( reinit_dt > 0. ) {
        reinit_dt_flag = true;
    }

    reinit_err = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, reinit_err, IFT_LSPCS_reinit_err, "rerr");

    nsd = 2;
    IR_GIVE_OPTIONAL_FIELD(ir, nsd, IFT_LSPCS_nsd, "nsd");

    return IRRT_OK;
}


void
LevelSetPCS :: updatePosition(TimeStep *atTime)
{
    int i, j, l, inodes, inode;
    int ndofman = domain->giveNumberOfDofManagers();
    bool twostage = true;

    double help, dt, volume, gfi_norm;

    FloatArray fs(ndofman), w(ndofman);
    FloatMatrix dN;
    FloatArray fi(4), gfi(nsd), n(nsd), k(4), dfii(4), alpha(4), un;
    LevelSetPCSElementInterface *interface;
    Element *ielem;
    ConnectivityTable *contable = domain->giveConnectivityTable();

    // needed for multistep update
    FloatArray ls_n;
    int __step = 0, __nstep = 10;

    levelSetValues = previousLevelSetValues;
    dt = atTime->giveTimeIncrement() / __nstep;

    do {
        ls_n = levelSetValues;

        this->pcs_stage1(levelSetValues, fs, w, atTime, PCS_levelSetUpdate);

        // update level set values
        for ( inode = 1; inode <= ndofman; inode++ ) {
            if ( fabs( w.at(inode) ) > 0.0 ) {
                // single stage integration
                levelSetValues.at(inode) = ls_n.at(inode) - dt *fs.at(inode) / w.at(inode);
            } else {
                // -------------------------
                // inflow into boundary node
                // -------------------------
                const IntArray *elems;
                IntArray mask(nsd);
                double v;
                // get velocity in inode
                if ( nsd == 2 ) {
                    mask.at(1) = V_u;
                    mask.at(2) = V_v;
                } else if ( nsd == 3 ) {
                    mask.at(1) = V_u;
                    mask.at(2) = V_v;
                    mask.at(3) = V_w;
                }

                domain->giveDofManager(inode)->giveUnknownVector( un, mask, EID_MomentumBalance, VM_Total, atTime->givePreviousStep() );
                elems = contable->giveDofManConnectivityArray(inode);
                // loop over shared elements
                volume = 0.0;
                help = 0.0;
                for ( l = 1; l <= elems->giveSize(); l++ ) {
                    // get element level set gradient
                    ielem = domain->giveElement( elems->at(l) );
                    inodes = ielem->giveNumberOfNodes();
                    interface = ( LevelSetPCSElementInterface * )
                                ielem->giveInterface(LevelSetPCSElementInterfaceType);

                    if ( interface ) {
                        interface->LS_PCS_computedN(dN);
                        // assemble element vector with  level set values
                        for ( i = 1; i <= inodes; i++ ) {
                            fi.at(i) = ls_n.at( ielem->giveDofManagerNumber(i) );
                        }

                        // compute gradient of level set
                        for ( j = 1; j <= nsd; j++ ) {
                            gfi.at(j) = 0.0;
                            for ( i = 1; i <= inodes; i++ ) {
                                gfi.at(j) += dN.at(i, j) * fi.at(i);
                            }
                        }

                        volume += ( v = interface->LS_PCS_computeVolume() );
                        gfi_norm = gfi.computeNorm();
                        if ( gfi_norm > 1.e-6 ) {
                            help += un.dotProduct(gfi) * v / gfi_norm;
                        }
                    }
                } // end loop over shared nodes

                levelSetValues.at(inode) = ls_n.at(inode) - dt * help / volume;
            }
        } // end loop over nodes

        if ( twostage ) {
            //this->pcs_stage1(levelSetValues, fs, w, atTime, PCS_levelSetUpdate); // ?

            for ( inode = 1; inode <= ndofman; inode++ ) {
                if ( w.at(inode) > 0.0 ) {
                    //two stage integration
                    // update
                    levelSetValues.at(inode) = 0.5 * ( ls_n.at(inode) + levelSetValues.at(inode) ) -
                                               0.5 *dt *fs.at(inode) / w.at(inode);
                }
            }
        }

        printf(".");
    } while ( ++__step <  __nstep );

    printf("\n");
    // print level set values to stdout (debug only)
    /*
     * printf ("Node: Level Set Value\n");
     * for (inode=1; inode<=ndofman; inode++) {
     * printf ("%5d %le\n",inode, levelSetValues.at(inode));
     * }
     */
    // redistance
    this->reinitialization(atTime);

    levelSetVersion++;
}


double
LevelSetPCS :: computeCriticalTimeStep(TimeStep *tStep)
{
    return 1.e6;
}


void
LevelSetPCS :: giveMaterialMixtureAt(FloatArray &answer, FloatArray &position)
{
    double ls;
    int i;
    FloatArray N(3);
    answer.resize(2);

    Element *elem = domain->giveSpatialLocalizer()->giveElementContainingPoint(position);
    LevelSetPCSElementInterface *interface = ( LevelSetPCSElementInterface * ) elem->giveInterface(LevelSetPCSElementInterfaceType);
    if ( interface ) {
        if ( elem->computeLocalCoordinates(N, position) ) {
            int inodes = elem->giveNumberOfNodes();
            for ( ls = 0.0, i = 1; i <= inodes; i++ ) {
                ls += N.at(i) * levelSetValues.at( elem->giveDofManagerNumber(i) );
            }

            if ( ls > 0.0 ) {
                answer.at(1) = 1.0;
                answer.at(2) = 0.0;
            } else {
                answer.at(1) = 0.0;
                answer.at(2) = 1.0;
            }
        } else {
            OOFEM_ERROR("LevelSetPCS::giveMaterialMixtureAt: computeLocalCoordinates failed");
        }
    } else {
        answer.at(1) = 1.0;
        answer.at(2) = 0.0;
    }
}

void
LevelSetPCS :: giveElementMaterialMixture(FloatArray &answer, int ie)
{
#ifdef LevelSetPCS_CACHE_ELEMENT_VOF
    if ( elemVofLevelSetVersion == levelSetVersion ) {
        answer = elemVof [ ie - 1 ];
        return;
    } else {
        Element *ielem;
        int i, _ie, inodes;
        LevelSetPCSElementInterface *interface;
        FloatArray fi;

        elemVof.resize( domain->giveNumberOfElements() );
        for ( _ie = 1; _ie <= domain->giveNumberOfElements(); _ie++ ) {
            ielem = domain->giveElement(_ie);
            inodes = ielem->giveNumberOfNodes();
            fi.resize(inodes);
            interface = ( LevelSetPCSElementInterface * ) ielem->giveInterface(LevelSetPCSElementInterfaceType);
            for ( i = 1; i <= inodes; i++ ) {
                fi.at(i) = levelSetValues.at( ielem->giveDofManagerNumber(i) );
            }

            interface->LS_PCS_computeVOFFractions(elemVof [ _ie - 1 ], fi);
        }

        elemVofLevelSetVersion = levelSetVersion;
        answer = elemVof [ ie - 1 ];
        return;
    }

#else

    Element *ielem = domain->giveElement(ie);
    int i, inodes = ielem->giveNumberOfNodes();
    LevelSetPCSElementInterface *interface = ( LevelSetPCSElementInterface * ) ielem->giveInterface(LevelSetPCSElementInterfaceType);
    FloatArray fi(inodes);

    for ( i = 1; i <= inodes; i++ ) {
        fi.at(i) = levelSetValues.at( ielem->giveDofManagerNumber(i) );
    }

    interface->LS_PCS_computeVOFFractions(answer, fi);
#endif
}

void
LevelSetPCS :: reinitialization(TimeStep *atTime)
{
     if ( reinit_alg == 0 ) {
        return;
     } else if ( reinit_alg == 1 ) {
        this->redistance(atTime);
     } else if ( reinit_alg == 2 ) {
        FloatArray ls1;
        this->FMMReinitialization(ls1);
        levelSetValues = ls1;
     } else {
        OOFEM_ERROR2("LevelSetPCS::reinitialization: unknown reinitialization scheme (%d)", reinit_alg);
     }
}


void
LevelSetPCS :: redistance(TimeStep *atTime)
{
    int nite = 0, inode, i, inodes;
    int ndofman = domain->giveNumberOfDofManagers();
    int nelem = domain->giveNumberOfElements();
    bool twostage = false;
    double dt, c, cm;

    FloatArray fs(ndofman), w(ndofman), d_old, d;
    FloatMatrix dN;
    FloatArray fi(4), gfi(nsd), n(nsd);
    //ConnectivityTable* contable = domain->giveConnectivityTable();
    //LevelSetPCSElementInterface* interface;
    Element *ielem;


    //return;
    if ( this->reinit_dt_flag ) {
        dt = this->reinit_dt;
    } else {
        dt = atTime->giveTimeIncrement();
    }

    IntArray _boundary(ndofman);
    int ie, pos, neg, _node;
    // check for boundary node
    for ( ie = 1; ie <= nelem; ie++ ) {
        ielem = domain->giveElement(ie);
        inodes = ielem->giveNumberOfNodes();
        pos = 0;
        neg = 0;
        for ( i = 1; i <= inodes; i++ ) {
            _node = ielem->giveDofManagerNumber(i);
            if ( levelSetValues.at(_node) >= 0 ) {
                pos++;
            } else {
                neg++;
            }
        }

        if ( pos && neg ) {
            for ( i = 1; i <= inodes; i++ ) {
                _node = ielem->giveDofManagerNumber(i);
                _boundary.at(_node) = 1;
            }
        }
    }


    d = levelSetValues;
    do {
        d_old = d;
        //levelSetValues = d;     // updated sign funtion
        pcs_stage1(levelSetValues, fs, w, atTime, PCS_levelSetRedistance);

        // update level set values
        // single stage integration
        cm = 0.0;
        for ( inode = 1; inode <= ndofman; inode++ ) {
            if ( _boundary.at(inode) ) {
                continue;
            }

            if ( fabs( w.at(inode) ) > 0.0 ) {
                c = dt * fs.at(inode) / w.at(inode);
                cm = max( cm, fabs(c / levelSetValues.at(inode)) );
                levelSetValues.at(inode) = levelSetValues.at(inode) - c;
            } else {
                //printf ("(%d) ", inode);
            }
        }

        if ( twostage ) {
            // this->pcs_stage1(d, fs, w, atTime, PCS_levelSetRedistance); //?
            cm = 0.0;
            for ( inode = 1; inode <= ndofman; inode++ ) {
                if ( _boundary.at(inode) ) {
                    continue;
                }

                if ( fabs( w.at(inode) ) > 0.0 ) {
                    //two stage integration
                    // update
                    d.at(inode) = 0.5 * ( d_old.at(inode) + d.at(inode) ) -
                        0.5 *dt *fs.at(inode) / w.at(inode);
                    cm = max( cm, fabs( ( d.at(inode) - d_old.at(inode) ) / d_old.at(inode) ) );
                }
            }
        }
    } while ( ( cm > this->reinit_err ) && ( ++nite < 2000 ) );

    //} while ((++nite < 200));
    printf("LS reinit: error %le in %d iterations\n", cm, nite);

    //levelSetValues = d;
}


void
LevelSetPCS :: pcs_stage1(FloatArray &ls, FloatArray &fs, FloatArray &w, TimeStep *atTime, PCSEqType t)
{
    int i, j, l, _ig, inodes;
    int ie, ndofman = domain->giveNumberOfDofManagers(), nelem   = domain->giveNumberOfElements();
    double alpha, dfi, help, sumkn, F, f, volume, gfi_norm;
    FloatMatrix dN;
    FloatArray gfi(nsd), fi(4), n(nsd), k(4), dfii(4);
    Element *ielem;
    LevelSetPCSElementInterface *interface;

    fs.resize(ndofman);
    w.resize(ndofman);
    fs.zero();
    w.zero();

    // loop over elements
    for ( ie = 1; ie <= nelem; ie++ ) {
        ielem = domain->giveElement(ie);
        inodes = ielem->giveNumberOfNodes();
        interface = ( LevelSetPCSElementInterface * )
                    ielem->giveInterface(LevelSetPCSElementInterfaceType);

        if ( interface ) {
            F = this->evalElemFContribution(t, ie, atTime);
            f = this->evalElemfContribution(t, ie, atTime);

            interface->LS_PCS_computedN(dN);
            volume = interface->LS_PCS_computeVolume();

            // assemble element vector with  level set values
            for ( i = 1; i <= inodes; i++ ) {
                fi.at(i) = ls.at( ielem->giveDofManagerNumber(i) );
            }

            // compute gradient of level set
            for ( j = 1; j <= nsd; j++ ) {
                gfi.at(j) = 0.0;
                for ( i = 1; i <= inodes; i++ ) {
                    gfi.at(j) += dN.at(i, j) * fi.at(i);
                }
            }

            // eval size of gfi
            gfi_norm = gfi.computeNorm();
            // compute ki
            for ( i = 1; i <= inodes; i++ ) {
                if ( gfi_norm > 1.e-12 ) {
                    // evaluate i-th normal (corresponding to side opposite to i-th vertex)
                    for ( j = 1; j <= nsd; j++ ) {
                        n.at(j) = nsd * dN.at(i, j) * volume; //?
                    }

                    k.at(i) = F * gfi.dotProduct(n) / ( nsd * gfi_norm );
                } else {
                    printf("zero gfi_norm for %d node\n", i);
                    k.at(i) = 0.0;
                }
            }
            dfi = fi.dotProduct(k);
            for ( i = 1; i <= inodes; i++ ) {
                help = 0.0;
                sumkn = 0.0;
                for ( l = 1; l <= inodes; l++ ) {
                    help += negbra( k.at(l) ) * ( fi.at(i) - fi.at(l) );
                    sumkn += negbra( k.at(l) );
                }

                if ( fabs(sumkn) > 1.e-12 ) {
                    dfii.at(i) = macbra( k.at(i) ) * help / sumkn;
                } else {
                    printf("zero sumkn for %d node\n", i);
                    dfii.at(i) = 0.0;
                }
            }

            //compute alpha_i
            for ( help = 0.0, l = 1; l <= inodes; l++ ) {
                help += max(0.0, dfii.at(l) / dfi);
            }

            for ( i = 1; i <= inodes; i++ ) {
                _ig = ielem->giveDofManagerNumber(i);
                if ( fabs(help) > 0.0 ) {
                    alpha = max(0.0, dfii.at(i) / dfi) / help;
                    fs.at(_ig) += alpha * ( dfi - f * volume );
                    w.at(_ig) += alpha * volume;
                }
            }
        } else {
            OOFEM_ERROR2("LevelSetPCS::updatePosition: element %d does not implement LevelSetPCSElementInterfaceType", ie);
        }
    } // end loop over elements

}


double
LevelSetPCS :: evalElemFContribution(PCSEqType t, int ie, TimeStep *atTime)
{
    LevelSetPCSElementInterface *interface = ( LevelSetPCSElementInterface * )
                                             domain->giveElement(ie)->giveInterface(LevelSetPCSElementInterfaceType);
    if ( t == PCS_levelSetUpdate ) {
        return interface->LS_PCS_computeF(this, atTime);
    } else if ( t == PCS_levelSetRedistance ) {
        return interface->LS_PCS_computeS(this, atTime);
    }

    return 0.0;
}


double
LevelSetPCS :: evalElemfContribution(PCSEqType t, int ie, TimeStep *atTime)
{
    LevelSetPCSElementInterface *interface = ( LevelSetPCSElementInterface * )
                                             domain->giveElement(ie)->giveInterface(LevelSetPCSElementInterfaceType);
    if ( t == PCS_levelSetUpdate ) {
        return 0.0;
    } else if ( t == PCS_levelSetRedistance ) {
        return interface->LS_PCS_computeS(this, atTime);
    }

    return 0.0;
}


void
LevelSetPCS :: FMMReinitialization(FloatArray &dmanValues)
{
    // tag points with boundary value as known
    // then tag as trial all points that are one grid point away
    // finally tag as far all other grid points
    int i, j, jnode, enodes, __pos, __neg, nelem = domain->giveNumberOfElements();
    double _lsval;
    Element *ie;
    std :: list< int >bcDofMans;
    std :: list< int > :: iterator it;

    dmanValues.resize( domain->giveNumberOfDofManagers() );
    // here we loop over elements and identify those, that have zero level set
    // then nodes belonging to these elements are boundary ones (with known distance)
    for ( i = 1; i <= nelem; i++ ) {
        ie = domain->giveElement(i);
        enodes = ie->giveNumberOfDofManagers();
        __pos = 0;
        __neg = 0;       // count positive and negative level set values in element nodes
        for ( j = 1; j <= enodes; j++ ) {
            _lsval = this->giveLevelSetDofManValue( ie->giveDofManagerNumber(j) );
            if ( _lsval > 0.0 ) {
                __pos++;
            }

            if ( _lsval < 0.0 ) {
                __neg++;
            }
        }

        if ( ( __pos && __neg ) || ( __pos + __neg < enodes ) ) {
            // zero level set within element
            // we have to tag element nodes as known and compute their boundary value
            for ( j = 1; j <= enodes; j++ ) {
                // simplified (here we use original level set values)
                jnode = ie->giveDofManagerNumber(j);
                if ( ( dmanValues.at(jnode) = this->giveLevelSetDofManValue(jnode) ) >= 0. ) {
                    bcDofMans.push_front(jnode);
                } else {
                    bcDofMans.push_front(-jnode);
                }
            }
        }
    }

    FastMarchingMethod fmm(domain);
    // fast marching for positive level set values
    fmm.solve(dmanValues, bcDofMans, 1.0);
    // revert bcDofMans signs
    for ( it = bcDofMans.begin(); it != bcDofMans.end(); ++it ) {
        * it = -* it;
    }

    // fast marching for negative level set values
    fmm.solve(dmanValues, bcDofMans, -1.0);
}


contextIOResultType
LevelSetPCS :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( !stream->write(& levelSetVersion, 1) ) {
      THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = levelSetValues.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
LevelSetPCS :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( !stream->read(&levelSetVersion, 1) ) {
       THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = levelSetValues.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    previousLevelSetValues=levelSetValues;
#ifdef LevelSetPCS_CACHE_ELEMENT_VOF
    elemVofLevelSetVersion = 0;
#endif

    return CIO_OK;
}
} // end namespace oofem
