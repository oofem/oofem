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

#include "skyline.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#include "element.h"
#include "mathfem.h"
#include "verbose.h"
#include "sparsemtrxtype.h"
#include "classfactory.h"
#include "activebc.h"
#include "contact/contactmanager.h"
#include "contact/contactdefinition.h"
#include "contact/contactelement.h"
#include "unknownnumberingscheme.h"


#include <climits>
#include <cstdlib>
#include <utility>

#ifdef TIME_REPORT
 #include "timer.h"
#endif

namespace oofem {
REGISTER_SparseMtrx(Skyline, SMT_Skyline);

Skyline :: Skyline(int n) : SparseMtrx(n, n),
    isFactorized(false)
{
}


Skyline :: Skyline(const Skyline &s) : SparseMtrx(s.giveNumberOfRows(), s.giveNumberOfColumns()),
    mtrx(s.mtrx),
    adr(s.adr),
    isFactorized(s.isFactorized)
{}


Skyline :: Skyline(int n, FloatArray mtrx1, IntArray adr1) : SparseMtrx(n, n),
    mtrx(std::move(mtrx1)),
    adr(std::move(adr1)),
    isFactorized(false)
{}


std::unique_ptr<SparseMtrx> Skyline :: clone() const
{
    return std::make_unique<Skyline>(*this);
}


double &
Skyline :: at(int i, int j)
{
#ifdef DEBUG
    // check size
    if ( ( i > this->giveNumberOfRows() ) || ( j > this->giveNumberOfRows() ) ) {
        OOFEM_ERROR("dimension mismatch - accessing value at (%d,%d)", i, j);
    }

#endif
    // only upper triangular part of skyline is stored
    if ( j < i ) {
        std::swap(i, j);
    }

    int d1 = this->adr.at(j);
    int ind = d1 + ( j - i );

    if ( ( adr.at(j + 1) - adr.at(j) ) <= ( j - i ) ) {
        OOFEM_ERROR("request for element which is not in sparse mtrx (%d,%d)", i, j);
        //
        // NOTE:
        //
        // don't return reference to some zero value; it is true, but possible change
        // of its value will require rebuilding internal storage structure
        // of sparse matrix
        //
    }

    this->version++;
    return mtrx [ ind ];
}

double
Skyline :: at(int i, int j) const
{
#ifdef DEBUG
    // check size
    if ( ( i > this->giveNumberOfRows() ) || ( j > this->giveNumberOfRows() ) ) {
        OOFEM_ERROR("dimension mismatch, when accessing value at (%d,%d)", i, j);
    }

#endif
    // only upper triangular part of skyline is stored
    if ( j < i ) {
        std::swap(i, j);
    }

    int d1 = this->adr.at(j);
    int ind = d1 + ( j - i );

    if ( ( adr.at(j + 1) - adr.at(j) ) <= ( j - i ) ) {
        return 0.;
    }

    return mtrx [ ind ];
}


bool
Skyline :: isAllocatedAt(int i, int j) const
{
    if ( j < i ) {
        std::swap(i, j);
    }

    if ( ( adr.at(j + 1) - adr.at(j) ) <= ( j - i ) ) {
        return false;
    }

    return true;
}


void
Skyline :: toFloatMatrix(FloatMatrix &answer) const
{
    int size = this->giveNumberOfColumns();

#  ifdef DEBUG
    if ( size != this->adr.giveSize() - 1 ) {
        OOFEM_ERROR("Internal error in skyline matrix: num columns != size(adr)-1: %d != %d", size, this->adr.giveSize() - 1);
    }
#  endif

    answer.resize(size, size);
    answer.zero();

    for ( int j = 1; j <= size; j++ ) {
        int d1 = adr.at(j);
        int d2 = adr.at(j + 1);
        int pk = j;
        for ( int i = d1; i < d2; i++ ) {
            answer.at(pk, j) = mtrx [ i ];
            pk--;
        }
    }
    answer.symmetrized();
}


int Skyline :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
#  ifdef DEBUG
    int dim = mat.giveNumberOfRows();
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("dimension of 'mat' and 'loc' mismatch");
    }

#  endif

    int ndofe = mat.giveNumberOfRows();

    for ( int i = 1; i <= ndofe; i++ ) {
        int ac1 = loc.at(i);
        if ( ac1 == 0 ) {
            continue;
        }

        for ( int j = 1; j <= ndofe; j++ ) {
            int ac2 = loc.at(j);
            if ( ac2 == 0 ) {
                continue;
            }

            if ( ac1 > ac2 ) {
                continue;
            }

            mtrx [ adr.at(ac2) + ac2 - ac1 ] += mat.at(i, j);
        }
    }

    this->version++;
    return 1;
}


int Skyline :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
    int dim1 = mat.giveNumberOfRows();
    int dim2 = mat.giveNumberOfColumns();
    for ( int i = 1; i <= dim1; i++ ) {
        int ii = rloc.at(i);
        if ( ii ) {
            for ( int j = 1; j <= dim2; j++ ) {
                int jj = cloc.at(j);
                if ( jj && ii <= jj ) {
                    this->at(ii, jj) += mat.at(i, j);
                }
            }
        }
    }

    this->version++;

    return 1;
}


FloatArray *Skyline :: backSubstitutionWith(FloatArray &y) const
{
    // allocation of answer
    FloatArray solution( y.giveSize() );
    int n = this->giveNumberOfRows();

    /************************************/
    /*  modification of right hand side */
    /************************************/
    for ( int k = 2; k <= n; k++ ) {
        int ack = adr.at(k);
        int ack1 = adr.at(k + 1);
        double s = 0.0;
        int acs = k - ( ack1 - ack ) + 1;
        for ( int i = ack1 - 1; i > ack; i-- ) {
            s += mtrx [ i ] * y.at(acs);
            acs++;
        }

        y.at(k) -= s;
    }

    /*****************/
    /*  zpetny chod  */
    /*****************/
    for ( int k = 1; k <= n; k++ ) {
        y.at(k) /= mtrx [ adr.at(k) ];
    }

    for ( int k = n; k > 0; k-- ) {
        int ack = adr.at(k);
        int ack1 = adr.at(k + 1);
        solution.at(k) = y.at(k);
        int acs = k - ( ack1 - ack ) + 1;
        for ( int i = ack1 - 1; i > ack; i-- ) {
            y.at(acs) -= mtrx [ i ] * solution.at(k);
            acs++;
        }
    }

    y = solution;
    return & y;
}

int Skyline :: setInternalStructure(IntArray a)
{
    adr = std::move(a);
    int n = adr.giveSize();
    int nwk = adr.at(n);
    this->mtrx.resize(nwk);

    nRows = nColumns = n - 1;

    this->version++;
    return true;
}


int Skyline :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s)
{
    // first create array of
    // maximal column height for assembled characteristics matrix
    //

    int maxle;
    int ac1;
    int neq;
    if ( s.isDefault() ) {
        neq = eModel->giveNumberOfDomainEquations(di, s);
    } else {
        neq = s.giveRequiredNumberOfDomainEquation();
    }
    if ( neq == 0 ) {
        mtrx.clear();
        adr.clear();
        return true;
    }

    IntArray loc;
    IntArray mht(neq);
    Domain *domain = eModel->giveDomain(di);

    for ( int j = 1; j <= neq; j++ ) {
        mht.at(j) = j; // initialize column height, maximum is line number (since it only stores upper triangular)
    }

    // loop over elements code numbers
    for ( auto &elem : domain->giveElements() ) {
        elem->giveLocationArray(loc, s);
        maxle = INT_MAX;
        for ( int ieq : loc ) {
            if ( ieq != 0 ) {
                maxle = min(maxle, ieq);
            }
        }

        for ( int ieq : loc ) {
            if ( ieq != 0 ) {
                mht.at(ieq) = min( maxle, mht.at(ieq) );
            }
        }
    }

    // loop over active boundary conditions (e.g. relative kinematic constraints)
    std :: vector< IntArray >r_locs;
    std :: vector< IntArray >c_locs;

    for ( auto &gbc : domain->giveBcs() ) {
        ActiveBoundaryCondition *bc = dynamic_cast< ActiveBoundaryCondition * >( gbc.get() );
        if ( bc != NULL ) {
            bc->giveLocationArrays(r_locs, c_locs, UnknownCharType, s, s);
            for ( std :: size_t k = 0; k < r_locs.size(); k++ ) {
                IntArray &krloc = r_locs [ k ];
                IntArray &kcloc = c_locs [ k ];
                maxle = INT_MAX;
                for ( int ii : krloc ) {
                    if ( ii > 0 ) {
                        maxle = min(maxle, ii);
                    }
                }
                for ( int jj : kcloc ) {
                    if ( jj > 0 ) {
                        mht.at(jj) = min( maxle, mht.at(jj) );
                    }
                }
            }
        }
    }


    if ( domain->hasContactManager() ) {
        ContactManager *cMan = domain->giveContactManager();

        for ( int i = 1; i <= cMan->giveNumberOfContactDefinitions(); i++ ) {
            ContactDefinition *cDef = cMan->giveContactDefinition(i);
            for ( int k = 1; k <= cDef->giveNumbertOfContactElements(); k++ ) {
                ContactElement *cEl = cDef->giveContactElement(k);
                cEl->giveLocationArray(loc, s);

                maxle = INT_MAX;
                for ( int ieq : loc ) {
                    if ( ieq != 0 ) {
                        maxle = min(maxle, ieq);
                    }
                }

                for ( int ieq : loc ) {
                    if ( ieq != 0 ) {
                        mht.at(ieq) = min( maxle, mht.at(ieq) );
                    }
                }
            }
        }
    }

    // NOTE
    // add there call to eModel if any possible additional equation added by
    // eModel
    // currently not supported

    // increases number of columns according to size of mht
    // mht is array containing minimal equation number per column
    // This method also increases column height.

    adr.resize(neq + 1);

    ac1 = 1;
    for ( int i = 1; i <= neq; i++ ) {
        adr.at(i) = ac1;
        ac1 += ( i - mht.at(i) + 1 );
    }

    adr.at(neq + 1) = ac1;
    nRows = nColumns = neq;

    mtrx.resize( ac1 );

    this->version++;
    return true;
}


SparseMtrx *Skyline :: factorized()
{
    // Returns the receiver in  U(transp).D.U  Crout factorization form.

    if ( isFactorized ) {
        return this;
    }

#ifdef TIME_REPORT
    Timer timer;
    timer.startTimer();
#endif

    int n = this->giveNumberOfRows();

    OOFEM_LOG_DEBUG("Skyline info: neq is %d, nwk is %d\n", n, this->giveNumberOfNonZeros());

    for ( int k = 2; k <= n; k++ ) {
        /*  smycka pres sloupce matice  */
        int ack = adr.at(k);
        int ack1 = adr.at(k + 1);
        int acrk = k - ( ack1 - ack ) + 1;
        for ( int i = acrk + 1; i < k; i++ ) {
            /*  smycka pres prvky jednoho sloupce matice  */
            int aci = adr.at(i);
            int aci1 = adr.at(i + 1);
            int acri = i - ( aci1 - aci ) + 1;
            int ac;
            if ( acri < acrk ) {
                ac = acrk;
            } else {
                ac = acri;
            }

            int acj = k - ac + ack;
            int acj1 = k - i + ack;
            int acs = i - ac + aci;
            double s = 0.0;
            for ( int j = acj; j > acj1; j-- ) {
                s += mtrx [ j ] * mtrx [ acs ];
                acs--;
            }

            mtrx [ acj1 ] -= s;
        }

        /*  uprava diagonalniho prvku  */
        double s = 0.0;
        for ( int i = ack1 - 1; i > ack; i-- ) {
            double g = mtrx [ i ];
            int acs = adr.at(acrk);
            acrk++;
            mtrx [ i ] /= mtrx [ acs ];
            s += mtrx [ i ] * g;
        }

        mtrx [ ack ] -= s;
    }

    isFactorized = true;

#ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_DEBUG( "Skyline info: user time consumed by factorization: %.2fs\n", timer.getUtime() );
#endif

    //this->version++;
    return this;
}


void Skyline :: times(const FloatArray &x, FloatArray &answer) const
{
    // Computes y, the results  of the  y = U.x, where U is
    // the receiver. Returns the result.

    int n = x.giveSize();

    //
    // first check sizes
    //
    if ( this->giveNumberOfRows() != n ) {
        OOFEM_ERROR("size mismatch");
    }

    answer.resize(n);
    answer.zero();

    int acc = 1;
    for ( int i = 1; i <= n; i++ ) {
        int aci = adr.at(i);
        int aci1 = adr.at(i + 1);
        int ac = i - ( aci1 - aci ) + 1;
        double s = 0.0;
        int acb = ac;
        for ( int k = aci1 - 1; k >= aci; k-- ) {
            s += mtrx [ k ] * x.at(acb);
            acb++;
        }

        answer.at(acc) = s;
        acc++;

        for ( int j = ac; j < i; j++ ) {
            aci1--;
            answer.at(j) += mtrx [ aci1 ] * x.at(i);
            aci++;
        }
    }
}


void Skyline :: times(double x)
{
    // Multiplies receiver by scalar value.
    mtrx.times(x);

    this->version++;
}


void Skyline :: add(double x, SparseMtrx &m)
{
    Skyline *M = dynamic_cast< Skyline* >( &m );

    mtrx.add(x, M->mtrx);

    this->version++;
}


void Skyline :: printYourself() const
{
    // Prints the receiver on screen.
    FloatMatrix copy;
    this->toFloatMatrix(copy);
    copy.printYourself();
}


void Skyline :: writeToFile(const char *fname) const
{
    FILE *file = fopen(fname, "w");
    FloatMatrix copy;
    this->toFloatMatrix(copy);
    for ( int i = 1; i <= nRows; ++i ) {
        for ( int j = 1; j <= nColumns; ++j ) {
            fprintf( file, "%10.3e  ", copy.at(i, j) );
        }
        fprintf(file, "\n");
    }
    fclose(file);
}


void Skyline :: zero()
{
    mtrx.zero();
    isFactorized = false;

    this->version++;
}


std::unique_ptr<SparseMtrx> Skyline :: giveSubMatrix(const IntArray &rows, const IntArray &cols) 
{
    IntArray positions( cols.giveSize() + 1 );

    FloatArray values( this->giveNumberOfNonZeros() ); //TODO choose a better initial size? 
    int diagPos = 1;
    int nnz = 0; // number of nonzeros 

    for ( int j = 1; j <= cols.giveSize(); j++ ) {
        for ( int i = rows.giveSize(); i >= 1; i-- ) { // start from the "bottom of the matrix"
            if( j < i ){
                continue;
            }

            bool hasValue = this->isAllocatedAt( rows.at(i), cols.at(j) );

            if ( hasValue  && i == j ) { // allocated diagonal element
                values.at(++nnz) = this->at( rows.at(i), cols.at(j) );
                positions.at(diagPos++) = nnz;
            } else if ( i == j  ) { // empty diagonal element
                values.at(++nnz) = 0.0;
                positions.at(diagPos++) = nnz;
            } else if ( hasValue  ) {
                values.at(++nnz) = this->at( rows.at(i), cols.at(j) );
            }
        }
    }

    positions.at(diagPos++) = ++nnz;

    FloatArray mtrxValues(nnz);

    for ( int i = 0; i < nnz-1; i++ ) {
        mtrxValues [ i + 1] = values [ i ];
    }
    int neq = rows.giveSize();

    return std::make_unique<Skyline>(neq, mtrxValues, positions);
}


void Skyline :: rbmodes(FloatMatrix &r, int &nse, IntArray &se,
                        double limit, int tc)
{
    /*
     * funkce rozlozi matici A na LDL tvar, pak vypocte
     * bazove vektory prostoru Ker A
     */
    int neq = this->giveNumberOfRows();
    IntArray adrb(7);
    FloatArray b(6 *neq);

    /**********************/
    /*  rozklad matice A  */
    /**********************/

    // report skyline statistics
    OOFEM_LOG_INFO("Skyline info: neq is %d, nwk is %d\n", neq, this->giveNumberOfNonZeros());

    if ( tc == 1 || tc == 3 ) {
        /*  pocitadlo singularnich rovnic  */
        int ise = 1;
        /*  pocitadlo v poli singularnich radku  */
        int ib = 1;
        adrb.at(1) = 1;

        /*  cyklus pres radky, ktere maji byt odkondezovany  */
        for ( int i = 2; i <= neq; i++ ) {
            int lj = adr.at(i);
            int uj = adr.at(i + 1) - 2;

            /*  minimalni radkovy index v i tem sloupci  */
            int mi = i - ( uj - lj ) - 1;
            int j = mi + 1;

            /*  cyklus pres mimodiagonalni prvky zpracovavaneho radku  */
            for ( int jj = uj; jj > lj; jj-- ) {
                int li = adr.at(j);
                int ui = adr.at(j + 1) - 1;
                int k = j - ( ui - li );

                /*  vyber nizsiho sloupce a tim urceni rozsahu cyklu  */
                int uk, ii;
                if ( k < mi ) {
                    uk = uj + 1;
                    ii = li + j - mi;
                } else {
                    uk = lj + i - k;
                    ii = ui;
                }

                /*  cyklus pres prvky nad zpracovavanym prvkem  */
                double s = 0.0;
                for ( int kk = uk; kk > jj; kk-- ) {
                    s += mtrx [ kk ] * mtrx [ ii ];
                    ii--;
                }

                mtrx [ jj ] -= s;
                j++;
            }

            /*  uprava diagonalniho prvku  */
            double s = 0.0;
            j = mi;
            for ( int jj = uj + 1; jj > lj; jj-- ) {
                double g = mtrx [ jj ];
                mtrx [ jj ] /= mtrx [ adr.at(j) ];
                s += mtrx [ jj ] * g;
                j++;
            }

            mtrx [ lj ] -= s;

            /*  kontrola diagonalniho prvku  */
            if ( fabs(mtrx [ lj ]) < limit ) {
                /*  pole cisel singularnich rovnic  */
                se.at(ise) = i;
                ise++;

                /*  vynulovani prvku sloupce v poli a a jejich uchovani v poli b  */
                for ( int jj = uj + 1; jj > lj; jj-- ) {
                    b.at(ib) = mtrx [ jj ];
                    ib++;
                    mtrx [ jj ] = 0.0;
                }

                mtrx [ lj ] = 1.0;
                /*  pole adres zacatku v poli obsahujicim singularni radky  */
                adrb.at(ise) = ib;

                /*  vynulovani prvku radku v poli a  */
                for ( int j = i + 1; j <= neq; j++ ) {
                    if ( j - ( adr.at(j + 1) - adr.at(j) ) < i ) {
                        mtrx [ adr.at(j) + j - i ] = 0.0;
                    }
                }
            }
        }

        nse = ise - 1;
    }

    if ( tc == 2 || tc == 3 ) {
        /*  navrat puvodne vynulovanych slozek  */
        int ise = nse;
        for ( int i = 1; i <= ise; i++ ) {
            int uj = adr.at(se.at(i) + 1) - 1;
            int lj = adr.at( se.at(i) );
            int ib = adrb.at(i);
            for ( int jj = uj; jj > lj; jj-- ) {
                mtrx [ jj ] = b.at(ib);
                ib++;
            }
        }

        /*  sestaveni baze ker A  */
        if ( ise ) {
            r.resize(neq, ise);
            r.zero();
        } else {
            r.clear();
        }

        for ( int i = 1; i <= ise; i++ ) {
            for ( int j = neq; j > 0; j-- ) {
                r.at(se.at(i), i) = 1.0;
                int uk = adr.at(j + 1) - 1;
                int lk = adr.at(j);
                int k = j - ( uk - lk );
                for ( int kk = uk; kk > lk; kk-- ) {
                    r.at(k, i) -= mtrx [ kk ] * r.at(j, i);
                    k++;
                }
            }
        }
    }

    //this->version++;
}


void Skyline :: ldl_feti_sky(FloatArray &x, FloatArray &y,
                             int nse, double limit, IntArray &se)
/*
 * funkce resi cast soustavy rovnic v metode FETI
 * matice soustavy a je ulozena ve skyline
 * matice soustavy muze byt singularni
 * matice soustavy musi byt jiz rozlozena na tvar LDL
 *
 * vstupy
 * y - vektor prave strany
 * adr - pole adres diagonalnich prvku
 * n - pocet neznamych
 * m - pocet neznamych redukovaneho problemu
 * nse - pocet linearne zavislych rovnic
 * limit - promenna urcujici linearni zavislost
 * se - pole cisel singularnich rovnic
 *
 * vystupy
 * x - vektor reseni
 */
{
    int neq = this->giveNumberOfRows();

    /*******************************************************/
    /*  vypocet pomocneho vektoru z                        */
    /*  vektor prave strany ma na prvnich m pozicich nuly  */
    /*  pro vektor z se nealokuje nove pole                */
    /*  slozky vektoru y se prepisuji na slozky vektoru z  */
    /*******************************************************/
    //k = 0;
    for ( int i = 1; i <= neq; i++ ) {
        int lj = adr.at(i);
        int uj = adr.at(i + 1) - 1;
        int ii = i - uj + lj;
        double s = 0.0;
        for ( int j = uj; j > lj; j-- ) {
            s += mtrx [ j ] * y.at(ii);
            ii++;
        }

        /*
         *  if (se[k]==i){
         *    y[i]=0.0;  k++;
         *  }
         *  else{
         */

        y.at(i) -= s;

        /*  }*/
    }

    /********************************************************/
    /*  kontrola zeliminovaneho vektoru prave strany        */
    /*  pokud na nektere pozici se[i] bude nenulove cislo,  */
    /*  je soustava neresitelna, protoze se lisi hodnosti   */
    /*  matice soustavy a rozsirene matice soustavy         */
    /********************************************************/

    /*
     * fprintf (s2,"\n");
     * for (i=1;i<=nse;i++){
     *  fprintf (s2,"\n eliminovana prava strana singularni rovnice  %ld     %le",se->at(i),y->at(se->at(i)));
     */
    /*
     *  if (fabs(y[se[i]])>limit){
     *    fprintf (s2,"\n\n v %ld. linearne zavisle rovnici je nenulovy prvek",se[i]);
     *    fprintf (s2,"\n ve vektoru zeliminovane prave strany.");
     *  }
     */
    /*
     * }
     */

    /**********************************************************/
    /*  deleni prvku vektoru prave strany diagonalnimi prvky  */
    /**********************************************************/
    for ( int i = 1; i <= neq; i++ ) {
        y.at(i) /= mtrx [ adr.at(i) ];
    }

    /*****************************************************/
    /*  vypocet vektoru x z vektoru z                    */
    /*  vypocet se provadi pro m redukovanych neznamych  */
    /*****************************************************/
    int k = nse;
    for ( int i = neq; i > 0; i-- ) {
        int lj = adr.at(i);
        int uj = adr.at(i + 1);
        if ( k > 0 ) {
            if ( se.at(k) == i ) {
                y.at(i) = 1.0;
                k--;
            }
        }

        x.at(i) = y.at(i);
        double s = x.at(i);
        int ii = i - 1;
        for ( int j = lj + 1; j < uj; j++ ) {
            y.at(ii) -= s * mtrx [ j ];
            ii--;
        }
    }

    //this->version++;
}
} // end namespace oofem
