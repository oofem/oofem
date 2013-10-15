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

#include "SparseConectivityMtx.h"
#include "BiSection.h"

DSS_NAMESPASE_BEGIN

SparseConectivityMtxII :: SparseConectivityMtxII(const SparseMatrixF &sm, char block_size)
{
    n = ( long ) sm.neq / block_size + ( ( ( sm.neq % block_size ) == 0 ) ? 0 : 1 );
    nonzeros = 0;
    ColumnsIndexes = new IntArrayList * [ n ];
    for ( long bj = 0; bj < n; bj++ ) {
        ColumnsIndexes [ bj ] = new IntArrayList();
    }

    long *pmask = new long [ n ];
    memset( pmask, 0, n * sizeof( long ) );
    IntArrayList *this_columnJ = NULL;

    // Fill the upper half of the matrix
    for ( unsigned long j = 0; j < sm.neq; j++ ) {
        long bj = j / block_size;

        if ( j % block_size == 0 ) {
            this_columnJ = ColumnsIndexes [ bj ];
        }

        long lastbi = -1;
        for ( unsigned long ad = sm.Adr(j); ad < sm.Adr(j + 1); ad++ ) {
            long i = sm.Ci(ad);
            long bi = i / block_size;
            if ( bi == bj || bi == lastbi ) {
                continue;
            }

            lastbi = bi;

            if ( pmask [ bi ] == 0 ) {
                nonzeros++;
                this_columnJ->SortInsert(bi);
                pmask [ bi ] = 1;
            }
        }

        if ( ( long ) ( j % block_size ) == ( long ) ( block_size - 1 ) ) {
            this_columnJ->SetIndexesTo(pmask, 0);
        }
    }

    // Fill lower half of the matrix
    for ( unsigned long i = 0; i < sm.neq; i++ ) {
        long bi = i / block_size;

        long lastbj = -1;
        for ( unsigned long ad = sm.Adr(i); ad < sm.Adr(i + 1); ad++ ) {
            long j = sm.Ci(ad);
            long bj = j / block_size;
            if ( bi == bj || bj == lastbj ) {
                continue;
            }

            lastbj = bj;

            if ( ColumnsIndexes [ bj ]->AddIfNotLast(bi) >= 0 ) {
                nonzeros++;
            }
        }
    }

    if ( pmask ) {
        delete [] pmask;
    }

    pmask = NULL;
}

SparseConectivityMtxII :: SparseConectivityMtxII(const SparseMatrixF &sm, Ordering *node_order, char block_size)
{
    long *nprm = NULL;
    if ( node_order ) {
        nprm = node_order->perm->Items;
        n = node_order->order->Count / block_size;
    } else {
        n = ( long ) sm.neq / block_size + ( ( ( sm.neq % block_size ) == 0 ) ? 0 : 1 );
    }

    nonzeros = 0;
    ColumnsIndexes = new IntArrayList * [ n ];
    for ( long bj = 0; bj < n; bj++ ) {
        ColumnsIndexes [ bj ] = new IntArrayList();
    }

    long *pmask = new long [ n ];
    memset( pmask, 0, n * sizeof( long ) );
    IntArrayList *this_columnJ = NULL;


    long old_bj = -1;

    // Fill the upper half of the matrix
    for ( unsigned long j = 0; j < sm.neq; j++ ) {
        long nj = node_order ? node_order->perm->Items [ j ] : j;
        long bj = nj / block_size;

        if ( bj != old_bj ) {
            if ( this_columnJ ) {
                this_columnJ->SetIndexesTo(pmask, 0);
            }

            this_columnJ = ColumnsIndexes [ bj ];
            old_bj = bj;
        }

        long lastbi = -1;
        for ( unsigned long ad = sm.Adr(j); ad < sm.Adr(j + 1); ad++ ) {
            long i = sm.Ci(ad);

            long ni = nprm ? nprm [ i ] : i;

            long bi = ni / block_size;
            if ( bi == bj || bi == lastbi ) {
                continue;
            }

            lastbi = bi;

            if ( pmask [ bi ] == 0 ) {
                nonzeros++;
                this_columnJ->SortInsert(bi);
                pmask [ bi ] = 1;
            }
        }
    }

    // Fill lower half of the matrix
    for ( unsigned long i = 0; i < sm.neq; i++ ) {
        long ni = nprm ? nprm [ i ] : i;

        long bi = ni / block_size;

        long lastbj = -1;
        for ( unsigned long ad = sm.Adr(i); ad < sm.Adr(i + 1); ad++ ) {
            long j = sm.Ci(ad);
            long nj = nprm ? nprm [ j ] : j;

            long bj = nj / block_size;
            if ( bi == bj || bj == lastbj ) {
                continue;
            }

            lastbj = bj;

            // I have to think about this better
            //if (ColumnsIndexes[bj]->AddIfNotLast(bi)>=0)
            //nonzeros++;
            if ( ColumnsIndexes [ bj ]->SortInsertIfNot(bi) >= 0 ) {
                nonzeros++;
            }
        }
    }

    if ( pmask ) {
        delete [] pmask;
    }

    pmask = NULL;
}

SparseConectivityMtxII :: SparseConectivityMtxII(SparseConectivityMtxII &mtx, Ordering *order)
{
    n = mtx.N();
    ColumnsIndexes = new IntArrayList * [ n ];
    nonzeros = mtx.Nonzeros();
    long nj, idx, i;
    long *lengths = new long [ n ];
    memset( lengths, 0, n * sizeof( long ) );

    // compute lengths
    for ( nj = 0; nj < n; nj++ ) {
        long j = order->order->Items [ nj ];
        IntArrayList &columnJ = * mtx.ColumnsIndexes [ j ];
        //foreach (long i in columnJ)
        for ( idx = columnJ.Count - 1; idx >= 0; idx-- ) {
            long i = columnJ.Items [ idx ];
            long ni = order->perm->Items [ i ];
            if ( ni > nj ) {
                lengths [ ni ]++;
            }
        }
    }

    // allocate
    for ( i = 0; i < n; i++ ) {
        ColumnsIndexes [ i ] = new IntArrayList(lengths [ i ]);
        ColumnsIndexes [ i ]->Count = ColumnsIndexes [ i ]->Capacity;
        lengths [ i ] = 0;
    }

    // fill
    for ( nj = 0; nj < n; nj++ ) {
        long j = order->order->Items [ nj ];
        IntArrayList &columnJ = * mtx.ColumnsIndexes [ j ];
        //foreach (long i in columnJ)
        long cnt = columnJ.Count;
        for ( long idx = 0; idx < cnt; idx++ ) {
            long i = columnJ.Items [ idx ];
            long ni = order->perm->Items [ i ];
            if ( ni > nj ) {
                ColumnsIndexes [ ni ]->Items [ lengths [ ni ]++ ] = nj;
            }
        }
    }

    delete [] lengths;
}

//Copy constructor
SparseConectivityMtxII :: SparseConectivityMtxII(const SparseConectivityMtxII &mtx)
{
    n = mtx.N();
    ColumnsIndexes = new IntArrayList * [ n ];
    nonzeros = mtx.Nonzeros();

    // copy
    for ( long i = 0; i < n; i++ ) {
        ColumnsIndexes [ i ] = new IntArrayList(* mtx.ColumnsIndexes [ i ]);
    }
}

SparseConectivityMtxII :: ~SparseConectivityMtxII()
{
    for ( long i = 0; i < n; i++ ) {
        if ( ColumnsIndexes [ i ] ) {
            delete ColumnsIndexes [ i ];
            ColumnsIndexes [ i ] = NULL;
        }
    }

    delete [] ColumnsIndexes;
    ColumnsIndexes = NULL;
}


void SparseConectivityMtxII :: AddGrow(long i, long j)
{
    if ( i < 0 || j < 0 || i >= n ) {
        //throw new IndexOutOfRangeException();
        return;
    }

    IntArrayList *ColumnIndexes = ColumnsIndexes [ i ];

    long myIndex = ColumnIndexes->BinarySearch(j);
    if ( myIndex < 0 ) {
        ColumnIndexes->Insert(~myIndex, j);
        nonzeros++;
    }
}


long SparseConectivityMtxII :: ColumnLength(long i)
{
    return ColumnsIndexes [ i ]->Count;
}

void SparseConectivityMtxII :: GetCmpRows(long * &adr, long * &ci)
{
    int n = N();
    adr = new long [ n + 1 ];

    int i, length = 0;
    adr [ 0 ] = 0;
    for ( i = 0; i < n; i++ ) {
        IntArrayList *plist = GetIndexesAboveDiagonalInColumn(i);
        length += plist->Count;
        adr [ i + 1 ] = length;
    }

    ci = new long [ length ];
    for ( i = 0; i < n; i++ ) {
        IntArrayList *plist = GetIndexesAboveDiagonalInColumn(i);
        Array :: Copy(plist->Items, ci + adr [ i ], plist->Count);
    }
}


IntArrayList *SparseConectivityMtxII :: GetOrder_Cuthill_McKee2()
{
    // the order
    IntArrayList *order = new IntArrayList(n);
    order->Alloc();
    long *p_order = order->Items;

    // r is the root..
    long *p_node_level = new long [ n ];
    memset( p_node_level, 0, n * sizeof( long ) );

    long last_level_start = 0;
    long next_level_start = 0;

    long last_level_count = 0;
    long next_level_count = 0;

    long l = 1;
    long k = 0;

    do {
        if ( last_level_count == 0 ) {
            long r = -1, cl, min = n + 1;
            for ( long i = 0; i < n; i++ ) {
                if ( ( p_node_level [ i ] == 0 ) && ( cl = ColumnLength(i) ) < min ) {
                    min = cl;
                    r = i;
                }
            }

            last_level_start = k;

            if ( r == -1 ) {
                //System.Diagnostics.Debug.Assert(false,"Cannot find any free node.");
                break;
            }

            p_node_level [ r ] = l;
            p_order [ k++ ] = r;
            last_level_count = 1;
        }

        next_level_start = k;
        next_level_count = 0;
        // mark new level
        for ( long ilr = last_level_start; ilr < next_level_start; ilr++ ) {
            long lr = p_order [ ilr ];
            //foreach (long neig in ColumnsIndexes[lr])
            IntArrayList &al = * ColumnsIndexes [ lr ];
            for ( long idx = al.Count - 1; idx >= 0; idx-- ) {
                long neig = al.Items [ idx ];
                if ( p_node_level [ neig ] != 0 ) {
                    continue;
                }

                p_node_level [ neig ] = l + 1;
                p_order [ k++ ] = neig;
                next_level_count++;
            }
        }

        l++;
        last_level_count = next_level_count;
        last_level_start = next_level_start;
    } while ( k < n );

    delete [] p_node_level;
    return order;
}

Ordering *SparseConectivityMtxII :: GenerateMD(IntArrayList *fixed)
{
    MD_Qqraph amd(this);
    return new Ordering( amd.GenerateMD(false, fixed) );
}

Ordering *SparseConectivityMtxII :: GenerateAMD(IntArrayList *fixed)
{
    MD_Qqraph amd(this);
    return new Ordering( amd.GenerateMD(true, fixed) );
}

Ordering *SparseConectivityMtxII :: GenerateAMD_AA(IntArrayList *fixed)
{
    MD_Qqraph amd(this);
    amd.aggressive_absorbtion = true;
    return new Ordering( amd.GenerateMD(true, fixed) );
}

Ordering *SparseConectivityMtxII :: Get_Cuthill_McKee()
{
    return new Ordering( this->GetOrder_Cuthill_McKee2() );
}

Ordering *SparseConectivityMtxII :: Get_Reverse_Cuthill_McKee()
{
    IntArrayList *order = this->GetOrder_Cuthill_McKee2();
    Array :: Reverse(order->Items, order->Count);
    return new Ordering(order);
}

Ordering *SparseConectivityMtxII :: Get_Unity()
{
    IntArrayList *order = new IntArrayList(n);
    order->InitIdentity();
    return new Ordering(order);
}

Ordering *SparseConectivityMtxII :: Get_RecursiveBiSection()
{
    IntArrayList *order = new IntArrayList(n);
    order->InitIdentity();

    CBiSection(this).RecurBiSectOrder(order);

    //RecurBiSect(order->Items,order->Count);

    return new Ordering(order);
}

#ifdef _LINK_METIS_
extern "C" void METIS_EdgeND(int *, int *, int *, int *, int *, int *, int *);
extern "C" void METIS_NodeND(int *, int *, int *, int *, int *, int *, int *);
#endif

Ordering *SparseConectivityMtxII :: Get_MetisDiSection()
{
#ifdef _LINK_METIS_
    IntArrayList *order = new IntArrayList(n);
    IntArrayList *perm = new IntArrayList(n);

    int n = N();
    long *adr;
    long *ci;

    GetCmpRows(adr, ci);

    int numflag = 0;
    /*int optionsNodeND[]=
     * {
     *      1,
     *      3,
     *      1,
     *      2,
     *      0,
     *      0,
     *      0,
     *      1
     * };
     *
     * METIS_NodeND(&n,adr,ci,&numflag,optionsNodeND,(int*)order->Items,(int*)perm->Items);
     */

    int optionsEdgeND[] =
    {
        0,
        3,
        1,
        2,
        0,
        0,
        0,
        1
    };

    METIS_NodeND(& n, ( int * ) adr, ( int * ) ci, & numflag, optionsEdgeND, ( int * ) order->Items, ( int * ) perm->Items);

    delete [] adr;
    delete [] ci;

    return new Ordering(perm, order);

#else
    Writeln(" The program was not linked with METIS library, use _LINK_METIS_ directive!");
    return NULL;

#endif
}

#ifdef _LINK_COLAMD_
 #include "colamd\colamd.h"
#endif

Ordering *SparseConectivityMtxII :: Get_ColAMD()
{
#ifdef _LINK_COLAMD_
    IntArrayList *order = new IntArrayList(n + 1);
    order->Count = n;

    int n = N();
    long *adr;
    long *ci;

    GetCmpRows(adr, ci);

    int stats [ COLAMD_STATS ];

    symamd                              /* return (1) if OK, (0) otherwise */
    (
        n,                                      /* number of rows and columns of A */
        ( int * ) ci,                                   /* row indices of A */
        ( int * ) adr,                                  /* column pointers of A */
        ( int * ) order->Items,                         /* output permutation, size n_col+1 */
        NULL,           /* parameters (uses defaults if NULL) */
        stats,                  /* output statistics and error codes */
        calloc,
        free
    );

    //colamd_report(stats);

    delete [] adr;
    delete [] ci;

    return new Ordering(order);

#else
    Writeln(" The program was not linked with COLAMD sources, use _LINK_COLAMD_ directive!");
    return NULL;

#endif
}

void SparseConectivityMtxII :: GenerateFillInPresorted(Ordering *ord)
{
    MD_Qqraph amd(this);
    amd.keep_sorted_order = true;
    amd.GenerateMD(false, ord->order);
}

Ordering *SparseConectivityMtxII :: GetOrdering(Ordering :: Type ord)
{
    if ( ord == Ordering :: None ) {
        return Get_Unity();
    } else
    if ( ord == Ordering :: ReverseCuthillMcKee ) {
        return Get_Reverse_Cuthill_McKee();
    } else
    if ( ord == Ordering :: CuthillMcKee ) {
        return Get_Cuthill_McKee();
    } else
    if ( ord == Ordering :: MinimumDegree ) {
        return GenerateMD();
    } else
    if ( ord == Ordering :: ApproxMinimumDegree ) {
        return GenerateAMD();
    } else
    if ( ord == Ordering :: MetisND ) {
        return Get_MetisDiSection();
    } else
    if ( ord == Ordering :: ColAMD ) {
        return Get_ColAMD();
    } else {
        return NULL;
    }
}

IntArrayList *SparseConectivityMtxII :: GetIndexesAboveDiagonalInColumn(long j)
{
    return ColumnsIndexes [ j ];
}

IntArrayList *SparseConectivityMtxII :: DetachIndexesAboveDiagonalInColumn(long j)
{
    IntArrayList *detach = ColumnsIndexes [ j ];
    ColumnsIndexes [ j ] = NULL;
    return detach;
}

Ordering *SparseConectivityMtxII :: GetPermutationAndPattern(Ordering :: Type ord, IntArrayList *fixed)
{
    Ordering *order = NULL;
    if ( ord == Ordering :: MinimumDegree ) {
        //CS();
        Writeln(" ordering            : MinimumDegree");
        Write("Symbolic QG factorization   : ");
        clock_t start = MT.ClockStart();
        order = GenerateMD(fixed);
        Write( MT.MeasureClock(start) );

        //Writeln(MC_());
        //CS();
        //Write("Permuting..");
        order->cm = new SparseConectivityMtxII(* this, order);
        //Writeln(MC_());
    } else if ( ord == Ordering :: ApproxMinimumDegree ) {
        //CS();Write("Symbolic factorization...");
        Writeln(" ordering            : ApproxMinimumDegree");
        Write("Symbolic QG factorization   : ");
        order = GenerateAMD(fixed);
        //Writeln(MC_());
        //CS();
        //Write("Permuting..");
        order->cm = new SparseConectivityMtxII(* this, order);
        //Writeln(MC_());
    } else if ( ord == Ordering :: ApproxMinimumDegreeAA ) {
        //CS();Write("Symbolic factorization...");
        Writeln(" ordering            : ApproxMinimumDegree (aggressive absorbtion)");
        Write("Symbolic QG factorization   : ");
        order = GenerateAMD(fixed);
        GenerateFillInPresorted(order);
        order->cm = new SparseConectivityMtxII(* this, order);
    } else if ( ord == Ordering :: ApproxMinimumDegreeIncomplete ) {
        //CS();Write("Symbolic factorization...");
        //backup this matrix
        Writeln(" ordering            : ApproxMinimumDegreeIncomplete");
        Write("Symbolic QG factorization   : ");
        SparseConectivityMtxII *mtxA = new SparseConectivityMtxII(* this);
        order = GenerateAMD(fixed);
        //Writeln(MC_());
        //CS();
        //Write("Permuting..");
        order->cm = new SparseConectivityMtxII(* mtxA, order);
        delete mtxA;
        //Writeln(MC_());
    } else if ( ord == Ordering :: None ) {
        //CS();Write("Symbolic factorization...");
        Writeln(" ordering            : None");
        Write("Symbolic QG factorization   : ");
        order = Get_Unity();
        GenerateFillInPresorted(order);

        //Writeln(MC_());
        //CS();
        //Write("Permuting..");
        order->cm = new SparseConectivityMtxII(* this, order);
        //Writeln(MC_());
    } else if ( ord == Ordering :: ReverseCuthillMcKee ) {
        Writeln(" ordering            : Reverse Cuthill-McKee");
        Write("Symbolic QG factorization   : ");
        clock_t start = MT.ClockStart();
        order = Get_Reverse_Cuthill_McKee();
        Write( MT.MeasureClock(start) );
        Write("...");
        GenerateFillInPresorted(order);
        order->cm = new SparseConectivityMtxII(* this, order);
    } else if ( ord == Ordering :: CuthillMcKee ) {
        Writeln(" ordering            : Cuthill-McKee");
        Write("Symbolic QG factorization   : ");
        clock_t start = MT.ClockStart();
        order = Get_Cuthill_McKee();
        Write( MT.MeasureClock(start) );
        Write("...");
        GenerateFillInPresorted(order);
        order->cm = new SparseConectivityMtxII(* this, order);
    } else if ( ord == Ordering :: NestedGraphBisection ) {
        Writeln(" ordering            : RecursiveBiSection");
        Write("Symbolic QG factorization   : ");
        clock_t start = MT.ClockStart();
        order = Get_RecursiveBiSection();
        Write( MT.MeasureClock(start) );
        Write("...");
        GenerateFillInPresorted(order);
        order->cm = new SparseConectivityMtxII(* this, order);
    } else if ( ord == Ordering :: MetisND ) {
        Writeln(" ordering            : Metis (G.Karypis, V.Kumar)");
        Write("Graph ordering optimization : ");
        clock_t start = MT.ClockStart();
        order = Get_MetisDiSection();
        Write( MT.MeasureClock(start) );
        Write("...");
        GenerateFillInPresorted(order);
        order->cm = new SparseConectivityMtxII(* this, order);
    } else if ( ord == Ordering :: ColAMD ) {
        Writeln(" ordering            : ColAMD (S.I.Larimore, T.A.Davis)");
        Write("Graph ordering optimization : ");
        clock_t start = MT.ClockStart();
        order = Get_ColAMD();
        Write( MT.MeasureClock(start) );
        Write("...");
        GenerateFillInPresorted(order);
        order->cm = new SparseConectivityMtxII(* this, order);
    } else {
        Writeln("!!!! Incorrect ordering selected !!!!");
    }

    return order;
}


/// <summary>
/// Computes the minimum degree or the approximate minimum degree ordering
/// using the quotient graph with support for mass elimination and supervariables.
/// </summary>
MD_Qqraph :: MD_Qqraph(SparseConectivityMtxII *mtx)
{
    degrees = NULL;
    amd_w = NULL;
    pfixed = NULL;

    MinDegA = MinDegB = MinDegC = vi_Min = 0;
    A = E = I = L = NULL;
    pLIdx = pPIdx = NULL;
    elements = NULL;
    no_elements = 0;

    approximate_degree = false;
    keep_sorted_order = false;
    aggressive_absorbtion = false;
    this->mtx = mtx;
    this->n = mtx->N();
    n_1 = n - 1;
}

MD_Qqraph :: ~MD_Qqraph()
{
    if ( degrees ) {
        delete [] degrees;
        degrees = NULL;
    }

    if ( elements ) {
        delete elements;
        elements = NULL;
    }

    if ( amd_w ) {
        delete [] amd_w;
        amd_w = NULL;
    }

    for ( long i = 0; i < n; i++ ) {
        //if (A[i]) {delete A[i];A[i] = NULL;}
        if ( E [ i ] ) {
            delete E [ i ];
            E [ i ] = NULL;
        }

        if ( L [ i ] ) {
            delete L [ i ];
            L [ i ] = NULL;
        }

        if ( I [ i ] ) {
            delete I [ i ];
            I [ i ] = NULL;
        }
    }

    //delete [] A; A = NULL;
    delete [] E;
    E = NULL;
    delete [] L;
    L = NULL;
    delete [] I;
    I = NULL;
}

void MD_Qqraph :: Insert(IntArrayList &al, long i)
{
    long myIndex = al.BinarySearch(i);
    if ( myIndex < 0 ) {
        al.Insert(~myIndex, i);
    }
}

void MD_Qqraph :: ComputeAMD_w_Le_Lp(long i)
{
    IntArrayList &el = * E [ i ];
    for ( long ie = el.Count - 1; ie >= 0; ie-- ) {
        long e = el [ ie ];
        if ( amd_w [ e ] < 0 ) {
            amd_w [ e ] = L [ e ]->Count;
        }

        amd_w [ e ] -= I [ i ]->Count;
        //aggresive absorbtion
        if ( aggressive_absorbtion && amd_w [ e ] == 0 ) {
            el.RemoveAt(ie);
            amd_w [ e ] = -1;
        }
    }
}

void MD_Qqraph :: ClearAMD_w(long i)
{
    IntArrayList &el = * E [ i ];
    for ( long ie = el.Count - 1; ie >= 0; ie-- ) {
        amd_w [ el [ ie ] ] = -1;
    }
}

long MD_Qqraph :: ApproximateMinumumDegree(long i, IntArrayList &Lp)
{
    IntArrayList &Ali = * A [ i ];
    IntArrayList &Eli = * E [ i ];

    long idx;
    //long na = Ali.Count;
    long ne = Eli.Count;
    long *pe = Eli.Items;
    {
        long Aic = 0;
        //foreach(long j in A[i]) if (pPIdx[j]==0) Aic++;
        Aic += Ali.CountWithoutMask(pPIdx);

        long Lpc = 0;
        //foreach(long j in Lp) if (pPIdx[j]==0) Lpc++;
        for ( idx = 0; idx < Lp.Count; idx++ ) {
            long j = Lp [ idx ];
            if ( pPIdx [ j ] == 0 ) {
                Lpc++;
            }
        }

        long d = min(n - this->no_elements, degrees [ i ] + Lpc);

        long Les = 0;
        //foreach(long e in E[i]) if (amd_w[e]<0) Les += L[e].Length;   else Les += amd_w[e];
        for ( idx = ne - 1; idx >= 0; idx-- ) {
            long e = pe [ idx ];
            if ( amd_w [ e ] < 0 ) {
                Les += L [ e ]->Count;
            } else {
                Les += amd_w [ e ];
            }
        }

        d = min(d, Aic + Lpc + Les);

        if ( MinDegB > d && IsFree(i) ) {
            MinDegB = d;
            vi_Min = i;
        }

        return d;
    }
}

/// <summary> This method computes the exact external degree </summary>
/// <param name="i">variable</param>
/// <returns>degree</returns>
long MD_Qqraph :: ExternalDegree(long i)
{
    long idx, idx2;
    long d = 0;
    long na = A [ i ]->Count;
    long ne = E [ i ]->Count;
    long *pa = A [ i ]->Items, *pe = E [ i ]->Items;
    {
        //foreach(long a in A[i]) if(pPIdx[a]==0)       {pLIdx[a] = 1;d++;}
        for ( idx = na - 1; idx >= 0; idx-- ) {
            long a = pa [ idx ];
            if ( pPIdx [ a ] == 0 ) {
                pLIdx [ a ] = 1;
                d++;
            }
        }

        //foreach(long e in E[i])       foreach(long j in L[e])if (pPIdx[j]==0 && pLIdx[j]!=1)  {pLIdx[j] = 1;d++;}
        for ( idx = ne - 1; idx >= 0; idx-- ) {
            long e = pe [ idx ];
            //foreach(long j in L[e])
            IntArrayList &Lp = * L [ e ];
            for ( idx2 = 0; idx2 < Lp.Count; idx2++ ) {
                long j = Lp [ idx2 ];
                if ( pPIdx [ j ] == 0 && pLIdx [ j ] != 1 ) {
                    pLIdx [ j ] = 1;
                    d++;
                }
            }
        }

        //foreach(long a in A[i]) pLIdx[a] = 0;
        A [ i ]->SetIndexesTo(pLIdx, 0);

        //foreach(long e in E[i]) foreach(long j in L[e]) if(pPIdx[j]==0) pLIdx[j] = 0;
        for ( idx = ne - 1; idx >= 0; idx-- ) {
            long e = pe [ idx ];
            //foreach(long j in L[e])
            IntArrayList &Lp = * L [ e ];
            for ( idx2 = 0; idx2 < Lp.Count; idx2++ ) {
                long j = Lp [ idx2 ];
                if ( pPIdx [ j ] == 0 ) {
                    pLIdx [ j ] = 0;
                }
            }
        }

        if ( MinDegB > d && IsFree(i) ) {
            MinDegB = d;
            vi_Min = i;
        }
    }
    return d;
}

/// <summary> Masselimination </summary>
/// <param name="p">pivot element</param>
void MD_Qqraph :: Eliminate(long p)
{
    if ( I [ p ] == NULL ) {
        return;
    }

    long ie, it, ii;

    // set I[p] marks
    IntArrayList &Ip = * I [ p ];

    Ip.SetIndexesTo(pLIdx, 1);

    // Compute L[p] = A[p] u (U(e e E[p])L[e]) / Vp
    IntArrayList &Lp = * A [ p ];
    A [ p ] = NULL;
    Lp.RemoveMarkByBattern(pLIdx);

    IntArrayList &Ep = * E [ p ];
    //foreach(long e in E[p])
    for ( ie = Ep.Count - 1; ie >= 0; ie-- ) {
        long e = Ep.Items [ ie ];
        pLIdx [ e ] = 1;
        //foreach (long i in L[e])
        IntArrayList &Le = * L [ e ];
        for ( it = 0; it < Le.Count; it++ ) {
            long i = Le [ it ];
            if ( pLIdx [ i ] == 0 ) {
                Lp.SortInsert(i);                //Lp.Add(i);
                pLIdx [ i ] = 1;
            }
        }
    }

    // pLIdx now contains L[p] and I[p] and E[p] marks
    L [ p ] = & Lp;
    //foreach(long i in L[p])
    for ( it = 0; it < Lp.Count; it++ ) {
        long i = Lp [ it ];
        if ( I [ i ] == NULL ) {
            continue;
        }

        // A[i] = (A[i] / L[p]) / I[p]
        A [ i ]->RemoveByBattern(pLIdx);
        // element absorbtion E[i] = (E[i] \ E[p]) u {p}
        E [ i ]->RemoveByBattern(pLIdx);
        Insert(* E [ i ], p);

        // This is for approximate minimum degree
        if ( !keep_sorted_order && approximate_degree ) {
            ComputeAMD_w_Le_Lp(i);
        }
    }

    // clear L[p] and I[p] and E[p] marks
    Lp.SetIndexesTo(pLIdx, 0);
    I [ p ]->SetIndexesTo(pLIdx, 0);
    E [ p ]->SetIndexesTo(pLIdx, 0);

    if ( !keep_sorted_order ) {
        // compute degree
        //foreach(long i in L[p])
        for ( it = 0; it < Lp.Count; it++ ) {
            long i = Lp [ it ];
            if ( I [ i ] == NULL ) {
                continue;
            }

            I [ i ]->SetIndexesTo(pPIdx, 1);            //set I[i] mask
            if ( approximate_degree ) {
                this->degrees [ i ] = ApproximateMinumumDegree(i, Lp);
            } else {
                this->degrees [ i ] = ExternalDegree(i);
            }

            I [ i ]->SetIndexesTo(pPIdx, 0);            //clear I[i] mask
        }

        if ( approximate_degree ) {     // null it again
            //foreach(long i in L[p])
            for ( it = 0; it < Lp.Count; it++ ) {
                long i = Lp [ it ];
                if ( I [ i ] != NULL ) {
                    ClearAMD_w(i);
                }
            }
        }

        //supervariable detection, pairs found via hash function
        SupervariablesDetection(Lp);
    }

    IntArrayList *Lpfull = new IntArrayList(Lp);

    //Make clique from the supervariable
    //Eliminate all in clique
    for ( ii = Ip.Count - 1; ii >= 0; ii-- ) {
        long k = Ip.Items [ ii ];
        Insert(* Lpfull, k);
    }

    Lpfull->TrimToSize();
    mtx->ColumnsIndexes [ p ] = Lpfull;

    for ( ii = Ip.Count - 1; ii >= 0; ii-- ) {
        long k = Ip.Items [ ii ];
        elements->Items [ no_elements++ ] = k;
        if ( k != p ) {
            mtx->ColumnsIndexes [ k ]->Fill(* Lpfull);
        }
    }

    variables.Remove(p);
    //if (A[p]) {delete A[p];A[p] = NULL}
    //A[p] = NULL;
    if ( E [ p ] ) {
        delete E [ p ];
        E [ p ] = NULL;
    }

    if ( I [ p ] ) {
        delete I [ p ];
        I [ p ] = NULL;
    }
}

bool MD_Qqraph :: IsIndistinguishable(long i, long j)
{
    if ( A [ i ]->Count != A [ j ]->Count ) {
        return false;
    }

    if ( E [ i ]->Count != E [ j ]->Count ) {
        return false;
    }

    bool ind = true;

    pPIdx [ j ] = pPIdx [ i ] = 1;
    A [ i ]->SetIndexesTo(pPIdx, 1);            //foreach(long k in A[i]) pPIdx[k] = 1;
    E [ i ]->SetIndexesTo(pPIdx, 2);            //foreach(long e in E[i]) pPIdx[e] = 2;

    A [ j ]->TestSetIndexesTo(pPIdx, 1, 0, ind); //foreach(long k in A[j]) {if (pPIdx[k] != 1) ind = false;pPIdx[k] = 0;}
    E [ j ]->TestSetIndexesTo(pPIdx, 2, 0, ind); //foreach(long e in E[j]) {if (pPIdx[e] != 2) ind = false;pPIdx[e] = 0;}

    pPIdx [ j ] = pPIdx [ i ] = 0;
    if ( ind ) {
        return true;
    }

    A [ i ]->SetIndexesTo(pPIdx, 0);            //foreach(long k in A[i]) pPIdx[k] = 0;
    E [ i ]->SetIndexesTo(pPIdx, 0);            //foreach(long e in E[i]) pPIdx[e] = 0;
    return false;
}

void MD_Qqraph :: SupervariablesDetection(IntArrayList &Lp)
{
    if ( keep_sorted_order ) {
        return;
    }

    long idx, ix;

    //supervariable detection, pairs found via hash function
    //foreach(long i in Lp)
    for ( idx = 0; idx < Lp.Count; idx++ ) {
        long i = Lp [ idx ];
        if ( I [ i ] == NULL || !IsFree(i) ) {
            continue;
        }

        long key = Hash(i);
        ht.AddValue(key, i);
    }

    for ( idx = ht.occupied_buckets.Count - 1; idx >= 0; idx-- ) {
        IntArrayList &al = * ht.buckets [ ht.occupied_buckets.Items [ idx ] ];
        if ( al.Count >= 2 ) {
            hash_parents.Alloc(al.Count);
            for ( ix = 0; ix < al.Count; ix++ ) {
                long i = al [ ix ];
                if ( hash_parents [ ix ] < 0 ) {
                    continue;
                }

                for ( long jx = ix + 1; jx < al.Count; jx++ ) {
                    if ( hash_parents [ jx ] < 0 ) {
                        continue;
                    }

                    long j = al [ jx ];
                    if ( IsIndistinguishable(i, j) ) {
                        hash_parents [ jx ] = ~i;
                    }
                }
            }

            for ( ix = 0; ix < al.Count; ix++ ) {
                if ( hash_parents [ ix ] < 0 ) {
                    long i = ~hash_parents [ ix ];
                    long j = al [ ix ];

                    //remove supervariable j
                    IntArrayList &Ij = * I [ j ];
                    for ( long ij = Ij.Count - 1; ij >= 0; ij-- ) {
                        I [ i ]->SortInsert(Ij.Items [ ij ]);
                    }

                    degrees [ i ] -= I [ j ]->Count;
                    if ( MinDegB > degrees [ i ] ) {
                        MinDegB = degrees [ i ];
                        vi_Min = i;
                    }

                    variables.Remove(j);
                    if ( I [ j ] ) {
                        delete I [ j ];
                        I [ j ] = NULL;
                    }

                    //if (A[j]) {delete A[j];A[j] = NULL;}
                    //A[j] = NULL;
                    if ( E [ j ] ) {
                        delete E [ j ];
                        E [ j ] = NULL;
                    }
                }
            }
        }
    }

    ht.Clear();
}


/// <summary>
/// Returns the permutation according to minimum degree algorithm.
/// </summary>
/// <param name="approximate_degree"> true=approximate degree, false=true external degree</param>
/// <returns>permutation vector</returns>
IntArrayList *MD_Qqraph :: GenerateMD(bool approximate_degree, IntArrayList *fixed)
{
    this->approximate_degree = approximate_degree;
    no_elements = 0;

    ht.Init(n);

    // temporarly used as the index array
    // eventually this array contains the permutation

    long *permutation = new long [ n ];

    // this array stores elements on position [0..no_elements-1] and then variables [no_elements..n]
    // Finally this array contains Order of eliminated elemtents
    long *vlasts = new long [ n ];
    long *vnexts = new long [ n ];

    // number of fixed nodes which must be ordered to the end of the list
    long no_fixed = 0;
    if ( fixed ) {
        pfixed = new long [ n ];
        Array :: Clear(pfixed, 0, n);
        fixed->SetIndexesTo(pfixed, 1);
        no_fixed = fixed->Count;
    }

    elements = new IntArrayList(n);
    elements->Alloc();
    long *pivot_pattern = new long [ n ];
    memset( pivot_pattern, 0, n * sizeof( long ) );
    degrees = new long [ n ];
    amd_w = new long [ n ];

    MinDegA = n;
    MinDegB = n;
    MinDegC = n;

    this->E = new IntArrayList * [ n ];
    this->L = new IntArrayList * [ n ];
    this->A = mtx->ColumnsIndexes;
    this->I = new IntArrayList * [ n ];
    // in the begining all DOFs are variables
    for ( long k = 0; k < n; k++ ) {
        I [ k ] = new IntArrayList(4);
        I [ k ]->Add(k);
        E [ k ] = new IntArrayList(4);
        //A[k] = mtx->ColumnsIndexes[k];
        L [ k ] = NULL;

        vlasts [ k ] = k - 1;
        vnexts [ k ] = k + 1;

        degrees [ k ] = A [ k ]->Count;
        if ( degrees [ k ] > n ) {
            mtx->Writeln("Connectivity matrix is corrupt !");
            return NULL;
        }

        if ( degrees [ k ] < MinDegA && IsFree(k) ) {
            MinDegA = degrees [ k ];
        }

        permutation [ k ] = 0;
        amd_w [ k ] = -1;
    }

    vnexts [ n - 1 ] = -1;

    // Draw the quotient graph
    //mtx.MT.DrawGraph(this);

    //long noDots = 0;
    //long denDots = n / 24; if (denDots==0) denDots = 1;

    long *pLIdx = permutation, *pPIdx = pivot_pattern, *vl = vlasts, *vn = vnexts;
    {
        this->pLIdx = pLIdx;
        this->pPIdx = pPIdx;
        variables.Init(vl, vn, n);

        if ( keep_sorted_order ) {
            for ( long i = 0; i < n; i++ ) {
                Eliminate(fixed->Items [ i ]);
            }
        } else {
            while ( no_elements < n ) {
                long MinimDeg = min(MinDegB, MinDegA);
                MinDegB = n;
                MinDegC = n;

                // find minimum degree
                for ( long vi = variables.first; vi >= 0; ) {
                    //long deg_vi = p_degrees[vi];
                    if ( IsFree(vi) && ( MinimDeg >= degrees [ vi ] ) ) {
                        // do elimination
                        Eliminate(vi);
                        if ( MinDegB < MinimDeg ) {
                            MinimDeg = MinDegB;
                            vi = vi_Min;
                            continue;
                        }
                    } else {
                        if ( IsFree(vi) && ( degrees [ vi ] < MinDegC ) ) {
                            MinDegC = degrees [ vi ];
                        }
                    }

                    vi = variables.next [ vi ];
                }

                MinDegA = MinDegC;

                if ( no_fixed && ( no_elements == n - no_fixed ) ) { // we can now unblock the fixed nodes
                    delete [] pfixed;
                    pfixed = NULL;
                    no_fixed = 0;
                }

                //for (long i=no_elements/denDots; i>noDots; noDots++)
                //{
                //      mtx->Write(".");
                //mtx.MT.NextStep();
                //}
            }
        }
    }
    // Permutation = order^-1
    delete [] pfixed;
    pfixed = NULL;

    delete [] permutation;
    delete [] vlasts;
    delete [] vnexts;

    delete [] pivot_pattern;

    IntArrayList *order = elements;
    elements = NULL;
    // return the order
    return order;
}


DSS_NAMESPASE_END
