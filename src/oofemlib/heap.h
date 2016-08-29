#ifndef heap_h
#define heap_h

namespace oofem {
/**
 * Class implementing a heap, which is an auxiliary data structure
 * used for efficient sorting and exploited e.g. by fast marching methods.
 * The main purpose is to sort large lists of real numbers efficiently,
 * at N*log(N) algorithmic complexity.
 */
  
class Heap
{
public:
    /// Constructor and destructor
    Heap(int N);
    ~Heap();

    /// Heap arrays:
    /// Keys contains certain real values that need to be sorted.
    /// The heap is organized according to these, using heap-sort mechanisms.
    /// H2T and T2H contain pointers to and from H and T.
    double *Keys;
    int *H2T;
    int *T2H;

    /// Interface with external algorithms (such as fast marching)
    bool isInHeap(int Ind);
    int nElems();
    void insert(double time, int Ind);
    void update(double time, int Ind);
    double getSmallest(int *Ind);

    /// Debugging tools
    int checkHeapProperty(int pInd);
    void print();
    void printTree();
    double *formMatrix(int m, int n);

private:
    /// Variables that control the memory allocation for the heap
    int Initial_Heap_Alloc_Size, allocatedSize;

    /// Keeps track of the number of elements in heap
    int heapCount;

    /// Elementary heap operations
    void upHeap(int Ind);
    void downHeap(int Ind);
    void swapElements(int Ind1, int Ind2);

    /// Index calculations
    int parentInd(int inInd);
    int leftChildInd(int inInd);
    int rightChildInd(int inInd);
    int lastParentInd();

    /// Used by printTree. Does the actual printing, top-down.
    void recurse(int row, int pad, int spacing, int S);
};
 
} // end namespace oofem
#endif // heap_h
