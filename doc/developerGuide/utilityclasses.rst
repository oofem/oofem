Utility classes
================

Vectors and Matrices
--------------------

The OOFEMlib provides the abstraction for integer vectors (`IntArray`)
and for real vectors and matrices (`FloatArray` and
`FloatMatrix` classes). The usual arithmetic operations (like vector and
matrix additions, subtractions, matrix vector multiplication, matrix
inverse, solution of linear system, finding eigenvalues and
eigenvectors) are provided. However, the usual math operators ('+','*')
are not overloaded, user has to call specific routines. 
The vector and matrices provide both 0-based and 1-based component
access, they allow for dynamic resize with optional chunk. In fact,
the current implementation of `resize` only grows the
receiver, the possible request for shrink does not cause the
reallocation of memory, since allocated memory is kept for future possible resize. If
resize operation wants to force allocation, then `hardResize` method
should be invoked instead. 
The preferred argument passing is by reference, even for function or
procedure return values.
The called function performs resize and fills up
the return parameter(s). This is motivated by aim to avoid memory
allocation/de-allocation problems. The programer should always avoid to
return pointers to newly allocated arrays or matrices. 
See section `Coding Standards`_ for
details.

Solution steps
-------------- 
The class representing solution step is `TimeStep`. It may represent either 
time step, load increment, or load case depending on current Engineering model used.

Solution step maintain the reference to corresponding Engineering model class instance.
It maintain also its "intrinsic time" and corresponding time increment. The meaning of these 
values is dependent on current Engineering model used. The time may represent either
current time, load increment number or load case number. See corresponding 
EngngModel reference for details.
	
Some components (typically integration points real stresses or integration points non-local values)
are computationally very demanding. Because in typical code, there are number of requests for same value 
during the computation process, it may be efficient to store these values and compute them only once.
The principal problem is to recognize, when is necessary to re-compute these stored values to reflect 
newly reached state. This cannot be determined form solution step "time", because solution step may 
represent for example load increment, inside which typically many iterations are needed to reach 
convergence. For this purpose, a concept of solution state counters is provided.
Whenever the solution state changes, the engineering model updates the solution state counter.
The solution state counter is guaranteed to grow up smoothly (it newer decreases) during solution process.
Other components of program (integration points) can then store their computationally expensive values
but have to store also corresponding solution state counter value valid when these were computed.
Then their can easily check for difference between freezed solution state counter for their value with 
current solution state requested from solution step and recompute the values if necessary.

Load Time Functions
-------------------
Abstract base class representing load time function. Classes derived from Load class typically 
describe load from spatial point of view. The purpose of introducing load time function is to express
variation of some components in time. Load time function typically belongs to domain and is 
attribute of one or more loads. Generally load time function is real function of time (:math:`y=f(t)`).

