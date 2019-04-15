# PerforAD
Automatic Differentiation for high-performance stencil loops

## Why would I need this?
Stencil computations are a common motif in computing. For example, a Finite Difference heat equation solver would look something like this:
```
for i in 2,n
  T_new[i] = (T_old[i-1] - 2*T_old[i] + T_old[i+1])/dx^2
end for
```
Note how the loop body writes to an array index that is also the loop counter, and gathers data from a neighbourhood of array indices that are constant offsets of the loop counter. A similar pattern emerges in many applications, including linear algebra, image processing, convolutional neural networks, etc.

The above loop is trivial to parallelise using shared-memory parallelism (OpenMP, CUDA, ...), since each iteration writes exclusively to its own index in the output array. When we use automatic differentiation to compute the adjoint of this loop, we instead get something like this:
```
for i in 2,n
  T_old_b[i-1] += T_new_b[i]/dx^2
  T_old_b[i] += T_new_b[i]/dx^2
  T_old_b[i+1] += T_new_b[i]/dx^2
end for
```
Note how each iteration now scatters data to neighbour indices. There are several options to parallelise such a loop, none of which is without problems. The operation could be implemented as a sum-reduction to `T_old_b`, requiring several times as much memory as the original code. Alternatively, some synchronisation (atomic updates, critical sections, barriers) could be used, but slow down execution. PerforAD is here to offer another, hopefully more attractive alternative.

## What does PerforAD do?
PerforAD splits the derivative computation of the loop body into several parts, each of which updates only one item. It then uses a combintation of loop splitting, index substitutions, and loop merging to recombine the computation into a core loop that only uses gather operations and therefore looks like a stencil computation, and a set of remainder loops that perform the correct treatment at boundaries.

## Why the name?
The name *stencil computation* is probably derived from the printing technique. From Wikipedia:
> Stencilling produces an image or pattern by applying pigment to a surface over an intermediate object with designed gaps in it which create the pattern or image by only allowing the pigment to reach some parts of the surface.

In the same way that stencilling in printing produces the same pattern wherever it is applied, a stencil computation has a fixed pattern in which it reads and writes data as it is applied throughout a domain.

When differentiating a stencil, PerforAD (*perforate*) creates a new stencil with the same shape as the old one. Think of it as a way to cut out parts of your adjoint computation that have the same shape as your primal computation.
Also, PerforAD results in really fast, high-performance Automatic Differentiation (AD).
