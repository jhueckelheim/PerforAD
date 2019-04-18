#   if defined(__cplusplus)
        extern "C" {
#   endif

/** Re-base waiting pointers *pp and *ppb
 * from their old base from the forward sweep
 * to their new base in the backward sweep.
 * When new base not available yet, schedules
 * this to be done when new base is allocated.
 * Usage: restoring a pointer pp (i.e. declared as T* pp ;)
 *  becomes/creates in the forward sweep:
 *  pushpointer8(pp) ;
 *  and becomes/creates in the backward sweep:
 *  poppointer8((void**)(&pp)) ;
 *  ADMM_Rebase((void**)(&pp), (void**)(&ppb)) ; */
#ifdef ADTOOL_TAPENADE
void ADMM_Rebase(void **pp, void **ppb);
#else
void ADMM_Rebase(void **pp);
#endif

/** Forward sweep correspondent of standard malloc().
 * Keeps track of base and size of the allocated chunk.
 * Usage: T* x = (T*)malloc(n*sizeof(T)) ;
 *  becomes/creates in the forward sweep:
 *  T* x = (T*)*FW_ADMM_Allocate(<if_dynamicMemory?0:x>, n*sizeof(T), <if_dynamicMemory?1:0>) ; */
void* FW_ADMM_Allocate(void *chunk, int size, int isDynamic);

/** Backward sweep correspondent of standard malloc().
 * Frees the chunk based at newbase, as well as its adjoint if presen.
 * Usage: T* x = (T*)malloc(n*sizeof(T)) ;
 *  becomes/creates in the backward sweep:
 *  BW_ADMM_Allocate((void*)x, (void*)xb, <if_dynamicMemory?1:0>) ; */
#ifdef ADTOOL_TAPENADE
void BW_ADMM_Allocate(void *newbase, void *newadjbase, int isDynamic);
#else
void BW_ADMM_Allocate(void *newbase, int isDynamic);
#endif

/** Forward sweep correspondent of standard free().
 * Pushes values from chunk before it is freed, if TBR.
 * Pushes base and size.
 * Usage: free(x) ;
 *  becomes/creates in the forward sweep:
 *  FW_ADMM_Deallocate((void*)x, <if_tbr?1:0>, <if_dynamicMemory?1:0>) ; */ 
 void FW_ADMM_Deallocate(void *base, int tbr, int isDynamic);


/** Backward sweep correspondent of standard free().
 * Pops old base and size.
 * Re-allocates a chunk, and an adjoint chunk if present.
 * If tbr is true, pops values into the new chunk.
 * If pointers is true, calls for rebase of the pointer in the chunk [better interface wanted !].
 * Remembers correspondence from old base to
 * new base and new adjoint base.
 * Re-base waiting pointers.
 * Usage: free(x) ;
 *  becomes/creates in the backward sweep:
 *  x = BW_ADMM_Deallocate(<if_dynamicMemory?0:x>, (void**)&xb, <if_tbr?1:0>,
 *                         <if_containsPointers?1:0>, <if_dynamicMemory?1:0>) ; */
#ifdef ADTOOL_TAPENADE
void* BW_ADMM_Deallocate(void *chunk, void **chunkb, int tbr, int pointers, int isDynamic);
#else
void* BW_ADMM_Deallocate(void *chunk, int tbr, int pointers, int isDynamic);
#endif

#ifndef ADTOOL_TAPENADE
#define allocateBytesForSizeOf_SgTypeInt(val, x) *(val)= x*sizeof(int);
#define allocateBytesForSizeOf_SgTypeDouble(val, x) *(val)= x*sizeof(double);
#endif

#   if defined(__cplusplus)
        }
#   endif
