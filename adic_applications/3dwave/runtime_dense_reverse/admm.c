#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "admm.h"
#ifdef ADTOOL_TAPENADE
#include "ADFirstAidKit/adBuffer.h"
#include "ADFirstAidKit/adStack.h"
#else
#include "ad_tape.h"
#endif

/** Cell of a chained list of void* */
typedef struct _ADMM_Cell {
  void* head ;
  struct _ADMM_Cell* tail ;
} ADMM_Cell ;

typedef struct {
  void* base ;
  int   size ;
  void* newbase ;
  void* newadjbase ;
} ADMM_ChunkInfo ;

typedef struct {
  void** pp ;
  void*  oldvalpp ;
  void** ppb ;
} ADMM_WaitingAddress ;

/** The first cell of the list of ChunkInfo's. It is just a hat */
ADMM_Cell firstChunkInfoCell = {NULL, NULL} ;
/** The current list of ChunkInfo's. */
ADMM_Cell* chunkInfoList = &firstChunkInfoCell ;

/** The first cell of the list of waiting pointers. It is just a hat */
ADMM_Cell firstWaitingRebaseCell = {NULL, NULL} ;
/** The current list of pointers waiting to be rebased */
ADMM_Cell* waitingRebaseList = &firstWaitingRebaseCell ;

/** TEMPORARY: FOR DEBUGGING. */
void dumpWaiting() ;
void dumpChunkInfos() ;

/** Registers into the chunkInfoList the data about a new chunk */
void ADMM_RegisterNewChunkInfo(void *base, int size, void* newbase, void* newadjbase) {
  ADMM_Cell* inChunkInfos = chunkInfoList ;
  ADMM_Cell* newCell = (ADMM_Cell*)malloc(sizeof(ADMM_Cell)) ;
  ADMM_ChunkInfo* newChunkInfo = (ADMM_ChunkInfo*)malloc(sizeof(ADMM_ChunkInfo)) ;
  newChunkInfo->base = base ;
  newChunkInfo->size = size ;
  newChunkInfo->newbase = newbase ;
  newChunkInfo->newadjbase = newadjbase ;
  newCell->head = (void*)newChunkInfo ;
  // Find the location to insert to, to keep chunkInfoList ordered by increasing base's:
  while (inChunkInfos->tail && ((ADMM_ChunkInfo*)inChunkInfos->tail->head)->base<base)
    inChunkInfos = inChunkInfos->tail ;
  newCell->tail = inChunkInfos->tail ;
  inChunkInfos->tail = newCell ;
/*   dumpChunkInfos() ; */
}

/** Frees the next Cell and its attached ChunkInfo */
void ADMM_FreeNextChunkInfo(ADMM_Cell *inChunkInfos) {
  free(inChunkInfos->tail->head) ;
  ADMM_Cell *newTail = inChunkInfos->tail->tail ;
  free(inChunkInfos->tail) ;
  inChunkInfos->tail = newTail ;
}

/** Registers the given pointer,adjoint-pointer pair as waiting for rebase. */
void ADMM_RegisterNewWaitingRebase(void** pp, void** ppb) {
  ADMM_Cell* inWaitingRebases = waitingRebaseList ;
  ADMM_Cell* newCell = (ADMM_Cell*)malloc(sizeof(ADMM_Cell)) ;
  ADMM_WaitingAddress* newWaiting = (ADMM_WaitingAddress*)malloc(sizeof(ADMM_WaitingAddress));
  newWaiting->pp = pp ;
  newWaiting->oldvalpp = *pp ;
  newWaiting->ppb = ppb ;
  newCell->head = (void*)newWaiting ;
  // Find the location to insert to, to keep waitingRebaseList ordered by increasing pp's:
  while (inWaitingRebases->tail && ((ADMM_WaitingAddress*)inWaitingRebases->tail->head)->pp<pp)
    inWaitingRebases = inWaitingRebases->tail ;
  newCell->tail = inWaitingRebases->tail ;
  inWaitingRebases->tail = newCell ;
}

/** Frees the next Cell and its attached pair of pointers. */
void ADMM_FreeNextWaitingRebase(ADMM_Cell *inWaitingRebases) {
  free(inWaitingRebases->tail->head) ;
  ADMM_Cell *newTail = inWaitingRebases->tail->tail ;
  free(inWaitingRebases->tail) ;
  inWaitingRebases->tail = newTail ;
}

/** Rebase the given pointer (and adjoint-pointer if present) coherently
 * with the base changed to newbase and its adjoint to newadjbase. */
void ADMM_DoRebase(void **pp, void *oldvalpp, void **ppb,
                    void* base, void* newbase, void* newadjbase) {
/*   printf("      now rebasing addresses in %x(%x) by adding offset (%x-%x=%i) to %x(%x)", */
/*          pp,ppb,oldvalpp,base,(oldvalpp-base),newbase,newadjbase) ; */
  if (newbase) {
    *pp = (char*)newbase+((char*)oldvalpp-(char*)base) ;
    if (ppb && newadjbase)
      *ppb = (char*)newadjbase+((char*)oldvalpp-(char*)base) ;
  }
/*   printf(", giving %x(%x) \n",*pp,(ppb?*ppb:0)) ; */
}





/** Re-base waiting pointers *pp and *ppb
 * from their old base from the forward sweep
 * to their new base in the backward sweep.
 * When new base not available yet, schedules
 * this to be done when new base is allocated. */
#ifdef ADTOOL_TAPENADE
void ADMM_Rebase(void **pp, void **ppb) {
#else
void ADMM_Rebase(void **pp) {
  void **ppb = NULL;
#endif
  // [Pb:] possible bug: the same pointer might be rebased twice if aliased!
  // Possible solution is to keep track of pointers already rebased.
  ADMM_ChunkInfo* chunkInfo ;
  ADMM_Cell* inChunkInfos = chunkInfoList->tail ;
/* dumpChunkInfos() ; */
/* dumpWaiting() ; */
/* printf("  Try to rebase pointers %x->%016x (%x->?)\n", pp, *pp, ppb) ; */
  while (inChunkInfos) {
    chunkInfo = (ADMM_ChunkInfo*)inChunkInfos->head ;
/* printf("  with chunk %x+>%i ==> %x+>? (%x+>?)\n",chunkInfo->base,chunkInfo->size, chunkInfo->newbase, chunkInfo->newadjbase) ; */
    if (*pp < (chunkInfo->base+chunkInfo->size)) {
      if (*pp >= chunkInfo->base) {
        ADMM_DoRebase(pp, *pp, ppb,
                       chunkInfo->base, chunkInfo->newbase, chunkInfo->newadjbase) ;
        pp = NULL ;
      }
      break ;
    }
    inChunkInfos = inChunkInfos->tail ;
  }
  // if the containing chunk is not found, it may appear later:
  if (pp) ADMM_RegisterNewWaitingRebase(pp, ppb) ;
}

/** Forward sweep correspondent of standard malloc().
 * Keeps track of base and size of the allocated chunk */
void* FW_ADMM_Allocate(void *chunk, int size, int isDynamic) {
  if(isDynamic!=0)
    chunk = malloc(size) ;
  ADMM_RegisterNewChunkInfo(chunk, size, NULL, NULL) ;
/* printf("Registered allocated address %x+>%i\n",chunk,size) ; */
  return chunk ;
}

/** Backward sweep correspondent of standard malloc().
 * Frees the chunk based at newbase, as well as its adjoint if present. */
#ifdef ADTOOL_TAPENADE
void BW_ADMM_Allocate(void *newbase, void *newadjbase, int isDynamic) {
#else
void BW_ADMM_Allocate(void *newbase, int isDynamic) {
  void *newadjbase = NULL;
#endif
  ADMM_ChunkInfo* chunkInfo ;
  ADMM_Cell* inChunkInfos = chunkInfoList ;
  int found = 0 ;
  while (!found && inChunkInfos->tail) {
    chunkInfo = (ADMM_ChunkInfo*)inChunkInfos->tail->head ;
    found = ((newbase && newbase==chunkInfo->newbase) ||
             (newadjbase && newadjbase==chunkInfo->newadjbase)) ;
    if (!found) inChunkInfos = inChunkInfos->tail ;
  }
  if (!found) {printf("Chunk not found in BW_ADMM_Allocate!\n"); exit(0);}
/* printf("re-de-Register allocated addresses %x(%x) +>%i\n",newbase,newadjbase,chunkInfo->size) ; */
  ADMM_FreeNextChunkInfo(inChunkInfos) ;
  if (newbase && isDynamic!=0) free(newbase) ;
  if (newadjbase&& isDynamic!=0) free(newadjbase) ;
}

/** Forward sweep correspondent of standard free().
 * Pushes values from chunk before it is freed, if TBR.
 * Pushes base and size. */ 
void FW_ADMM_Deallocate(void *base, int tbr, int isDynamic) {
  ADMM_ChunkInfo* chunkInfo ;
  ADMM_Cell* inChunkInfos = chunkInfoList ;
  while (inChunkInfos->tail) {
    chunkInfo = (ADMM_ChunkInfo*)inChunkInfos->tail->head ;
    if (chunkInfo->base >= base) {
      if (chunkInfo->base!=base) inChunkInfos = NULL ;
      break ;
    }
    inChunkInfos = inChunkInfos->tail ;
  }
  if (!inChunkInfos) {printf("Chunk not found in FW_ADMM_Deallocate!\n"); exit(0);}
#ifdef ADTOOL_TAPENADE
  if (tbr) pushNarray(base, (unsigned int)chunkInfo->size) ;
#endif
/* printf("de-Register and push allocated address %x size %i\n",chunkInfo->base,chunkInfo->size) ; */
#ifdef ADTOOL_TAPENADE
  pushinteger4(chunkInfo->size) ;
  pushpointer8((char*)base) ;
#else
  push_i_s0(chunkInfo->size) ;
  push_p_s0((char*)base) ;
#endif
  if(isDynamic!=0)
    free(base);
  ADMM_FreeNextChunkInfo(inChunkInfos) ;
}

/** Backward sweep correspondent of standard free().
 * Pops old base and size.
 * Re-allocates a chunk, and an adjoint chunk if present.
 * If tbr is true, pops values into the new chunk.
 * If pointers is true, calls for rebase of the pointer in the chunk.
 * Remembers correspondence from old base to
 * new base and new adjoint base.
 * Re-base waiting pointers. */
#ifdef ADTOOL_TAPENADE
void* BW_ADMM_Deallocate(void *chunk, void **chunkb, int tbr, int pointers, int isDynamic) {
#else
void* BW_ADMM_Deallocate(void *chunk, int tbr, int pointers, int isDynamic) {
  void **chunkb = NULL;
#endif
  void* oldBase ;
  int size ;
#ifdef ADTOOL_TAPENADE
  poppointer8((void**)&oldBase) ;
  popinteger4(&size) ;
#else
  pop_p_s0(oldBase) ;
  pop_i_s0(size) ;
#endif
  if (isDynamic!=0 && chunkb!=NULL) {
    *chunkb = malloc(size) ;
    // [Pb:] here we should initialize to zero the derivative *chunkb, but how?
    // int i ;
    // for (i=0 ; i<size ; ++i) (*chunkb)(i) = 0.0 ;
  }
  if (isDynamic!=0) {
    chunk = malloc(size) ;
  }
  #ifdef ADTOOL_TAPENADE   
    if (tbr) popNarray((void *)chunk, (unsigned int)size) ;
  #endif    
  if (pointers) {
    // If the popped *chunk contains pointer addresses, they should be rebased:
    int i ;
    for (i=size/8-1 ; i>=0 ; --i) {
  #ifdef ADTOOL_TAPENADE   
      ADMM_Rebase(((void **)(chunk))+i, (chunkb?((void **)(*chunkb))+i:NULL)) ;    
  #else
      ADMM_Rebase(((void **)(chunk))+i) ;    
  #endif
    }
  }
  // [Pb:] bizarre that the behavior seems the same if chunkb==NULL and if *chunk==NULL ?
  ADMM_RegisterNewChunkInfo(oldBase, size, (chunk?chunk:0), (chunkb?*chunkb:0)) ;
/* printf("Re-Registered allocated addresses %x(%x) for old base %x size %i\n",(chunk?chunk:0),(chunkb?*chunkb:0),oldBase,size) ; */
  // now try to re-base waiting pointers:
  ADMM_WaitingAddress* waitingRebase ;
  ADMM_Cell* inWaitingRebases = waitingRebaseList ;
/* dumpChunkInfos() ; */
/* dumpWaiting() ; */
/* printf("  Try to apply chunk %x+>%i ==> %x+>? (%x+>?)\n", oldBase, size, (chunk?*chunk:0), (chunkb?*chunkb:0)) ; */
  while (inWaitingRebases->tail) {
    waitingRebase = (ADMM_WaitingAddress*)inWaitingRebases->tail->head ;
/* printf("  to rebase waiting pointers %x->%x (%x->?)\n", waitingRebase->pp, waitingRebase->oldvalpp, waitingRebase->ppb) ; */
    if (waitingRebase->oldvalpp>=(oldBase+size))
      break ; //because inWaitingRebases is sorted, no more rebasing is possible.
    else if (waitingRebase->oldvalpp>=oldBase) {
      // rebase and then remove from waitingRebaseList:
      ADMM_DoRebase(waitingRebase->pp, waitingRebase->oldvalpp, waitingRebase->ppb,
                     oldBase, (chunk?chunk:0), (chunkb?*chunkb:0)) ;
      ADMM_FreeNextWaitingRebase(inWaitingRebases) ;
    } else {
      // keep in waitingRebaseList for later rebase:
      inWaitingRebases = inWaitingRebases->tail ;
    }
  }
  //SHK: This should be ignored by the caller if isDynamic==0
  return chunk;
}





/** TEMPORARY: FOR DEBUGGING. */
void dumpWaiting() {
  ADMM_WaitingAddress* waitingRebase ;
  ADMM_Cell* inWaitingRebases = waitingRebaseList->tail ;
  printf(" Waiting pointers:") ;
  while (inWaitingRebases) {
    waitingRebase = (ADMM_WaitingAddress*)inWaitingRebases->head ;
    if (waitingRebase)
      printf("%x->%x (%x->?) ; ",
             waitingRebase->pp,
             waitingRebase->oldvalpp,
             waitingRebase->ppb) ;
    else
      printf("??? ; ") ;
    inWaitingRebases = inWaitingRebases->tail ;
  }
  printf("end\n") ;
}

void dumpChunkInfos() {
  printf(" Registered chunks: ") ;
  ADMM_ChunkInfo *chunkInfo ;
  ADMM_Cell* inChunkInfos = chunkInfoList->tail ;
  while (inChunkInfos) {
    chunkInfo = (ADMM_ChunkInfo*)inChunkInfos->head ;
    if (chunkInfo)
      printf("%x+>%i ==> %x+>? (%x+>?) ; ",chunkInfo->base, chunkInfo->size, chunkInfo->newbase, chunkInfo->newadjbase) ;
    else
      printf("??? ; ") ;
    inChunkInfos = inChunkInfos->tail ;
  }
  printf("end\n") ;
}
