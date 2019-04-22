#include "mystack.h"

// The stack is used by Tapenade to store intermediate results
// in the original computation that are then needed in the
// reverse accumulation of the derivatives.
//
// The original implementation in the Tapenade library seems
// to be something that CIVL does not like, so I implemented a
// simpler version with fixed maximum size for verification
// purposes.
//
// Of course this will mean that we can not verify if the
// stack library is doing the correct thing, I just hope that
// the Tapenade developers got this right (it has been in
// productive use for a decade now, so should be OK).

#define maxSize 10000000
double stack[maxSize];
int stacki[maxSize/100];
int stackc[maxSize/100];
int head = 0;
int headi = 0;
int headc = 0;

void pushreal8array(double *x, int n) {
  for(int i=0; i<n; i++) pushreal8(x[i]);
}

void popreal8array(double *x, int n) {
  for(int i=0; i<n; i++) popreal8(&(x[i]));
}

void pushreal8(double x) {
  stack[head] = x;
  head++;
}

void popreal8(double *x) {
  head--;
  *x = stack[head];
}

void pushinteger4(int x) {
  stacki[headi] = x;
  headi++;
}

void popinteger4(int *x) {
  headi--;
  *x = stacki[headi];
}

void pushcontrol1b(int x) {
  stackc[headc] = x;
  headc++;
}

void popcontrol1b(int *x) {
  headc--;
  *x = stackc[headc];
}
