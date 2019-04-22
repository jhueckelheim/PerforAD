/*        Generated by TAPENADE     (INRIA, Ecuador team)
    Tapenade 3.14 (r7079) -  5 Oct 2018 09:55
*/
#include <adBuffer.h>

/*
  Differentiation of fmax in reverse (adjoint) mode:
   gradient     of useful results: fmax b
   with respect to varying inputs: b
*/
void fmax_b(double a, double b, double *bb, double fmaxb) {
    double fmax;
    if (a <= b)
        *bb = *bb + fmaxb;
}

double fmax_nodiff(double a, double b) {
    if (a > b)
        return a;
    else
        return b;
}

/*
  Differentiation of fmin in reverse (adjoint) mode:
   gradient     of useful results: fmin b
   with respect to varying inputs: b
*/
void fmin_b(double a, double b, double *bb, double fminb) {
    double fmin;
    if (a >= b)
        *bb = *bb + fminb;
}

double fmin_nodiff(double a, double b) {
    if (a < b)
        return a;
    else
        return b;
}

/*
  Differentiation of burgers1d in reverse (adjoint) mode:
   gradient     of useful results: *u *u_1
   with respect to varying inputs: *u *u_1
   RW status of diff variables: *u:in-out *u_1:incr
   Plus diff mem management of: u:in u_1:in
*/
void burgers1d_b(double *u, double *ub, double *u_1, double *u_1b, double D, 
        double C, int n) {
    int i;
    double result1;
    double result1b;
    double result2;
    double result2b;
    double tempb;
    double tempb0;
#pragma omp parallel for private(i, result1, result2, tempb, tempb0, result1b, result2b)
    for (i = n-2; i > 0; --i) {
        result1 = fmin_nodiff(0, u_1[i]);
        result2 = fmax_nodiff(0, u_1[i]);
        tempb = -(C*ub[i]);
        tempb0 = D*ub[i];
	#pragma omp atomic
        u_1b[i + 1] = u_1b[i + 1] + tempb0 + result1*tempb;
	#pragma omp atomic
        u_1b[i] = u_1b[i] + ub[i] - 2.0*tempb0 + (result2-result1)*tempb;
        result1b = (u_1[i+1]-u_1[i])*tempb;
	#pragma omp atomic
        u_1b[i - 1] = u_1b[i - 1] + tempb0 - result2*tempb;
        result2b = (u_1[i]-u_1[i-1])*tempb;
	tempb = 0;
        fmax_b(0, u_1[i], &(tempb), result2b);
	#pragma omp atomic
	u_1b[i] += tempb;
	tempb = 0;
        fmin_b(0, u_1[i], &(tempb), result1b);
	#pragma omp atomic
	u_1b[i] += tempb;
    }
}