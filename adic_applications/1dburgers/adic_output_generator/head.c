void head (double *u, double *u_1, double C, double D,int n){
int i;
for ( i=1; i<=n - 2; i++ ) {
    u[i] += -C*((-u_1[i] + u_1[i + 1])*((0< u_1[i])?0: u_1[i]) + (u_1[i] - u_1[i - 1])*((0> u_1[i])?0: u_1[i])) + D*(-2.0*u_1[i] + u_1[i - 1] + u_1[i + 1]) + u_1[i];
}
}
