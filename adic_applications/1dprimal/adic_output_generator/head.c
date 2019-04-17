//#define N 10
//void head (double outv[N], double invec[N],double  vel[N], int n){
void head (double *outv, double *invec,double *vel, int n){
int i;
for ( i=1; i<=n - 2; i++ ) {
    outv[i] += (-2.0*invec[i] + invec[i - 1] + invec[i + 1])*vel[i];
}
}
