void head (double ***u, double ***u_1, double  ***u_2, double ***c, double D, int n){
int i,j,k;

for ( k=1; k<=n - 2; k++ ) {
    for ( j=1; j<=n - 2; j++ ) {
        for ( i=1; i<=n - 2; i++ ) {
            u[i][j][k] += D*(-6*u_1[i][ j][ k] + u_1[i][ j][ k - 1] + u_1[i][ j][ k + 1] + u_1[i][ j - 1][ k] + u_1[i][ j + 1][ k] + u_1[i - 1][ j][ k] + u_1[i + 1][ j][ k])*c[i][ j][ k] + 2.0*u_1[i][ j][ k] - u_2[i][ j][ k];
        }
    }
}
}
