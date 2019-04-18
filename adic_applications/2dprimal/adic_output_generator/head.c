#define N 10
void head (double **outv, double **inv, double  **vel, int n){
int i,j;

for ( j=1; j<=n - 2; j++ ) {
    for ( i=1; i<=n - 2; i++ ) {
        outv[i][j] += (-4.0*inv[i][j] + inv[i][j - 1] + inv[i][ j + 1] + inv[i - 1][ j] + inv[i + 1][ j])*vel[i][ j];
    }
}
}
