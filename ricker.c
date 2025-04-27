/*
RICKER: Ricker wavelet of central frequency f.
IN
     fm : central freq. in Hz (fm <<1/(2dt) )
     dt: sampling interval in sec
     nt: number of sample
OUT  wt:  the Ricker wavelet
*/
#include "cstd.h"
int main(int argc, char **argv) {
    initargs(argc, argv);
    float fm, dt;
    int nt;
    char *s;
    if (!getparfloat("fm", &fm))
        err("need fm");
    if (!getparfloat("dt", &dt))
        err("need dt");
    if (!getparint("nt", &nt))
        err("need nt");
    if (!getparstring("wt", &s))
        err("need wt=");
    float nw = 2.2 / fm / dt;
    int nnw = 2 * floor(nw / 2) + 1;
    int nc = floor(nnw / 2);
    float *wt = alloc1float(nt);
    float alpha, beta;
    if (nt < nnw)
        err("nt is smaller than condition!");
    for (int k = 0; k < nnw; k++) {
        alpha = (nc - k) * fm * dt * PI;
        beta = SQUARE(alpha);
        wt[k] = (1 - SQUARE(beta)) * exp(-beta);
    }
    for (int k = nnw; k < nt; k++) {
        wt[k] = 0;
    }
    FILE *fp = efopen(s, "wb");
    efwrite(wt, sizeof(float), nt, fp);
    efclose(fp);
}