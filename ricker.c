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
    float fm, dt, t0;
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
    if (!getparfloat("t0", &t0))
        t0 = 1. / fm;
    float *wt = alloc1float(nt);
    float alpha;
    for (int k = 0; k < nt; k++) {
        alpha = (k * dt - t0) * fm * PI;
        alpha = SQUARE(alpha);
        wt[k] = (1 - 2 * alpha) * exp(-alpha);
    }
    FILE *fp = efopen(s, "wb");
    efwrite(wt, sizeof(float), nt, fp);
    efclose(fp);
    free(wt);
}