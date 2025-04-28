/*
 *有限差分法2-8阶声波方程模拟
 *吸收边界采用有效吸收边界
 */
#include "cstd.h"
void read_data(char *filename, float *v, int n);
void write_data(char *filename, float *v, int n);
void a2d_mod_eal28(float **rcd, sim *s, float ***wfd, int snapint);
int main(int argc, char **argv) {
    int nz, nx, nbt, nbb, nbl, nbr;
    int nt;
    int flag;
    int rz;
    float dt, dx, dz;
    int sx, sz, jsx, ns, sx0;
    float **rcd = NULL, ***wfd = NULL;
    int snap;
    initargs(argc, argv);
    if (!getparint("flag", &flag))
        flag = 1;
    if (!getparint("snap", &snap))
        flag = 0;
    if (!getparint("nz", &nz))
        err("need nz");
    if (!getparint("nx", &nx))
        err("need nx");
    if (!getparint("nbt", &nbt))
        nbt = 30;
    if (!getparint("nbb", &nbb))
        nbb = 30;
    if (!getparint("nbl", &nbl))
        nbl = 30;
    if (!getparint("nbr", &nbr))
        nbr = 30;
    if (!getparint("nt", &nt))
        err("need nt");
    if (!getparint("sx", &sx0))
        err("need sx");
    if (nbt < 4)
        nbt = 4;
    if (nbb < 4)
        nbb = 4;
    if (nbl < 4)
        nbl = 4;
    if (nbr < 4)
        nbr = 4;
    if (!getparint("sz", &sz)) {
        sz = 0;
        warn("sz=0");
    }
    if (!getparint("rz", &rz)) {
        rz = 0;
        warn("rz=0");
    }
    if (!getparint("jsx", &jsx))
        err("need jsx");
    if (!getparint("ns", &ns)) {
        ns = 1;
        warn("ns=1");
    }
    if (!getparfloat("dt", &dt))
        err("need dt");
    if (!getparfloat("dx", &dx))
        err("need dx");
    if (!getparfloat("dz", &dz))
        err("need dz");
    char *vpf, *wtf, *wfdf;
    if (!getparstring("vp", &vpf))
        err("need vp");
    if (!getparstring("wt", &wtf))
        err("need wt");
    if (!getparstring("wfd", &wfdf))
        err("need wfd");
    sim *s = sim_init(nz, nx, dx, dz, nt, dt);
    sim_set_bnd(s, nbt, nbb, nbl, nbr);
    read_data(vpf, s->v[0], nz * nx);
    read_data(wtf, s->wt, nt);
    sim_set_rev(s, 0, rz, 1);
    sim_pad(s);
    s->flag = flag;
    if (snap) {
        wfd = alloc3float(s->nzb, s->nxb, nt);
    } else {
        rcd = alloc2float(nt, nx);
    }
    for (int is = 0; is < ns; is++) {
        sx = sx0 + jsx * is;
        sim_set_src(s, sx, sz);
        if (snap) {
            a2d_mod_eal28(NULL, s, wfd, 1);
        } else {
            a2d_mod_eal28(rcd, s, NULL, 1);
        }
        if (snap) {
            write_data(wfdf, &wfd[0][0][0], nt * s->nzb * s->nxb);
        } else {
            write_data(wfdf, rcd[0], nt * nx);
        }
    }
    sim_close(s);
    if (snap) {
        free3(wfd);
    } else {
        free2(rcd);
    }
}
void read_data(char *filename, float *v, int n) {
    FILE *fp = fopen(filename, "rb");
    fread(v, sizeof(float), n, fp);
    fclose(fp);
}
void write_data(char *filename, float *v, int n) {
    FILE *fp = efopen(filename, "wb");
    efwrite(v, sizeof(float), n, fp);
    efclose(fp);
}