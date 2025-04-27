/*
 *有限差分法2-8阶声波方程模拟
 *吸收边界采用有效吸收边界
 */
#include "cstd.h"
sim *sim_init(int nbt, int nbb, int nbl, int nbr, int nz, int nx, float dx,
              float dz, int nt, float dt, int rz);
void free_sim(sim *s);
void read_data(char *filename, float *v, int n);
void write_data(char *filename, float *v, int n);

void eal2d(sim *s);
void abc2dm(sim *s);
void pad2d(sim *s);
void a2d_mod_eal28(float **rcd, sim *s, float ***wfd, int snapint);
void a2d_mod_sponge28(float **rcd, sim *s, float ***wfd, int snapint);
int main(int argc, char **argv)
{
    int nz, nx, nbt, nbb, nbl, nbr;
    int nt;
    int mode;
    int rz;
    float dt, dx, dz;
    int sx, sz, jsx, ns, sx0;
    initargs(argc, argv);
    if (!getparint("mode", &mode))
        mode = 1;
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
    if (!getparint("sz", &sz))
    {
        sz = 0;
        warn("sz=0");
    }
    if (!getparint("rz", &rz))
    {
        rz = 0;
        warn("sz=0");
    }
    if (!getparint("jsx", &jsx))
        err("need jsx");
    if (!getparint("ns", &ns))
    {
        ns = 1;
        warn("ns=0");
    }
    if (!getparfloat("dt", &dt))
        err("need dt");
    if (!getparfloat("dx", &dx))
        err("need dx");
    if (!getparfloat("dz", &dz))
        err("need dz");
    char *vpf, *wtf;
    if (!getparstring("vp", &vpf))
        err("need vp");
    if (!getparstring("wt", &wtf))
        err("need wt");
    sim *s = sim_init(nbt, nbb, nbl, nbr, nz, nx, dx, dz, nt, dt, rz);
    read_data(vpf, s->v[0], nz * nx);
    read_data(wtf, s->wt, nt);
    pad2d(s);
    float **rcd = alloc2float(nt, nx);
    for (int is = 0; is < ns; is++)
    {
        char wfd[200];
        sprintf(wfd, "wavefield_%d.dat", is + 1);
        s->sx = sx0 + jsx * is;
        s->sz = sz;
        warn("Forward is=%d/%d", is + 1, ns);
        a2d_mod_abc28(rcd, s, NULL, 1);
        write_data(wfd, rcd[0], nt * nx);
    }
    free_sim(s);
}
sim *sim_init(int nbt, int nbb, int nbl, int nbr, int nz, int nx, float dx,
              float dz, int nt, float dt, int rz)
{
    sim *s = alloc1(1, sizeof(sim));
    s->nbt = nbt, s->nbb = nbb, s->nbl = nbl, s->nbr = nbr;
    s->nz = nz, s->nx = nx;
    int nxb = nx + nbl + nbr, nzb = nz + nbt + nbb;
    s->nzb = nzb, s->nxb = nxb;
    s->v = alloc2float(nz, nx);
    s->wt = alloc1float(nt);
    s->dx = dx, s->dz = dz;
    s->nt = nt, s->dt = dt;
    s->rz = rz;
}

void free_sim(sim *s)
{
    free(s->wt);
    free2(s->v);
    free2(s->vv);
    free2(s->d);
    free(s);
}
void read_data(char *filename, float *v, int n)
{
    FILE *fp = fopen(filename, "rb");
    fread(v, sizeof(float), n, fp);
    fclose(fp);
}
void write_data(char *filename, float *v, int n)
{
    FILE *fp = fopen(filename, "wb");
    fwrite(v, sizeof(float), n, fp);
    fclose(fp);
}