#include "sim.h"
sim *sim_init(int nz, int nx, float dx, float dz, int nt, float dt)
{
    sim *s = alloc1(1, sizeof(sim));
    s->paded = 0;
    s->nz = nz, s->nx = nx;
    s->v = alloc2float(nz, nx);
    s->wt = alloc1float(nt);
    s->dx = dx, s->dz = dz;
    s->nt = nt, s->dt = dt;
    s->nxb = nx, s->nzb = nz;
    s->nbt = 0, s->nbb = 0, s->nbl = 0, s->nbr = 0;
    return s;
}
void sim_set_bnd(sim *s, int nbt, int nbb, int nbl, int nbr)
{
    s->nbt = nbt;
    s->nbb = nbb;
    s->nbl = nbl;
    s->nbr = nbr;
    s->nzb = s->nz + nbb + nbt;
    s->nxb = s->nx + nbl + nbr;
    if (s->paded)
    {
        free2(s->vv);
        s->paded = 0;
    }
}
void sim_set_src(sim *s, int sx, int sz)
{
    s->sx = sx;
    s->sz = sz;
}
void sim_set_rev(sim *s, int rx, int rz, int jrx)
{
    int nr, nx, k;
    s->rx = rx;
    s->rz = rz;
    s->jrx = jrx;
    nx = s->nx;
    k = (nx - rx) / jrx;
    if ((nx - rx) % jrx == 0)
    {
        s->nr = k;
    }
    else
    {
        s->nr = k + 1;
    }
}
void sim_pad(sim *s)
{
    int nbt = s->nbt, nbl = s->nbl, nbb = s->nbb, nbr = s->nbr;
    int nx = s->nx, nz = s->nz;
    int nxb = s->nxb, nzb = s->nzb;
    if (!s->paded)
    {
        s->vv = alloc2float(nzb, nxb);
        s->paded = 1;
    }
    float **vv = s->vv;
    float **v = s->v;
    /* inner */
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iz = 0; iz < nz; iz++)
        {
            vv[ix + nbl][iz + nbt] = v[ix][iz];
        }
    }

    for (int ix = nbl; ix < nx + nbl; ix++)
    {
        for (int iz = 0; iz < nbt; iz++)
        {
            vv[ix][iz] = vv[ix][nbt];
        }
        for (int iz = nz + nbt; iz < nzb; iz++)
        {
            vv[ix][iz] = vv[ix][nz + nbt - 1];
        }
    }
    for (int iz = 0; iz < nzb; iz++)
    {
        for (int ix = 0; ix < nbl; ix++)
        {
            vv[ix][iz] = vv[nbl][iz];
        }
        for (int ix = nx + nbl; ix < nxb; ix++)
        {
            vv[ix][iz] = vv[nbl + nx - 1][iz];
        }
    }
}
void sim_close(sim *s)
{
    free(s->wt);
    free2(s->v);
    if (s->paded)
    {
        free2(s->vv);
    }
    free(s);
}