#include "cstd.h"
void pad2d(sim *s)
{
    int nbt = s->nbt, nbl = s->nbl, nbb = s->nbb, nbr = s->nbr;
    int nx = s->nx, nz = s->nz;
    int nxb = s->nxb, nzb = s->nzb;
    float **vv = alloc2float(nzb, nxb);
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
    s->vv = vv;
}