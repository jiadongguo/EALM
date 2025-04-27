/* generate the absorption boundary coefficient */
#include "cstd.h"
static float findMax(float *p, int n)
{
    float ret = p[0];
    for (int i = 0; i < n; i++)
    {
        ret = MAX(ret, p[i]);
    }
    return ret;
}
static float findMin(float *p, int n)
{
    float ret = p[0];
    for (int i = 0; i < n; i++)
    {
        ret = MIN(ret, p[i]);
    }
    return ret;
}
float **eal2d(sim *s)
{
    int nbt = s->nbt, nbl = s->nbl, nbb = s->nbb, nbr = s->nbr;
    int nx = s->nx, nz = s->nz;
    int nxb = s->nxb, nzb = s->nzb;
    float dx = s->dx, dz = s->dz;
    float **d = alloc2float(nzb, nxb);
    float **v = s->vv;
    float ref = 1e-6, alpha = 1.5 * log(1 / ref);
    float width, dis;
    float vmax = findMax(v[0], nzb * nxb);
    /* inner */
    for (int ix = nbl; ix < nbl + nx; ix++)
    {
        for (int iz = nbt; iz < nbt + nz; iz++)
        {
            d[ix][iz] = 0;
        }
    }
    /* top */
    width = nbt * dz;
    for (int iz = 0; iz < nbt; iz++)
    {
        dis = (nbt - iz) * dz;
        for (int ix = nbl; ix < nbl + nx; ix++)
        {
            d[ix][iz] = alpha / width * SQUARE(dis / width) * v[ix][iz];
        }
    }
    /*bottom*/
    width = nbb * dz;
    for (int iz = nbt + nz; iz < nzb; iz++)
    {
        dis = (iz + 1 - nbt - nz) * dz;
        for (int ix = nbl; ix < nbl + nx; ix++)
        {
            d[ix][iz] = alpha / width * SQUARE(dis / width) * v[ix][iz];
        }
    }
    /* left */
    width = nbl * dx;
    for (int ix = 0; ix < nbl; ix++)
    {
        dis = (nbl - ix) * dx;
        for (int iz = nbt; iz < nbt + nz; iz++)
        {
            d[ix][iz] = alpha / width * SQUARE(dis / width) * v[ix][iz];
        }
    }
    /* right */
    width = nbr * dx;
    for (int ix = nbl + nx; ix < nxb; ix++)
    {
        dis = (ix + 1 - nbl - nx) * dx;
        for (int iz = nbt; iz < nbt + nz; iz++)
        {
            d[ix][iz] = alpha / width * SQUARE(dis / width) * v[ix][iz];
        }
    }
    /* left-top */
    width = hypot(nbt * dz, nbl * dx);
    for (int ix = 0; ix < nbl; ix++)
    {
        for (int iz = 0; iz < nbt; iz++)
        {
            dis = hypot((nbl - ix) * dx, (nbt - iz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * v[ix][iz];
        }
    }
    /* left-bottom */
    width = hypot(nbb * dz, nbl * dx);
    for (int ix = 0; ix < nbl; ix++)
    {
        for (int iz = nbt + nz; iz < nzb; iz++)
        {
            dis = hypot((nbl - ix) * dx, (iz + 1 - nbt - nz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * v[ix][iz];
        }
    }
    /* right-top */
    width = hypot(nbt * dz, nbr * dx);
    for (int ix = nbl + nx; ix < nxb; ix++)
    {
        for (int iz = 0; iz < nbt; iz++)
        {
            dis = hypot((ix + 1 - nbl - nx) * dx, (nbt - iz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * v[ix][iz];
        }
    }
    /* right-bottom */
    width = hypot(nbb * dz, nbr * dx);
    for (int ix = nbl + nx; ix < nxb; ix++)
    {
        for (int iz = nbt + nz; iz < nzb; iz++)
        {
            dis = hypot((ix + 1 - nbl - nx) * dx, (iz + 1 - nbt - nz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * v[ix][iz];
        }
    }
    return d;
}
void sponge2d(sim *s)
{
    float alpha = 0.15;
    int nbt = s->nbt, nbl = s->nbl, nbb = s->nbb, nbr = s->nbr;
    int nx = s->nx, nz = s->nz;
    int nxb = s->nxb, nzb = s->nzb;
    float dx = s->dx, dz = s->dz;
    float **d = alloc2float(nzb, nxb);
    /* inner */
    for (int ix = 0; ix < nxb; ix++)
    {
        for (int iz = 0; iz < nzb; iz++)
        {
            d[ix][iz] = 1;
        }
    }
    /* top */
    for (int ix = nbl; ix < nx + nbl; ix++)
    {
        for (int iz = 0; iz < nbt; iz++)
        {
            d[ix][iz] = exp(-SQUARE(alpha * (nbt - iz)));
        }
    }
    /* bottom */
    for (int ix = nbl; ix < nx + nbl; ix++)
    {
        for (int iz = nbt + nz; iz < nzb; iz++)
        {
            d[ix][iz] = exp(-SQUARE(alpha * (iz + 1 - nz - nbt)));
        }
    }
    /* left */
    for (int ix = 0; ix < nbl; ix++)
    {
        for (int iz = nbt; iz < nbt + nz; iz++)
        {
            d[ix][iz] = exp(-SQUARE(alpha * (nbl - ix)));
        }
    }
    /* right */
    for (int ix = nbl + nx; ix < nxb; ix++)
    {
        for (int iz = nbt; iz < nbt + nz; iz++)
        {
            d[ix][iz] = exp(-SQUARE(alpha * (ix + 1 - nbl - nx)));
        }
    }
    /* left-top */

    for (int ix = 0; ix < nbl; ix++)
    {
        for (int iz = 0; iz < nbt; iz++)
        {
            d[ix][iz] = exp(-SQUARE(alpha * hypot(nbl - ix, nbt - iz)));
        }
    }
    /* left-bottom */

    for (int ix = 0; ix < nbl; ix++)
    {
        for (int iz = nbt + nz; iz < nzb; iz++)
        {

            d[ix][iz] = exp(-SQUARE(alpha * hypot(nbl - ix, iz + 1 - nbt - nz)));
        }
    }
    /* right-top */

    for (int ix = nbl + nx; ix < nxb; ix++)
    {
        for (int iz = 0; iz < nbt; iz++)
        {

            d[ix][iz] = exp(-SQUARE(alpha * hypot(ix + 1 - nx - nbl, nbt - iz)));
        }
    }
    /* right-bottom */

    for (int ix = nbl + nx; ix < nxb; ix++)
    {
        for (int iz = nbt + nz; iz < nzb; iz++)
        {
            d[ix][iz] = exp(-SQUARE(alpha * hypot(ix + 1 - nx - nbl, iz + 1 - nz - nbt)));
        }
    }
    return d;
}
void abc2dm(sim *s)
{
    int nbt = s->nbt, nbl = s->nbl, nbb = s->nbb, nbr = s->nbr;
    int nx = s->nx, nz = s->nz;
    int nxb = s->nxb, nzb = s->nzb;
    float dx = s->dx, dz = s->dz;
    float **d = alloc2float(nzb, nxb);
    float **v = s->vv;
    float ref = 1e-6, alpha = 1.5 * log(1 / ref);
    float width, dis;
    float vmax = findMax(v[0], nzb * nxb);
    /* inner */
    for (int ix = nbl; ix < nbl + nx; ix++)
    {
        for (int iz = nbt; iz < nbt + nz; iz++)
        {
            d[ix][iz] = 0;
        }
    }
    /* top */
    width = nbt * dz;
    for (int iz = 0; iz < nbt; iz++)
    {
        dis = (nbt - iz) * dz;
        for (int ix = nbl; ix < nbl + nx; ix++)
        {
            d[ix][iz] = alpha / width * SQUARE(dis / width) * vmax;
        }
    }
    /*bottom*/
    width = nbb * dz;
    for (int iz = nbt + nz; iz < nzb; iz++)
    {
        dis = (iz + 1 - nbt - nz) * dz;
        for (int ix = nbl; ix < nbl + nx; ix++)
        {
            d[ix][iz] = alpha / width * SQUARE(dis / width) * vmax;
        }
    }
    /* left */
    width = nbl * dx;
    for (int ix = 0; ix < nbl; ix++)
    {
        dis = (nbl - ix) * dx;
        for (int iz = nbt; iz < nbt + nz; iz++)
        {
            d[ix][iz] = alpha / width * SQUARE(dis / width) * vmax;
        }
    }
    /* right */
    width = nbr * dx;
    for (int ix = nbl + nx; ix < nxb; ix++)
    {
        dis = (ix + 1 - nbl - nx) * dx;
        for (int iz = nbt; iz < nbt + nz; iz++)
        {
            d[ix][iz] = alpha / width * SQUARE(dis / width) * vmax;
        }
    }
    /* left-top */
    width = hypot(nbt * dz, nbl * dx);
    for (int ix = 0; ix < nbl; ix++)
    {
        for (int iz = 0; iz < nbt; iz++)
        {
            dis = hypot((nbl - ix) * dx, (nbt - iz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * vmax;
        }
    }
    /* left-bottom */
    width = hypot(nbb * dz, nbl * dx);
    for (int ix = 0; ix < nbl; ix++)
    {
        for (int iz = nbt + nz; iz < nzb; iz++)
        {
            dis = hypot((nbl - ix) * dx, (iz + 1 - nbt - nz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * vmax;
        }
    }
    /* right-top */
    width = hypot(nbt * dz, nbr * dx);
    for (int ix = nbl + nx; ix < nxb; ix++)
    {
        for (int iz = 0; iz < nbt; iz++)
        {
            dis = hypot((ix + 1 - nbl - nx) * dx, (nbt - iz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * vmax;
        }
    }
    /* right-bottom */
    width = hypot(nbb * dz, nbr * dx);
    for (int ix = nbl + nx; ix < nxb; ix++)
    {
        for (int iz = nbt + nz; iz < nzb; iz++)
        {
            dis = hypot((ix + 1 - nbl - nx) * dx, (iz + 1 - nbt - nz) * dz);
            d[ix][iz] = alpha / width * SQUARE(dis / width) * vmax;
        }
    }
    return d;
}