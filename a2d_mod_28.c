#include "cstd.h"
static float c[] = {-1 / 560, 8 / 315, -1 / 5,  8 / 5,   -205 / 72,
                    8 / 5,    -1 / 5,  8 / 315, -1 / 560};
float **eal2d(sim *s);
float **sponge2d(sim *s);
void a2d_mod_eal28(float **rcd, sim *s, float ***wfd, int snapint) {
    int snap = 0;
    int flag = s->flag;
    int nx = s->nx, nz = s->nz;
    int nxb = s->nxb, nzb = s->nzb;
    float dx = s->dx, dz = s->dz;
    int nbt = s->nbt, nbl = s->nbl;
    int sx = s->sx + nbl, sz = s->sz + nbt;
    int rz = s->rz + s->nbt;
    int nt = s->nt;
    float dt = s->dt;
    float *wt = s->wt;
    float **d = eal2d(s);
    float **v = s->vv;
    float **p0 = alloc2float(nzb, nxb);
    float **p1 = alloc2float(nzb, nxb);
    float **p2 = alloc2float(nzb, nxb);
    for (int ix = 0; ix < nxb; ix++) {
        for (int iz = 0; iz < nzb; iz++) {
            p0[ix][iz] = 0;
            p1[ix][iz] = 0;
            p2[ix][iz] = 0;
        }
    }

    float **tmp, lap_x, lap_z, lap;
    float alpha;
    for (int it = 0; it < nt; it++) {
        if (flag)
            warn("Forward it=%d/%d;", it + 1, nt);
        for (int ix = 4; ix < nxb - 4; ix++) {
            for (int iz = 4; iz < nzb - 4; iz++) {
                alpha = 1 + d[ix][iz] * dt;
                lap_x = 0, lap_z = 0;
                for (int i = -4; i <= 4; i++) {
                    lap_x += c[i + 4] * p1[ix + i][iz];
                    lap_z += c[i + 4] * p1[ix][iz + i];
                }
                lap = lap_x / SQUARE(dx) + lap_z / SQUARE(dz);
                p2[ix][iz] = (2 - SQUARE(d[ix][iz] * dt)) / alpha * p1[ix][iz] +
                             (d[ix][iz] * dt - 1) / alpha * p0[ix][iz] +
                             SQUARE(dt * v[ix][iz]) / alpha * lap;
            }
        }
        tmp = p0, p0 = p1, p1 = p2, p2 = tmp;
        p1[sx][sz] += wt[it];
        if (rcd != NULL) {
            for (int ix = 0; ix < nx; ix++) {
                rcd[ix][it] = p1[ix + nbl][rz];
            }
        }
        if (wfd != NULL && it % snapint == 0) {
            // memcpy(wfd[snap++][0], p1[0], sizeof(float) * nzb * nxb);
            for (int ix = 0; ix < nxb; ix++) {
                for (int iz = 0; iz < nzb; iz++) {
                    wfd[snap][ix][iz] = p1[ix][iz];
                }
            }
        }
        snap++;
    }
    free2(p0);
    free2(p1);
    free2(p2);
    free2(d);
}
void a2d_mod_sponge28(float **rcd, sim *s, float ***wfd, int snapint) {
    int snap = 0;
    int nx = s->nx, nz = s->nz;
    int nxb = s->nxb, nzb = s->nzb;
    float dx = s->dx, dz = s->dz;
    int nbt = s->nbt, nbl = s->nbl;
    int sx = s->sx + nbl, sz = s->sz + nbt;
    int rz = s->rz + s->nbt;
    int nt = s->nt;
    float dt = s->dt;
    float *wt = s->wt;
    float **d = sponge2d(s);
    float **v = s->vv;
    float **p0 = alloc2float(nzb, nxb);
    float **p1 = alloc2float(nzb, nxb);
    float **p2 = alloc2float(nzb, nxb);
    for (int ix = 0; ix < nxb; ix++) {
        for (int iz = 0; iz < nzb; iz++) {
            p0[ix][iz] = 0;
            p1[ix][iz] = 0;
            p2[ix][iz] = 0;
        }
    }

    float **tmp, lap_x, lap_z, lap;
    float alpha;
    for (int it = 0; it < nt; it++) {
        for (int ix = 4; ix < nxb - 4; ix++) {
            for (int iz = 4; iz < nzb - 4; iz++) {
                alpha = 1 + d[ix][iz] * dt;
                lap_x = 0, lap_z = 0;
                for (int i = -4; i <= 4; i++) {
                    lap_x += c[i + 4] * p1[ix + i][iz];
                    lap_z += c[i + 4] * p1[ix][iz + i];
                }
                lap = lap_x / SQUARE(dx) + lap_z / SQUARE(dz);
                p2[ix][iz] =
                    2 * p1[ix][iz] - p0[ix][iz] + SQUARE(v[ix][iz] * dt) * lap;
                p2[ix][iz] *= d[ix][iz];
            }
        }
        tmp = p0, p0 = p1, p1 = p2, p2 = tmp;
        p1[sx][sz] += wt[it];
        if (rcd != NULL) {
            for (int ix = 0; ix < nx; ix++) {
                rcd[ix][it] = p1[ix + nbl][rz];
            }
        }
        if (wfd != NULL && it % snapint == 0) {
            memcpy(wfd[snap++][0], p1[0], sizeof(float) * nzb * nxb);
        }
    }
    free2(p0);
    free2(p1);
    free2(p2);
    free2(d);
}