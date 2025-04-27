#ifndef __SIM_H__
#define __SIM_H__
#include "cstd.h"
typedef struct sim {
    int paded;
    int flag;
    int nz, nx;
    int nbt, nbb, nbl, nbr;
    int nzb, nxb;
    int nt;
    float dt, dx, dz;
    float **v, **vv;     /* velocity */
    float **wt;          /* wavelet */
    int sx, sz;          /* source parameters*/
    int rx, rz, jrx, nr; /* receiver parameters*/
} sim;
sim *sim_init(int nz, int nx, float dx, float dz, int nt,
              float dt); /* 初始化波场模拟工具 */
void sim_set_bnd(sim *s, int nbt, int nbb, int nbl, int nbr); /* 设置扩展边界 */
void sim_set_src(sim *s, int sx, int sz);          /* 设置单个炮点位置 */
void sim_set_rev(sim *s, int rx, int rz, int jrx); /* 设置检波器的分布 */
void sim_pad(sim *s);                              /* 根据边界对网格进行扩展 */
void sim_close(sim *s);
#endif