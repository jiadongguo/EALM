# Seismic Modelling with An Effective Absorbing Boundary

![](https://img.shields.io/badge/License-GPLv3-blue)  ![](https://img.shields.io/badge/Author-Jiadong_Guo-blue)  ![](https://img.shields.io/badge/Email-jdongguo@126.com-blue)  ![](https://img.shields.io/badge/Language-C_Shell_Python-blue)  ![](https://img.shields.io/badge/System-Linux-blue)  ![](https://img.shields.io/badge/Dependencies-MPI_OpenBlas-blue)


Solution method: **High-order finite-dfference time-domain (FDTD) for modelling on regular grid**

Governing equation: **2nd order acoustic wave equation**

## Credit

- Yao G, Da Silva N V, Wu D. An effective absorbing layer for the boundary condition in acoustic seismic wave simulation[J]. Journal of Geophysics and Engineering, 2018, 15(2): 495–511.

- 张平民, 姚刚. 基于改进的 EAL 边界条件的地震波场数值模拟[J].

## Code structure

- src: the source code in .c

- include: the header files in .h

- doc: documents for theoretic background

- bin: the folder to store executable after compilation

- layer_model: a quick modeling example in layered medium

## Instructions to run

1. gcc main_model.c cstd.c pad2d.c a2d_mod_28.c -o main_model -lm
2. go to running template and test: bash run.sh
