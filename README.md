# Seismic Modelling with An Effective Absorbing Boundary

Seismic Modelling with An Effective Absorbing Boundary

Author: Jiadong Guo, China University of Petroleum-Beijing, China

Email: [jdongguo@126.com](mailto:jdongguo@126.com)

Programming language: C, Shell, Python

Operating System: Linux

Software dependencies: MPI

Solution method: High-order finite-dfference time-domain (FDTD) for modelling on regular grid

Governing equation: 2nd order acoustic wave equation

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

1. go to /src and compile: cd /src;make
2. configure input parameters in run.sh
3. go to running template and test: bash run.sh
