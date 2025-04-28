echo "flag=1 
snap=0
nz=81 nx=201
dx=5 dz=5
nt=2001 dt=0.0005
sx=100 sz=5 jsx=10 ns=1 
rz=0
vp=vel.dat
wt=wt.dat
nbt=50 nbb=50 nbl=50 nbr=50
wfd=snap.dat" > par.txt

./jd_model $(cat par.txt)
