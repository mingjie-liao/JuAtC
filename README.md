# Julia-based Atomistic-to-Continuum coupling

This is a implementations of BGFC (Blended Ghost Force Correction) method for atomistic-to-continuum coupling in 3D.

## Prerequisites

The 3D mesh toolkit is supported by Mr. Kejie Fu, a PhD candidate of Zhejiang University. See [Mesher3D](https://github.com/kjfu/Mesher3DForSJTU.git/) for more details.

The atomistic computations involved are heavily depends on the pure Julia package [JuLIP](https://github.com/JuliaMolSim/JuLIP.jl) (Julia Library for Interatomic Potentials).

Julia Packages involved: JuLIP, DelimitedFiles, Printf, NeighbourLists, QHull, Optim, LineSearches, SparseArrays, nsoli.

## FIO

The module named "ACFIO" is used to cope with the geometrical operations.

The functionals of FIO are to read/write:   

- .mesh
- .remesh   
- .value   
and write .dump files for visualization.

## AtC

Major struct contains geometrical and computational information.

To construct an atc objective invokes:   
```
function AtC(Ra::Int64, bw::Int64, Lmsh, h; Rbuf=2, sp=:W, r0=rnn(:W), defects=:SingVac, meshpath="/home/mliao/Program/Mesh/Mesher3DForSJTU/build/mesher3d")
```
Note that `meshpath` is the path of mesher toolkit that may differ from each devices.

##Solve
we use `minimise!` which mimics `JuLIP.Solve` and also `nsoli`.