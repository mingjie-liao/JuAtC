include("AtC.jl")
include("Solve.jl")
include("utils.jl")
include("../Adaptive/adaptive.jl")
# include("FIO.jl")
using SparseArrays
using Isaac
using Printf

Ra, bw, Lmsh = 3, 2, 15
meshpath="/Users/mliao/Program/Mesh/Mesher3DForSJTU/build/mesher3d"

atc = AtC(Ra, bw, Lmsh, 15; meshpath=meshpath);
# TODO: mesh generation is re-built. 
dataPath = joinpath(pathof(JuLIP)[1:end-13], "data/")
calc = JuLIP.Potentials.FinnisSinclair(dataPath*"W-pair-Wang-2014.plt", dataPath*"W-e-dens-Wang-2014.plt");

U = zero(atc.X)
update!(atc, U, Val{:U}())

# Solve 
println("--------------------- Solving -----------------------------")
P = laplace_matrix(atc)

obj_g = x -> nsoli_condition_gradient(atc, x, calc, P)
@time x, it_hist, ierr, x_hist = nsoli(get_x(atc), obj_g; atol=1e-3, rtol=1e-3, debug = 1, lmaxit=30, maxit=30)
@show size(it_hist, 1)

# Estimate
println("--------------------- Estimating -----------------------------")
# ∇U = atc.∇U
∇U, _, J = gradient(atc)
Est = [norm(∇U[:,:,i] - I,2) for i in 1:size(∇U, 3)] # ML: This `I` is necessary, or Est will be ranged betwwen (1.7176407571699979, 1.740646586716024)

# Mark
println("--------------------- Marking -----------------------------")
X, TIdx = Adaptive.mark(atc, Est; Rbuf=Ra)

# Refine
println("--------------------- Refining -----------------------------")
atcnew = Adaptive.refine!(atc, X, TIdx; meshpath=meshpath); # TODO: check if atcnew == atc
