"""
## module Adaptive

## Description:
Adaptive functionals
	1. mark
	2. refine
""" 
# include("../BGFC/AtC.jl")
module Adaptive
using JuLIP
using LinearAlgebra
# include("../BGFC/AtC.jl")
# include("../BGFC/utils.jl")
# include("../BGFC/FIO.jl")

function mark(atc::Main.AtC{Float64}, Est::Vector{Float64}; Rbuf = 2)
	v = sortperm(Est, rev=true)
	TType = atc.data["TType"]
	Tet = atc.T
	v = v[findall(x -> x ≠ 0.0, TType[v])]
	R = Int((atc.Ra + atc.bw)/rnn(:W)) + Rbuf

	at = bulk(:W, cubic=true)*(R+1)  # Rbuf = 2 by default and expand one layer
	XIdx = Main.get_at_boundary(at)
	X = positions(at) |> mat
	X = X[:, XIdx]
	XType = atc.data["XType"]

	vpop = []
	println("Interface Tetrahedra:")
	for i in v
		# XType[T[:,vrefine[1]]]
		if any(XType[Tet[:,i]].==2.0)
			# println(i)
			append!(vpop, i)
		end
	end

	setdiff!(v, vpop)

	vs = cumsum(Est[v])

	vtol = 0.5 * vs[end]
	idx = findfirst(x -> x >= vtol, vs)
	TIdx = v[1:idx]
	return X, TIdx
end

function refine!(atc::Main.AtC{Float64}, X::Array{Float64,2}, TIdx::Vector{Int64}; filename = "out3d", meshpath="/Users/mliao/Program/Mesh/Mesher3DForSJTU/build/mesher3d", Rbuf = 2)
	mfn = "../FIO/adaptive/$(filename).mesh"
	run(`cp ../FIO/out3d.mesh $mfn`)

	rfn = "../FIO/adaptive/$(filename).remesh"
	Main.ACFIO.write_remesh(rfn, X, TIdx)

	vfn = "../FIO/adaptive/$(filename).value"
	Main.ACFIO.write_value(vfn, atc.U)
	rfn = "../FIO/adaptive/$(filename)"
	run(`$meshpath -r $rfn`)
	ufn = "../FIO/adaptive/$(filename)_out.value"
	U = Main.ACFIO.read_value(ufn)
	# update!(atc, u)
	
	ofn = "../FIO/adaptive/$(filename)_out.mesh"
	X, T = Main.ACFIO.read_mesh(ofn)
	iBdry = findall(x->x==1.0, X[4,:])

	## alternative data
	data = Dict{String, Real}()
	∇U, volT, J = Main.gradient(T, X, U)
	wat = bulk(:W, cubic = true)
	V0=det(cell(wat))/2

	r0 = rnn(:W)
	R = Int((atc.Ra + atc.bw)/r0) + Rbuf
	# RaI = Int(atc.Ra/r0)
	# bwI = Int(atc.bw/r0)
	Ra = atc.Ra + r0;
	bw = atc.bw

	at = bulk(:W, cubic=true)*(R+1)  # Rbuf = 2 by default and expand one layer
	atc = Main.AtC{eltype(X)}(at, V0, X[1:3, :], U, ∇U, T[1:4, :], Ra, bw, J, wat, iBdry, data)
	atc.data["xc"] = [0.0, 0.0, 0.0]
	atc.data["volT"] = copy(volT)
	atc.data["XType"] = copy(X[4,:])
	atc.data["TType"] = copy(T[5,:])
	return atc
end


end # module Adaptive
