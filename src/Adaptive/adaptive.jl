"""
## module Adaptive

## Description:
Adaptive functionals
	1. mark
	2. refine
""" 
module Adaptive
using JuLIP
using LinearAlgebra

using JuAtC
import JuAtC: ACFIO
import JuAtC: AtC, get_at_boundary

function mark(atc::AtC{Float64}, Est::Vector{Float64}; Rbuf = 2)
	v = sortperm(Est, rev=true)
	TType = atc.data["TType"]
	Tet = atc.T
	v = v[findall(x -> x ≠ 0.0, TType[v])]
	R = Int((atc.Ra + atc.bw)/rnn(:W)) + Rbuf

	at = bulk(:W, cubic=true)*(R+1)  # Rbuf = 2 by default and expand one layer
	XIdx = get_at_boundary(at)
	X = positions(at) |> mat
	X = X[:, XIdx]
	XType = atc.data["XType"]

	vpop = []
	println("Interface Tetrahedra:")
	for i in v
		if any(XType[Tet[:,i]].==2.0)
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

function refine!(atc::AtC{Float64}, X::Array{Float64,2}, TIdx::Vector{Int64}; filename = "out3d", meshpath="/Users/mliao/Program/Mesh/Mesher3DForSJTU/build/mesher3d", Rbuf = 2)
	FPath = joinpath(pathof(JuAtC)[1:end-8], "FIO/")
	mfn = FPath*"adaptive/$(filename).mesh"
	run(`cp $(FPath)out3d.mesh $mfn`)

	rfn = FPath*"adaptive/$(filename).remesh"
	ACFIO.write_remesh(rfn, X, TIdx)

	vfn = FPath*"adaptive/$(filename).value"
	ACFIO.write_value(vfn, atc.U)
	rfn = FPath*"adaptive/$(filename)"
	run(`$meshpath -r $rfn`)
	ufn = FPath*"adaptive/$(filename)_out.value"
	U = ACFIO.read_value(ufn)
	# update!(atc, u)
	
	ofn = FPath*"adaptive/$(filename)_out.mesh"
	X, T = ACFIO.read_mesh(ofn)
	iBdry = findall(x->x==1.0, X[4,:])

	## alternative data
	data = Dict{String, Real}()
	∇U, volT, J = gradient(T, X, U)
	wat = bulk(:W, cubic = true)
	V0=det(cell(wat))/2

	r0 = rnn(:W)
	R = Int((atc.Ra + atc.bw)/r0) + Rbuf
	# RaI = Int(atc.Ra/r0)
	# bwI = Int(atc.bw/r0)
	Ra = atc.Ra + r0;
	bw = atc.bw

	at = bulk(:W, cubic=true)*(R+1)  # Rbuf = 2 by default and expand one layer
	atc = AtC{eltype(X)}(at, V0, X[1:3, :], U, ∇U, T[1:4, :], Ra, bw, J, wat, iBdry, atc.calc, data)
	atc.data["xc"] = [0.0, 0.0, 0.0]
	atc.data["volT"] = copy(volT)
	atc.data["XType"] = copy(X[4,:])
	atc.data["TType"] = copy(T[5,:])
	return atc
end


end # module Adaptive
