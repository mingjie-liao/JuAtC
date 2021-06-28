using JuLIP
using NeighbourLists
using Printf
using QHull
using LinearAlgebra

import JuLIP: energy, forces, gradient

export AtC, update!

abstract type AbstractAtC{T} end

mutable struct AtC{T} <: AbstractAtC{T} #where T # <: AbstractFloat
	at::AbstractAtoms{T}
	V0::T
	X::Array{T, 2}			# nodes' positions
	U::Array{T, 2}			# nodes' displacements
	∇U::Array{T, 3} 
	T::Array{Int64, 2}		# elements
	Ra::Float64				# radius of atomistic region
	bw::Float64				# additional width of blending region
	J::Array{T}				# Jacobian?
	Wat::AbstractAtoms{T}			# unitCell with the same type of atomistic
	iBdry::Vector{Int64}
	calc::AbstractCalculator
	data::Dict{String, Array{T}}
end 

function AtC(Ra::Int64, bw::Int64, Lmsh, h, calc; 
						Rbuf=2, sp=:W, r0=rnn(:W), defects=:SingVac, meshpath="/home/mliao/Program/Mesh/Mesher3DForSJTU/build/mesher3d")
	## construct atomistic region
	at = bulk(sp, cubic=true)*(Ra+bw+Rbuf)
	set_pbc!(at, false)
	# locate the 'center' of atomistic lattice
	Xat = positions(at)
	xcell = diag(cell(at))/2
	xcidx = findall(x->isapprox(x, xcell), Xat)
	if isnothing(xcidx)
		r = [ norm(x - xcell) for x in Xat ]
		_, xcidx = findmin(r)
		Xcidx = xcidx.I[2]
		@assert !isnothing(xcidx)
	end
	xc = Xat[xcidx]
	Xat .-= xc
	set_positions!(at, Xat)

	# defect case, need to be optimized
	if defects == :SingVac
		deleteat!(at, xcidx)
		Xat = positions(at)
	elseif defects == :MultiVac
		error("Haven't implemented yet!")
	elseif defects == :MicroCrack
		error("Haven't implemented yet!")
	elseif defects == :Disloc
		error("Haven't implemented yet!")
	end

	iBdry_at = get_at_boundary(at; xc=[0.0, 0.0, 0.0])

	## construct continuum information
	Ra *= r0
	bw *= r0
	Rbuf *= r0
	Lmsh *= r0
	cb = [-1 1 1 -1 -1 1 1 -1; -1 -1 1 1 -1 -1 1 1; -1 -1 -1 -1 1 1 1 1]
	Xcb = Lmsh .* cb
	Xcb = Xcb
	X = hcat(mat(Xat), Xcb)
	Xtype = zeros(Int64, length(at))
	Xtype[iBdry_at] .= 2
	append!(Xtype, ones(Int64, 8))  # TODO: 1 or 3?
	FPath = joinpath(pathof(JuAtC)[1:end-8], "FIO/")
	fn = FPath*"cp.mesh"
	ACFIO.write_mesh(fn, X, Xtype)
	# call `mesher3d` to build coupled mesh
	ofn = FPath*"out3d.mesh"
	run(`$meshpath -s $h $fn -o $ofn`)
	X, T = ACFIO.read_mesh(ofn)
	iBdry = findall(x->x==1.0, X[4,:])

	## alternative data
	data = Dict{String, Real}()
	U = zeros(3, size(X,2))
	∇U, volT, J = gradient(T, X, U)
	wat = bulk(:W, cubic = true)
	V0=det(cell(wat))/2

	atc = AtC{eltype(X)}(at, V0, X[1:3, :], U, ∇U, T[1:4, :], Ra, bw, J, wat, iBdry, calc, data)
	atc.data["xc"] = [0.0, 0.0, 0.0]
	atc.data["volT"] = copy(volT)
	atc.data["XType"] = copy(X[4,:])
	atc.data["TType"] = copy(T[5,:])
	return atc
end

gradient(config::AtC{T};domain=1:size(config.T,2)) where T = gradient(config.T, config.X, config.U;domain=domain)

"""
gradient: (T::Tetrahedron, X::positions, U::displacements) ⟶ (∇U::gradients, volT::volumes)   
"""
function gradient(T::Array{Int64,2}, X::Array{Float64,2}, U::Array{Float64,2}; domain=1:size(T,2))
	nT = size(T,2)
	J = zeros(Float64, 3, 3)
	Du = zeros(Float64, 3, 3)
	∇U = zeros(Float64, 3, 3, nT)
	Jc = zeros(Float64, 3, 3, nT)
	VolT = zeros(Float64, nT)
	for i in domain
		t = T[1:4,i]
		xx = X[1:3,t]
		uu = U[1:3,t]

		for j = 1:3
			J[:,j] = xx[:,j+1] - xx[:,1]
			Du[:,j] = uu[:,j+1] - uu[:,1]
		end

		∇U[:,:,i] = (Du + J) / J # check Du or Du+J would be suitable for minimise!
		Jc[:,:,i] = J
		VolT[i] = abs(det(J)/6)
	end
	return ∇U, VolT, Jc
end

# forces(config::AtC{Float64}, U::Array{Float64,2}, calc::AbstractCalculator) = 
		# forces(update!(config, U, Val{:U}()), calc, Val{:BGFC}())

function gradient(config::AtC{Float64}, U::Array{Float64,1})
	if "xfree" in keys(config.data)
		xfree = convert(Array{Int64,1}, config.data["xfree"])
	else 
		xfree = get_free!(config)
	end
	Ux = zero(config.U)
	Ux[xfree] = U
	Frc = forces(update!(config, Ux, Val{:U}()), Val{:BGFC}())
	return rmul!(Frc[xfree], -1.0)
end

"""
Wcb (Cauchy-Born energy density): (∇U::gradients, Wat::JuLIP.Atoms, Tidx::Tetrahedron index) ⟶ (W::energy density, dW::wrt F)
"""
function Wcb(∇U::Array{Float64}, Wat::Atoms{Float64}, Tidx::Int64, J::Array{Float64}, calc::AbstractCalculator; V0=det(cell(Wat))/2)
	∇u = ∇U[:,:,Tidx]

	at = deepcopy(Wat)
	if isnothing(calculator(at))
		set_calculator!(at, calc)
	end
	nlist = neighbourlist(at, cutoff(at))
	tmp = JuLIP.alloc_temp_d(calc, at) # tmp.R and tmp.dV
	_, R = neigs!(tmp.R, nlist, 1) # fillin tmp.R
	Rref = mat(tmp.R)[:,:]
	W0 = Potentials.evaluate!(tmp, calc, R)# /V0

	X0 = positions(at)
	C0 = cell(at)

	apply_defm!(at, ∇u)
	Rdef = ∇u * Rref
	Ridx = findall(x->x≤cutoff(at), vec(mapslices(x->norm(x,2), Rdef, dims=1)))
	Rref = Rref[:,Ridx]

	nlist = neighbourlist(at, cutoff(at))
	tmp = JuLIP.alloc_temp_d(calc, at) # tmp.R and tmp.dV #ML: will this be a problem?
	_, R = neigs!(tmp.R, nlist, 1)

	W = Potentials.evaluate!(tmp, calc, R) - W0
	Potentials.evaluate_d!(tmp.dV, tmp, calc, R) # fillin tmp.dV 
	dV = mat(tmp.dV)[:,:]

	@assert size(dV,2) == size(Rref,2)
	dW = dV*adjoint(Rref)/adjoint(J)
	return W, dW
end

"""
continuum energy: (W::energy density, volT::volumes) ⟶ (Ecs::continuum elements' energy)
"""
# function energy(W::Array{Float64,1}, volT::{Float64,1}) ML: seems useless...
function energy(config::AtC{Float64}, model::Val{:continuum})
	if "W" in keys(config.data)
		W = config.data["W"]
		volT = config.data["volT"]
		return sum(volT .*W)
	end
	∇U, _, J = gradient(config)
	nT = size(config.T, 2)
	W = zeros(Float64, nT, 1)
	dW = zeros(Float64, 3, 3, nT)
	volT = config.data["volT"]
	for i = 1:nT
		if isapprox(volT[i], 0, atol=1e-10)
			continue
		end
		w, dw = Wcb(∇U, config.Wat, i, J[:,:,i], config.calc)
		W[i] = w
		dW[:,:,i] = dw
	end
	config.data["W"] = W
	config.data["dW"] = dW
	return sum(volT.*W)
end

"""
continuum forces: 
"""
function forces(config::AtC{Float64}, model::Val{:continuum})
	Frc = zeros(Float64, 3, size(config.X,2))
	nT = size(config.T,2)
	if "dW" in keys(config.data)
		dW = config.data["dW"]
		volT = config.data["volT"]
		Tet = config.T
		for j = 1:nT
			t = Tet[1:4,j]
			dw = dW[:,:,j]
			for n = 1:length(t)-1
				Frc[:,t[n+1]] -= volT[j] * dw[:,n]
				Frc[:,t[1]] += volT[j] *dw[:,n]
			end
		end
		Frc[:, config.iBdry] .= 0
		return Frc
	end
	∇U, _, J = gradient(config)
	volT = config.data["volT"]
	W = zeros(Float64, nT, 1)
	dW = zeros(Float64, 3, 3, nT)
	Tet = config.T
	for i = 1:nT
		if isapprox(volT[i], 0, atol=1e-10)
			continue
		end
		t = Tet[1:4,i]
		w, dw = Wcb(∇U, config.Wat, i, J[:,:,i], config.calc)
		W[i] = w
		dW[:,:,i] = dw
		for n = 1:length(t)-1
			Frc[:,t[n+1]] -= volT[i] *dw[:,n]
			Frc[:,t[1]] += volT[i] *dw[:,n]
		end
	end
	config.data["W"] = W
	config.data["dW"] = dW
	Frc[:, config.iBdry] .= 0
	return Frc
end

function bqce_prep_config!(config::AtC{Float64}; bfcn=:affine)
	X = config.X
	xc = config.data["xc"]
	r = [norm(X[:,i] - xc,2) for i=1:size(X,2)]
	β = blendfcn(r, config.Ra, config.Ra+config.bw, bfcn)

	if "volT" in keys(config.data)
		volT = Vector{Float64}(undef, size(config.T,2))
		copyto!(volT, config.data["volT"])
	else
		_, volT, _ = gradient(config)
	end
	Tet = config.T
	@assert β[config.iBdry[1]] == 1.0
	for i = 1:size(Tet,2)
		t = Tet[1:4,i]
		if all(x->x==1.0, β[t])
			continue
		elseif all(x->x==0.0, β[t])
			volT[i] = 0.0
		end
		xt = X[1:3,t]
		mpt, mpIdx = Tet_midpoint(xt)
		tfcc, tfccIdx = Tet_circumcenter_face(xt)
		tc = Tet_circumcenter(xt)
		for n = 1:4
			p = [xt[:,n] mpt[:,mpIdx[n]] tfcc[:,tfccIdx[n]] tc]
			ph = chull(p'[:,:])
			volT[i] -= β[t[n]]*ph.volume
		end
		volT[i] = max(volT[i], 0.0)
	end
	config.data["volT"] = volT ./ config.V0
	config.data["β"] = 1 .- β
	return 1 .- β
end

"""
Atomistic forces for blending schemds
"""
function energy(at::AbstractAtoms{T}, calc::AbstractCalculator, β::Vector{T}; domain=1:length(at)) where {T}
	tmp = JuLIP.alloc_temp(calc, at)
	Ea = zero(T)
	nlist = neighbourlist(at, cutoff(calc))
	for i in domain
		_, R = neigs!(tmp.R, nlist, i)
		Ea += β[i] * JuLIP.Potentials.evaluate!(tmp, calc, R)
	end
	return Ea
end

"""
Atomistic forces for blending schemes, modified from `JuLIP.forces!`
"""
function forces(at::AbstractAtoms{T}, calc::AbstractCalculator, β::Vector{T}; domain=1:length(at)) where {T}#, reset=true)
	frc = zeros(JVec{T}, length(at))
	tmp = JuLIP.alloc_temp_d(calc, at) 
	nlist = neighbourlist(at, cutoff(calc))
	for i in domain
		j, R = neigs!(tmp.R, nlist, i)
		JuLIP.Potentials.evaluate_d!(tmp.dV, tmp, calc, R)
		for a = 1:length(j)
			frc[j[a]] -= β[i]*tmp.dV[a]
			frc[i]    += β[i]*tmp.dV[a]
		end
	end
	return frc
end

"""
BQCE energy:
"""
function energy(config::AtC{Float64}, model::Val{:BQCE}; bfcn=:affine)
	if "β" in keys(config.data)
		β = config.data["β"]
		volT = config.data["volT"]
	else
		β = bqce_prep_config!(config; bfcn=bfcn)
	end
	@assert β[config.iBdry[1]] == 0.0
	Ea = energy(config.at, config.calc, β)
	Ec = energy(config, Val{:continuum}()) 
	return Ea + Ec
end

"""
BQCE forces:
"""
function forces(config::AtC{Float64}, model::Val{:BQCE}; bfcn=:affine)
	if "β" in keys(config.data)
		β = config.data["β"]
	else
		β = bqce_prep_config!(config; bfcn=bfcn)
	end
	@assert β[config.iBdry[1]] == 0.0
	Fa = zero(config.X)
	Fa[:, 1:length(config.at)] = forces(config.at, config.calc, β) |> mat
	Fc = forces(config, Val{:continuum}())
	return Fc + Fa
end

function bgfc_prep_config!(config::AtC{Float64}; bfcn=:affine)
	U = deepcopy(config.U)
	update!(config, zero(config.U), Val{:U}())
	E0 = energy(config, Val{:BQCE}(); bfcn=bfcn)
	β = config.data["β"]
	@assert β[config.iBdry[1]] == 0.0

	F0 = forces(config, Val{:BQCE}(); bfcn=bfcn)

	xc = config.data["xc"]
	r = [norm(config.X[:,i] - xc,2) for i=1:size(config.X,2)]
	IBC = findall(x-> x < max(config.Ra-2*config.bw, config.bw), r)
	F0[:, IBC] .= 0
	update!(config, U, Val{:U}())
	config.data["E0"] = [E0]
	config.data["F0"] = F0
	return E0, F0
end

# energy(config::AtC{Float64}, U::Array{Float64,2}, calc::AbstractCalculator) = 
# 		energy(update!(config, U, Val{:U}()), calc, Val{:BGFC}())
function energy(config::AtC{Float64}, U::Array{Float64,1})
	if "xfree" in keys(config.data)
		xfree = convert(Array{Int64,1}, config.data["xfree"])
	else 
		xfree = get_free!(config)
	end
	Ux = zero(config.U)
	Ux[xfree] = U
	E = energy(update!(config, Ux, Val{:U}()), Val{:BGFC}())
	return E
end

"""
BGFC energy
"""
function energy(config::AtC{Float64}, model::Val{:BGFC}; bfcn=:affine)
	if "E0" in keys(config.data)
		β = config.data["β"]
		Eqce = config.data["E0"][1]
		Fqce = config.data["F0"]
	else
		Eqce, Fqce = bgfc_prep_config!(config; bfcn=bfcn)
		β = config.data["β"]
	end

	@assert β[config.iBdry[1]] == 0.0
	Ea = energy(config.at, config.calc, β)
	Ec = energy(config, Val{:continuum}()) 
	return Ea + Ec + sum(Fqce .* config.U)
end

# forces(config::AtC{Float64}, U::Array{Float64,2}, calc::AbstractCalculator) = 
# 		forces(update!(config, U, Val{:U}()), calc, Val{:BGFC}())

"""
BGFC forces:
"""
function forces(config::AtC{Float64}, model::Val{:BGFC}; bfcn=:affine)
	if "F0" in keys(config.data)
		β = config.data["β"]
		Eqce = config.data["E0"]
		Fqce = config.data["F0"]
	else
		Eqce, Fqce = bgfc_prep_config!(config; bfcn=bfcn)
		β = config.data["β"]
	end
	@assert β[config.iBdry[1]] == 0.0

	Fa = zero(config.X)
	Fa[:, 1:length(config.at)] = forces(config.at, config.calc, β) |> mat
	Fc = forces(config, Val{:continuum}())
	return Fc + Fa - Fqce
end

"""
blending function: (r::radii, R1::inner, R2::outer, flag::type) ⟶ (β::weights on lattices)
"""
#TODO: update to optimal (s3) blending function
function blendfcn(r::Array{Float64,1}, R1::Float64, R2::Float64, flag::Symbol)
    β = zeros(length(r))
    if flag == :none
		return β
	elseif flag == :const
		fill!(β, 1.0)
    elseif flag == :shock
        β[findall(r .> R2)] .= 1.0
    elseif flag == :affine
        β[findall(r .> R2)] .= 1.0
        ib = findall(R1 .< r .< R2)
        β[ib] = (r[ib].-R1)/(R2-R1)
    elseif flag == :cos
        β[findall(r .> R2)] .= 1.0
        ib = findall(R1 .< r .< R2)
		β[ib] = cos.(2*(r[ib].-R1)/(π*(R2-R1)))
    else
        @info(`Not implemented yet.`)   
    end
    return β
end

function update!(config::AtC{Float64}, U::Array{Float64,1})
	if "xfree" in keys(config.data)
		xfree = convert(Array{Int64,1}, config.data["xfree"])
	else 
		xfree = get_free!(config)
	end
	Ux = zero(config.U)
	Ux[xfree] = U
	update!(config, Ux, Val{:U}())
end

function update!(config::AtC{Float64}, U::Array{Float64,2}, attr::Val{:U})
	@assert size(config.U) == size(U)
	Uprv = config.U
	dU = Uprv - U
	iV = findall(x->abs(x)>1e-10, sum(dU,dims=1)[:])
	Tet = config.T
	
	iTet = findall(x -> x in iV, Tet[:])
	iTet = ceil.(Int, iTet./4) |> unique

	copyto!(config.U, U)

	nat = length(config.at)
	Uat = config.X[:,1:nat] + U[:, 1:nat]
	set_positions!(config.at, Uat)

	nT = size(config.T,2)
	∇U, = gradient(config;domain=iTet) 
	config.∇U[:,:,iTet] = ∇U[:,:,iTet]
	if "W" in keys(config.data)
		volT = config.data["volT"]
		for i in iTet
			if isapprox(volT[i], 0, atol=1e-10)
				continue
			end
			w, dw = Wcb(∇U, config.Wat, i, config.J[:,:,i], config.calc)
			config.data["W"][i] = w
			config.data["dW"][:,:,i] = dw
		end
	end
	return config
end

function get_free!(config::AtC{Float64})
	nX = size(config.X, 2)
	free = setdiff(1:nX, config.iBdry)
	mask = fill(false, 3, nX)
	mask[:, free] .= true
	xfree = findall(mask[:])
	config.data["xfree"] = xfree
	return xfree
end

function get_x(config::AtC{Float64})
	U = deepcopy(config.U)
	if "xfree" in keys(config.data)
		xfree = convert(Array{Int64,1}, config.data["xfree"])
	else
		xfree = get_free!(config)
	end
	return U[xfree]
end

# function get_interface(X::Array{Float64,2}, Tet::Array{Int64,2})
function get_interface(config::AtC{Float64})
    X = config.X
    Tet = config.T
    xc = config.data["xc"]
    XType = config.data["XType"]
    Xat = X[:, XType .== 0.0]
    r = [norm(Xat[:,i] - xc,Inf) for i=1:size(Xat,2)]
	rin = config.Ra + config.bw
	Idx = findall(x-> x > rin, r)
    Iinterface = []

    for idx in Idx
        iTet = findall(x -> x == idx, Tet[:]) 
        iTet = ceil.(Int, iTet./4) |> unique
        neigs = Tet[1:4, iTet] |> unique
        if any(x -> x != 0.0, XType[neigs])
            push!(Iinterface, idx)
        end
    end
    return Iinterface
end
