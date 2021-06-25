using LinearAlgebra
using ScatteredInterpolation

# for 3 dimentional interpolation
function interp(Xlast::Array{Float64,2}, Xnext::Array{Float64,2}, Ulast::Array{Float64,2})
    Unext = zero(Xnext)
    for i = 1:3
        F = interpolate(Multiquadratic(), Xlast, Ulast[i,:])
        Unext[i, :] = evaluate(F, Xnext)
    end
    return Unext
end

function meshgrid(xx::Array{Float64}, yy::Array{Float64}, zz::Array{Float64})
	m, n, o = length(xx), length(yy), length(zz)
	xx = reshape(xx, 1, n, 1)
    yy = reshape(yy, m, 1, 1)
    zz = reshape(zz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
	return adjoint([vec(xx[om, :, oo]) vec(yy[:, on, oo]) vec(zz[om, on, :])])[:,:]
end

function meshgrid(xx::Array{Int}, yy::Array{Int}, zz::Array{Int})
	m, n, o = length(xx), length(yy), length(zz)
	xx = reshape(xx, 1, n, 1)
    yy = reshape(yy, m, 1, 1)
    zz = reshape(zz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
	return adjoint([vec(xx[om, :, oo]) vec(yy[:, on, oo]) vec(zz[om, on, :])])[:,:]
end

function GausWeights(ww::Array{Float64}, ii::Array{Int})
    nW = size(AA,2)
    W = zeros(nW, 1)
    for i = 1:nW
        W[i] = prod(ww[ii[:,i]])
    end
    return W
end 

function Tet_midpoint(X::Array{Float64,2})
    mpt = Array{Float64,2}(undef, 3, 6)
    for (i, j, k) ∈ zip([1,1,1,2,2,3],[2,3,4,3,4,4],[1,2,3,4,5,6])
        mpt[:,k] = 0.5*(X[:,i]+X[:,j])
    end
    mpIdx = [[1,2,3], [1,4,5], [2,4,6], [3,5,6]]
    return mpt, mpIdx
end

function Tet_circumcenter_face(X::Array{Float64,2})
    tfcc = Array{Float64,2}(undef, 3, 4)
    for (i,j,k,l) ∈ zip([1,1,1,2], [2,2,3,3], [3,4,4,4], [1,2,3,4])
        tfcc[:,l] = circumcenter_tri(X[:,[i,j,k]])
    end
    tfccIdx = [[1,2,3], [1,2,4], [1,3,4], [2,3,4]]
    return tfcc, tfccIdx
end

function circumcenter_tri(X::Array{Float64,2})
    x1, y1, z1 = X[:,1]
    x2, y2, z2 = X[:,2]
    x3, y3, z3 = X[:,3]
    a1 = y1*z2 - y1*z3 - z1*y2 + z1*y3 + y2*z3 - y3*z2
    b1 = -x1*z2 + x1*z3 + z1*x2 - z1*x3 -x2*z3 + x3*z2
    c1 = x1*y2 - x1*y3 - y1*x2 + y1*x3 + x2*y3 -x3*y2
    d1 = -x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x3*y1*z2 - x2*y3*z1 + x3*y2*z1
    a2 = 2*(x2 - x1)
    b2 = 2*(y2 - y1)
    c2 = 2*(z2 - z1)
    d2 = x1^2 + y1^2 + z1^2 - x2^2 -y2^2 - z2^2
    a3 = 2*(x3 - x1)
    b3 = 2*(y3 - y1)
    c3 = 2*(z3 - z1)
    d3 = x1^2 + y1^2 + z1^2 - x3^2 -y3^2 - z3^2

    A = [a1 b1 c1; a2 b2 c2; a3 b3 c3]
    b = [d1,d2,d3]
    return -A\b
end

function Tet_circumcenter(X::Array{Float64,2})
    r1 = sum(abs2.(X[:,1]))/2
    r2 = sum(abs2.(X[:,2]))/2
    r3 = sum(abs2.(X[:,3]))/2
    r4 = sum(abs2.(X[:,4]))/2
    A = X[:,[2,3,4]] - repeat(X[:,1], 1, 3)
    b = [r2-r1,r3-r1, r4-r1]
    return A'\b
end

function classify_tet(Tet::Array{Int64,2}, β::Vector{Float64})
    nT = size(Tet,2)
    Ttype = Vector{Int}(undef, nT)
    for i = 1:nT
        t = Tet[1:4, i]
        if all(x->x==1.0, β[t])
            Ttype[i] = 1
        elseif all(x->x==0.0, β[t])
            Ttype[i] = 0
        else
            Ttype[i] = -1
        end
    end
    return Ttype
end

function get_at_boundary(at::Atoms{Float64}; xc = nothing)
    set_pbc!(at, false)
	# locate the 'center' of atomistic lattice
	Xat = positions(at)
    if isnothing(xc)
        xcell = diag(cell(at))/2
        xcidx = findall(x->isapprox(x, xcell), Xat)
        if isnothing(xcidx)
            r = [ norm(x - xcell) for x in Xat ]
            _, xcidx = findmin(r)
            Xcidx = xcidx.I[2]
            @assert !isnothing(xcidx)
        end
        xc = Xat[xcidx]
        Xat .-= xc[:]
        set_positions!(at, Xat)
        
        deleteat!(at, xcidx)
        Xat = positions(at)
        Xm = mat(Xat)
    else
        Xm = mat(Xat)
        Xm .-= xc
    end

    Xm = mat(Xat)
    xmin, xmax = extrema(Xm)
    @assert abs(xmin) > abs(xmax)
    # xc = mat(xc)
    r = [norm(Xm[:,i],Inf) for i=1:size(Xm,2)]

    Idx = findall(x-> x > (xmax-0.1), r)
    Idx1 = findall(x->x > (abs(xmin)-0.1), r)

    Idx2 = setdiff(Idx, Idx1)

    for idx in Idx2
        x = Xm[:, idx]
        if any(i -> i>(xmax-0.1), x)
            push!(Idx1, idx)
        end
    end

    return Idx1
end