"""
## module FIO

## Description:
File I/O stream controler
1. read geometrical contents X, T from .msh file exported by dealII
""" 

module ACFIO
using DelimitedFiles
using Printf
# read mesh information in 3D
function read_dealii_3D(fn::AbstractString)
    open(fn, "r") do io
        read_msh_dealii_3D(io)
    end
end
    
function read_msh_dealii_3D(io)
    dim = 3
    thisLine = io |> readline |> strip
    thisLine = io |> readline |> strip
    NV = parse(Int, thisLine)
    X = zeros(dim+1, NV)
    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Float64)
        X[:,i] = d
    end

    thisLine = io |> readline |> strip
    thisLine = io |> readline |> strip
    thisLine = io |> readline |> strip
    NV = parse(Int, thisLine)
    T = Array{Int64, 2}(zeros(dim+10, NV))
    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Int64)
        T[:,i] = d
    end
    return X[2:4, :], T[6:end, :]
end

"""
    read_mesh(filename) -> mesh::Mesh
"""
function read_mesh(fn::AbstractString)
    open(fn, "r") do io
        read_mesh(io)
    end
end

"""
    read_mesh(iostream) -> mesh::Mesh

Reads the mesh nodes, edges and elements stored in the input .mesh file

Returns an object `mesh` of type `Mesh`, comprising both vector arrays.
"""
function read_mesh(io)

    thisLine = io |> readline |> strip
    while thisLine != "Dimension"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    dim = parse(Int, thisLine)
#   println("Dimension = ", dim)
    if dim != 3
        error("Only 3 dimensional problem are considered so far!")
    end

	thisLine = io |> readline |> strip
    while thisLine != "Vertices"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    NV = parse(Int, thisLine)

#    P = SVector{dim+1,Float64}
    Vex = Array{Float64}(undef,dim+1,NV)
    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Float64)
        Vex[:,i] = d
    end

    thisLine = io |> readline |> strip
    while thisLine != "Tetrahedra"
        thisLine = io |> readline |> strip
    end

	thisLine = io |> readline |> strip
    NT = parse(Int, thisLine)

#    P = SVector{dim+2,Int64}
    Tet = Array{Int64}(undef,dim+2,NT)
    for i in 1:NT
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Int64)
        Tet[:,i] = d
    end
    return (Vex,Tet)
end

"""
	write_mesh(filename, particle positions, particle types; Edges) -> filename |> mesher3d
	rule:
		0: atoms
		1: boundary
		2: a/c interface
		3: continuum
"""
function write_mesh(filename, X::Array{Float64,2}, Xtype::Vector{Int64}; Edges = nothing)
	@assert size(X, 2) == length(Xtype)
    f = open(filename, "w")
	@printf(f, "MeshVersionFormatted 2\n\n\nDimension 3\n\n\n");
	NV = size(X, 2)
	@printf(f, "Vertices\n%d\n", NV);
	for i = 1:NV
		@printf(f, "%f %f %f %d\n", X[1,i], X[2,i], X[3,i], Xtype[i])
	end

    if !isnothing(Edges)
        NE = size(Edges, 3)
        Ne = size(Edges, 1)
        @printf(f, "# Edges\n# %d, %d\n", NE, Ne);
        for i = 1:NE
            edge = Edges[:,:,i]
            @printf(f, "# %d\n", i)
            for j = 1:Ne
                @printf(f, "# %d %d\n", edge[j,1], edge[j,2])
            end
        end
    end

	@printf(f, "\n\nEnd\n")
    close(f)
end

function write_remesh(filename, X::Array{Float64,2}, TetIdx::Vector{Int64})
    f = open(filename, "w")
    NV = size(X,2)
	@printf(f, "Append_points\n%d\n", NV);
	for i = 1:NV
		@printf(f, "%f %f %f\n", X[1,i], X[2,i], X[3,i])
    end
    @printf(f, "\n")
    NT = length(TetIdx)
    @printf(f, "Refine_elements\n%d\n", NT)
    for j = 1:NT
        @printf(f, "%d\n", TetIdx[j])
    end
	# @printf(f, "\n\nEnd\n")
    close(f)
end

function write_value(filename, U::Array{Float64,2})
    f = open(filename, "w")
    NU = size(U,2)
	@printf(f, "vector displacement\n");
    @printf(f, "%d\n", NU);
	for i = 1:NU
		@printf(f, "%f %f %f\n", U[1,i], U[2,i], U[3,i])
    end
    close(f)
end

function read_value(fn::AbstractString)
    open(fn, "r") do io
        read_value(io)
    end
end

function read_value(io)

    thisLine = io |> readline |> strip
    while thisLine != "vector displacement"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    NV = parse(Int, thisLine)
    U = Array{Float64}(undef, 3, NV)

    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Float64)
        U[:,i] = d
    end

    return U
end

"""
return filename.dump file for visualization with X, the vertices of the coupled mesh.  
X should be moved with: X[1:3,:] .-= minimum(X), to put particles inside the bounding box.
"""
function write_dump(filename, X::Array{Float64,2})
    dl = " "
    n = size(X, 2)
    L = ceil(maximum(X))
    open(filename, "w") do f
        println(f, "ITEM: TIMESTEP")
        println(f, 1)
        println(f, "ITEM: NUMBER OF ATOMS")
        println(f, n)
        println(f, "ITEM: BOX BOUNDS pp pp pp")
        println(f, 0.00, dl, L)
        println(f, 0.00, dl, L)
        println(f, 0.00, dl, L)
        println(f, "ITEM: ATOMS id type x y z fx fy fx radius")
        for i = 1:n
            u = X[1:3,i]
            # r = config.M[i]
            z = X[4, i]
            r = 1.5
            if z != 0
                z = 1
                r = 1.0
            end
            println(f, i, dl, z, dl, u[1], dl, u[2], dl, u[3], dl, 0.00, dl, 0.00, dl, 0.00, dl, r)
        end
    end
end

end
