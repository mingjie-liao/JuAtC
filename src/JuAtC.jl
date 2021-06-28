module JuAtC

using JuLIP
using LinearAlgebra

include("BGFC/AtC.jl")
include("BGFC/FIO.jl")
include("BGFC/Solve.jl")
include("BGFC/utils.jl")
include("BGFC/beta.jl")

include("Adaptive/adaptive.jl")
end