using Test, LinearAlgebra, TOML
include(joinpath(@__DIR__,"..","src","MolecularResponseBounds.jl"))
using .MolecularResponseBounds
using .MolecularResponseBounds.BoundsCommon
using .MolecularResponseBounds.ResponseBounds
using .MolecularResponseBounds.ExactResponse
using .MolecularResponseBounds.IndependentCheck
using .MolecularResponseBounds.StudyAPI
include("bounds_tests.jl")
include("exact_response_tests.jl")
include("driver_tests.jl")
include("caps_tests.jl")
