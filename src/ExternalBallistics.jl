module ExternalBallistics
using Distributions, LinearAlgebra
using DifferentialEquations, Distances
using Interpolations, DataFrames
#using Meshes

# Write your package code here.
include("Types.jl")
#using .Types
include("Scenario.jl")
#using .Scenario
include("Utils.jl")
#using .Utils
include("OneD.jl")
include("Aero.jl")
#using .OneD
include("Mpmm.jl")
include("fragmentation.jl")

export Projectile, AbstractTarget, AbstractPenetrator, Wind, createTarget, createGun, Air, trajectoryMPMM, Gun,
QEfinderMPMM!,Zone,Fragment,FragShapes,frag

end
