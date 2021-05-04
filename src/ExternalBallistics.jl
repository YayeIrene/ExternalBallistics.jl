module ExternalBallistics
using Distributions, LinearAlgebra

# Write your package code here.
include("Types.jl")
#using .Types
include("Scenario.jl")
#using .Scenario
include("Utils.jl")
#using .Utils
include("OneD.jl")
#using .OneD

export Projectile1D, createProjectile, TargetRect, TargetCirc, createTargetRect, createTargetCirc, AbstractTarget

end
