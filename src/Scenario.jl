module Scenario
using  ..Types
using Distributions, LinearAlgebra
#import generic functions
#import .Ballistics: InternalBallistics!, ExternalBallistics!, Vulnerability!

export createProjectile, createFragment, createTarget

function createProjectile(calibre::Float64,mass::Float64, velocity::Float64,  position::Float64 ;ev = nothing, spin = nothing, bf = nothing, type=nothing, mat = nothing, yaw=nothing, cd = 0.0,disp=0.1)
            Projectile1D(calibre, type, mat, mass, velocity, ev, position, spin, yaw, bf,cd, disp)
end

function createTarget(size::Float64, position::Float64)
    Target1D(size,position)
end



end
