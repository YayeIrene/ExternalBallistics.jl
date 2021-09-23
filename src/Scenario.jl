#module Scenario
#using  ..Types

#import generic functions
#import .Ballistics: InternalBallistics!, ExternalBallistics!, Vulnerability!

#export createProjectile, createFragment, createTarget



function createTarget(size::Float64, position::Float64)
    Target1D(size,position)
end

"""
    creates a rectangular target by specifying the length a, the height b and the position (range)
"""
function createTargetRect(a::Float64,b::Float64, position::Float64)
    TargetRect(a,b,position)
end

"""
    creates a circular target by specifying the radius ρ and the position (range)
"""
function createTargetCirc(size::Float64, position::Float64)
    TargetCirc(size,position)
end

"""
    Creates a target. Optional options are the size and the position
"""
function createTarget(;ρ=nothing,position=nothing)
    Target(position,ρ)
end

"""
    Creates a Gun by specifying the muzzle velocity its latitude, the elevation angle, azimuth angle and the twist. Optional options are the gun length and the height of the gun
"""
function createGun(u₀::Float64,lat::Float64,QE::Float64,AZ::Float64,tc::Float64; lw=nothing, X2w=nothing)
    Gun(u₀,lat,tc,lw,X2w,QE,AZ)
end



#end
