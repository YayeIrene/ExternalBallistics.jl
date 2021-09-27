#module Scenario
#using  ..Types

#import generic functions
#import .Ballistics: InternalBallistics!, ExternalBallistics!, Vulnerability!

#export createProjectile, createFragment, createTarget



function createTarget(size::Float64, position::Float64)
    Target1D(size,position)
end

"""
    createTargetRect(a,b, position)
Return a rectangular target by specifying the length 'a', the height 'b' and the 'position' (range)
"""
function createTargetRect(a::Float64,b::Float64, position::Float64)
    TargetRect(a,b,position)
end

"""
    createTargetCirc(size, position)
Returns a circular target by specifying the radius (size) and the position (range)
"""
function createTargetCirc(size::Float64, position::Float64)
    TargetCirc(size,position)
end

"""
    createTarget()
Retursns an object of type 'Target'.
Optional arguments are:
* ρ the size of the target
* position, the range to the target
"""
function createTarget(;ρ=nothing,position=nothing)
    Target(position,ρ)
end

"""
    createGun(u₀,lat,QE,AZ,tc; args ...)
Returns an object of type 'Gun' by specifying the muzzle velocity (u₀) its latitude (lat), the elevation angle (QE),
azimuth angle (AZ) and the twist (tc).
Optional arguments are
* lw, the gun length
* X2w, the height of the gun
"""
function createGun(u₀::Float64,lat::Float64,QE::Float64,AZ::Float64,tc::Float64; lw=nothing, X2w=nothing)
    Gun(u₀,lat,tc,lw,X2w,QE,AZ)
end



#end
