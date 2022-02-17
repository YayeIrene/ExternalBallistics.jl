#module Scenario
#using  ..Types

#import generic functions
#import .Ballistics: InternalBallistics!, ExternalBallistics!, Vulnerability!

#export createProjectile, createFragment, createTarget




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

"""
    generateFragment(mass;optional arguments)

Creates an object fragment with mass (kg).
optional arguments are: position, velocity, time of flight, moment of inertia, distance to the center of gravity.
"""
function generateFragment(mass::Float64; rad=nothing, ax=nothing, density=nothing, cd_sub=nothing, Aₚ=nothing, position=nothing, velocity =nothing, tof =nothing)
    Fragment( rad,ax,mass,density,cd_sub,Aₚ, position,velocity,tof)
        
end


#end
