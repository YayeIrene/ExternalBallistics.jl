#module Types

#export   Target1D, Projectile1D
"""
    abstract type under which projectiles types are defined.
"""
abstract type AbstractPenetrator end

"""
    abstract type under which target types are defined.
"""
abstract type AbstractTarget end

abstract type AbstractMaterial end

abstract type AbstractZone end 


mutable struct Target<:AbstractTarget
    position::Union{Array{Float64,1},Nothing}
    ρ::Union{Float64,Nothing}
end

"""
    Zone(angle, μ_number, σ_number, μ_mass, σ_mass, μ_velocity, σ_velocity, density)
defines a framentation zone:
* angleBegin: begin angle of the segment (degrees)
* angleEnd: end angle of the segment (degrees)
* μ_number: average number of fragments
* σ_number: standard deviation of the distribution of number of fragments
* μ_mass: average mass of the fragments
* σ_mass: standard deviation of the distribution of the mass of fragments
* μ_velocity: average velocity of the fragments
* σ_velocity: standard deviation of the distribution of the velocity of fragments
* density: density of the fragments
"""
struct Zone<:AbstractZone
    angle::Tuple{Float64, Float64}
    number::Tuple{Float64, Float64}
    mass::Tuple{Float64, Float64}
    velocity::Tuple{Float64, Float64}
    density::Float64
end

"""
    ZoneZdata(angle, number, mass, velocity, density)
defines a framentation zone:
* angle (degrees): vector with 3 numbers: 
    - angle[1] is the lower point angle 
    - angle[2] is the midpoint angle
    - angle[3] is the upper point angle
* number: Vector containing the number of fragments in each mass group
* mass (kg): average mass of the fragments in each mass group
!!! warning
    The vectors number and mass have equal length
* velocity (m/s): vector with 3 numbers:
    - velocity[1] initial fragment speed for lower point angle
    - velocity[2] initial fragment speed for midpoint angle
    - velocity[3] initial fragment speed for upper point angle
* density: density of the fragments
"""
struct ZoneZdata<:AbstractZone
    angle::Vector{Float64}
    number::Vector{Float64}
    mass::Vector{Float64}
    velocity::Vector{Float64}
    density::Float64
end 

struct ZoneList{T} <: AbstractZone
    list::Vector{T}
end

"""
    Fragment(position, position, density)
defines the fragment by specifying:
* position
* position
* density
"""
mutable struct Fragment <:AbstractPenetrator
    rad::Union{Float64,Nothing}
    ax::Union{Float64,Nothing}
    mass::Union{Float64,Nothing}
    density::Union{Float64,Nothing}
    cd_sub::Union{Float64,Nothing}
    Aₚ::Union{Float64,Nothing}
    position::Union{Array{Float64,1},Nothing}
    velocity::Union{Array{Float64,1},Nothing}
    tof::Union{Float64,Nothing}

end

"""
    Projectile(mass,calibre,position,tof,Ix,Iy,Xcg)
Projectile is defined by the following parameters:
* mass (kg)
* calibre (m)
* positionsion: coordinate of the 'aimpoint'
* velocity: velocity vector
* tof: time of flight
* Ix: axial moment of inertia
* Iy: transversal moment of inertia
* Xcg: distance to the center of gravity
"""
mutable struct Projectile <:AbstractPenetrator
    mass::Float64
    calibre::Float64
    position::Union{Array{Float64,1}, Nothing}
    velocity::Union{Array{Float64,1}, Nothing}
    tof::Union{Float64,Nothing}
    Ix::Union{Float64,Nothing}
    Iy::Union{Float64,Nothing}
    Xcg::Union{Float64,Nothing}
    spin::Union{Float64,Nothing}
    αₑ::Union{Array{Float64,1},Nothing}


end




"""
    Gun(u₀,lat,tc,lw,X2w,QE,AZ)
Gun is defined by the following parameters:
* u₀: the muzzle velocity
* lat: the latitude
* tc: the twist
* lw: the gun length
* X2w: the height
* QE: the the elevation angle
* AZ: the azimuth
"""
mutable struct Gun
    u₀::Float64
    lat::Union{Float64,Nothing}
    tc::Union{Float64,Nothing}
    lw::Union{Float64,Nothing}
    X2w::Union{Float64,Nothing}
    QE::Union{Float64,Nothing}
    AZ::Union{Float64,Nothing}
end

"""
    Wind(cross, range)
Wind is defined by the:
* cross wind magnitude
* range wind magnitude
"""
mutable struct Wind
    cross::Float64
    range::Float64
end

"""
    Air(p,t,rh)
defines the atmosphere by specifying:
* p: the pressure
* t: the temperature
* rh: the relative humidity
"""
mutable struct Air
    p::Float64
    t::Float64
    rh::Float64
end


"""
    FragShapes(cd, cf)
defines the fragment by specifying:
* cd : drag coefficient
* cf : shape factor
"""
mutable struct FragShapes
    cd::Tuple{Float64, Float64}
    cf::Tuple{Float64, Float64}
end
