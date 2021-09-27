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


mutable struct Target1D <:AbstractTarget
    size::Float64
    position::Float64
end

"""
    TargetRect(a,b,position)
defines a rectangular target by specifying the lenght a, the height b and the position (range)
"""
mutable struct TargetRect <:AbstractTarget
    a::Float64
    b::Float64
    position::Float64
end

"""
    TargetCirc(ρ,position)
defines a circular target by specifying the radius ρ and the position (range)
"""
mutable struct TargetCirc <:AbstractTarget
    ρ::Float64 #diameter
    position::Float64
end

mutable struct Target<:AbstractTarget
    position::Union{Array{Float64,1},Nothing}
    ρ::Union{Float64,Nothing}
end

struct Zone
    angle::Tuple{Float64, Float64}
    number::Tuple{Float64, Float64}
    mass::Tuple{Float64, Float64}
    velocity::Tuple{Float64, Float64}
    density::Float64
end

mutable struct Fragment <:AbstractPenetrator
    rad::Float64
    ax::Float64
    mass::Float64
    density::Float64
    cd_sub::Float64
    Aₚ::Float64
    position::Array{Float64,1}
    velocity::Array{Float64,1}
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

end


mutable struct FragShapes
    cd::Tuple{Float64, Float64}
    cf::Tuple{Float64, Float64}
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


#end
