#module Types

#export   Target1D, Projectile1D

abstract type AbstractPenetrator end

abstract type AbstractTarget end

abstract type AbstractMaterial end


mutable struct Target1D <:AbstractTarget
    size::Float64
    position::Float64
end

mutable struct TargetRect <:AbstractTarget
    a::Float64
    b::Float64
    position::Float64
end

mutable struct TargetCirc <:AbstractTarget
    ρ::Float64 #diameter
    position::Float64
end

mutable struct Target<:AbstractTarget
    position::Array{Float64,1}
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

mutable struct Gun
    u₀::Float64
    lat::Union{Float64,Nothing}
    tc::Union{Float64,Nothing}
    lw::Union{Float64,Nothing}
    X2w::Union{Float64,Nothing}
    QE::Union{Float64,Nothing}
    AZ::Union{Float64,Nothing}
end

mutable struct Wind
    cross::Float64
    range::Float64
end

mutable struct Air
    p::Float64
    t::Float64
    rh::Float64
end


#end
