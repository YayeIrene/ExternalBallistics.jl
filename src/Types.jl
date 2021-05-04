#module Types

#export   Target1D, Projectile1D

abstract type AbstractPenetrator end

abstract type AbstractTarget end

abstract type AbstractMaterial end

mutable struct Projectile1D{T<: Union{AbstractPenetrator,Nothing}, R<:Union{AbstractMaterial,Nothing}} <: AbstractPenetrator
    calibre::Float64
    type::T
    mat::R
    mass::Float64
    velocity_m::Float64 #magnitude
    e_v::Union{Array{Float64,1},Nothing} #velocity_vec #direction cosines of the unit vector directed along the initial shotline
    position::Float64
    spin::Union{Float64,Nothing}
    yaw::Union{Float64,Nothing}
    b_f::Union{Array{Float64,1},Nothing} #body_vec#b_1f, b_2f, b_3f
    coeffAero
end


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

#end
