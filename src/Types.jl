module Types

export   Target1D, Projectile1D
"""
AbstractPenetrator is the super type which contains the different types of penetrators
"""
abstract type AbstractPenetrator end
"""
AbstractTarget is the super type which contains all the different types of targets
"""
abstract type AbstractTarget end

"""
AbstractMaterial is the super type which contains all the different types of materials
"""
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
    σ::Float64
end


mutable struct Target1D <:AbstractTarget
    size::Float64
    position::Float64
end

end
