#module Scenario
#using  ..Types

#import generic functions
#import .Ballistics: InternalBallistics!, ExternalBallistics!, Vulnerability!

#export createProjectile, createFragment, createTarget



function createTarget(size::Float64, position::Float64)
    Target1D(size,position)
end

function createTargetRect(a::Float64,b::Float64, position::Float64)
    TargetRect(a,b,position)
end

function createTargetCirc(size::Float64, position::Float64)
    TargetCirc(size,position)
end

function createTarget(position::AbstractVector{Float64})
    Target(position)
end

function createGun(u₀::Float64,lat::Float64,QE::Float64,AZ::Float64,tc::Float64; lw=nothing, X2w=nothing)
    Gun(u₀,lat,tc,lw,X2w,QE,AZ)
end



#end
