#module Utils
#using ..Types
#using LinearAlgebra

#export machNumber, density


function machNumber(v::Float64)
    ma = v / 340.0
    return ma
end
function machNumber(v::Float64, y::Float64)
    t = temperature(y)
    a = speedOfSound(t)
    ma = v / a
    return ma
end

function temperature(y::Float64;y0 = 0.0, t0 = 288.15, β = -6.5)
    t = (t0+β * (y-y0)*1e-3)#-273.15 #deg Celsius
    return t
end

function speedOfSound(t::Float64; γ = 1.4, R = 287.05287)
    a = sqrt(γ*R*t)
    return a
end
function gravity(y::Float64; g0 = 9.80665, r = 6356766.0)
    g = g0*(r/(r+y))^2
    return g
end
function pressure(y::Float64;p0 = 101325.0, y0 = 0.0, t0 = 288.15, β = -6.5, R = 287.05287)
    p = p0*(1+β/t0*(y-y0)*1e-3)^(-gravity(y)/(β*R))
    return p
end
function density(y::Float64; R = 287.05287)
    t = temperature(y)
    p = pressure(y)
    ρ = p/(R*t)
    return ρ
end

#end
