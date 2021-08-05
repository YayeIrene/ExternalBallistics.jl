function inter(mach::Vector{Float64},coeff::Vector{Float64})
    coeff_inter = LinearInterpolation(mach, coeff, extrapolation_bc=Flat())
    return coeff_inter

end

function inter_matrix(mach::Vector{Float64},coeff::Array{Float64, 2},α::Vector{Float64})
    coeff_inter = LinearInterpolation((mach,α), coeff,extrapolation_bc=Flat())
    return coeff_inter

end

function drag(CD0,CDδ2,αₜ,mach)
    CD = CD0(mach) + CDδ2(mach) * (sin(αₜ))^2
    return CD
end
function drag(CX0,CX2,αₜ,mach,CNa)
    CX = CX0(mach) + CX2(mach) * (sin(αₜ))^2
    CD = CNa(mach) *sin(αₜ) - CX*cos(αₜ)
    return -CD
end

function lift(CLα,CLα³,αₜ,mach)
    CL = CLα + CLα³ * (sin(αₜ))^2
    return CL
end
function lift(CX0,CX2,αₜ,mach,CNa)
    CX = CX0(mach) + CX2(mach) * (sin(αₜ))^2
    CL = CNa(mach) *cos(αₜ) + CX
    return CL
end

function matrix(m,args...)

for arg in args
    m=hcat(m,arg)
end
return m
end
