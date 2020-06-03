
xR2vars(x) = (cosθ1=atan(x[1])*2/π, cosθ2=atan(x[2])*2/π, ϕ=atan(x[3])*2)
function vars2xR(vars)
    @unpack cosθ1,cosθ2,ϕ = vars
    return [tan(cosθ1*π/2), tan(cosθ2*π/2), tan(ϕ/2)]
end

function I4μ(vars; H=error("helicity coupling"))
    @unpack cosθ1,cosθ2,ϕ = vars
    return real(sum(
        (λ1-λ2 == λ1′-λ2′ ? 1 : 0) *
            wignerd(1,λ1 ,ξ1,cosθ1) * wignerd(1,λ2 ,ξ2,cosθ2) * (isodd(1-λ2 ) ? -1 : 1) *
            wignerd(1,λ1′,ξ1,cosθ1) * wignerd(1,λ2′,ξ2,cosθ2) * (isodd(1-λ2′) ? -1 : 1) *
            cis((λ1-λ1′)*ϕ) *
                 H[ λ1+2,λ2+2] *
            conj(H[λ1′+2,λ2′+2])
            for ξ1 in [-1,1], ξ2 in [-1,1],
                λ1  in -1:1, λ2  in -1:1,
                λ1′ in -1:1, λ2′ in -1:1))
end

function I4K(vars; H=error("helicity coupling"))
    @unpack cosθ1,cosθ2,ϕ = vars
    return real(sum(
        (λ1-λ2 == λ1′-λ2′ ? 1 : 0) *
            wignerd(1,λ1 ,0,cosθ1) * wignerd(1,λ2 ,0,cosθ2) * (isodd(1-λ2 ) ? -1 : 1) *
            wignerd(1,λ1′,0,cosθ1) * wignerd(1,λ2′,0,cosθ2) * (isodd(1-λ2′) ? -1 : 1) *
            cis((λ1-λ1′)*ϕ1) *
                H[ λ1+2, λ2+2] *
           conj(H[λ1′+2,λ2′+2])
            for λ1  in -1:1, λ2  in -1:1,
                λ1′ in -1:1, λ2′ in -1:1))
end

# I4μ_LS(vars; LS, J) = I4μ(vars; H=[Hλλ(λ1, λ2;  LS=LS, J=J) for λ1 in -1:1, λ2  in -1:1])
# I4K_LS(vars; LS, J) = I4K(vars; H=[Hλλ(λ1, λ2;  LS=LS, J=J) for λ1 in -1:1, λ2  in -1:1])
#
function I4μ_intϕ(vars; H=error("helicity coupling"))
    @unpack cosθ1,cosθ2 = vars
    return real(sum(
            wignerd(1,λ1,ξ1,cosθ1) * wignerd(1,λ2,ξ2,cosθ2) *
            wignerd(1,λ1,ξ1,cosθ1) * wignerd(1,λ2,ξ2,cosθ2) *
                abs2(H[ λ1+2, λ2+2])
            for ξ1 in [-1,1], ξ2 in [-1,1],
                λ1  in -1:1, λ2  in -1:1))
end
