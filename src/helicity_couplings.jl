#
LS_for_J0 = [(0,0), (2,2)];
LS_for_J1 = [(0,1), (2,1), (2,2)];
LS_for_J2 = [(0,2), (2,0), (2,1), (2,2), (4,2)];
#
function Hλλ(λ1,λ2; LS, J, j1=1, j2=1)
    L, S = LS
    Δλ = λ1-λ2
    (isodd(j2-λ2) ? -1 : 1) *
        ClGd(2*j1,2*λ1,2*j2,-2*λ2,2*S,2*Δλ) *
        ClGd(2*L,0,2*S,2*Δλ,2*J,2*Δλ)
end
#
function A4μ(vars; LS, J)
    @unpack cosθ1,ϕ1,cosθ2,ϕ2 = vars
    return real(sum(
        (λ1-λ2 == λ1′-λ2′ ? 1 : 0) *
            wignerd(1,λ1 ,ξ1,cosθ1) * wignerd(1,λ2 ,ξ2,cosθ2) *
            wignerd(1,λ1′,ξ1,cosθ1) * wignerd(1,λ2′,ξ2,cosθ2) *
            cis((λ1-λ1′)*ϕ1)*cis((λ2-λ2′)*ϕ2) *
                 Hλλ(λ1, λ2;  LS=LS, J=J) *
            conj(Hλλ(λ1′,λ2′; LS=LS, J=J))
            for ξ1 in [-1,1], ξ2 in [-1,1],
                λ1  in -1:1, λ2  in -1:1,
                λ1′ in -1:1, λ2′ in -1:1))
end
# 
function A4μ_intϕ(vars; LS, J)
    @unpack cosθ1,cosθ2 = vars
    return real(sum(
            wignerd(1,λ1,ξ1,cosθ1) * wignerd(1,λ2 ,ξ2,cosθ2) *
            wignerd(1,λ1,ξ1,cosθ1) * wignerd(1,λ2,ξ2,cosθ2) *
                 Hλλ(λ1,λ2; LS=LS, J=J) *
            conj(Hλλ(λ1,λ2; LS=LS, J=J))
            for ξ1 in [-1,1], ξ2 in [-1,1],
                λ1  in -1:1, λ2  in -1:1))
end
