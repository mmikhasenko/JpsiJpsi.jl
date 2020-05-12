
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

function A4K(vars; LS, J)
    @unpack cosθ1,ϕ1,cosθ2,ϕ2 = vars
    return real(sum(
        (λ1-λ2 == λ1′-λ2′ ? 1 : 0) *
            wignerd(1,λ1 ,0,cosθ1) * wignerd(1,λ2 ,0,cosθ2) *
            wignerd(1,λ1′,0,cosθ1) * wignerd(1,λ2′,0,cosθ2) *
            cis((λ1-λ1′)*ϕ1)*cis((λ2-λ2′)*ϕ2) *
                 Hλλ(λ1, λ2;  LS=LS, J=J) *
            conj(Hλλ(λ1′,λ2′; LS=LS, J=J))
            for λ1  in -1:1, λ2  in -1:1,
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
