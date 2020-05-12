#
function Hλλ(λ1,λ2; LS, J, j1=1, j2=1)
    L, S = LS
    Δλ = λ1-λ2
    (isodd(j2-λ2) ? -1 : 1) *
        CG(j1,λ1,j2,-λ2,S,Δλ) *
        CG(L,0,S,Δλ,J,Δλ)
end
#
intensity(H) = sum(abs2, H)
#
# LS couplings
const LS_for_J0 = [(0,0), (2,2)];
const LS_for_J1 = [(0,1), (2,1), (2,2)];
const LS_for_J2 = [(0,2), (2,0), (2,1), (2,2), (4,2)];
#
const JHs = [[[Hλλ(λ1,λ2; LS=LS, J=i-1) for λ1=-1:1, λ2=-1:1] for LS in LSs]
    for (i,LSs) in enumerate([LS_for_J0,LS_for_J1,LS_for_J2])]
const nJHs = [[h ./ sqrt(intensity(h)) for h in Hs] for Hs in JHs]
#
