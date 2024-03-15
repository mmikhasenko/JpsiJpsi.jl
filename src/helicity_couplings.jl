#
function Hλλ(λ1, λ2; Pm1J, m1J, hvalues)
    # @unpack h11, h10, h1m1, h00 = hvalues
    (λ1, λ2) == (-λ1, -λ2) && Pm1J == '-' && return 0.0
    (λ1, λ2) == (λ2, λ1) && m1J == '-' && return 0.0
    (λ1, λ2) == (1, -1) && m1J != Pm1J && return 0.0
    λ1 < λ2 && return (m1J == '+' ? 1.0 : -1.0) * Hλλ(λ2, λ1; Pm1J=Pm1J, m1J=m1J, hvalues=hvalues)
    (λ2 < 0) && (λ1 != 1) && return (Pm1J == '+' ? 1.0 : -1.0) * Hλλ(-λ1, -λ2; Pm1J=Pm1J, m1J=m1J, hvalues=hvalues)
    # return 1.0
    (λ1, λ2) == (1, 1) && return hvalues.h11
    (λ1, λ2) == (1, 0) && return hvalues.h10
    (λ1, λ2) == (0, 0) && return hvalues.h00
    (λ1, λ2) == (1, -1) && return hvalues.h1m1
    error("(λ1,λ2)=($λ1,$λ2) does not in the set of expected values")
end
#
# [Hλλ(λ1,λ2; Pm1J='+', m1J='+', hvalues = (h11=1, h10=2, h1m1=3, h00=4)) for λ1=-1:1, λ2=-1:1]
# [Hλλ(λ1,λ2; Pm1J='-', m1J='+', hvalues = (h11=1, h10=2)) for λ1=-1:1, λ2=-1:1]
# [Hλλ(λ1,λ2; Pm1J='+', m1J='-', hvalues = (h10=2,)) for λ1=-1:1, λ2=-1:1]
# [Hλλ(λ1,λ2; Pm1J='-', m1J='-', hvalues = (h10=2, h1m1=3)) for λ1=-1:1, λ2=-1:1]
#
# general fro groups
@with_kw struct group
    Pm1J::Char
    m1J::Char
    hs::Vector{NamedTuple}
end
#
groupI = group(Pm1J='+', m1J='+',
    hs=[(h11=1.0, h10=0.0, h1m1=0.0, h00=0.0),
        (h11=0.0, h10=1.0, h1m1=0.0, h00=0.0),
        (h11=0.0, h10=0.0, h1m1=1.0, h00=0.0),
        (h11=0.0, h10=0.0, h1m1=0.0, h00=1.0)]);
groupII = group(Pm1J='-', m1J='+',
    hs=[(h11=1.0, h10=0.0),
        (h11=0.0, h10=1.0)]);
groupIII = group(Pm1J='+', m1J='-', hs=[(h10=1.0,)]);
groupIV = group(Pm1J='-', m1J='-',
    hs=[(h10=1.0, h1m1=0.0),
        (h10=0.0, h1m1=1.0)]);
#
groups = [groupI, groupII, groupIII, groupIV]
#
spec0p = group(Pm1J='+', m1J='+',
    hs=[(h11=1.0, h10=0.0, h1m1=0.0, h00=0.0),
        (h11=0.0, h10=0.0, h1m1=0.0, h00=1.0)]);
spec0m = group(Pm1J='-', m1J='+',
    hs=[(h11=1.0, h10=0.0)]);
spec1p = group(Pm1J='-', m1J='-',
    hs=[(h10=1.0, h1m1=0.0)]);
specs = [spec0p, spec0m, spec1p]
#
#
const gHs = [[[Hλλ(λ1, λ2; Pm1J=g.Pm1J, m1J=g.m1J, hvalues=h) for λ1 = -1:1, λ2 = -1:1] for h in g.hs]
             for g in vcat(groups, specs)]
#
intensity(H) = sum(abs2, H)
const ngHs = [[h ./ sqrt(intensity(h)) for h in Hs] for Hs in gHs]
#
#