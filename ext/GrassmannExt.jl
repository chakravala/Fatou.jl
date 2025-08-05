module GrassmannExt

#   This file is part of Fatou.jl.
#   It is licensed under the MIT license
#   Copyright (C) 2017 Michael Reed
#       _           _                         _
#      | |         | |                       | |
#   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _
#  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` |
# | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| |
#  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_|
#
#   https://github.com/chakravala
#   https://crucialflow.com

using Fatou
isdefined(Fatou, :Requires) ? (import Fatou: Grassmann) : (using Grassmann)

function Fatou.orbit(K::Fatou.Define{FT,QT,CT,M,N,P,D,B},Z0::Complex{Float64}) where {FT,QT,CT,M,N,P,D,B}
    V = Grassmann.Manifold(B)
    z0 = Grassmann.Couple{V,B}(Z0)
    M ? (z = Grassmann.Couple{V,B}(K.seed)) : (z = P ? Grassmann.Couple{V,B}(Fatou.plane(Z0)) : z0)
    zn = 0x0000
    while (N ? (Grassmann.value(K.Q(z,z0))::Float64>K.ϵ)::Bool : (Grassmann.value(K.Q(z,z0))::Float64<K.ϵ))::Bool && K.N>zn
        z = K.F(z,z0)::Grassmann.Couple{V,B,Float64}
        zn+=0x0001
    end; #end
    # return the normalized argument of z or iteration count
    return (zn::UInt16,(D ? Fatou.disk(Grassmann.value(z)) : Grassmann.value(z))::Complex{Float64})
end

end # module
