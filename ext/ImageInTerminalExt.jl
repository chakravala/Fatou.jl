module ImageInTerminalExt

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
isdefined(Fatou, :Requires) ? (import Fatou: ImageInTerminal) : (using ImageInTerminal)

import Fatou: ColorSchemes
function Base.show(io::IO,K::Fatou.FilledSet;c::String="",bare::Bool=false)
    isempty(c) && (c = K.meta.cmap)
    display(getproperty(ColorSchemes, isempty(c) ? :balance : Symbol(c))(K))
    !bare && print(io,String(K))
end

end # module
