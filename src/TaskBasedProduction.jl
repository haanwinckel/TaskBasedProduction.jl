module TaskBasedProduction

include("common.jl")
include("unitInputDemand.jl")
include("component_positive_ups.jl")
include("component_negative_ups.jl")
include("margProdLabor.jl")

export unitInputDemand, component_positive_ups, component_negative_ups, margProdLabor

end # module