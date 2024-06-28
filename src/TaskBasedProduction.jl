module TaskBasedProduction

include("common.jl")
include("unitInputDemand.jl")
include("component_positive_ups.jl")
include("component_negative_ups.jl")
include("margProdLabor.jl")
include("prod_fun.jl")
include("elasticity_sub_comp.jl")

export unitInputDemand, margProdLabor, prod_fun, elasticity_sub_comp

end # module