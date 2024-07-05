module TaskBasedProduction

include("common.jl")
include("unitInputDemand.jl")
include("component_positive_ups.jl")
include("component_negative_ups.jl")
include("margProdLabor.jl")
include("prod_fun.jl")
include("elasticity_sub_comp.jl")
include("is_density_function.jl")
include("unitInputDemand_general.jl")
include("prod_fun_general.jl")
include("margProdLabor_general.jl")
include("elasticity_sub_comp_general.jl")
export unitInputDemand, margProdLabor,margProdLabor_general, prod_fun, elasticity_sub_comp, unitInputDemand_general, prod_fun_general, elasticity_sub_comp_general

end # module