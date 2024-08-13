module TaskBasedProduction

include("common.jl")
include("unitInputDemand.jl")
include("component_positive_ups.jl")
include("component_negative_ups.jl")
include("margProdLabor.jl")
include("prodFun.jl")
include("elasticitySubComp.jl")
include("is_density_function.jl")
include("unitInputDemandGeneral.jl")
include("prodFunGeneral.jl")
include("margProdLaborGeneral.jl")
include("elasticitySubCompGeneral.jl")
include("numerical_derivative.jl")
include("getStartGuess_xT.jl")
include("getStartGuessGen_xT.jl")
export getStartGuess_xT, getStartGuessGen_xT, unitInputDemand, margProdLabor,margProdLaborGeneral, prodFun, elasticitySubComp, unitInputDemandGeneral, prodFunGeneral, elasticitySubCompGeneral

end # module