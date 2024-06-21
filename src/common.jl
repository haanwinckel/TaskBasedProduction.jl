using SpecialFunctions: gamma_inc

low_norm_gamma_inc(a, x) = gamma_inc(a, x, 0)[1]
