using Test
using Unitful: Quantity, @u_str
using Random: randstring
using RNAstructure
using RNAstructure: run_AllSub, run_ct2dot, run_CycleFold, run_draw,
    run_dot2ct, run_EDcalculator, run_efn2, run_EnsembleEnergy,
    run_Fold, run_MaxExpect, run_partition!, run_ProbabilityPlot,
    run_RemovePseudoknots, run_stochastic

@testset verbose=true "RNAstructure" begin
    include("ct-format.jl")
    include("plot.jl")
    include("RNAstructure.jl")
end
