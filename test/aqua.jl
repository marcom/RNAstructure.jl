import Aqua
using RNAstructure

@testset "Aqua.test_all" begin
    showtestset()
    Aqua.test_all(RNAstructure)
end
