using Test
using Aqua
import RegisterMismatchCommon

# Tests are in RegisterMismatch and RegisterMismatchCuda

@testset "RegisterMismatchCommon" begin
    Aqua.test_all(RegisterMismatchCommon)
end
