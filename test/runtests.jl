using Test
using Aqua
using ExplicitImports
import RegisterMismatchCommon

# Tests are in RegisterMismatch and RegisterMismatchCuda

@testset "RegisterMismatchCommon" begin
    Aqua.test_all(RegisterMismatchCommon)
    @testset "ExplicitImports" begin
        test_explicit_imports(RegisterMismatchCommon)
    end
end
