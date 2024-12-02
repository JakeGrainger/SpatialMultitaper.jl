@testset "resampling" begin
    region = Box(Point(0,0), Point(100,100))
    shift = ToroidalShift(region)
    @test all(shift.shift.min .≈ (-50.0,-50.0))
    @test all(shift.shift.max .≈ (50.0,50.0))
    shift_0 = rand(shift.shift)
    @test all(x->-50≤x≤50, shift_0)
end