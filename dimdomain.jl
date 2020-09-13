
# Need to use `] add DimensionalData#dataset` to get the latest things used here

# we will extend behavior from:
import GeoStatsBase
import DimensionalData

# and users will use features from:
using GeoStats
using DimensionalData
using Plots, Test

# first we create a wrapper type for the underlying domain of a
# dimensional array, similar to the wrapper DimTable for the values
"""
    DimDomain{C} <: AbstractDomain

A wrapper for the underlying domain of an AbstractDimArray.

See also DimTable.
"""
struct DimDomain{C} <: AbstractDomain
    dimcolumns::C
end
GeoStatsBase.nelms(dom::DimDomain) = length(dom.dimcolumns)
GeoStatsBase.ncoords(dom::DimDomain) = length(dims(dom.dimcolumns))

# this should return the type of the lat/lon pairs, probably Float64
# I will just hard code it for now
GeoStatsBase.coordtype(dom::DimDomain) = Float64

# this should set the buf to contain the coordinates of the i-th cell
# in the dimensional array
GeoStatsBase.coordinates!(buf, dom::DimDomain, i) = @error "not implemented"

# now that we have both DimDomain and DimTable, let's just forward
# them whenever we get a dimensional array
GeoStatsBase.domain(A::AbstractDimArray) = DimDomain(A)
GeoStatsBase.values(A::AbstractDimArray) = DimTable(A)

"""
    interpolate(da::AbstractDimArray, solver::AbstractSolver; dims=dims(da))

Either interpolate to fill missings in the current array,
or interpolate to a new domain specified by `dims`.
"""
interpolate(da::AbstractDimArray, solver::AbstractSolver; dims=dims(da)) = begin
    # it seems that an abstract array stores a single layer?
    # in this case, get the layer name, and interpolate on
    # a possibly finer grid. how do we construct such grid?
    # ideally we should be able to easily construct a DimDomain
    # from the dims keyword argument. because I am not familiar
    # with the DD.jl api, I think it is safer if you do it.
    # for now I will just interpolate on the same original grid,
    # effectively doing nothing
    prob = EstimationProblem(da, domain(da), DimensionalData.name(da))
    sol  = solve(prob, solver)
end

#########################################################################
# Tests

@testset "Test simulation on DimDomain" begin
    da = DimArray(rand(20, 30), (X(LinRange(-19.0, 19.0, 20)), Y(LinRange(3.0, 90.0, 30))))

    P = SimulationProblem(domain(da), :Z => Float64, 2)
    S = LUGaussSim(:Z=>(variogram=GaussianVariogram(range=30.0),))
    sol = solve(P, S)

    plot(sol)
end

@testset "Test DimDomain's GeoStatsBase.jl interface" begin
    da = DimArray(zeros(20, 30), (X([-19.0:2.0:19.0...]), Y([3.0:3.0:90.0...])))

    D = domain(da)
    @test nelms(D) == 20 * 30
    @test ncoords(D) == 2
    @test coordtype(D) == Float64
    @test coordinates(D, 4) == [-13.0, 3.0]
end

@testset "interpolate current array size to remove missing" begin
    # Define a regular index DimArray that can hold `missing`
    A = convert(Array{Union{Float64,Missing}}, rand(20, 30))
    # Make some values missing
    for i in 1:4:20, j in 1:5:20
        A[i, j] = missing
    end
    # Wrap as a DimArray
    da = DimArray(A, (X(LinRange(-19.0, 19.0, 20)), Y(LinRange(3.0, 90.0, 30))), :data)

    # interpolate the array with tne Kriging solver
    res = interpolate(da, Kriging())

    # @test res isa DimArray
    # @test dims(res) == dims(da)
    # @test size(res) == size(da)
    # @test any(ismissing, da) == true
    # @test any(ismissing, res) == false

    heatmap(da)
    heatmap(res)
end

@testset "interpolate to a larger array" begin
    # Define a regular index DimArray
    A = rand(20, 30)
    # Wrap as a DimArray
    da = DimArray(A, (X(LinRange(-19.0, 19.0, 20)), Y(LinRange(3.0, 90.0, 30))), :data)
    newdims = (X(LinRange(-19.0, 19.0, 45))), Y(LinRange(3.0, 90.0, 53))

    # interpolate the array with tne Kriging solver and newdims for the domain
    result = interpolate(da, Kriging(); dims=newdims)

    # @test result isa DimArray
    # # `dims` for the returned array wont be identical to `newdims` as they are formatted
    # # to match the array. But theyre basically the same type:
    # @test dims(result) isa Tuple{<:X{<:LinRange,<:Sampled,Nothing},<:Y{<:LinRange,<:Sampled,Nothing}}
    # # Array is the right size
    # @test size(result) == map(length, newdims)
    # # Dims have the same index as newdims
    # @test index(result) == index(newdims)

    # `heatmap` the interpolated array
    heatmap(result)
    # And compare the original
    heatmap(da)
end

# @testset "interpolate all arrays in dataset to larger arrays" begin
#     # Define a regular index DimArray
#     layers = (a=rand(20, 30), b=rand(Float32, 20, 30))
#     # Wrap as a DimArray
#     dimz = (X(LinRange(-19.0, 19.0, 20)), Y(LinRange(3.0, 90.0, 30)))
#     ds = DimDataset(layers, dimz)
#     newdims = (X(LinRange(-19.0, 19.0, 45))), Y(LinRange(3.0, 90.0, 53))

#     # interpolate the array with tne Kriging solver and newdims for the domain
#     result = interpolate(ds, Kriging(); dims=newdims)

#     # The return value is a DimArray
#     @test result isa DimDataset
#     # `dims` for the returned array wont be identical to `newdims` as they are formatted
#     # to match the array. But theyre basically the same type:
#     @test dims(result) isa Tuple{<:X{<:LinRange,<:Sampled,Nothing},<:Y{<:LinRange,<:Sampled,Nothing}}
#     # Array is the right size
#     @test size(first(values(result))) == map(length, newdims)
#     # Dims have the same index as newdims
#     @test index(result) == index(newdims)

#     # `heatmap` the interpolated arrays
#     heatmap(result[:a])
#     heatmap(result[:b])
# end

# TODO: test more things



#####################################################################################


#= Discussion point: the domain is nearly identical to the table:

struct DimDomain{C} <: AbstractDomain
    dimcolumns::C
end

struct DimTable{Keys,A,C} <: Tables.AbstractColumns
    array::A
    dimcolumns::C
end

This begs the question: why is the domain even separate?
The key thing that needs to be specified is which columns are domain,
which are values.

That could be done with extra interface methods for AbstractTable

- `dimcolnames`/`domaincolnames`
- `valcolnames`

Or something similar.

Then you can implement coordinates! generically
on any table that defines `dimcolnames` or whatever it is.

=#
