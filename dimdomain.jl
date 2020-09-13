
# Need to use `] add DimensionalData#dataset` to get the latest things used here

using DimensionalData, GeoStatsBase, GeoStats, Plots, BenchmarkTools,
      DirectGaussianSimulation, Test, StaticArrays, LinearAlgebra, KrigingEstimators

using DimensionalData: DimColumn, formatdims, name


"""
    DimDomain{C} <: AbstractDomain

A domain object that can uses a tuple of `Dimension`s to generate coordinates.
"""
struct DimDomain{C} <: AbstractDomain
    dimcolumns::C
end
dimcolumns(dom::DimDomain) = dom.dimcolumns
# This syntax should be fixed in DD for DimColumn
DimensionalData.dims(dom::DimDomain) = map(DimensionalData.dim, dimcolumns(dom))  

# Domain interface for DimDomain

GeoStatsBase.nelms(dom::DimDomain) = length(first(dimcolumns(dom)))
GeoStatsBase.ncoords(dom::DimDomain) = length(dimcolumns(dom))
GeoStatsBase.coordtype(dom::DimDomain) = eltype(first(dimcolumns(dom)))

# Why do we need both of these methods?
GeoStatsBase.coordinates!(buf::AbstractVector, dom::DimDomain, i::Int) = 
    buf[i] = SA[map(c -> c[i], dimcolumns(dom))...]
GeoStatsBase.coordinates!(buf::Union{SubArray,MArray}, dom::DimDomain, i::Int) =
    buf .= map(c -> c[i], dimcolumns(dom))



# Data interface

"""
    domain(dims::Tuple{Vararg{<:Dimension}})

Get the domain from a tuple of Dimension
"""
GeoStatsBase.domain(dims::Tuple{Vararg{<:Dimension}}) = begin
    # Here we check if this these dimensions are all regular
    isregular = map(dims) do d
        mode(d) isa Sampled && span(d) isa Regular
    end

    # If they are regular, use the RegularGrid Domain
    if all(isregular)
        lengths = map(length, dims)
        origins = map(first, dims)
        steps = map(d -> val(span(d)), dims)
        # RegularGrid could just accept AbstractRange?
        RegularGrid(lengths, origins, steps)

    # Otherwise, generate a DimDomain from the dimensions
    else 
        dimcolumns = map(d -> DimColumn(d, dims), dims)
        DimDomain(dimcolumns) 
    end
end
GeoStatsBase.domain(x::Union{AbstractDimArray,AbstractDimDataset}) = 
    GeoStatsBase.domain(dims(A))

"""
    domain(A::AbstractDimArray)

Get the domain from any AbstractDimArray
"""
GeoStatsBase.domain(A::AbstractDimArray) = GeoStatsBase.domain(dims(A))

# This table also has the dimension columns in it, because its a general
# tables interface. It should still work? 
GeoStatsBase.values(A::AbstractDimArray) = DimTable(A)


"""
    interpolate(da::AbstractDimArray, solver::AbstractSolver; dims=dims(da))
    interpolate(ds::AbstractDimDataset, solver::AbstractSolver; dims=dims(ds))

Either interpolate to fill missings in the current array, 
or interpolate to a new domain specified by `dims`.

If the passed in object is a `AbstractDimDataset`, all layers will be interpolated.
"""
interpolate(da::AbstractDimArray, solver::AbstractSolver; dims=dims(da)) = begin
    key = name(da)
    data = SpatialData(domain(da), DimTable(da)) 
    problem = EstimationProblem(data, domain(dims), key, mapper=CopyMapper())
    sol = solve(problem, solver)
    newA = reshape(sol[key][:mean], map(length, dims))
    rebuild(da, newA, formatdims(newA, dims))
end
interpolate(ds::AbstractDimDataset, solver::AbstractSolver; dims=dims(ds)) =
    map(da -> interpolate(da, solver; dims=dims), ds)


#########################################################################
# Tests


@testset "RegularGrid is detected where possible" begin
    # Define a regular index DimArray
    da_r = DimArray(rand(20, 30), (X(LinRange(-19.0, 19.0, 20)), Y(LinRange(3.0, 90.0, 30))))

    # Get the domain from it
    rg = GeoStatsBase.domain(da_r)

    # It's a RegularGrid domain
    @test rg isa RegularGrid

    # Use it for something
    P = SimulationProblem(rg, :Z => Float64, 2)
    S = DirectGaussSim(:Z => (variogram=GaussianVariogram(range=30.0),)) 
    sol = solve(P, S)

    # @btime 17.760 ms (266 allocations: 11.06 MiB)
    plot(sol)
end


@testset "DimDomain is used if dims are irregular" begin
    # Define an identical but "Irregular" DimArray (the default for a vector index)
    da_i = DimArray(zeros(20, 30), (X([-19.0:2.0:19.0...]), Y([3.0:3.0:90.0...])))

    # Get the domain from it
    D = GeoStatsBase.domain(da_i)
    # It returns a DimDomain
    @test D isa DimDomain

    # Check the GeoStatsBase interface works
    @test GeoStatsBase.nelms(D) == 20 * 30
    @test GeoStatsBase.ncoords(D) == 2
    @test GeoStatsBase.coordtype(D) == Float64
    buffer = [SA[0.0, 0.0] for i in 1:20 * 30]
    @test GeoStatsBase.coordinates!(buffer, D, 4) == [-13.0, 3.0]

    P = SimulationProblem(da_i, :Z => Float64, 2)
    S  = DirectGaussSim(:Z=>(variogram=GaussianVariogram(range=30.0),))
    sol = solve(P, S)
    # @btime 19.930 ms (265 allocations: 11.06 MiB) 
    # - a little slower but we are indexing into a vector
    
    # But it doesn't plot. This could "just work" using the interface?
    # plot(sol)
end


@testset "We can just use any `AbstractDimArray` directly as the domain" begin
    da = DimArray(zeros(20, 30), (X([-19.0:2.0:19.0...]), Y([3.0:3.0:90.0...])))

    # Check the geostatsbase interface works
    @test GeoStatsBase.domain(da) isa DimDomain
    @test GeoStatsBase.nelms(da) == 20 * 30
    @test GeoStatsBase.ncoords(da) == 2
    @test GeoStatsBase.coordtype(da) == Float64
    buffer = [SA[0.0, 0.0] for i in 1:20 * 30]
    @test GeoStatsBase.coordinates!(buffer, da, 4) == [-13.0, 3.0]

    P = SimulationProblem(da, :Z => Float64, 2)
    S  = DirectGaussSim(:Z=>(variogram=GaussianVariogram(range=30.0),))
    sol = solve(P, S)
end

@testset "interpolate current array size to remove missing" begin
    # Define a regular index DimArray that can hold `missing`
    A = convert(Array{Union{Float64,Missing}}, rand(20, 30)) 
    # Make some values missing
    for i in 1:4:20, j in 1:5:20
        A[i, j] = missing
    end
    A

    @testset "Regular" begin
        # Wrap as a DimArray
        da = DimArray(A, (X(LinRange(-19.0, 19.0, 20)), Y(LinRange(3.0, 90.0, 30))), :data)

        # interpolate the array missing values with tne Kriging solver
        kr = Kriging(:data => (degree=1, variogram=SphericalVariogram(range=20.)))
        result = interpolate(da, kr) 


        @test result isa DimArray
        @test dims(result) == dims(da)
        @test size(result) == size(da)
        # We have intepolated the missing values
        @test any(ismissing, da) == true
        @test any(ismissing, result) == false

        # `heatmap` the interpolated array
        heatmap(da)
        heatmap(result)
        # Compare the original
    end

    @testset "Irregular" begin
        # Wrap as a DimArray
        da = DimArray(A, (X([-19.0:2.0:19.0...]), Y([3.0:3.0:90.0...])), :data)

        # interpolate the array with tne Kriging solver
        kr = Kriging(:data => (degree=1, variogram=SphericalVariogram(range=20.)))

        @test result isa DimArray
        @test dims(result) == dims(da)
        @test size(result) == size(da)
        # We have intepolated the missing values
        @test any(ismissing, da) == true
        @test any(ismissing, result) == false

        # `heatmap` the interpolated array
        heatmap(result)
        # Compare the original
        heatmap(da)
    end
end

@testset "interpolate to a larger array" begin
    # Define a regular index DimArray
    A = rand(20, 30) 
    # Wrap as a DimArray
    da = DimArray(A, (X(LinRange(-19.0, 19.0, 20)), Y(LinRange(3.0, 90.0, 30))), :data)
    newdims = (X(LinRange(-19.0, 19.0, 75))), Y(LinRange(3.0, 90.0, 73))
    kr = Kriging(:data => (degree=1, variogram=SphericalVariogram(range=20.)))

    # interpolate the array with tne Kriging solver and newdims for the domain
    result = interpolate(da, kr; dims=newdims) 

    # The return value is a DimArray
    @test result isa DimArray
    # `dims` for the returned array wont be identical to `newdims` as they are formatted 
    # to match the array. But theyre basically the same type:
    @test dims(result) isa Tuple{<:X{<:LinRange,<:Sampled,Nothing},<:Y{<:LinRange,<:Sampled,Nothing}}
    # Array is the right size
    @test size(result) == map(length, newdims)
    # Dims have the same index as newdims
    @test index(result) == index(newdims)

    # `heatmap` the interpolated array
    heatmap(da)
    savefig("original.png")
    heatmap(result)
    savefig("interpolated.png")
    # And compare the original
end

@testset "interpolate all arrays in dataset to larger arrays" begin
    # Define a regular index DimArray
    layerz = (a=rand(20, 30), b=rand(Float32, 20, 30))
    # Wrap as a DimArray


    @testset "Regular" begin
        dimz = (X(LinRange(-19.0, 19.0, 20)), Y(LinRange(3.0, 90.0, 30)))
        ds = DimDataset(layerz, dimz)
        newdims = (X(LinRange(-19.0, 19.0, 45))), Y(LinRange(3.0, 90.0, 53))
        # interpolate the array with tne Kriging solver and newdims for the domain
        result = interpolate(ds, Kriging(); dims=newdims) 

        @test result isa DimDataset
        # `dims` for the returned array wont be identical to `newdims` as they are formatted 
        # to match the array. But theyre basically the same type:
        @test dims(result) isa Tuple{<:X{<:LinRange,<:Sampled,Nothing},<:Y{<:LinRange,<:Sampled,Nothing}}
        # Array is the right size
        @test size(first(values(result))) == map(length, newdims)
        # result dims have the same index as newdims
        @test index(result) == index(newdims)

        # `heatmap` the interpolated arrays
        heatmap(result[:a])
        heatmap(result[:b])
    end

    @testset "Irregular" begin
        dimz = X([-19.0:2.0:19.0...]), Y([3.0:3.0:90.0...])
        ds = DimDataset(layerz, dimz)
        newdims = X([-19.0:1.0:19.0...]), Y([3.0:2.0:90.0...])
        # interpolate the array with tne Kriging solver and newdims for the domain
        result = interpolate(ds, Kriging(); dims=newdims) 

        @test result isa DimDataset
        # `dims` for the returned array wont be identical to `newdims` as they are formatted 
        # to match the array. But theyre basically the same type:
        @test dims(result) isa Tuple{<:X{<:Vector,<:Sampled,Nothing},<:Y{<:Vector,<:Sampled,Nothing}}
        # Array is the right size
        @test size(first(values(result))) == map(length, newdims)
        # result dims have the same index as newdims
        @test index(result) == index(newdims)

        # `heatmap` the interpolated arrays
        heatmap(result[:a])
        heatmap(result[:b])
    end
end

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
