
using DimensionalData, GeoStatsBase, GeoStats, Plots, BenchmarkTools
      DirectGaussianSimulation, Test, StaticArrays, LinearAlgebra, 

using DimensionalData: DimColumn

# Define a domain object that can use Array Dimensions
# to generate coordinates
struct DimDomain{C} <: AbstractDomain
    dimcolumns::C
end

getdomain(dims::Tuple{Vararg{<:Dimension}}) = begin

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

    # Otherwise, generate the irregular domain from the index
    else 
        dimcolumns = map(d -> DimColumn(d, dims), dims)
        DimDomain(dimcolumns) 
    end
end

# Allow getting the domain from any AbstractDimArray
getdomain(A::AbstractDimArray) = DimDomain(dims(A))


dimcolumns(dom::DimDomain) = dom.dimcolumns
# This syntax should be fixed in DD for DimColumn
DimensionalData.dims(dom::DimDomain) = map(DimensionalData.dim, dimcolumns(dom))  


# GeoStatsBase interface
GeoStatsBase.nelms(dom::DimDomain) = length(first(dimcolumns(dom)))
GeoStatsBase.ncoords(dom::DimDomain) = length(dimcolumns(dom))
GeoStatsBase.coordtype(dom::DimDomain) = eltype(first(dimcolumns(dom)))

# Why do we need both of these?
GeoStatsBase.coordinates!(buf::AbstractVector, dom::DimDomain, i::Int) = 
    buf[i] = SA[map(c -> c[i], dimcolumns(dom))...]
GeoStatsBase.coordinates!(buf::MArray, dom::DimDomain, i::Int) =
    buf .= map(c -> c[i], dimcolumns(dom))


@testset "RegularGrid is detected where possible" begin
    # Define a regular index DimArray
    da_r = DimArray(rand(20, 30), (X(LinRange(-19.0, 19.0, 20)), Y(LinRange(3.0, 90.0, 30))))

    # Get the domain from it
    rg = getdomain(da_r)

    # It's a RegularGrid domain
    @test rg isa RegularGrid

    # Use it for something
    P = SimulationProblem(rg, :Z => Float64, 2)
    S = DirectGaussSim(:Z=>(variogram=GaussianVariogram(range=30.0),))
    sol_rg2 = solve(P, S)

    # @btime 17.760 ms (266 allocations: 11.06 MiB)
    plot(sol)
end


@testset "DimDomain is used otherwise" begin

    # Define an identical but "Irregular" DimArray (the default for a vector index 
    da_i = DimArray(zeros(20, 30), (X([-19.0:2.0:19.0...]), Y([3.0:3.0:90.0...])))

    # Get the domain from it
    D = getdomain(da_i)
    # It returns a DimDomain
    @test D isa DimDomain

    # Check the geostatsbase interface works
    @test GeoStatsBase.nelms(D) == 20 * 30
    @test GeoStatsBase.ncoords(D) == 2
    @test GeoStatsBase.coordtype(D) == Float64
    buffer = [SA[0.0, 0.0] for i in 1:20 * 30]
    @test GeoStatsBase.coordinates!(buffer, D, 4) == [-13.0, 3.0]

    # Run the same problem, now with the 
    P = SimulationProblem(D, :Z => Float64, 2)
    S  = DirectGaussSim(:Z=>(variogram=GaussianVariogram(range=30.0),))
    sol_D = solve(P, S)
    # @btime 19.930 ms (265 allocations: 11.06 MiB) 
    # - a little slower but we are indexing into a vector
    
    # But it doesn't plot. This could "just work" using the interface?
    # plot(sol)
end



# Also be good to test something non-randomised






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
