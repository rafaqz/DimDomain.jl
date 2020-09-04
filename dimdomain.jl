
using DimensionalData, GeoStatsBase, GeoStats, Plots, BenchmarkTools
      DirectGaussianSimulation, Test, StaticArrays, LinearAlgebra, 

using DimensionalData: DimColumn


struct DimDomain{C} <: AbstractDomain
    dimcolumns::C
end
DimDomain(dims::Tuple{Vararg{<:Dimension}}) = begin
    isregular = map(dims) do d
        mode(d) isa Sampled && span(d) isa Regular
    end
    println(isregular)
    if all(isregular) # Use the RegularGrid Domain
        lengths = map(length, dims)
        origins = map(first, dims)
        steps = map(d -> val(span(d)), dims)
        # RegularGrid could just accept AbstractRange?
        RegularGrid(lengths, origins, steps)
    else 
        # Generate the irregular domain from the index
        dimcolumns = map(d -> DimColumn(d, dims), dims)
        DimDomain(dimcolumns) 
    end
end
DimDomain(A::AbstractDimArray) = DimDomain(dims(A))

dimcolumns(dom::DimDomain) = dom.dimcolumns
# This syntax should be fixed in DD for DimColumn
DimensionalData.dims(dom::DimDomain) = map(DimensionalData.dim, dimcolumns(dom))  


# GeoStatsBase interface
GeoStatsBase.nelms(dom::DimDomain) = length(first(dimcolumns(dom)))
GeoStatsBase.ncoords(dom::DimDomain) = length(dimcolumns(dom))
GeoStatsBase.coordtype(dom::DimDomain) = eltype(first(dimcolumns(dom)))
GeoStatsBase.coordinates!(buf::AbstractVector, dom::DimDomain, i::Int) = 
    buf[i] = SA[map(c -> c[i], dimcolumns(dom))...]
GeoStatsBase.coordinates!(buf::MArray, dom::DimDomain, i::Int) =
    buf .= map(c -> c[i], dimcolumns(dom))


# Define a regular index DimArray
da = DimArray(rand(20, 30), (X(-19.0:2.0:19.0), Y(3.0:3.0:90.0)))
rg = DimDomain(da)
# Get a RegularGrid domain
@test rg isa RegularGrid

P = SimulationProblem(rg, :Z => Float64, 2)
S  = DirectGaussSim(:Z=>(variogram=GaussianVariogram(range=30.0),))
@btime sol = solve(P, S)
# 17.760 ms (266 allocations: 11.06 MiB)

plot(sol)

# Define an identical but "Irregular" marked DimArray 
da = DimArray(zeros(20, 30), (X([-19.0:2.0:19.0...]), Y([3.0:3.0:90.0...])))
D = DimDomain(da)
# Get a DimDomain
@test D isa DimDomain
@test GeoStatsBase.nelms(D) == 20 * 30
@test GeoStatsBase.ncoords(D) == 2
@test GeoStatsBase.coordtype(D) == Float64

# It works jus the same
P = SimulationProblem(D, :Z => Float64, 2)
S  = DirectGaussSim(:Z=>(variogram=GaussianVariogram(range=30.0),))
@btime sol = solve(P, S)
# 19.930 ms (265 allocations: 11.06 MiB) - a little slower
# but we are indexing into a vector

# But it doesn't plot. This could "just work" using the interface?
# plot(sol)




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
