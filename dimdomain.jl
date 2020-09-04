
using DimensionalData, GeoStatsBase, GeoStats, Plots, 
      DirectGaussianSimulation, Test

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

GeoStatsBase.nelms(dom::DimDomain) = length(first(dimcolumns(dom)))

GeoStatsBase.ncoords(dom::DimDomain) = length(dimcolumns(dom))

GeoStatsBase.coordtype(dom::DimDomain) = eltype(first(dimcolumns(dom)))

GeoStatsBase.coordinates!(buf, dom::DimDomain, i::Int) = 
    map(c -> c[i], dimcolumns(dom))


# Define a regular index DimArray
da = DimArray(rand(20, 30), (X(-19.0:2.0:19.0), X(1.0:3.0:90.0)))
D = DimDomain(da)
# Get a RegularGrid domain
@test D isa RegularGrid

P = SimulationProblem(D, :Z => Float64, 2)
S  = DirectGaussSim(:Z=>(variogram=GaussianVariogram(range=30.0),))
sol = solve(P, S)
plot(sol)

# Define an Irregular index DimArray
da = DimArray(rand(9, 20), (X([1, 1, 2, 3, 5, 8, 13, 21, 34]), X(1.0:3.0:60.0)))
D = DimDomain(da)
# Get a DimDomain
@test D isa DimDomain

P = SimulationProblem(D, :Z => Float64, 2)
S  = DirectGaussSim(:Z=>(variogram=GaussianVariogram(range=30.0),))
# This errors 1 in 10 times with:
# ERROR: PosDefException: matrix is not positive definite; Cholesky factorization failed.
sol = solve(P, S)
# How to make plotting work? Could the interface make this easier?
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
