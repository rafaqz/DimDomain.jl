# ] add DimensionalData#tables

using DimensionalData, GeoStatsBase


struct DimDomain{C} <: AbstractDomain
    dimcolumns::C
end
DimDomain(dims::Tuple{Vararg{<:Dimension}}) = begin
    dimcolumns = map(d -> DimColumn(d, dims), dims)
    DimDomain(dimcolumns)
end

columns(dom::DimDomain) = dom.columns

GeoStatsBase.nelms(dom::DimDomain) = prod(map(length, dims(dom)))
GeoStatsBase.ncoords(dom::DimDomain{Tuple{Vararg{<Dimension,N}}}) where N = N
GeoStatsBase.coordtype(dom::DimDomain) = eltype(first(dims(dom)))
GeoStatsBase.coordinates!(buf, dom::DimDomain, i::Int) = map(c -> c[i], columns(dom))


#= Discussion point: the domain is nearly identical to the table: 

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
