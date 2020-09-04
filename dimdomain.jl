# ] add DimensionalData#tables

using DimensionalData, GeoStatsBase

using DimensionalData: DimColumn


struct DimDomain{C} <: AbstractDomain
    dimcolumns::C
end
DimDomain(dom::Tuple{Vararg{<:Dimension}}) = begin
    dimcolumns = map(d -> DimColumn(d, dims), dims)
    DimDomain(dimcolumns) 
end
dimcolumns(dom::DimDomain) = dom.dimcolumns
# This syntax should be fixed in DD for DimColumn
dims(dom::DimDomain) = map(DimensionalData.dim, dimcolumns(dom))  

GeoStatsBase.nelms(dom::DimDomain) = length(first(dimcolumns))

GeoStatsBase.ncoords(dom::DimDomain{Tuple{Vararg{<:DimColumn,N}}}) where N = N

GeoStatsBase.coordtype(dom::DimDomain) = eltype(first(dimcolumns(dom)))

GeoStatsBase.coordinates!(buf, dom::DimDomain, i::Int) = 
    map(c -> c[i], dimcolumns(dom))


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
