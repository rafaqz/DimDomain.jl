using DimensionalData, GeoStatsBase

struct DimDomain{D} <: AbstractDomain
    dims::D
end
DimDomain(dims::Tuple{Vararg{<:Dimension}}) = begin
    # Do some type conversions / checks here
    DimDomain{typeof(dims)}(dims)
end

GeoStatsBase.nelms(dom::DimDomain) = product(map(length, dims(dom)))
GeoStatsBase.ncoords(dom::DimDomain{Tuple{Vararg{<Dimension,N}}}) where N = N
GeoStatsBase.coordtype(dom::DimDomain) = eltype(first(dims(dom)))
GeoStatsBase.coordinates!(buf, dom::DimDomain, ind) = nothing # compute the coords from dims