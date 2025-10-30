module Utils

function boltzmann(arr::Array{Float64}, alpha::Float64, dims::Union{Int64, Tuple})
    return sum(arr .* exp.(alpha .* arr), dims=dims) ./ sum(exp.(alpha .* arr), dims=dims)
end

end