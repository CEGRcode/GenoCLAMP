module Engine

using SpecialFunctions: loggamma
using Base.Threads
include("./utils.jl")
using .Utils: boltzmann

struct GreedyItem
    idx::Int64
    pfm::Matrix{Float64}
    revcomp::Matrix{Float64}
    width::Int64
    source::Tuple{String, Int64, Float64}
    sites::Set{Tuple{String, Int64, Int64, String}}
end
function GreedyItem(idx::Int64, pfm::Matrix{Float64},
        source::Tuple{String, Int64, Float64},
        sites::Set{Tuple{String, Int64, Int64, String}})::GreedyItem
    return GreedyItem(idx, pfm, reverse(pfm), size(pfm, 1), source, sites)
end

struct GreedyCluster
    idx::Int64
    items::Vector{GreedyItem}
    aligned_pfms::Array{Float64, 3}
    width::Int64
    llr::Float64
    min_bits::Vector{Float64}
    bits::Vector{Float64}
    sites::Dict{Tuple{String, Int64, Int64, String}, Set{String}}
    merged_from::Union{Nothing, Tuple{Int64, Int64}}
end
function GreedyCluster(idx::Int64, items::Vector{GreedyItem},
        aligned_pfms::Array{Float64, 3}, llr::Float64,
        sites::Dict{Tuple{String, Int64, Int64, String}, Set{String}};
        merged_from::Union{Nothing, Tuple{Int64, Int64}} = nothing)::GreedyCluster
    _, width, _ = size(aligned_pfms)

    aligned_pfms_eps = aligned_pfms .+ 1e-20
    aligned_posterior_pwms = aligned_pfms_eps ./ sum(aligned_pfms_eps, dims=3)
    min_bits = dropdims(boltzmann(dropdims(sum(aligned_posterior_pwms .* log2.(aligned_posterior_pwms), dims=3), dims=3) .+ 2., -2., 1), dims=1)
    
    consensus_pwm = dropdims(sum(aligned_pfms_eps, dims=1) ./ sum(aligned_pfms_eps, dims=(1, 3)), dims=1)
    bits = dropdims(sum(consensus_pwm .* log2.(consensus_pwm), dims=2), dims=2) .+ 2
    
    return GreedyCluster(idx, items, aligned_pfms, width, llr,
        min_bits, bits, sites, merged_from)
end

struct GreedyEngine
    items::Vector{GreedyItem}
    clusters::Vector{GreedyCluster}
    alphabet_length::Int64
    pc::Vector{Float64}
    lgpc::Vector{Float64}
    pc_sum::Float64
    lga::Float64
    min_base_overlap::Int64
    min_information_overlap::Float64
    max_information_overhang::Float64
    concentration::Float64
    clusters_trace::Vector{Vector{Int64}}
    llr_trace::Vector{Float64}
    cache::Dict{Tuple{Int64, Int64}, Tuple{Array{Float64, 3}, Float64, Float64, Float64, Float64, Int64, Int64, Int64, Int64, Bool}}
end
function GreedyEngine(items::Vector{GreedyItem}; alphabet_length::Int64 = 4,
        pc::Vector{Float64} = [2., 2., 2., 2.], min_base_overlap::Int64 = 4,
        min_information_overlap::Float64 = 8., max_information_overhang::Float64 = 12.,
        concentration::Float64 = .5)::GreedyEngine
    clusters = Vector{GreedyCluster}(undef, length(items))
    for (idx, item) in enumerate(items)
        if size(item.pfm, 2) != alphabet_length
            throw(ArgumentError("alphabet_length does not match size of pfm(s)"))
        end
        sites = Dict{Tuple{String, Int64, Int64, String}, Set{String}}()
        for site in item.sites
            sites[site] = Set([item.source[1]])
        end
        clusters[idx] = GreedyCluster(idx, [item], reshape(item.pfm, 1, size(item.pfm)...), 0., sites)
    end
    lgpc = loggamma.(pc)
    pc_sum = sum(pc)
    lga = loggamma(pc_sum)
    return GreedyEngine(items, clusters, alphabet_length, pc, lgpc, pc_sum, lga, min_base_overlap,
        min_information_overlap, max_information_overhang, concentration, [Vector{Int64}(1:length(items))], [0.],
        Dict{Tuple{Int64, Int64}, Tuple{Array{Float64, 3}, Float64, Float64, Float64, Float64, Int64, Int64, Int64, Int64, Bool}}())
end

function compute_llr(aligned_pfms::Array{Float64, 3}, pc::Vector{Float64},
        lgpc::Vector{Float64}, pc_sum::Float64, lga::Float64)::Float64
    npfms, width, alphabet_length = size(aligned_pfms)
    complement_pfm = dropdims(sum(aligned_pfms, dims=1), dims=1)
    previous_pfm = zeros(Float64, width, alphabet_length)
    llr = 0.
    for i in 1:npfms
        current_pfm = aligned_pfms[i, :, :]
        complement_pfm .+= previous_pfm .- current_pfm
        previous_pfm = current_pfm
        for j in 1:width
            current_sum = sum(current_pfm[j, :])
            complement_sum = sum(complement_pfm[j, :])
            llr += loggamma(current_sum + pc_sum) + loggamma(complement_sum + pc_sum) - loggamma(current_sum + complement_sum + pc_sum) - lga
            for b in 1:alphabet_length
                current_count = current_pfm[j, b]
                complement_count = complement_pfm[j, b]
                llr += loggamma(current_count + complement_count + pc[b]) + lgpc[b] - loggamma(current_count + pc[b]) - loggamma(complement_count + pc[b])
            end
        end
    end
    return llr
end

function compute_maximal_llr(engine::GreedyEngine, c1::Int64, c2::Int64)::Tuple{Array{Float64, 3}, Float64, Float64, Float64, Float64, Int64, Int64, Int64, Int64, Bool}
    cluster1 = engine.clusters[c1]
    width1 = cluster1.width
    n1 = length(cluster1.items)
    bits1 = cluster1.bits
    min_bits1 = cluster1.min_bits

    cluster2 = engine.clusters[c2]
    width2 = cluster2.width
    n2 = length(cluster2.items)
    bits2 = cluster2.bits
    min_bits2 = cluster2.min_bits
    bits2_reverse = reverse(cluster2.bits)
    min_bits2_reverse = reverse(cluster2.min_bits)
    cluster2_reverse_pfms = reverse(cluster2.aligned_pfms, dims=(2, 3))
    
    min_width = min(width1, width2)

    maximal_llr::Float64 = -Inf
    aligned_pfms::Array{Float64, 3} = zeros(Float64, 0, 0, 0)
    left_offset1::Int64 = 0
    right_offset1::Int64 = 0
    left_offset2::Int64 = 0
    right_offset2::Int64 = 0
    rc::Bool = false
    for i in 1:width1 + width2 - 1
        start1 = max(i - width1 + 1, 1)
        start2 = max(width1 - i + 1, 1)
        overlap = min(i, width1 + width2 - i, min_width)

        info_overlap_forward = sum(min_bits1[start2:start2 + overlap - 1] .*
            min_bits2[start1:start1 + overlap - 1])
        info_overhang_forward = sum(bits1[1:start2 - 1]) + sum(bits1[start2 + overlap:end]) +
            sum(bits2[1:start1 - 1]) + sum(bits2[start1 + overlap:end]) +
            sum(abs.(bits1[start2:start2 + overlap - 1] - bits2[start1:start1 + overlap - 1]))
        
        combined_width = width1 + max(i, width2) - min(i, width1)

        if info_overlap_forward >= engine.min_information_overlap &&
            info_overhang_forward <= engine.max_information_overhang
            combined_pfms = zeros(Float64, n1 + n2, combined_width, engine.alphabet_length)
            combined_pfms[1:n1, start1:start1 + width1 - 1, :] = cluster1.aligned_pfms
            combined_pfms[n1 + 1:n1 + n2, start2:start2 + width2 - 1, :] = cluster2.aligned_pfms
            potential_llr = compute_llr(combined_pfms, engine.pc, engine.lgpc, engine.pc_sum, engine.lga)
            if potential_llr > maximal_llr
                maximal_llr = potential_llr
                aligned_pfms = combined_pfms
                left_offset1 = start1 - 1
                right_offset1 = combined_width - width1 - start1 + 1
                left_offset2 = start2 - 1
                right_offset2 = combined_width - width2 - start2 + 1
                rc = false
            end
        end

        info_overlap_reverse = sum(min_bits1[start2:start2 + overlap - 1] .*
            min_bits2_reverse[start1:start1 + overlap - 1])
        info_overhang_reverse = sum(bits1[1:start2 - 1]) + sum(bits1[start2 + overlap:end]) +
            sum(bits2_reverse[1:start1 - 1]) + sum(bits2_reverse[start1 + overlap:end]) +
            sum(abs.(bits1[start2:start2 + overlap - 1] - bits2_reverse[start1:start1 + overlap - 1]))
        
        if info_overlap_reverse >= engine.min_information_overlap &&
            info_overhang_reverse <= engine.max_information_overhang
            combined_pfms = zeros(Float64, n1 + n2, combined_width, engine.alphabet_length)
            combined_pfms[1:n1, start1:start1 + width1 - 1, :] = cluster1.aligned_pfms
            combined_pfms[n1 + 1:n1 + n2, start2:start2 + width2 - 1, :] = cluster2_reverse_pfms
            potential_llr = compute_llr(combined_pfms, engine.pc, engine.lgpc, engine.pc_sum, engine.lga)
            if potential_llr > maximal_llr
                maximal_llr = potential_llr
                aligned_pfms = combined_pfms
                left_offset1 = start1 - 1
                right_offset1 = combined_width - width1 - start1 + 1
                left_offset2 = start2 - 1
                right_offset2 = combined_width - width2 - start2 + 1
                rc = true
            end
        end
    end

    scaled_llr = maximal_llr * (n1 + n2) ^ engine.concentration
    scaled_llr1 = cluster1.llr * n1 ^ engine.concentration
    scaled_llr2 = cluster2.llr * n2 ^ engine.concentration
    return aligned_pfms, maximal_llr, maximal_llr - cluster1.llr - cluster2.llr,
        scaled_llr, scaled_llr - scaled_llr1 - scaled_llr2, left_offset1,
        right_offset1, left_offset2, right_offset2, rc
end

function one_iteration!(engine::GreedyEngine, lk::ReentrantLock)
    current_clusters = copy(engine.clusters_trace[end])
    nclusters = length(current_clusters)
    uncached_combos = Vector{Tuple{Int64, Int64}}(undef, div(nclusters * (nclusters - 1), 2) - length(engine.cache))
    cnt = 1
    for i in 1:nclusters - 1
        for ii in i + 1:nclusters
            combo = (current_clusters[i], current_clusters[ii])
            if !(combo in keys(engine.cache))
                uncached_combos[cnt] = combo
                cnt += 1
            end
        end
    end

    Threads.@threads for j in 1:length(uncached_combos)
        res = compute_maximal_llr(engine, uncached_combos[j]...)
        lock(lk) do
            engine.cache[uncached_combos[j]] = res
        end
    end

    function compare_keys(a::Tuple{Int64, Int64}, c::Tuple{Int64, Int64})::Tuple{Int64, Int64}
        if engine.cache[c][3] < 0
            return a
        end
        return engine.cache[c][5] > engine.cache[a][5] ? c : a
    end

    c1, c2 = foldl(compare_keys, keys(engine.cache))
    aligned_pfms, llr, _, _, _, left_offset1, right_offset1, left_offset2, right_offset2, rc = engine.cache[(c1, c2)]
    if isinf(llr)
        return false
    end

    popidx1 = searchsortedfirst(current_clusters, c1)
    deleteat!(current_clusters, popidx1)
    deleteat!(current_clusters, findnext(c -> c == c2, current_clusters, popidx1))
    delete!(engine.cache, (c1, c2))
    for c in current_clusters
        delete!(engine.cache, (c, c1))
        delete!(engine.cache, (c1, c))
        delete!(engine.cache, (c, c2))
        delete!(engine.cache, (c2, c))
    end
    
    cluster1 = engine.clusters[c1]
    cluster2 = engine.clusters[c2]
    cluster_idx = length(engine.clusters) + 1
    sites = Dict{Tuple{String, Int64, Int64, String}, Set{String}}()
    for site in keys(cluster1.sites)
        chrom, start, stop, strand = site
        if strand == '+'
            expanded_site = (chrom, start - left_offset1, stop + right_offset1, '+')
        else
            expanded_site = (chrom, start - right_offset1, stop + left_offset1, '-')
        end
        if haskey(sites, expanded_site)
            union!(sites[expanded_site], cluster1.sites[site])
        else
            sites[expanded_site] = copy(cluster1.sites[site])
        end
    end
    for site in keys(cluster2.sites)
        chrom, start, stop, strand = site
        if (strand == '+' && rc) || (strand == '-' && !rc)
            expanded_site = (chrom, start - right_offset2, stop + left_offset2, '-')
        else
            expanded_site = (chrom, start - left_offset2, stop + right_offset2, '+')
        end
        if haskey(sites, expanded_site)
            union!(sites[expanded_site], cluster2.sites[site])
        else
            sites[expanded_site] = copy(cluster2.sites[site])
        end
    end

    cluster = GreedyCluster(cluster_idx, vcat(cluster1.items, cluster2.items), aligned_pfms, llr, sites, merged_from=(c1, c2))
    push!(engine.clusters, cluster)
    push!(current_clusters, cluster_idx)
    push!(engine.clusters_trace, current_clusters)
    push!(engine.llr_trace, sum(ck -> engine.clusters[ck].llr * length(engine.clusters[ck].items) ^ engine.concentration, current_clusters))

    return true
end

function cluster_motifs!(engine::GreedyEngine)
    n_iter = length(engine.clusters_trace[end]) - 1
    lk = ReentrantLock()
    for i in 1:n_iter
        print("\r", i, "/", n_iter, "\t")
        if !(one_iteration!(engine, lk))
            println("\nNo more valid merges... done")
            break
        end
    end
    sizehint!(engine.cache, 0)
    println()
end

end