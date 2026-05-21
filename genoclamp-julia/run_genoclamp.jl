using Printf
import JSON
using XML
include("./engine.jl")
using .Engine: GreedyEngine, GreedyCluster, GreedyItem, cluster_motifs!
include("./utils.jl")
using .Utils: highest_n_info_sum, check_periodicity, trim_motif, xlog2x
using ArgParse

### Input parsing
function parse_meme_files(meme_files::Vector{String}; get_sites::Bool = true,
    nsites_pattern::Regex = r"letter-probability matrix: alength= \d+ w= (\d+) nsites= (\d+) E= ([\d.+e-]+)",
    weights_pattern::Regex = r"\s*([\d\.]+)\s*([\d\.]+)\s*([\d\.]+)\s*([\d\.]+)",
    id_pattern::Regex = r"MOTIF\s+(.+?)\s")::Vector{GreedyItem}

    items = Vector{GreedyItem}()
    sites = Dict{Int64, Set{Tuple{String, Int64, Int64, Char}}}()
    n = 0
    for fn in meme_files
        file_n::Int64 = 1
        width::Int64 = 0
        open(fn, "r") do io
            while !eof(io)
                line = readline(io)
                id_match = match(id_pattern, line)
                if id_match !== nothing
                    motif_id = id_match[1]
                end
                if endswith(strip(line), "sites sorted by position p-value") && get_sites
                    readline(io)
                    readline(io)
                    readline(io)
                    line = readline(io)
                    sites[n] = Set{Tuple{String, Int64, Int64, Char}}()
                    while startswith(line, "chr")
                        region = split(strip(line), "(")[1]
                        offset = parse(Int64, line[31:37])
                        chrom, coords = split(region, ":")
                        start = parse(Int64, split(coords, "-")[1]) + offset - 1
                        stop = start + width
                        push!(sites[n], (chrom, start, stop, line[30]))
                        line = readline(io)
                    end
                end

                nsites_match = match(nsites_pattern, line)
                if nsites_match !== nothing
                    width = parse(Int64, nsites_match[1])
                    nsites = parse(Int64, nsites_match[2])
                    evalue = parse(Float64, nsites_match[3])
                    pfm = Matrix{Float64}(undef, width, 4)
                    line = readline(io)
                    weights_match = match(weights_pattern, line)
                    i = 1
                    while weights_match !== nothing
                        for j in 1:4
                            pfm[i, j] = parse(Float64, weights_match[j]) * nsites
                        end
                        line = readline(io)
                        weights_match = match(weights_pattern, line)
                        i += 1
                    end
                    source = (fn * "-motif" * string(file_n), nsites, evalue)
                    push!(items, GreedyItem(n, pfm, source, get_sites ? sites[n] : Set{Tuple{String, Int64, Int64, Char}}()))
                    file_n += 1
                    n += 1
                end
            end
        end
    end
    
    return items
end

### Output files
genoclamp_dir = dirname(@__DIR__)
glyph_data = open("$genoclamp_dir/logo_symbols/glyphs.json", "r") do io
    JSON.parse(io)
end
symbol_library = open("$genoclamp_dir/logo_symbols/symbol_library.json", "r") do io
    JSON.parse(io)
end

function write_aligned_transfac(cluster::GreedyCluster, filename::String)
    open(filename, "w") do io
        for i in eachindex(cluster.items)
            motif_id = cluster.items[i].source[1]
            @printf io "AC\t%s\n" motif_id
            @printf io "XX\n"
            @printf io "ID\t%s\n" motif_id
            @printf io "P0\tA\tC\tG\tT\n"
            pfm = cluster.aligned_pfms[i, :, :]
            for j in axes(pfm, 1)
                @printf io "%02d\t%06f\t%06f\t%06f\t%06f\n" j pfm[j, 1] pfm[j, 2] pfm[j, 3] pfm[j, 4]
            end
            @printf io "XX\n//\n"
        end
    end
end

function write_consensus_transfac(cluster::GreedyCluster, filename::String; info_thresh::Float64 = .5)
    c = cluster.idx
    trimmed_pfm, _, _, _ = trim_motif(cluster.aligned_pfms; info_thresh=info_thresh)
    open(filename, "w") do io
        @printf io "AC\t%s%d\n" "cluster" c
        @printf io "XX\n"
        @printf io "ID\t%s%d\n" "cluster" c
        @printf io "P0\tA\tC\tG\tT\n"
        for j in axes(trimmed_pfm, 1)
            @printf io "%02d\t%06f\t%06f\t%06f\t%06f\n" j trimmed_pfm[j, 1] trimmed_pfm[j, 2] trimmed_pfm[j, 3] trimmed_pfm[j, 4]
        end
        @printf io "XX\n//\n"
    end
end

function write_consensus_meme(cluster::GreedyCluster, filename::String; info_thresh::Float64 = .5)
    c = cluster.idx
    trimmed_pfm, _, _, _ = trim_motif(cluster.aligned_pfms; info_thresh=info_thresh)
    open(filename, "w") do io
        @printf io "MEME version 5\n"
        @printf io "ALPHABET= ACGT\n"
        @printf io "MOTIF %s%d %s%d\n" "cluster" c "cluster" c
        @printf io "letter-probability matrix: alength= 4 w= %d nsites= %d E= 0\n" size(trimmed_pfm, 1) sum(item.source[2] for item in cluster.items)
        for j in axes(trimmed_pfm, 1)
            @printf io "%f %f %f %f\n" trimmed_pfm[j, 1] trimmed_pfm[j, 2] trimmed_pfm[j, 3] trimmed_pfm[j, 4]
        end
    end
end

symbol_sets = Dict{Symbol, Tuple{Float64, Vector{Tuple{String, String}}}}()
symbol_sets[:DNA] = (
    symbol_library["DNA"]["max_bits"],
    [(glyph_data[s["name"]]["path"], s["color"]) for s in symbol_library["DNA"]["symbols"]]
)
symbol_sets[:RNA] = (
    symbol_library["RNA"]["max_bits"],
    [(glyph_data[s["name"]]["path"], s["color"]) for s in symbol_library["RNA"]["symbols"]]
)
symbol_sets[:AA] = (
    symbol_library["AA"]["max_bits"],
    [(glyph_data[s["name"]]["path"], s["color"]) for s in symbol_library["AA"]["symbols"]]
)

function plot_logo_stack(aligned_pfms::Array{Float64, 3}; symbol::Symbol = :DNA, glyph_width::Float64 = 100., stack_height::Float64 = 200.)
    height, width, _ = size(aligned_pfms)
    svg = XML.Element("svg"; baseProfile="full", version="1.1", xmlns="http://www.w3.org/2000/svg", viewBox="0 0 $(width * glyph_width) $(height * stack_height)")

    logo_stack = XML.Element("g")
    push!(svg, logo_stack)

    for y in 1:height
        pfm = aligned_pfms[y, :, :]
        row = XML.Element("g"; transform="translate(0 $((y - 1) * stack_height))")
        push!(logo_stack, row)
        for i in 1:width
            pwv = pfm[i, :]
            n = sum(pwv)
            if n == 0.
                continue
            end
            vec = xlog2x.(pwv ./ n)
            bits = sum(vec) + symbol_sets[symbol][1]
            heights = pwv .* bits .* (stack_height / (n * symbol_sets[symbol][1]))
            idx = sortperm(pwv)

            stack = XML.Element("g"; transform="translate($((i - 1) * glyph_width) 0)")
            push!(row, stack)

            y_offset = 0.
            for j in idx
                if heights[j] > 0.
                    base = symbol_sets[symbol][2][j]
                    y_offset += heights[j]
                    glyph = XML.Element("path"; d=base[1], fill=base[2], transform="matrix($(glyph_width / 100.) 0 0 $(heights[j] / 100.) 0 $(stack_height - y_offset))")
                    push!(stack, glyph)
                end
            end
        end
    end
    
    return svg
end

function write_bed_file(cluster::GreedyCluster, filename::String)
    open(filename, "w") do io
        for ((chrom, start, stop, strand), sources) in cluster.sites
            @printf io "%s\t%d\t%d\t%s\t.\t%s\n" chrom start stop join(sources, ";") strand
        end
    end
end

### Motif filtering
function filter_motifs(items::Vector{GreedyItem}; nsites_thresh::Int64 = 10,
        evalue_thresh::Float64 = .01, info_score_thresh::Float64 = 5.,
        periodicity1_thresh::Float64 = .6, periodicity2_thresh::Float64 = .75,
        periodicity3_thresh::Float64 = .75)::Vector{GreedyItem}
    function passes_filters(item::GreedyItem)::Bool
        nsites = item.source[2]
        evalue = item.source[3]
        pfm = item.pfm
        info_score = highest_n_info_sum(pfm)
        periodicity1 = check_periodicity(pfm, p=1)
        periodicity2 = check_periodicity(pfm, p=2)
        periodicity3 = check_periodicity(pfm, p=3)
        return nsites >= nsites_thresh && evalue <= evalue_thresh &&
            info_score >= info_score_thresh && periodicity1 <= periodicity1_thresh &&
            periodicity2 <= periodicity2_thresh && periodicity3 <= periodicity3_thresh
    end
    return filter(passes_filters, items)
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--meme"
        arg_type = String
        nargs = '+'
        help = "MEME files to process"
        "--meme-list"
        arg_type = String
        help = "File containing list of MEME files to process"
        "--nsites-thresh"
        arg_type = Int64
        default = 10
        help = "Minimum number of sites to consider a motif"
        "--evalue-thresh"
        arg_type = Float64
        default = .01
        help = "Maximum E-value to consider a motif"
        "--info-score-thresh"
        arg_type = Float64
        default = 5.
        help = "Minimum information score to consider a motif"
        "--periodicity1-thresh"
        arg_type = Float64
        default = .6
        help = "Maximum periodicity (p=1) to consider a motif"
        "--periodicity2-thresh"
        arg_type = Float64
        default = .75
        help = "Maximum periodicity (p=2) to consider a motif"
        "--periodicity3-thresh"
        arg_type = Float64
        default = .75
        help = "Maximum periodicity (p=3) to consider a motif"
        "--pc"
        arg_type = Float64
        nargs = '+'
        default = [2., 2., 2., 2.]
        help = "Alpha parameters for Dirichlet prior"
        "--min-base-overlap"
        arg_type = Int64
        default = 4
        help = "Minimum number of bases to overlap for merging"
        "--min-information-overlap"
        arg_type = Float64
        default = 0.
        help = "Minimum information overlap for merging"
        "--max-information-overhang"
        arg_type = Float64
        default = 12.
        help = "Maximum information overhang for merging"
        "--concentration"
        arg_type = Float64
        default = .5
        help = "Concentration parameter for merging"
        "--trim-thresh"
        arg_type = Float64
        default = .5
        help = "Information threshold for trimming motifs"
        "--get-sites"
        action = :store_true
        help = "Whether to extract binding sites from MEME files"
        "--output-dest", "-o"
        arg_type = String
        default = "genoclamp_out"
        help = "Folder to save results, will be created if it does not exist"
    end
    args = parse_args(s)

    if length(args["meme"]) > 0
        if args["meme-list"] != nothing
            println("Warning: both --meme and --meme-list specified, ignoring --meme-list")
        end
        meme_files = args["meme"]
    elseif args["meme-list"] != nothing
        meme_files = readlines(args["meme-list"])
    else
        error("Must specify either --meme or --meme-list")
    end

    items = parse_meme_files(meme_files; get_sites=args["get-sites"])
    items = filter_motifs(items; nsites_thresh=args["nsites-thresh"],
        evalue_thresh=args["evalue-thresh"], info_score_thresh=args["info-score-thresh"],
        periodicity1_thresh=args["periodicity1-thresh"],
        periodicity2_thresh=args["periodicity2-thresh"],
        periodicity3_thresh=args["periodicity3-thresh"])
    engine = GreedyEngine(items; pc=args["pc"], min_base_overlap=args["min-base-overlap"],
        min_information_overlap=args["min-information-overlap"],
        max_information_overhang=args["max-information-overhang"],
        concentration=args["concentration"])
    cluster_motifs!(engine)

    if !ispath(args["output-dest"])
        mkdir(args["output-dest"])
    end

    maximal_clusters = engine.clusters_trace[argmax(engine.llr_trace)]

    for c in maximal_clusters
        clust = engine.clusters[c]

        cluster_output_dir = args["output-dest"] * "/cluster" * string(c - 1)
        if !ispath(cluster_output_dir)
            mkdir(cluster_output_dir)
        end

        write_aligned_transfac(clust, cluster_output_dir * "/cluster" * string(c - 1) * "_aligned-motifs.transfac")
        write_consensus_transfac(clust, cluster_output_dir * "/cluster" * string(c - 1) * "_consensus-motif.transfac"; info_thresh=args["trim-thresh"])
        write_consensus_meme(clust, cluster_output_dir * "/cluster" * string(c - 1) * "_consensus-motif.meme"; info_thresh=args["trim-thresh"])

        svg = plot_logo_stack(clust.aligned_pfms)
        open(cluster_output_dir * "/cluster" * string(c - 1) * "_aligned-motifs.svg", "w") do io
            XML.write(io, svg)
        end

        trimmed_pfm, _, _, _ = trim_motif(clust.aligned_pfms; info_thresh=args["trim-thresh"])
        svg = plot_logo_stack(reshape(trimmed_pfm, (1, size(trimmed_pfm)...)))
        open(cluster_output_dir * "/cluster" * string(c - 1) * "_consensus-motif.svg", "w") do io
            XML.write(io, svg)
        end

        if args["get-sites"]
            write_bed_file(clust, cluster_output_dir * "/cluster" * string(c - 1) * "_binding-sites.bed")
        end
    end
end

main()