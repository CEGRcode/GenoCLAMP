include("./engine.jl")
using .Engine: GreedyEngine, GreedyItem
include("./utils.jl")
using .Utils: highest_n_info_sum, check_periodicity, trim_motif
include("./input.jl")
using .Input: parse_meme_files
include("./output.jl")
using .Output
using ArgParse

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
    @add_arg_table s begin
        "--meme"
        arg_type = String
        nargs = "+"
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
        nargs = "+"
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
        "--n-workers"
        arg_type = Int64
        help = "Number of workers to use for parallelization"
        "--trim-thresh"
        arg_type = Float64
        default = .5
        help = "Information threshold for trimming motifs"
        "--get-sites"
        action = :store_true
        help = "Whether to extract binding sites from MEME files"
    end
end