# Simulate motifs

## Getting started
First simulate the base pwms:
```
python simulate_base_pwms.py --entr-dist cisbp_entropy.txt --seed {seed} --out {base_pwm_filename}
```

Then simulate the motifs:
```
python simulate_motifs.py --out {output_pfm_filename} --seed {seed} {base_pwm_filename}
```

## Flags
### simulate_base_pwms.py
 - `--entr-dist`: Text file with nucleotide entropy distribution to sample from (required)
 - `--n`: Number of base PWMs to simulate
 - `--a-default`: Relative number of default (non-palindromic and non-periodic) motifs to simulate
 - `--a-palindromic`: Relative number of palindromic motifs to simulate
 - `--a-periodic`: Relative number of periodic motifs to simulate
   - Proportion of each category of motif is sampled from a Dirichlet Distribution using `--a-default`, `--a-palindromic`, and `--a-periodic` as parameters
 - `--min-width`: Minimum motif width
 - `--max-width`: Maximum motif width
 - `--min-logits`: Minimum total information for each PWM in logits
 - `--bg-alpha`: Pseudocounts for background nucleotide frequencies
 - `--min-period-width`: Minimum period width for periodic motfs
 - `--max-period-width`: Maximum period width for periodic motfs
 - `--min-period-stability`: Minimum period stability for periodic motifs
 - `--max-period-stability`: Maxmimum period stability for periodic motifs
   - Period stability is the inverse of the rate of information degeneration near the peripheries
 - `--seed`: Random seed
 - `--out`: Output file path

### simulate_motifs.py
 - The base PWMs output from `simulate_base_pwms.py` is a required positional output
 - `--exp-n`: Expected number of PFMs for each base PWM (Poisson)
 - `--min-width`: Minimum width for each PFM
 - `--min-sites`: Minimum number of sites for each PFM
 - `--max-sites`: Maximum number of sites for each PFM
 - `--padding-p`: Geometric probability of padding nucleotides for PFMs of periodic motifs
 - `--bg-alpha`: Pseudocounts for background nucleotide frequencies
 - `--seed`: Random seed
 - `--out`: Output file path