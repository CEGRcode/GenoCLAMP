import re
import numpy as np
from engine import GreedyItem

def parse_meme_files(meme_files, get_sites=False,
    nsites_pattern = re.compile('letter-probability matrix: alength= \d+ w= (\d+) nsites= (\d+) E= ([\d.+e-]+)'),
    weights_pattern = re.compile('\s*([\d\.]+)\s*([\d\.]+)\s*([\d\.]+)\s*([\d\.]+)'),
    id_pattern = re.compile('MOTIF\s+(.+?)\s')):

    items = []
    sites = {}
    n = 0
    for fn in meme_files:
        with open(fn, 'r') as f:
            for line in f:
                id_match = id_pattern.match(line)
                if id_match:
                    motif_id = id_match.groups()[0]

                if line.strip().endswith('sites sorted by position p-value') and get_sites:
                    f.readline()
                    f.readline()
                    f.readline()
                    line = f.readline()
                    sites[n] = {'chrom': [], 'start': [], 'strand': []}
                    while line.startswith('chr'):
                        region = line.strip().split('(')[0]
                        offset = int(line[30:37])
                        chrom, coords = region.split(':')
                        start = int(coords.split('-')[0]) + offset - 1
                        sites[n]['chrom'].append(chrom)
                        sites[n]['start'].append(start)
                        sites[n]['strand'].append(line[29])
                        line = f.readline()

                nsites_match = nsites_pattern.match(line)
                if nsites_match:
                    width = int(nsites_match.groups()[0])
                    nsites = int(nsites_match.groups()[1])
                    evalue = float(nsites_match.groups()[2])
                    rows = []
                    line = f.readline()
                    weights_match = weights_pattern.match(line)
                    while weights_match:
                        rows.append([float(w) * nsites for w in weights_match.groups()])
                        line = f.readline()
                        weights_match = weights_pattern.match(line)
                    source = ('{}-{}'.format(fn, motif_id), nsites, evalue)
                    motif_sites = {(chrom, start, start + width, strand) for chrom, start, strand in zip(
                        sites[n]['chrom'], sites[n]['start'], sites[n]['strand'])} if get_sites else set()
                    items.append(GreedyItem(n, np.array(rows), source, motif_sites))
                    n += 1
    
    return items
