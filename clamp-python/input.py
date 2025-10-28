import re
import numpy as np
from engine import GreedyItem

def parse_meme_files(meme_files, get_sites=False,
    nsites_pattern = re.compile('letter-probability matrix: alength= \d+ w= \d+ nsites= (\d+) E= ([\d.+e-]+)'),
    weights_pattern = re.compile('\s*([\d\.]+)\s*([\d\.]+)\s*([\d\.]+)\s*([\d\.]+)'),
    width_pattern = re.compile('MOTIF.+width =\s+(\d+)')):

    items = []
    sources = []
    sites = {}
    n = 0
    for fn in meme_files:
        file_n = 1
        with open(fn, 'r') as f:
            for line in f:
                nsites_match = nsites_pattern.match(line)
                width_match = width_pattern.match(line)
                if width_match:
                    width = int(width_match.groups()[0])
                if line.strip().endswith('sites sorted by position p-value') and get_sites:
                    f.readline()
                    f.readline()
                    f.readline()
                    line = f.readline()
                    sites[n] = set()
                    while line.startswith('chr'):
                        region = line.strip().split('(')[0]
                        offset = int(line[30:37])
                        chrom, coords = region.split(':')
                        start = int(coords.split('-')[0]) + offset - 1
                        stop = start + width
                        sites[n].add((chrom, start, stop, line[29]))
                        line = f.readline()

                if nsites_match:
                    nsites = int(nsites_match.groups()[0])
                    evalue = float(nsites_match.groups()[1])
                    rows = []
                    line = f.readline()
                    weights_match = weights_pattern.match(line)
                    while weights_match:
                        rows.append([float(w) * nsites for w in weights_match.groups()])
                        line = f.readline()
                        weights_match = weights_pattern.match(line)
                    source = ('{}-motif{}'.format(fn, file_n), nsites, evalue)
                    items.append(GreedyItem(n, np.array(rows), source, sites[n] if get_sites else set()))
                    sources.append(source)
                    file_n += 1
                    n += 1
    
    return items, sources, sites