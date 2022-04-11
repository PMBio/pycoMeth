from typing import IO
from io import StringIO

import numpy as np
import pandas as pd

from pycoMeth.meth_seg.math import bs_from_llrs, compute_differential_methylation


class SegmentsWriterBED:
    def __init__(self, outfile: IO, chrom: str):
        self.outfile = outfile
        self.chrom = chrom
        self.first = True
    
    def write_segments_llr(
        self,
        llrs: np.ndarray,
        segments: np.ndarray,
        genomic_starts: np.ndarray,
        genomic_ends: np.ndarray,
        samples: np.ndarray,
        compute_diffmet: bool = False,
        **kwargs,
    ):
        df_rowvals = []
        
        for seg in sorted(list(set(segments))):
            seg_pos = np.arange(llrs.shape[1])[segments == seg]
            start = genomic_starts[seg_pos[0]]
            end = genomic_ends[seg_pos[-1]]
            
            rowval = {"chrom": self.chrom, "start": start, "end": end, "num_sites": (segments == seg).sum()}
            if samples is not None and compute_diffmet:
                samples_unique = list(set(samples))
                seg_llrs = llrs[:, segments == seg]
                sample_llr = {sample: seg_llrs[samples == sample] for sample in samples_unique}
                for sample in samples_unique:
                    rowval[f"met_rate_{sample}"] = bs_from_llrs((sample_llr[sample]))
                    
                    for i, s_a in enumerate(samples_unique):
                        for s_b in samples_unique[i + 1 :]:
                            up, pp = compute_differential_methylation(sample_llr[s_a], sample_llr[s_b])
                            rowval["unpaired_pval_%s_vs_%s" % (s_a, s_b)] = up
                            rowval["paired_pval_%s_vs_%s" % (s_a, s_b)] = pp
            
            df_rowvals.append(rowval)
        
        df = pd.DataFrame(df_rowvals, columns=["chrom", "start", "end", "num_sites"])
        df = df.astype({"chrom": str, "start": int, "end": int, "num_sites": int})
        
        df.to_csv(self.outfile, sep="\t", header=False, index=False, mode="w" if self.first else "a")
        self.first = False


class SegmentsWriterBedGraph:
    def __init__(self, outfile_base, chrom: str):
        self.outfiles = {}
        self.outfile_base = outfile_base
        self.chrom = chrom
    
    def get_outfile(self, sample):
        if sample not in self.outfiles:
            self.outfiles[sample] = f"{self.outfile_base}.{sample}.bedGraph"
            with open(self.outfiles[sample], "w") as f:
                f.write(
                    f"track type=bedGraph name={sample} description=center_label visibility=display_mode color=252,127,44 altColor=25,4,248 graphType=heatmap viewLimits=0:1 midRange=0.50:0.50 midColor=255,255,255\n"
                )
        return self.outfiles[sample]
    
    def write_segments_llr(
        self,
        llrs: np.ndarray,
        segments: np.ndarray,
        genomic_starts: np.ndarray,
        genomic_ends: np.ndarray,
        samples: np.ndarray,
        **kwargs,
    ):
        samples_unique = list(set(samples))
        for sample in samples_unique:
            with open(self.get_outfile(sample), "a") as f:
                for seg in sorted(list(set(segments))):
                    seg_pos = np.arange(llrs.shape[1])[segments == seg]
                    start = genomic_starts[seg_pos[0]]
                    end = genomic_ends[seg_pos[-1]]
                    seg_llrs = llrs[:, segments == seg]
                    sample_llr = seg_llrs[samples == sample]
                    sample_metrate = (sample_llr > 2).sum() / (np.abs(sample_llr) > 2).sum()
                    row = f"{self.chrom} {start} {end} {sample_metrate}\n"
                    f.write(row)


class SegmentsReaderCSV:
    """
    Reads a concatenation of CSV files (with potentially differenct columns) and merges
    them into one dataframe, filling in null for missing values. Accepts multiple header
    lines, but requires them to start with the word "chrom" in order to identify the
    header.

    This is used to read a CSV file that has been created in chunks by multiple
    concurrent worker threads.
    """
    
    def __init__(self, *argc):
        self.inputfiles = argc
    
    def read(self):
        ret = []
        for inputfile in self.inputfiles:
            with open(inputfile) as f:
                cur_string = None
                while True:
                    line = f.readline().strip()
                    if line.startswith("chrom") or not line:
                        # new header line
                        if cur_string is not None:
                            cur_pd = pd.read_csv(StringIO("\n".join(cur_string)), sep="\t")
                            ret.append(cur_pd)
                        cur_string = []
                    if not line:
                        break
                    cur_string.append(line)
        ret = pd.concat(ret, sort=False)
        ret = ret.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
        return ret
