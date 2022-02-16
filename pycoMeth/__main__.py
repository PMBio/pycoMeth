#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import argparse
import sys

# Local imports
import pycoMeth as pkg
from pycoMeth.common import *
from pycoMeth.Meth_Comp import Meth_Comp
from pycoMeth.CGI_Finder import CGI_Finder
from pycoMeth.Comp_Report import Comp_Report
from pycoMeth.Meth_Seg import Meth_Seg

# ~~~~~~~~~~~~~~TOP LEVEL ENTRY POINT~~~~~~~~~~~~~~#
def main(args=None):
    """ Main entry point for pycoMeth command line interface"""
    
    # Parser and subparsers for command
    parser = argparse.ArgumentParser(description=pkg.__description__)
    parser.add_argument("--version", action="version", version="{} v{}".format(pkg.__name__, pkg.__version__))
    subparsers = parser.add_subparsers(description="%(prog)s implements the following subcommands", dest="subcommands")
    subparsers.required = True

    # Meth_Comp subparser
    f = Meth_Comp
    sp_met = subparsers.add_parser("Meth_Comp", description=doc_func(f))
    sp_met.set_defaults(func=f)
    sp_met_io = sp_met.add_argument_group("Input/Output options")
    arg_from_docstr(sp_met_io, f, "h5_file_list", "i")
    arg_from_docstr(sp_met_io, f, "ref_fasta_fn", "f")
    arg_from_docstr(sp_met_io, f, "read_groups_key", "r")
    arg_from_docstr(sp_met_io, f, "interval_bed_fn", "a")
    arg_from_docstr(sp_met_io, f, "output_bed_fn", "b")
    arg_from_docstr(sp_met_io, f, "output_tsv_fn", "t")
    arg_from_docstr(sp_met_io, f, "interval_size", "n")
    arg_from_docstr(sp_met_io, f, "min_num_reads_per_interval", "c")
    sp_met_ms = sp_met.add_argument_group("Misc options")
    arg_from_docstr(sp_met_ms, f, "max_missing", "m")
    arg_from_docstr(sp_met_ms, f, "worker_processes", "w")
    arg_from_docstr(sp_met_ms, f, "min_diff_llr", "l")
    arg_from_docstr(sp_met_ms, f, "sample_id_list", "s")
    arg_from_docstr(sp_met_ms, f, "pvalue_adj_method")
    arg_from_docstr(sp_met_ms, f, "pvalue_threshold")
    arg_from_docstr(sp_met_ms, f, "only_tested_sites")
    arg_from_docstr(sp_met_ms, f, "hypothesis")
    arg_from_docstr(sp_met_ms, f, "do_independent_hypothesis_weighting")
    
    # Comp_Report subparser
    f = Comp_Report
    sp_cr = subparsers.add_parser("Comp_Report", description=doc_func(f))
    sp_cr.set_defaults(func=f)
    sp_cr_io = sp_cr.add_argument_group("Input/Output options")
    arg_from_docstr(sp_cr_io, f, "h5_file_list", "i")
    arg_from_docstr(sp_cr_io, f, "ref_fasta_fn", "f")
    arg_from_docstr(sp_cr_io, f, "read_groups_key", "r")
    arg_from_docstr(sp_cr_io, f, "methcomp_fn", "c")
    arg_from_docstr(sp_cr_io, f, "gff3_fn", "g")
    arg_from_docstr(sp_cr_io, f, "outdir", "o")
    sp_cr_ms = sp_cr.add_argument_group("Misc options")
    arg_from_docstr(sp_cr_ms, f, "sample_id_list", "s")
    arg_from_docstr(sp_cr_ms, f, "n_top", "n")
    arg_from_docstr(sp_cr_ms, f, "max_tss_distance", "d")
    arg_from_docstr(sp_cr_ms, f, "pvalue_threshold")
    arg_from_docstr(sp_cr_ms, f, "min_diff_llr")
    arg_from_docstr(sp_cr_ms, f, "min_diff_bs")
    arg_from_docstr(sp_cr_ms, f, "n_len_bin")
    arg_from_docstr(sp_cr_ms, f, "export_static_plots")
    arg_from_docstr(sp_cr_ms, f, "report_non_significant")
    
    # CGI_Finder subparser
    f = CGI_Finder
    sp_cgi = subparsers.add_parser("CGI_Finder", description=doc_func(f))
    sp_cgi.set_defaults(func=f)
    sp_cgi_io = sp_cgi.add_argument_group("Input/Output options")
    arg_from_docstr(sp_cgi_io, f, "ref_fasta_fn", "f")
    arg_from_docstr(sp_cgi_io, f, "output_bed_fn", "b")
    arg_from_docstr(sp_cgi_io, f, "output_tsv_fn", "t")
    sp_cgi_ms = sp_cgi.add_argument_group("Misc options")
    arg_from_docstr(sp_cgi_ms, f, "merge_gap", "m")
    arg_from_docstr(sp_cgi_ms, f, "min_win_len", "w")
    arg_from_docstr(sp_cgi_ms, f, "min_CG_freq", "c")
    arg_from_docstr(sp_cgi_ms, f, "min_obs_CG_ratio", "r")

    # Meth_Seg subparser
    f = Meth_Seg
    sp_seg = subparsers.add_parser("Meth_Seg", description=doc_func(f))
    sp_seg.set_defaults(func=f)
    sp_seg_io = sp_seg.add_argument_group("Input/Output options")
    arg_from_docstr(sp_seg_io, f, "h5_file_list", "i")
    arg_from_docstr(sp_seg_io, f, "chromosome", "c")
    arg_from_docstr(sp_seg_io, f, "chunks", "n")
    arg_from_docstr(sp_seg_io, f, "output_tsv_fn", "t")
    arg_from_docstr(sp_seg_io, f, "read_groups_keys", "r")
    arg_from_docstr(sp_seg_io, f, "chunk_size")
    sp_seg_ms = sp_seg.add_argument_group("Misc options")
    arg_from_docstr(sp_seg_ms, f, "max_segments_per_window", "m")
    arg_from_docstr(sp_seg_ms, f, "workers", "p")
    arg_from_docstr(sp_seg_ms, f, "reader_workers")
    arg_from_docstr(sp_seg_ms, f, "window_size", "w")
    arg_from_docstr(sp_seg_ms, f, "print_diff_met")

    # Add common group parsers
    for sp in [sp_met, sp_cr, sp_cgi]:
        sp_vb = sp.add_argument_group("Verbosity options")
        sp_vb.add_argument("-v", "--verbose", action="store_true", default=False, help="Increase verbosity")
        sp_vb.add_argument("-q", "--quiet", action="store_true", default=False, help="Reduce verbosity")
        sp_vb.add_argument("-p", "--progress", action="store_true", default=False, help="Display a progress bar")
    
    # Parse args and call subfunction
    args = parser.parse_args()
    
    args.func(**vars(args))


if __name__ == "__main__":
    main()
