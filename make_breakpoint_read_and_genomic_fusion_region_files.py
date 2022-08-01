from genomic_bead_capture_translocation import SplitReadFinder, BreakpointReadSetListFinder
import argparse

## add the arguments
parser = argparse.ArgumentParser(description='Get a file of breakpoint reads and a file of genomic fusion regions for each breakpoint in a bead capture experiment.')
parser.add_argument("outfile_prefix", help="prefix for the output files",type=str)
parser.add_argument("samfile_path", help="path to the samfile of the alignment of the bead capture reads to the human genome",type=str)
parser.add_argument("fastq_file_path", help="path to the fastq file containing all the reads of the bead capture experiment",type=str)
parser.add_argument("bedfile_path", help="path to the bedfile specifying the probed regions, type=str")
parser.add_argument("-mro", "--max_read_overlap",default=500, help="maximum overlap in the read segments that are aligned to different translocation partners to be considered a split read",type=int)
parser.add_argument("-brt", "--breakpoint_read_thresh",default=5, help="lower bound for the number of reads at a specific breakpoint for it to be considered a potential translocation",type=int)
parser.add_argument("-mmd","--max_merge_dist",default=20, help="the maximum distance away from a breakpoint with at least thresh reads to consider merging other breakpoints",type=int)
parser.add_argument("-smt","--shared_read_merge_thresh",default=0.05, help="upper bound for the fraction of shared reads between two breakpoint read sets for them to be merged",type=float)
parser.add_argument("-pmt","--post_merge_thresh",default=8, help="lower bound for the number of reads at a specific post-merged breakpoint for it to be considered a potential translocation",type=int)
parser.add_argument("-mft","--max_brs_read_fraction_thresh",default=0.1, help="",type=float)
args = parser.parse_args()

srf = SplitReadFinder(args.samfile_path, args.max_read_overlap)
brf = BreakpointReadSetListFinder(srf.SplitReadAlignment_list, args.breakpoint_read_thresh, args.max_merge_dist, args.shared_read_merge_thresh, args.post_merge_thresh, args.max_brs_read_fraction_thresh, args.bedfile_path)
brf.make_all_BRS_fasta_files(1000,outfile_path_prefix=args.outfile_prefix,id_prefix=args.outfile_prefix)
brf.make_all_BRS_read_fastq_files(args.fastq_file_path,args.outfile_prefix)
