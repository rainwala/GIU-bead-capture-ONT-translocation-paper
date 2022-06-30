from genomic_bead_capture_translocation import SplitReadFinder, BreakpointReadSetListFinder
import argparse

## add the arguments
parser = argparse.ArgumentParser(description='Get a file of breakpoint reads and a file of genomic fusion regions for each breakpoint in a bead capture experiment.')
parser.add_argument("outfile_prefix", help="prefix for the output files",type=str)
parser.add_argument("samfile_path", help="path to the samfile of the alignment of the bead capture reads to the human genome",type=str)
parser.add_argument("fastq_file_path", help="path to the fastq file containing all the reads of the bead capture experiment",type=str)
parser.add_argument("bedfile_path", help="path to the bedfile specifying the probed regions, type=str")
parser.add_argument("thresh", help="lower bound for the number of reads at a specific breakpoint for it to be considered a potential translocation",type=int)
parser.add_argument("max_merge_dist", help="the maximum distance away from a breakpoint with at least thresh reads to consider merging other breakpoints",type=int)
parser.add_argument("post_merge_thresh", help="lower bound for the number of reads at a specific post-merged breakpoint for it to be considered a potential translocation",type=int)
args = parser.parse_args()

srf = SplitReadFinder(args.samfile_path)
brf = BreakpointReadSetListFinder(srf.SplitReadAlignment_list, args.thresh, args.max_merge_dist, args.post_merge_thresh, args.bedfile_path)
for brs in sorted(brf.region_BreakpointReadSet_dict,key=lambda x: len(brf.region_BreakpointReadSet_dict[x].SplitReadAlignment_set)):
	print(brs,len(brf.region_BreakpointReadSet_dict[brs].SplitReadAlignment_set))
	reads_out_filepath = f"{args.outfile_prefix}{brs}_breakpoint_reads.fq"
	brf.region_BreakpointReadSet_dict[brs].make_selected_read_fastq_file(args.fastq_file_path,reads_out_filepath)
	brf.region_BreakpointReadSet_dict[brs].make_breakpoint_region_fasta_file(1000, outfile_path_prefix=args.outfile_prefix,id_prefix=args.outfile_prefix)
