import pysam
import sys
import re
from dataclasses import dataclass
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyfaidx

class SplitReadAlignment:
	""" class to hold a pair of pysam alignments corresponding to a single split read 
	indexed in order of chr first, then lb """
	
	def __init__(self,pysam_aln1,pysam_aln2):
			if (pysam_aln1.query_name != pysam_aln2.query_name) or (pysam_aln1 == pysam_aln2):
				return None
			self.alignments = self._get_alignment_tuple_from_input_alignments(pysam_aln1,pysam_aln2)
			self.read_id = self.alignments[0].query_name
			self.chr1 = self.alignments[0].reference_name
			self.chr2 = self.alignments[1].reference_name
			self.lb1 = self.alignments[0].reference_start
			self.lb2 = self.alignments[1].reference_start
			self.ub1 = self.alignments[0].reference_end
			self.ub2 = self.alignments[1].reference_end
			self.MAX_INDEL_DIST = 1000 # parameter by which we judge if the "translocation" is an indel or not

	def _get_alignment_tuple_from_input_alignments(self,pysam_aln1,pysam_aln2):
		""" using pysam_aln1, pysam_aln2 return a tuple ordered by the chr then by the lb """
		if pysam_aln1.reference_name < pysam_aln2.reference_name:
			return (pysam_aln1,pysam_aln2)
		if pysam_aln1.reference_name > pysam_aln2.reference_name:
			return (pysam_aln2,pysam_aln1)
		if pysam_aln1.reference_start < pysam_aln2.reference_start:
			return (pysam_aln1,pysam_aln2)
		return (pysam_aln2,pysam_aln1)

	def _get_alignment_tuple_from_input_alignments_ordered_by_probe_region_bedfile(self,pysam_aln1,pysam_aln2):
		pass

	def _are_alignments_small_indel(self):
		""" return True if the alignments in self.alignments are on the same chromosome and any of their co-ords are separated by less than self.MAX_INDEL_DIST, otherwise return False """
		if self.alignments[0].reference_name != self.alignments[1].reference_name:
			return False
		if (abs(self.alignments[0].reference_start - self.alignments[1].reference_start) < self.MAX_INDEL_DIST) or (abs(self.alignments[0].reference_end - self.alignments[1].reference_end) < self.MAX_INDEL_DIST):
			return True
		if (abs(self.alignments[0].reference_end - self.alignments[1].reference_start) < self.MAX_INDEL_DIST) or (abs(self.alignments[0].reference_end - self.alignments[1].reference_start) < self.MAX_INDEL_DIST):
			return True
		return False

	def __repr__(self):
		return ",".join(map(str,[self.read_id,self.chr1,self.lb1,self.ub1,self.chr2,self.lb2,self.ub2]))

class SplitReadFinder:
	""" methods to get all split reads as a list of SplitReadAlignment objects from a sam alignment file """

	def __init__(self,samfile_path,max_read_overlap):
		self.samfile_path = samfile_path
		self.MAX_READ_OVERLAP = max_read_overlap # the maximum allowed overlap amount in bases for the same part of a read to the main and supplemental alignments to consider this a split read 
		self.supp_alignments = self._get_supplementary_alignments()
		self.SplitReadAlignment_list = self._get_SplitReadAlignments()

	def _get_supplementary_alignments(self):
		""" return all the supplementary reads in the samfile """
		supp_alignments = defaultdict(list)
		for read in pysam.AlignmentFile(self.samfile_path).fetch():
			if read.is_supplementary:
				supp_alignments[read.query_name].append(read)	
		return supp_alignments

	def _is_pysam_alignment_on_canonical_chromosome(self,pysam_read):
		""" return True if the alignment if on one of the canonical chrs 1-22, or X, or Y """
		match = re.search(r'chr[A-Za-z0-9]+$',pysam_read.reference_name)
		if match is not None:	
			return True
		return False 

	def _is_read_alignment_region_overlapping_on_main_and_supplemental(self,pysam_aln,sup_aln):
		""" return True if the amount of overlap on the main and supplemental alignments exceeds self.MAX_READ_OVERLAP, False otherwise"""
		overlap = min(pysam_aln.query_alignment_end,sup_aln.query_alignment_end) - max(pysam_aln.query_alignment_start,sup_aln.query_alignment_start)
		#print(overlap,pysam_aln.reference_name,pysam_aln.reference_start,pysam_aln.reference_end,sup_aln.reference_name,sup_aln.reference_start,sup_aln.reference_end,pysam_aln.query_alignment_start,pysam_aln.query_alignment_end,sup_aln.query_alignment_start,sup_aln.query_alignment_end)
		return overlap > self.MAX_READ_OVERLAP

	def _get_SplitReadAlignments(self):
		""" return a list of SplitReadAlignment objects from the samfile """
		SplitReadAlignment_list = []
		for aln in pysam.AlignmentFile(self.samfile_path).fetch():
			if (not aln.is_supplementary) and (aln.query_name in self.supp_alignments) and (len(self.supp_alignments[aln.query_name]) >= 1):
				for sup_aln in self.supp_alignments[aln.query_name]:
					if (self._is_pysam_alignment_on_canonical_chromosome(aln)) and (self._is_pysam_alignment_on_canonical_chromosome(sup_aln)): 
						if not self._is_read_alignment_region_overlapping_on_main_and_supplemental(aln,sup_aln):
							SplitReadAlignment_list.append( SplitReadAlignment(aln,sup_aln) )
		return SplitReadAlignment_list


class Region:
	""" methods and data relating to a genomic region """
	
	def __init__(self,chr,lb,ub,strand):
		self.chr = chr
		self.lb = lb
		self.ub = ub
		self.strand = strand
		self.chr_fasta_filepath = "/home/aawan/SHARED/REFERENCE/hg38/Homo_sapiens_assembly38.fasta"

	@property
	def	seq(self):
		if self.strand is None:
			return None
		faidx_fasta_chrs = pyfaidx.Fasta(self.chr_fasta_filepath)
		if self.strand == '-':
			return faidx_fasta_chrs[self.chr][self.lb:self.ub].reverse.complement.seq	
		return faidx_fasta_chrs[self.chr][self.lb:self.ub].seq

	def overlaps(self,other):
		""" return True if self overlaps other Region, else False """
		if (self.chr == other.chr) and ((min(self.ub,other.ub) - max(self.lb,other.lb)) > 0):
			return True
		return False

	def __str__(self):
		rep = f"{self.chr}-{self.lb}-{self.ub}"
		if self.strand == '-':
			rep += '-rc'
		return rep 

class BreakpointReadSet:
	""" class to hold a set of SplitReadAlignment objects that map to a common breakpoint region between two different regions on one or two chromosomes 
	Also contains information about the modal chromosomal co-ordinates for the breakpoint """

	def __init__(self,SplitReadAlignment_set,chr1,chr2,chr1_coord,chr2_coord,modal_read_count):
		self.SplitReadAlignment_set = SplitReadAlignment_set
		self.chr1 = chr1
		self.chr2 = chr2
		self.chr1_coord = chr1_coord
		self.chr2_coord = chr2_coord
		self.chr_fasta_filepath = "/home/aawan/SHARED/REFERENCE/hg38/Homo_sapiens_assembly38.fasta"
		self.modal_read_count = modal_read_count # number of reads at the exact breakpoint with the most number of reads in this set

	@property
	def read_id_set(self):
		return set([SRA.read_id for SRA in self.SplitReadAlignment_set])

	@property
	def length(self):
		return len(self.SplitReadAlignment_set)

	def __repr__(self):
		return f"{self.chr1}-{self.chr1_coord}-{self.chr2}-{self.chr2_coord}-{len(self.SplitReadAlignment_set)}"

	def __eq__(self,other):
		if (self.read_id_set == other.read_id_set) and (self.chr1 == other.chr1) and (self.chr2 == other.chr2) and (self.chr1_coord == other.chr1_coord) and (self.chr2_coord == other.chr2_coord):
			return True
		return False

	def _does_SRA_match_coords(self,SRA):
		""" return True if the given SRA's co-ordinates match the co-ordinates of this BreakpointReadSet, otherwise False """
		if (SRA.chr1 == self.chr1) and (SRA.chr2 == self.chr2) and ((SRA.lb1 == self.chr1_coord) or (SRA.ub1 == self.chr1_coord)) and ((SRA.lb2 == self.chr2_coord) or (SRA.ub2 == self.chr2_coord)):
			return True
		return False

	def _get_single_SRA_single_direction_given_lb_ub_breakpoint_coord(self,lb,ub,breakpoint_coord):
		""" given one lb and ub from one chromosome and the corresponding breakpoint_coord from the same chromosome, return the scalar direction of the SRA away from the breakpoint """
		if ((ub - breakpoint_coord) + (lb - breakpoint_coord)) > 0:
			return '+'
		return '-'

	def _get_direction_of_single_SRA_from_breakpoints(self,SRA):
		""" return a  2-tuple representing the direction of the given SplitReadAlignment away from this object's breakpoint coords """
		direction1 = self._get_single_SRA_single_direction_given_lb_ub_breakpoint_coord(SRA.lb1,SRA.ub1,self.chr1_coord)
		direction2 = self._get_single_SRA_single_direction_given_lb_ub_breakpoint_coord(SRA.lb2,SRA.ub2,self.chr2_coord)
		return (direction1, direction2)

	def _get_direction_of_SRAs_from_breakpoints(self):
		""" iterate through this objects SplitReadAlignment_set, and return a 2-tuple representing the majority direction of the SplitReadAlignments away from the breakpoint coords """
		directions = [defaultdict(int), defaultdict(int)]
		for SRA in self.SplitReadAlignment_set:
			SRA_directions = self._get_direction_of_single_SRA_from_breakpoints(SRA)
			directions[0][SRA_directions[0]] += 1
			directions[1][SRA_directions[1]] += 1
		return (max(directions[0], key=directions[0].get),  max(directions[1], key=directions[1].get)) # this line returns the majority vote for each of the two directions over all SRAs

	def _get_breakpoint_Regions(self,lb_ub_distance):
		""" given the desired distance between lb and ub, returns a 2-tuple of a dict of Region objects representing both sides of the breakpoint region 
		The strand of each region is chosen to change the directions obtained from self._get_direction_of_SRAs_from_breakpoints() to the desired directions: ('-','+')"""
		directions = self._get_direction_of_SRAs_from_breakpoints()
		Region1 = Region(self.chr1,self.chr1_coord,self.chr1_coord+lb_ub_distance,'-')
		Region2 = Region(self.chr2,self.chr2_coord,self.chr2_coord+lb_ub_distance,'+')
		if directions[0] == '-':
			Region1 = Region(self.chr1,self.chr1_coord-lb_ub_distance,self.chr1_coord,'+')		
		if directions[1] == '-':
			Region2 = Region(self.chr2,self.chr2_coord-lb_ub_distance,self.chr2_coord,'-')
		return (Region1,Region2)

	def _is_Region_in_Region_list(self,Region,Region_list):
		""" return True if Region overlaps with any Region in Region_list, else False """
		for r in Region_list:
			if Region.overlaps(r):
				return True
		return False

	def _get_breakpoint_Region_order_using_probed_Region_list(self,breakpoint_Regions,probed_Region_list):
		""" 
		given a 2-tuple of breakpoint Regions and a list of probed Regions, return the order of breakpoint Regions
		that puts the probed region first
		"""	
		order = (0,1)
		if self._is_Region_in_Region_list(breakpoint_Regions[1],probed_Region_list):
			order = (1,0)
		return order

	def make_breakpoint_region_fasta_file(self,lb_ub_distance,probed_Region_list,outfile_path_prefix="",id_prefix=""):
		""" given an lb_ub_distance and an outfile_path, make a fasta file of the breakpoint genomic fusion region """
		Regions = self._get_breakpoint_Regions(lb_ub_distance)
		order = self._get_breakpoint_Region_order_using_probed_Region_list(Regions,probed_Region_list)
		outfile_path = f"{outfile_path_prefix}{Regions[order[0]]}_{Regions[order[1]]}_breakpoint_genomic_fusion.fa"
		rec_id = f"{Regions[order[0]]};{Regions[order[1]]}"
		seq = Regions[0].seq + Regions[1].seq
		if order == (1,0):
			seq = str(Seq(seq).reverse_complement())
		SeqIO.write([ SeqRecord( Seq(seq),id=f"{id_prefix}{rec_id}",name="",description="" ) ], outfile_path, "fasta-2line")

	def make_selected_read_fastq_file(self,all_read_fastq_filepath,out_filepath):
		""" given the path to a fastq file of all reads, and an out_filepath, make a fastq file of just this object's read_id_set """
		seq_records = []
		for seq_rec in SeqIO.parse(all_read_fastq_filepath, "fastq"):
			if seq_rec.id in self.read_id_set:
				seq_records.append(seq_rec)
		SeqIO.write(seq_records, out_filepath, "fastq")

	def add_SplitReadAlignment_to_SplitReadAlignment_set(self,SRA):
		""" add a single SplitReadAlignment to self.SplitReadAlignment_set"""
		if self._does_SRA_match_coords(SRA):
			self.modal_read_count += 1
		self.SplitReadAlignment_set.add(SRA)
		
	def add_SplitReadAlignment_set_to_SplitReadAlignment_set(self,SplitReadAlignment_set):
		""" add a SplitReadAlignment_set to self.SplitReadAlignment_set"""
		if len(self.SplitReadAlignment_set.intersection(SplitReadAlignment_set)) > 0:
			for SRA in SplitReadAlignment_set:
				if self._does_SRA_match_coords(SRA):
					self.modal_read_count += 1 
		self.SplitReadAlignment_set = self.SplitReadAlignment_set.union(SplitReadAlignment_set)


class BreakpointReadSetListFinder:
	""" methods to obtain a list of BreakpointReadSet objects from a list of SplitReadAlignment objects """
	
	def __init__(self, SplitReadAlignment_list, breakpoint_read_thresh, max_merge_dist, shared_read_merge_thresh, post_merge_thresh, max_brs_read_fraction_thresh, probed_region_bedfile_path):
		self.SplitReadAlignment_list = SplitReadAlignment_list
		self.BREAKPOINT_READ_THRESH = breakpoint_read_thresh
		self.MAX_MERGE_DIST = max_merge_dist # the maximum distance two breakpoint co-ordinates on the same chromosome are away from each other to be considered the same
		self.POST_MERGE_THRESH = post_merge_thresh # after merging, the minimum number of reads needed to say a breakpoint is real
		self.SHARED_READ_MERGE_THRESH = shared_read_merge_thresh # percentage of reads shared by two breakpoint read sets above which to consider merging
		self.MAX_BRS_READ_FRACTION_THRESH = max_brs_read_fraction_thresh # percentage of read count from breakpoint with most reads for a breakpoint to be kept
		self.probed_Regions = self._get_Region_list_from_bedfile_path(probed_region_bedfile_path)
		self.exact_coord_BreakpointReadSet_dict = self._get_exact_coord_BreakpointReadSet_dict()
		self.region_BreakpointReadSet_dict = self._get_region_seed_BreakpointReadSet_dict()
		self._expand_region_seed_BreakpointReadSet_dict()
		self._make_merged_region_BreakpointReadSet_dict()
		self._prune_merged_region_BreakpointReadSet_dict()
		self._prune_BreakpointReadSet_dict_by_probed_regions()
		self._split_all_reciprocal_translocation_BreakpointReadSets()
		self._prune_BreakpointReadSet_dict_by_biggest_breakpoint_read_fraction()

	def _get_Region_list_from_bedfile_path(self,bedfile_path):
		""" given a bedfile path, return a list of Region objects """
		Region_list = []
		with open(bedfile_path) as f:
			for line in f:
				tabs = line.rstrip('\n').split('\t')
				chr = tabs[0]
				lb = int(tabs[1])
				ub = int(tabs[2])
				Region_list.append(Region(chr,lb,ub,None))
		return Region_list

	def _get_breakpoint_key_from_chrs_and_coords(self,chr1,chr2,chr1_coord,chr2_coord):
		""" given chr1 and chr2, and one coord each from chr1 and chr2, return a hashable key string """
		return f"{chr1}-{chr1_coord}-{chr2}-{chr2_coord}"

	def _get_constituents_from_breakpoint_key(self,breakpoint_key):
		""" inverse method to _get_breakpoint_key_from_chrs_and_coords """
		chr1, chr1_coord, chr2, chr2_coord =  breakpoint_key.split("-")
		return [chr1,chr2,int(chr1_coord),int(chr2_coord)]

	def _get_all_potential_breakpoint_keys_from_SplitReadAlignment(self,SRA):
		""" return a list of strings in the format {lower chr}-{coord from lower chr}-{upper chr}-{coord from upper chr}, where both coords are considered """
		keys = []
		for chr1_coord in [SRA.lb1, SRA.ub1]:
			for chr2_coord in [SRA.lb2, SRA.ub2]:
				keys.append(self._get_breakpoint_key_from_chrs_and_coords(SRA.chr1,SRA.chr2,chr1_coord,chr2_coord))
		return keys

	def _get_new_BreakpointReadSet_from_breakpoint_key(self,breakpoint_key):
		""" given a breakpoint key in the format generated by the _get_breakpoint_key_from_chrs_and_coords method, return the a new 
		BreakpointReadSet object with an empty SplitReadAlignment_set and the chrs and coords from the breakpoint_key """
		chr1, chr2, chr1_coord, chr2_coord = self._get_constituents_from_breakpoint_key(breakpoint_key)
		return BreakpointReadSet(set(),chr1, chr2, chr1_coord, chr2_coord,0)

	def _get_distance_between_two_breakpoint_keys(self,bk1,bk2):
		""" given two breakpoint keys, return a two-tuple of the scalar distances between the chr1_coords and chr2 coords of both keys,
		or None if they are on different chromosomes """
		bk1_chr1, bk1_chr2, bk1_chr1_coord, bk1_chr2_coord = self._get_constituents_from_breakpoint_key(bk1)
		bk2_chr1, bk2_chr2, bk2_chr1_coord, bk2_chr2_coord = self._get_constituents_from_breakpoint_key(bk2)
		if (bk1_chr1 != bk2_chr1) or (bk1_chr2 != bk2_chr2):
			return None
		return (abs(bk1_chr1_coord-bk2_chr1_coord), abs(bk1_chr2_coord-bk2_chr2_coord))

	def _is_breakpoint_distance_mergable(self, breakpoint_distance):
		""" return True if both values in the two-tuple of breakpoint_distance are <= self.MAX_MERGE_DIST and > 0, otherwise return False """
		if (breakpoint_distance is None) or (breakpoint_distance[0] > self.MAX_MERGE_DIST) or (breakpoint_distance[1] > self.MAX_MERGE_DIST) or (breakpoint_distance == (0,0)):
			return False
		return True

	def _get_exact_coord_BreakpointReadSet_dict(self):
		""" use the SplitReadAlignment_list to return a dict of all exact-coordinate BreakpointReadSets that do not search a range of co-ordinates 
		where the key is {lower chr}-{bound from lower chr}-{upper chr}-{bound from upper chr}, where both bounds are considered. So each read gets indexed 
		in four different BreakpointReadSets"""
		BreakpointReadSet_dict = {}
		for SRA in self.SplitReadAlignment_list:
			if SRA._are_alignments_small_indel(): # don't consider breakpoinbts that are small indels
				continue
			for SRA_key in self._get_all_potential_breakpoint_keys_from_SplitReadAlignment(SRA):
				if SRA_key not in BreakpointReadSet_dict:
					BreakpointReadSet_dict[SRA_key] = self._get_new_BreakpointReadSet_from_breakpoint_key(SRA_key)
				BreakpointReadSet_dict[SRA_key].add_SplitReadAlignment_to_SplitReadAlignment_set(SRA)
		return BreakpointReadSet_dict

	def _get_region_seed_BreakpointReadSet_dict(self):
		""" using self.exact_coord_BreakpointReadSet_dict, return a dict of key: BreakpointReadSet for all BreakpointReadSets where the number of reads is > self.BREAKPOINT_READ_THRESH """	
		BreakpointReadSet_dict = {}
		for brs in self.exact_coord_BreakpointReadSet_dict:
			if len(self.exact_coord_BreakpointReadSet_dict[brs].SplitReadAlignment_set) >= self.BREAKPOINT_READ_THRESH:	
				BreakpointReadSet_dict[brs] = self.exact_coord_BreakpointReadSet_dict[brs]
		return BreakpointReadSet_dict	

	def _expand_region_seed_BreakpointReadSet_dict(self):
		""" compare each key in self.exact_coord_BreakpointReadSet_dict to each key in self.region_BreakpointReadSet_dict, and if the distance between the 
		keys is <= self.MAX_MERGE_DIST, add the SplitReadAlignment_set from the former to the SplitReadAlignment_set of the latter """
		for breakpoint_key1 in self.exact_coord_BreakpointReadSet_dict:
			for breakpoint_key2 in self.region_BreakpointReadSet_dict:
				breakpoint_distance = self._get_distance_between_two_breakpoint_keys(breakpoint_key1,breakpoint_key2)
				if self._is_breakpoint_distance_mergable(breakpoint_distance):
					self.region_BreakpointReadSet_dict[breakpoint_key2].add_SplitReadAlignment_set_to_SplitReadAlignment_set( self.exact_coord_BreakpointReadSet_dict[breakpoint_key1].SplitReadAlignment_set )

	def _get_merged_BreakpointReadSet_from_pair(self,BRS1,BRS2,shared_read_thresh):
		""" given two input BreakpointReadSet objects, if they share any SplitReadAlignments, return a merged BreakpointReadSet, where the new coordinates are those of the 
		single co-ordinate with the most reads otherwise return None """
		if (BRS1 == BRS2) or (len( BRS1.read_id_set.union(BRS2.read_id_set)) == 0) or ( len( BRS1.read_id_set.intersection(BRS2.read_id_set) ) / len( BRS1.read_id_set.union(BRS2.read_id_set) ) < shared_read_thresh):
			return None
		old_BRS = BRS1
		other_BRS = BRS2
		if BRS2.modal_read_count > BRS1.modal_read_count:
			old_BRS = BRS2
			other_BRS = BRS1
		new_BRS = BreakpointReadSet(old_BRS.SplitReadAlignment_set,old_BRS.chr1,old_BRS.chr2,old_BRS.chr1_coord,old_BRS.chr2_coord,old_BRS.modal_read_count)
		new_BRS.add_SplitReadAlignment_set_to_SplitReadAlignment_set(other_BRS.SplitReadAlignment_set)
		return new_BRS

	def _do_single_region_BreakpointReadSet_dict_merge_iteration(self):	
		""" loop through every breakpoint_key pair in the self.region_BreakpointReadSet_dict, updating the dict as required, and returning True if there's a merge, False otherwise """
		for breakpoint_key1 in self.region_BreakpointReadSet_dict:
			for breakpoint_key2 in self.region_BreakpointReadSet_dict:
				if breakpoint_key1 == breakpoint_key2:
					continue
				if len(self.region_BreakpointReadSet_dict[breakpoint_key1].read_id_set.intersection( self.region_BreakpointReadSet_dict[breakpoint_key2].read_id_set )) > 0:
					new_BRS = self._get_merged_BreakpointReadSet_from_pair(self.region_BreakpointReadSet_dict[breakpoint_key1], self.region_BreakpointReadSet_dict[breakpoint_key2],self.SHARED_READ_MERGE_THRESH)
					if new_BRS is not None:
						new_breakpoint_key = self._get_breakpoint_key_from_chrs_and_coords(new_BRS.chr1,new_BRS.chr2,new_BRS.chr1_coord,new_BRS.chr2_coord)
						del self.region_BreakpointReadSet_dict[breakpoint_key1]
						del self.region_BreakpointReadSet_dict[breakpoint_key2]
						self.region_BreakpointReadSet_dict[new_breakpoint_key] = new_BRS
						return True
		return False

	def _make_merged_region_BreakpointReadSet_dict(self):
		""" using self.region_BreakpointReadSet_dict, return an updated version where any two BreakpointReadSets that share enough SplitReadAlignments are merged """
		merged = True
		while(merged):
			merged = self._do_single_region_BreakpointReadSet_dict_merge_iteration()

	def _prune_merged_region_BreakpointReadSet_dict(self):
		""" get rid of any merged BreakpointReadSets where the number of reads is lower than self.POST_MERGE_THRESH """
		pruned_region_BreakpointReadSet_dict = {}
		for breakpoint_key in self.region_BreakpointReadSet_dict:
			if self.region_BreakpointReadSet_dict[breakpoint_key].length >= self.POST_MERGE_THRESH:
				pruned_region_BreakpointReadSet_dict[breakpoint_key] = self.region_BreakpointReadSet_dict[breakpoint_key]
		self.region_BreakpointReadSet_dict = pruned_region_BreakpointReadSet_dict

	def _is_BreakpointReadSet_in_Region(self,BRS,Region):
		""" Return true if the given BreakpointReadSet is within the given Region, else False """
		if (BRS.chr1 == Region.chr) and (int(BRS.chr1_coord) >= int(Region.lb)) and (int(BRS.chr1_coord) <= int(Region.ub)):
			return True
		if (BRS.chr2 == Region.chr) and (int(BRS.chr2_coord) >= int(Region.lb)) and (int(BRS.chr2_coord) <= int(Region.ub)):
			return True
		return False

	def _prune_BreakpointReadSet_dict_by_probed_regions(self):
		pruned_region_BreakpointReadSet_dict = {}
		for breakpoint_key in self.region_BreakpointReadSet_dict:
			for region in self.probed_Regions:
				if self._is_BreakpointReadSet_in_Region(self.region_BreakpointReadSet_dict[breakpoint_key], region):
					pruned_region_BreakpointReadSet_dict[breakpoint_key] = self.region_BreakpointReadSet_dict[breakpoint_key]
					break
		self.region_BreakpointReadSet_dict = pruned_region_BreakpointReadSet_dict

	def _split_reciprocal_translocation_BreakpointReadSet(self,BRS):
		""" If this BreakpointReadSet consists of two reciprocal translocations sharing a breakpoint, split the BRS into two """
		left = set()
		right = set()
		other = set()
		for SRA in BRS.SplitReadAlignment_set:
			if (abs(SRA.lb1-BRS.chr1_coord) < abs(SRA.ub1-BRS.chr1_coord)) and (abs(SRA.lb2-BRS.chr2_coord) < abs(SRA.ub2-BRS.chr2_coord)):
				left.add(SRA)
			elif (abs(SRA.lb1-BRS.chr1_coord) > abs(SRA.ub1-BRS.chr1_coord)) and (abs(SRA.lb2-BRS.chr2_coord) > abs(SRA.ub2-BRS.chr2_coord)):
				right.add(SRA)
			else:
				other.add(SRA)
		if (len(left) > 0) and (len(right) > 0) and (len(other) == 0):
			left_BRS = BreakpointReadSet(left,BRS.chr1,BRS.chr2,BRS.chr1_coord,BRS.chr2_coord,sum([1 for x in left if BRS._does_SRA_match_coords(x)]))		
			right_BRS = BreakpointReadSet(right,BRS.chr1,BRS.chr2,BRS.chr1_coord,BRS.chr2_coord,sum([1 for x in right if BRS._does_SRA_match_coords(x)]))
			return [left_BRS, right_BRS]
		return [BRS]

	def _split_all_reciprocal_translocation_BreakpointReadSets(self):
		""" perform _split_reciprocal_translocation_BreakpointReadSet on all BreakpointReadSets in self.region_BreakpointReadSet_dict """
		new_region_BreakpointReadSet_dict = {}
		split_BreakpointReadSet_dict = {}
		for breakpoint_key in self.region_BreakpointReadSet_dict:
			## do the split, and add the first result to new_region_BreakpointReadSet_dict, and the second if it exists to split_BreakpointReadSet_dict
			split = self._split_reciprocal_translocation_BreakpointReadSet(self.region_BreakpointReadSet_dict[breakpoint_key])
			key = self._get_breakpoint_key_from_chrs_and_coords(split[0].chr1,split[0].chr2,split[0].chr1_coord,split[0].chr2_coord)
			new_region_BreakpointReadSet_dict[key] = split[0]
			if len(split) > 1:
				split_BreakpointReadSet_dict[key] = split[1]
		self.region_BreakpointReadSet_dict = new_region_BreakpointReadSet_dict
		self.split_BreakpointReadSet_dict	= split_BreakpointReadSet_dict		
			
	def _prune_BreakpointReadSet_dict_by_biggest_breakpoint_read_fraction(self):
		""" remove those BreakpointReadSets whose read count as a fraction of the biggest BreakpointReadSet is not greater than self.MAX_BRS_READ_FRACTION_THRESH """
		pruned_region_BreakpointReadSet_dict = {}
		pruned_split_BreakpointReadSet_dict = {}
		max_read_count = max([len(self.region_BreakpointReadSet_dict[BRS].SplitReadAlignment_set) for BRS in self.region_BreakpointReadSet_dict])
		if len(self.split_BreakpointReadSet_dict) > 0:
			mrc2 = max([len(self.split_BreakpointReadSet_dict[BRS].SplitReadAlignment_set) for BRS in self.split_BreakpointReadSet_dict])
			max_read_count = max(max_read_count,mrc2)
		for breakpoint_key in self.region_BreakpointReadSet_dict:
			if len(self.region_BreakpointReadSet_dict[breakpoint_key].SplitReadAlignment_set)/max_read_count > self.MAX_BRS_READ_FRACTION_THRESH:
				pruned_region_BreakpointReadSet_dict[breakpoint_key] = self.region_BreakpointReadSet_dict[breakpoint_key]
		self.region_BreakpointReadSet_dict = pruned_region_BreakpointReadSet_dict
		for breakpoint_key in self.split_BreakpointReadSet_dict:
			if len(self.split_BreakpointReadSet_dict[breakpoint_key].SplitReadAlignment_set)/max_read_count > self.MAX_BRS_READ_FRACTION_THRESH:
				pruned_split_BreakpointReadSet_dict[breakpoint_key] = self.split_BreakpointReadSet_dict[breakpoint_key]
		self.split_BreakpointReadSet_dict = pruned_split_BreakpointReadSet_dict

	def make_all_BRS_read_fastq_files(self,all_read_fastq_filepath,outfile_prefix):
		""" make the read fastq files for all BRSs """
		for breakpoint_key in self.region_BreakpointReadSet_dict:
			self.region_BreakpointReadSet_dict[breakpoint_key].make_selected_read_fastq_file(all_read_fastq_filepath,f"{outfile_prefix}{breakpoint_key}_breakpoint_reads.fq")
		for breakpoint_key in self.split_BreakpointReadSet_dict:
			self.split_BreakpointReadSet_dict[breakpoint_key].make_selected_read_fastq_file(all_read_fastq_filepath,f"{outfile_prefix}reciprocal_{breakpoint_key}_breakpoint_reads.fq")

	def make_all_BRS_fasta_files(self,lb_ub_distance,outfile_path_prefix="",id_prefix=""):
		""" make the genomic fusion fasta files for all BRSs """
		for breakpoint_key in sorted(self.region_BreakpointReadSet_dict,key=lambda x: len(self.region_BreakpointReadSet_dict[x].SplitReadAlignment_set)): 
			print(breakpoint_key,len(self.region_BreakpointReadSet_dict[breakpoint_key].SplitReadAlignment_set))
			self.region_BreakpointReadSet_dict[breakpoint_key].make_breakpoint_region_fasta_file(lb_ub_distance,self.probed_Regions,outfile_path_prefix,id_prefix)
		for breakpoint_key in sorted(self.split_BreakpointReadSet_dict,key=lambda x: len(self.split_BreakpointReadSet_dict[x].SplitReadAlignment_set)):
			print(breakpoint_key,len(self.split_BreakpointReadSet_dict[breakpoint_key].SplitReadAlignment_set))
			self.split_BreakpointReadSet_dict[breakpoint_key].make_breakpoint_region_fasta_file(lb_ub_distance,self.probed_Regions,f"{outfile_path_prefix}reciprocal_",f"{id_prefix}reciprocal_")

if __name__ == "__main__":
	pass
