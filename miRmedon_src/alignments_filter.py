###################################################################################################################

# This script includes 3 three functions for filtering bam file after alignment to e-miRbase with STAR. The script
# filter alignment accriding to these 3 creiteria:
#   1. Read wasn't soft-clipped more than 3 cases v
#   2. Mis matches aren't located at the putatve editing position
#   3. For read with multiple alignmnet, only alignment with minimum editing positions would be kept. 

###################################################################################################################
import pickle
import re
import pysam
import argparse 


def get_mismatch_pos(alignment):
    reference_sequence = alignment.get_reference_sequence()
    try:
        pos = re.search(r'[a-z]', reference_sequence).start()
    except AttributeError:
        pos = -1
    return pos

    
def cigar_inspect(alignment, soft_clip_thresh, total_mod_thresh):
    """This fucntion will take an alignment as an input,
    and returns False if STAR soft clipped more than 3bp at
    the ends, and True if not."""
    CIGAR = alignment.cigarstring
    search = re.search(r'[IDN]+', CIGAR)
    if search:
        return False
    else:
        split_cigar = list(filter(lambda x: x != '', re.split(r'(\d+)', CIGAR)))
        try:
            softclipped_count = sum(map(int, list(zip(*list(filter(lambda x: x[1] == 'S', 
                                                                   list(zip(split_cigar[::2], split_cigar[1::2]))))))[0]))
            if softclipped_count > soft_clip_thresh:
                return False
        except IndexError:
            softclipped_count = 0
    if get_mismatch_pos(alignment) != -1:
        if softclipped_count+1 > total_mod_thresh:
            return False
    return True

 
def mismatch_pos_inspect(alignment, e_miRbase_dict):
    """This function takes an alignment as input,
    and returns Flase if the mismatch position isn't 
    intersect with the putative editing sites of the specific 
    edited miRNA"""
    reference_name = alignment.reference_name
    if reference_name.find('_e') != -1:
        reference_sequence = alignment.get_reference_sequence()
        mismatch_query_index = get_mismatch_pos(alignment)
        if mismatch_query_index != -1:
            reference_name = alignment.reference_name
            basename = reference_name[:reference_name.find('_')]
            reference_full_suquence = e_miRbase_dict[basename][reference_name]['sequence']
            mismatch_reference_index = reference_full_suquence.find(reference_sequence.upper()) + mismatch_query_index
            editingSites = e_miRbase_dict[basename][reference_name]['editingSites']
            if mismatch_reference_index in editingSites:
                return False
            else:
                return True
        else:
            return True
    return True

    
def min_editing_filter(list_of_alignments, e_miRbase_dict):
    """This function takes a list of alignments,
    and retruns a filtered list with only alignments that were mapped to
    the miRNA with the minumum number of editing sites"""
    editing_sites_numbers = []
    for alignment in list_of_alignments:
        reference_name = alignment.reference_name
        basename = reference_name[:reference_name.find('_')]
        editing_sites = e_miRbase_dict[basename][reference_name]['editingSites']
        editing_sites_numbers.append(len(editing_sites))
    minimum_editing_sites = min(editing_sites_numbers)
    filtered_alignmnets_list = []
    for alignment in list_of_alignments:
        reference_name = alignment.reference_name
        basename = reference_name[:reference_name.find('_')]
        editing_sites = e_miRbase_dict[basename][reference_name]['editingSites']
        if len(editing_sites) == minimum_editing_sites:
            filtered_alignmnets_list.append(alignment)
    return filtered_alignmnets_list


def single_alignment_inspect(alignment, e_miRbase_dict, soft_clip_thresh, total_mod_thresh):
    if cigar_inspect(alignment, soft_clip_thresh, total_mod_thresh) == True and \
        mismatch_pos_inspect(alignment, e_miRbase_dict=e_miRbase_dict) == True:
            return True
    return False



def main_alignments_filter(bam_path, path_to_emiRbase, soft_clip_thresh, total_mod_thresh):
    filt_bam_path = bam_path[:bam_path.find('.bam')] + '.filt.st.bam'
    e_miRbase_dict = pickle.load(open(path_to_emiRbase, 'rb'))
    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    filtered_bam_file = pysam.AlignmentFile(filt_bam_path, 'wb', template=bam_file)
    alignments_list = []
    for alignment in bam_file.fetch(until_eof=True):
        if len(alignments_list) == 0:
            alignments_list.append(alignment)
        else:
            if alignment.query_name == alignments_list[-1].query_name:
                alignments_list.append(alignment)
            else:
                valid_alignments = [sub_alignment for sub_alignment in alignments_list 
                                    if single_alignment_inspect(alignment=sub_alignment, e_miRbase_dict=e_miRbase_dict,
                                                                soft_clip_thresh=soft_clip_thresh, total_mod_thresh=total_mod_thresh) == True]
                if len(valid_alignments) > 0:
                    final_filt_alignment_list = min_editing_filter(valid_alignments, e_miRbase_dict=e_miRbase_dict)
                    for sub_alignment in final_filt_alignment_list:
                        filtered_bam_file.write(sub_alignment)
                alignments_list = [alignment]
    valid_alignments = [sub_alignment for sub_alignment in alignments_list
                        if single_alignment_inspect(alignment=sub_alignment, e_miRbase_dict=e_miRbase_dict,
                                                    soft_clip_thresh=soft_clip_thresh, total_mod_thresh=total_mod_thresh) == True]
    if len(valid_alignments) > 0:
        final_filt_alignment_list = min_editing_filter(valid_alignments, e_miRbase_dict=e_miRbase_dict)
        for sub_alignment in final_filt_alignment_list:
            filtered_bam_file.write(sub_alignment)
    filtered_bam_file.close()
    return None
