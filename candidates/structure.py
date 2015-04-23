
TYPE_UNDECIDED = -1
TYPE_CANDIDATE = 0
TYPE_DEAD = 1
TYPE_HIGH_CONF = 2
TYPE_LOW_CONF = 3


class Candidate:
    ''' structure to store a miRNA candidate '''
    
    all_mapped_sequences = set()
    
    def __init__(self, chromosome, strand_dir, hairpin_start, hairpin_end,
                 begin_5p, end_5p, begin_3p, end_3p, mapped_sequences):

        self.chromosome = chromosome
        self.chromosome_direction = strand_dir
        
        self.hairpin_start = hairpin_start
        self.hairpin_end = hairpin_end
        self.pos_5p_begin = begin_5p # in the genome... large value
        self.pos_5p_end = end_5p
        self.pos_3p_begin = begin_3p
        self.pos_3p_end = end_3p
        
        self.hairpin = None
        self.mapped_sequences = set()
        self.small_subs = 0
        self.small_subs_5p = 0
        self.small_subs_3p = 0
        
        self.padding_size = 40
        self.candidate_type = TYPE_UNDECIDED
        self.miRNAid = None
        
        self.bulge_factor = -1
        self.peak_5b = -1
        self.peak_5e = -1
        self.peak_3b = -1
        self.peak_3e = -1
        
        self.has_overhang = False
        self.has_5p = not (begin_5p == end_5p == -1)
        self.has_3p = not (begin_3p == end_3p == -1)
        self.has_hairpin_struct = False
        self.has_hairpin_struct_5p = False
        self.has_hairpin_struct_3p = False
        
        self.estimate_5pb = -1 # only if there is no 5p or 3p seq
        self.estimate_5pe = -1
        self.estimate_3pb = -1
        self.estimate_3pe = -1

        self.overhang_inner = 15
        self.overhang_outer = 15
        self.overhang_inner_conf = 0.0
        self.overhang_outer_conf = 0.0
        self.loop_size = 0
        self.folds_5p = 0
        self.folds_3p = 0
        self.folds_before = 0
        self.folds_after = 0
        
        self.stops_before_5p = 0.0
        self.starts_after_3p = 0.0
        
        self.ratio_short_long = 1.0
        self.ratio_short_long_logval = 1.0
#         self.ratio_short_long_5p = 1.0 # merge
#         self.ratio_short_long_3p = 1.0 # merge
        self.has_short_seqs_5p = False
        self.has_short_seqs_3p = False
        self.short_seq_5p_stdev = -1 # not used
        self.short_seq_3p_stdev = -1 # not used
        self.short_seq_5p_offset = -1 # not used
        self.short_seq_3p_offset = -1 # not used
        
        
        self.short_seq_align = 0.0
        self.short_seq_align_10_13 = 0.0
        self.short_seq_align_10_14 = 0.0
        self.short_seq_align_10_15 = 0.0
        self.short_seq_align_10_16 = 0.0
        self.short_seq_align_10_17 = 0.0
        self.short_seq_align_10_18 = 0.0 # too long ? 
        self.short_seq_align_10_19 = 0.0 # too long ? 

        self.short_seq_align_8_17 = 0.0 # too long ? 
        self.short_seq_align_9_17 = 0.0 # too long ? 
        self.short_seq_align_10_17 = 0.0 # too long ? 
        self.short_seq_align_11_17 = 0.0 # too long ? 
        self.short_seq_align_12_17 = 0.0 # too long ? 
        self.short_seq_align_13_17 = 0.0 # too long ? 
        self.short_seq_align_14_17 = 0.0 # too long ? 
        self.short_seq_align_15_17 = 0.0 # too long ? 
        self.short_seq_align_16_17 = 0.0 # too long ? 


        
        self.leading_au = 0
        self.tailing_au = 0
        
        self.mirBase_matures = None
        
        if mapped_sequences is not None:
            self.mapped_sequences = set(mapped_sequences)
            self.all_mapped_sequences.update(mapped_sequences)

    def set_seq_outside(self, sequences_before, sequences_after):
        self.sequences_before = sequences_before
        self.sequences_after = sequences_after


    def set_hairpin_padding(self, hairpin, padded_40):
        self.hairpin = hairpin
        self.hairpin_padded_40 = padded_40
    
    def set_fold_hairpin(self, fold, en):
        self.hairpin_fold = fold
        self.hairpin_energy = en

    
    def set_fold_10(self, fold_10, en_10):        
        self.hairpin_energy_10 = en_10
        self.hairpin_fold_10 = fold_10

    
    def set_fold_40(self, fold_40, en_40):
        self.hairpin_fold_40 = fold_40
        self.hairpin_energy_40 = en_40

    def set_hairpin_pos(self, begin, end):
        self.pos_hairpin_begin = begin
        self.pos_hairpin_end = end
        
    def add_sequence(self, sequences):
        self.all_mapped_sequences.update(sequences)
        self.mapped_sequences.update(sequences)
        
    def set_entropy(self, entropy_nucleotides, entropy_structure):
        self.entropy_nucleotides = entropy_nucleotides
        self.entropy_structure = entropy_structure
        
    def set_heterogenity(self, h_5_begin, h_5_end, h_3_begin, h_3_end):
        self.heterogenity_5_begin = h_5_begin
        self.heterogenity_5_end = h_5_end
        self.heterogenity_3_begin = h_3_begin
        self.heterogenity_3_end = h_3_end
        
    def set_reads(self, reads_5p_begin, reads_5p_end, reads_3p_begin, reads_3p_end):
        self.reads_5p_begin = reads_5p_begin
        self.reads_5p_end = reads_5p_end
        self.reads_3p_begin = reads_3p_begin
        self.reads_3p_end = reads_3p_end
        
    def set_quality(self, quality):
        self.quality = quality
        
    def set_tailing(self,tailing):
        self.tailing = tailing
    
    def set_bitpair_entropy(self, bitpair_dict):
        self.bitpair_entropy_dict = bitpair_dict
        
    def set_bitpair_probs(self, bitpair_p):
        self.bitpair_probabilities = bitpair_p
    
    def set_junction_pos(self, pos_5, pos_3):
        self.junction_pos_5 = pos_5
        self.junction_pos_3 = pos_3
        
    def set_alignment_10(self, max_10, lev_10_out, oh_10_out, lev_10_in, oh_10_in):
        self.bindings_max_10 = max_10
        self.overhang_level_outer_10 = lev_10_out
        self.overhang_outer_10 = oh_10_out
        self.overhang_level_inner_10 = lev_10_in
        self.overhang_inner_10 = oh_10_in

    def set_alignment_40(self, max_40, lev_40_out, oh_40_out, lev_40_in, oh_40_in):
        self.overhang_outer_40 = oh_40_out
        self.overhang_inner_40 = oh_40_in
        self.bindings_max_40 = max_40
        self.overhang_level_outer_40 = lev_40_out
        self.overhang_level_inner_40 = lev_40_in


class Sequence:
    ''' contains a sequence from a sequencing file, and its duplicate number '''
    
    def __init__(self, number_id, duplicates, nucleotides):
        self.candidates = set([])
        self.number_id = number_id
        self.duplicates = duplicates
        self.nucleotides = nucleotides
#         self.number = int(name.split("-")[0])
#         self.duplicates = int(name.split("-")[1])
        
    def add_candidates(self, candidates):
        self.candidates.update(candidates)
    
    def add_candidate(self, candidate):
        self.candidates.add(candidate)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    