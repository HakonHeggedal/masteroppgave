class Candidate:
    ''' structure to store a miRNA candidate '''
    
    all_mapped_sequences = set()
    
    def __init__(self, chromosome, strand_dir, begin_5p, end_5p,
                 begin_3p, end_3p, mapped_sequences):

        self.chromosome = chromosome
        self.chromosome_direction = strand_dir
        self.pos_5p_begin = begin_5p
        self.pos_5p_end = end_5p
        self.pos_3p_begin = begin_3p
        self.pos_3p_end = end_3p
        self.hairpin = None
        self.mapped_sequences = set()
        self.small_subs = 0
        self.small_subs_5p = 0
        self.small_subs_3p = 0
        self.padding_size = 40
        
        if mapped_sequences is not None:
            self.mapped_sequences = set(mapped_sequences)
            self.all_mapped_sequences.update(mapped_sequences)
    
    def set_hairpin_padding(self, hairpin, padded_40):
        self.hairpin = hairpin
        self.hairpin_padded_40 = padded_40

    
    def set_viennafold(self, fold, en, fold_10, en_10, fold_40, en_40):        
        self.hairpin_fold = fold 
        self.hairpin_fold_10 = fold_10
        self.hairpin_fold_40 = fold_40
        self.hairpin_energy = en
        self.hairpin_energy_10 = en_10
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
    
    def set_small_seq(self, index, copies):
        self.small_subs += copies
        if index == self.padding_size:
            self.small_subs_5p += copies
        elif index == self.padding_size + self.pos_3p_begin - self.pos_5p_begin:
            self.small_subs_3p += copies

        
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    