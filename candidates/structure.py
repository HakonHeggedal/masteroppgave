


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
        self.hairpin_padded = None
        self.hairpin_fold = None
        self.hairpin_energy = None
        self.mapped_sequences = set()
        
        if mapped_sequences is not None:
            self.mapped_sequences = set(mapped_sequences)
            self.all_mapped_sequences.update(mapped_sequences)
    
    def set_hairpin_padding(self, hairpin, padded, padding_size):
        self.hairpin = hairpin
        self.hairpin_padded = padded
        self.padding_size = padding_size
    
    def set_viennafold(self, hairpin_fold, hairpin_energy):
        self.hairpin_fold = hairpin_fold 
        self.hairpin_energy = hairpin_energy
        
    def set_hairpin_pos(self, begin, end):
        self.pos_hairpin_begin = begin
        self.pos_hairpin_end = end
        
    def add_sequence(self, sequences):
        self.all_mapped_sequences.update(sequences)
        self.mapped_sequences.update(sequences)
        
    def set_entropy(self, entropy_nucleotides, entropy_structure):
        self.entropy_nucleotides = entropy_nucleotides
        self.entropy_structure = entropy_structure
        
    def set_heterogenity(self, heterogenity_5, heterogenity_3):
        self.heterogenity_5 = heterogenity_5
        self.heterogenity_3 = heterogenity_3
        
    def set_quality(self, quality):
        self.quality = quality
        
    def set_tailing(self,tailing):
        self.tailing = tailing
        
        
#     def __hash__self(self):
#         return hash(self)
#     
#     def __eq__(self, other):
#         return self == other
        
#     def __hash__(self):
#         return hash( (self.chromosome, self.chromosome_direction, self.pos_5p_begin) )
# 
#     def __eq__(self, other):
#         return ((self.chromosome, self.chromosome_direction, self.pos_5p_begin) ==
#             (other.chromosome, other.chromosome_direction, other.pos_5p_begin) )




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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    