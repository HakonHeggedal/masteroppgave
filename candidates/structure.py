


class Candidate:
    ''' structure to store a miRNA candidate '''
    
    all_mapped_sequences = set()
    
    

#     seq_hairpin = ""
#     seq_hairpin_padded = ""
#     seq_prime5 = ""
#     seq_prime3 = ""
#     strand = "" # + is forward, - is backward
#     id_prime5 = ""
#     id_prime3 = ""
#     energy_hairpin = 0.0
#     energy_padding = 0.0
#     heterogenity_prime5 = 0
#     heterogenity_prime3 = 0
#      
#     entropy_nucleotide = 0
#     entropy_structure = 0
#      
#     overhang_prime5 = 0
#     overhang_prime3 = 0
#      
#     tailing_au = 0
#      
#     loki_quality = 0.0 # 0 to 1

    def __init__(self, chromosome, strand_dir, begin_5, end_5,
                 begin_3, end_3, mapped_sequences):
        
        self.chromosome = chromosome
        self.chromosome_direction = strand_dir
        self.pos_5_begin = begin_5
        self.pos_5_end = end_5
        self.pos_3_begin = begin_3
        self.pos_3_end = end_3
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

#     def __init__(self, interval_data, chromosome):
# #         [strand, 5'name, 5'sequence, 3'name, 3'sequence]
#         self.strand = interval_data[0]
#         self.id_prime5 = interval_data[1]        
#         self.seq_prime5 = interval_data[2]
#         self.id_prime3 = interval_data[3]
#         self.seq_prime3 = interval_data[4]
#         self.chromosome = chromosome

#     def get_hairpin(self):
#         return self.seq_hairpin
#     
#     def get_padded_hairpin(self):
#         return self.seq_hairpin_padded








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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    