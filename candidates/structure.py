

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
        
        self.all_mapped_sequences.update(mapped_sequences)
        
        self.chromosome = chromosome
        self.chromosome_direction = strand_dir
        self.pos_5_begin = begin_5
        self.pos_5_end = end_5
        self.pos_3_begin = begin_3
        self.pos_3_end = end_3
        self.mapped_sequences = set(mapped_sequences)
    
    
    def set_hairpin_padding(self, hairpin, padded):
        self.hairpin = hairpin
        self.hairpin_padded = padded
    
    def set_viennafold(self, hairpin_fold, hairpin_energy):
        self.hairpin_fold = hairpin_fold 
        self.hairpin_energy = hairpin_energy
        
    def add_sequence(self, sequences):
        
        self.all_mapped_sequences.update(sequences)
        self.mapped_sequences.update(sequences)

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