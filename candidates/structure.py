

class Candidate:
    ''' structure to store a miRNA candidate '''
    
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
    
#     def __init__(self, hairpin, padded, prime5, prime3, id5, id3, strand):
#         self.seq_hairpin = hairpin
#         self.seq_hairpin_padded = padded
#         self.seq_prime5 = prime5
#         self.seq_prime3 = prime3
#         self.id_prime5 = id5        
#         self.id_prime3 = id3
#         self.strand = strand
        
    def __init__(self, hairpin, padded, interval_data):
#         [strand, 5'name, 5'sequence, 3'name, 3'sequence]
        self.seq_hairpin = hairpin
        self.seq_hairpin_padded = padded
        self.strand = interval_data[0]
        self.id_prime5 = interval_data[1]        
        self.seq_prime5 = interval_data[2]
        self.id_prime3 = interval_data[3]
        self.seq_prime3 = interval_data[4]
     
    def set_viennafold(self, hairpin_fold, hairpin_energy):
        self.hairpin_fold = hairpin_fold 
        self.hairpin_energy = hairpin_energy
#         self.padded_fold = padded_fold
#         self.padded_energy = padded_energy
#     def get_hairpin(self):
#         return self.seq_hairpin
#     
#     def get_padded_hairpin(self):
#         return self.seq_hairpin_padded