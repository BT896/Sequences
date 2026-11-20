# dict (codon)
bases = "tcag".upper()
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

# Superclass
class Sequence:
    # constructor
    def __init__(self, id, seq):
        self._id = id
        self._seq = seq.upper()
        self._seq_len = len(seq)
    # getters
    @property
    def id(self):
        return self._id 
    @property
    def seq(self):
        return self._seq
    def __len__(self):
        return self._seq_len
    
    # fasta format method
    def fasta(self):
        return(f'>{self.id}\n{self.seq}')

# DNASequence Subclass
class DNASequence(Sequence):
    def __str__(self):
        return f'ID: {self.id}\nSequence: {self.seq}'

    # gc content method
    def calc_gc_content(self, dp=2):
        c_count = self.seq.count('C')
        g_count = self.seq.count('G')
        gc_content = (c_count + g_count) / len(self.seq)
        return round(gc_content, dp)

    # translation method
    def translate(self):
        protein = []
        for i in range(0, len(self.seq), 3):
            codon = self.seq[i:i+3]
            amino_acid = codon_table.get(codon, 'X')
            protein.append(amino_acid)

            if amino_acid == '*':
                break

        return ''.join(protein)

# Protein sequence
class ProteinSequence(Sequence):
    def __init__(self, id, seq, descr=None):
        super().__init__(id, seq)
        if descr is None:
            descr = []
        self._descr = descr
    @property
    def descr(self):
        return self._descr
    def __str__(self):
        return f'ID: {self.id}\nSequence: {self.seq}\nDescription: {self.descr}'

    # gc content method
    def calc_hydrophobic(self):
        a_count = self.seq.count('A')
        i_count = self.seq.count('I')
        l_count = self.seq.count('L')
        m_count = self.seq.count('M')
        f_count = self.seq.count('F')
        w_count = self.seq.count('W')
        y_count = self.seq.count('Y')
        v_count = self.seq.count('V')
        hydrophobic_ratio = (a_count + i_count + l_count + m_count + f_count + w_count + y_count + v_count) / len(self.seq)
        hydrophobic_percent = hydrophobic_ratio * 100
        return round(hydrophobic_percent)