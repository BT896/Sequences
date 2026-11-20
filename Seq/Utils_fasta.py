# Reading a fasta file
class FastaFile:
    # Constructor
    def __init__(self, file, type):
        self._file = file
        self._type = type
        self._records = []
    @property
    def file(self):
        return self._file
    @property
    def records(self):
        return self._records
    @property
    def type(self):
        return self._type

    def parse(self):
        with open(self.file) as file:
            gene = None
            seq = ""

            for line in file:
                line = line.strip()

                if line.startswith(">"):
                    gene = line.rstrip().lstrip('>')
                    record = DNASequence(gene, seq)
                    self.records.append(record)

                    seq = ""
                else:
                    seq += line

            # Save the last sequence
            if gene is not None:
                record = DNASequence(gene, seq)
                self._records.append(record)

        return self._records