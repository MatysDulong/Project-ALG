import mmh3
from bitarray import bitarray

class BloomFilter:

    def __init__(self, size, num_hashes, k ):
        self.size = size
        self.num_hashes = num_hashes
        self.k = k 
        self.bits = bitarray(size)
        self.bits.setall(0)

    def _positions(self, kmer):
        return [mmh3.hash(kmer, seed)%self.size for seed in range(self.num_hashes)]


    def add(self, kmer):
        """Insère un k-mer dans le filtre en mettant à 1 ses positions."""
        for position in self._positions(kmer):
            self.bits[position] = 1

    def contains(self, kmer):
        return all(self.bits[i] for i in self._positions(kmer))
     
    def add_sequence(self, seq):
        for i in range(len(seq) - self.k + 1):
            self.add(seq[i:i+self.k])



bf = BloomFilter(1000, 5, 21)
positions = bf._positions("ACGTACGTACGTACGTACGTA")
