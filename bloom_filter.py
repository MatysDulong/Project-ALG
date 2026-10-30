import mmh3
from bitarray import bitarray


class BloomFilter:
    """Filtre de Bloom des k-mers d'un génome (bitarray + mmh3)."""

    def __init__(self, size, num_hashes, k):
        """Crée un filtre vide de size bits, avec num_hashes hachages."""
        self.size = size
        self.num_hashes = num_hashes
        self.k = k
        self.bits = bitarray(size)
        self.bits.setall(0)

    def _positions(self, kmer):
        """Renvoie les num_hashes positions du k-mer (une graine par hash)."""
        return [mmh3.hash(kmer, seed) % self.size
                for seed in range(self.num_hashes)]

    def add(self, kmer):
        """Insère un k-mer dans le filtre en mettant à 1 ses positions."""
        for position in self._positions(kmer):
            self.bits[position] = 1

    def contains(self, kmer):
        """Renvoie True si toutes les positions du k-mer sont à 1.

        True peut être un faux positif, False est toujours sûr.
        """
        return all(self.bits[i] for i in self._positions(kmer))

    def add_sequence(self, seq):
        """Insère tous les k-mers de seq (aucun si len(seq) < k)."""
        for i in range(len(seq) - self.k + 1):
            self.add(seq[i:i+self.k])
