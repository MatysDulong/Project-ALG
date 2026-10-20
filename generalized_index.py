import io_utils
import fm_index
from pydivsufsort import divsufsort


class GeneralizedIndex:
    """Approche 1 : un FM-index de G0$G1$...$ et un color array."""

    def __init__(self, genomes):
        """Construit l'index de genomes : liste de listes de contigs.

        L'ordre des génomes (celui du FOF) donne leur couleur, de 0 à N-1.
        """
        contig_list = []
        position = []
        for color, genome in enumerate(genomes):
            for contig in genome:
                contig_list.append(contig+"$")
                position.extend([color] * (len(contig)+1))
        concatenated_text = "".join(contig_list)
        suffix_array = divsufsort(concatenated_text.encode())
        self.fm = fm_index.FMIndex(concatenated_text, suffix_array)
        self.color_array = [position[p] for p in suffix_array]
        self.nb_genomes = len(genomes)

    def query(self, pattern):
        """Renvoie N valeurs : 1 si pattern est dans le génome i, 0 sinon."""
        presence = [0 for _ in range(self.nb_genomes)]
        interval = self.fm.backward_search(pattern)
        if interval is None:
            return presence
        top, bottom = interval
        for i in range(top, bottom+1):
            presence[self.color_array[i]] = 1
        return presence


if __name__ == "__main__":
    genomes = [["ACG"], ["CGT"], ["TAACG"]]
    test = GeneralizedIndex(genomes)
    print(test.fm.bwt)
    print(test.color_array)
    print(test.query("CG"))
    print(test.query("GT"))
    print(test.query("AA"))
    print(test.query("TT"))
    print(test.query("TTT"))
    print(test.query("T$T"))
    genome_paths = io_utils.read_fof("data/test/fof_coli_10k.txt")
    genomes = [io_utils.read_genome(path) for path in genome_paths]
    index = GeneralizedIndex(genomes)
    print(index.nb_genomes, "génomes indexés")