import io_utils
import fm_index
from pydivsufsort import divsufsort

#On  ́etend la notion de FM-index aux collections de s ́equences en construisant le FM-index (BW T, Count, Rank)1
#de la chaˆıne de caract`eres S obtenue en concat ́enant les N g ́enomes, s ́epar ́es par le caract`ere sp ́ecial $
#( ́egalement pr ́esent en fin de chaˆıne). On associe `a ce FM-index, un color array. C’est une table, similaire
#au suffix array, mais au lieu de stocker la position de chaque suffixe, on stocke l’identifiant (la couleur) du
#g ́enome (un entier dans [0, N −1]) o`u apparaˆıt chaque suffixe.

class GeneralizedIndex:
    def __init__(self, genomes):
        contig_list = []
        position = []
        for color, genome in enumerate(genomes):
            for contig in genome:
                contig_list.append(contig+"$")
                position.extend([color] * (len(contig)+1))
        S = "".join(contig_list)
        suffix_array = divsufsort(S)
        self.fm = fm_index.FMIndex(S, suffix_array)
        self.color_array = [position[p] for p in suffix_array]
        



if __name__ == "__main__":
    genomes = [["ACG"], ["CGT"], ["TAACG"]]
    test = GeneralizedIndex(genomes)
    print(test.fm.bwt)
    print(test.color_array)