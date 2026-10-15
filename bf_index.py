from bloom_filter import BloomFilter
from fm_index import FMIndex
from pydivsufsort import divsufsort


class BloomFMIndex:
    """Index de N génomes : un filtre de Bloom et un FM-index par génome.

    Les deux listes sont alignées : l'élément i de chacune décrit le
    génome i, dans l'ordre du FOF.
    """

    def __init__(self, k, bloom_size, num_hashes):
        """Enregistre les paramètres communs à tous les filtres.

        Args:
            k (int) : longueur des k-mers insérés dans les filtres.
            bloom_size (int) : nombre de bits de chaque filtre.
            num_hashes (int) : nombre de fonctions de hachage par filtre.
        """
        self.k = k
        self.bloom_size = bloom_size
        self.num_hashes = num_hashes
        self.bloom_filters = []
        self.fm_indexes = []

    def add_genome(self, genome):
        """Indexe un génome : ajoute son filtre de Bloom et son FM-index.

        Les deux objets sont ajoutés en fin de liste, donc dans l'ordre
        d'appel (celui du FOF).

        Args:
            genome (str) : séquence du génome, en majuscules, sur
                l'alphabet A, C, G, T et sans '$'.
        """
        # Le filtre ne reçoit que des k-mers du génome : pas de '$'.
        bloom_filter = BloomFilter(self.bloom_size, self.num_hashes, self.k)
        bloom_filter.add_sequence(genome)

        # Le FM-index exige un texte terminé par '$'.
        text = genome + "$"
        fm_index = FMIndex(text, divsufsort(text))

        # Ajouts groupés en dernier pour garder les deux listes alignées.
        self.bloom_filters.append(bloom_filter)
        self.fm_indexes.append(fm_index)


    def query(self, sequence):
        """Indique dans quels génomes la séquence apparaît en entier.

        Le filtre de Bloom écarte vite les génomes qui ne contiennent
        sûrement pas le premier k-mer de la requête. Le FM-index confirme
        pour les autres : le résultat est donc exact, même si le filtre
        donne des faux positifs.

        Une requête plus courte que k n'a pas de premier k-mer : le filtre
        est alors sauté et seul le FM-index est interrogé.

        Args:
            sequence (str) : requête en majuscules, sur l'alphabet ACGT.

        Returns:
            list[int] : une valeur par génome, dans l'ordre du FOF,
            1 si la séquence est présente, 0 sinon.
        """
        results = []
        is_short = len(sequence) < self.k
        first_kmer = sequence[:self.k]

        for bloom_filter, fm_index in zip(self.bloom_filters, self.fm_indexes):
            if is_short or bloom_filter.contains(first_kmer):
                is_found = fm_index.backward_search(sequence) is not None
                results.append(int(is_found))
            else:
                results.append(0)
        return results