from pydivsufsort import divsufsort
import io_utils

ALPHABET = "$ACGT"


class FMIndex:
    """FM-index d'un texte terminé par '$' (sans suffix array conservé).

    Attributs :
        bwt (str) : transformée de Burrows-Wheeler du texte.
        occ (dict[str, list[int]]) : occ[lettre][i] = nombre de lettre
            dans bwt[0:i] (structure Rank), len(bwt) + 1 cases par lettre.
        first_occurrence (dict[str, int]) : N[lettre] = nombre de lettres
            strictement plus petites dans le texte (structure Count).
    """

    def __init__(self, text, suffix_array):
        """Construit la BWT, N (Count) et occ (Rank) du texte.

        Args:
            text (str) : texte terminé par '$', sur l'alphabet ALPHABET.
            suffix_array : suffix array de text (calculé par l'appelant).

        Ni text ni suffix_array ne sont conservés : ils ne partent donc
        pas dans le fichier sérialisé.
        """
        self.bwt = "".join(text[i - 1] for i in suffix_array)
        self.occ, letter_counts = self.compute_occ()
        self.first_occurrence = self.compute_first_occurrence(letter_counts)

    def compute_first_occurrence(self, counts):
        """Renvoie N : N[lettre] = nombre de lettres strictement plus petites.

        C'est aussi la ligne de F où apparaît la première occurrence de lettre.
        """
        first_occurrence = {}
        total = 0
        for letter in ALPHABET:
            first_occurrence[letter] = total
            total += counts[letter]
        return first_occurrence

    def compute_occ(self):
        """Renvoie (occ, counts) en un seul passage sur la BWT.

        occ[lettre][i] est le nombre de lettre dans bwt[0:i], c'est-à-dire
        strictement avant la ligne i ; chaque liste a len(bwt) + 1 cases.
        counts[lettre] est le nombre total de cette lettre dans la BWT.
        """
        occ = {letter: [0] for letter in ALPHABET}
        counts = {letter: 0 for letter in ALPHABET}
        for bwt_letter in self.bwt:
            counts[bwt_letter] += 1
            for letter in ALPHABET:
                occ[letter].append(counts[letter])
        return occ, counts

    def backward_search(self, pattern):
        """Cherche pattern par backward search, en O(len(pattern)).

        Args:
            pattern (str) : motif à chercher, sur l'alphabet ALPHABET.

        Returns:
            tuple[int, int] | None : intervalle (top, bottom), bornes
            incluses, des lignes de la BWT dont le suffixe commence par
            pattern ; None si pattern est absent. Le motif vide renvoie
            toutes les lignes, (0, len(bwt) - 1).
        """
        top = 0
        bottom = len(self.bwt) - 1
        for c in reversed(pattern):
            top = self.first_occurrence[c] + self.occ[c][top]
            bottom = self.first_occurrence[c] + self.occ[c][bottom + 1] - 1
            if top > bottom:
                return None
        return (top, bottom)


if __name__ == "__main__":
    ### En attendant d'avoir io_utils
    s = 'ACAACG$'
    sa = divsufsort(s)
    test = FMIndex(s, sa)
    print(test.bwt)
    print(test.first_occurrence)
    print(test.occ["C"])

    for motif in ["AC", "CA", "GG", "ACAACG", "T", ""]:
        print(motif, test.backward_search(motif))


    #source ~/envs/alg_projet/Scripts/activate