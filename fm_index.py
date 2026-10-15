from pydivsufsort import divsufsort

ALPHABET = "$ACGT"


class FMIndex:
    """FM-index d'un texte terminé par '$' : BWT, Count (N) et Rank (occ)."""

    def __init__(self, text, suffix_array):
        """Construit l'index à partir du texte et de son suffix array.

        Ni le texte ni le suffix array ne sont conservés.
        """
        self.bwt = "".join(text[i - 1] for i in suffix_array)
        self.occ, letter_counts = self.compute_occ()
        self.first_occurrence = self.compute_first_occurrence(letter_counts)

    def compute_first_occurrence(self, counts):
        """Renvoie N : N[lettre] = nombre de lettres plus petites."""
        first_occurrence = {}
        total = 0
        for letter in ALPHABET:
            first_occurrence[letter] = total
            total += counts[letter]
        return first_occurrence

    def compute_occ(self):
        """Renvoie (occ, counts) : occ[lettre][i] = nb de lettre dans bwt[0:i].

        counts[lettre] est le nombre total de lettre dans la BWT.
        """
        occ = {letter: [0] for letter in ALPHABET}
        counts = {letter: 0 for letter in ALPHABET}
        for bwt_letter in self.bwt:
            counts[bwt_letter] += 1
            for letter in ALPHABET:
                occ[letter].append(counts[letter])
        return occ, counts

    def backward_search(self, pattern):
        """Renvoie (top, bottom) : lignes dont le suffixe commence par pattern.

        Renvoie None si pattern est absent. Complexité : O(len(pattern)).
        """
        top = 0
        bottom = len(self.bwt) - 1
        for letter in reversed(pattern):
            top = self.first_occurrence[letter] + self.occ[letter][top]
            bottom = self.first_occurrence[letter] + self.occ[letter][bottom + 1] - 1
            if top > bottom:
                return None
        return (top, bottom)


if __name__ == "__main__":
    s = 'ACAACG$'
    sa = divsufsort(s)
    test = FMIndex(s, sa)
    print(test.bwt)
    print(test.first_occurrence)
    print(test.occ["C"])

    for motif in ["AC", "CA", "GG", "ACAACG", "T", ""]:
        print(motif, test.backward_search(motif))