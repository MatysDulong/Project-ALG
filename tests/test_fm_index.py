import random

from pydivsufsort import divsufsort

from fm_index import FMIndex


def build(text):
    """Construit le FM-index d'un texte terminé par '$'."""
    return FMIndex(text, divsufsort(text.encode()))


def naive_bwt(text):
    """BWT par tri de toutes les rotations (lent, mais sûr)."""
    rotations = sorted(text[i:] + text[:i] for i in range(len(text)))
    return "".join(rotation[-1] for rotation in rotations)


def naive_count(text, pattern):
    """Nombre d'occurrences de pattern dans text, chevauchements compris."""
    return sum(1 for i in range(len(text) - len(pattern) + 1)
               if text[i:i + len(pattern)] == pattern)


def test_bwt_example():
    # Exemple du cours, BWT vérifiée à la main
    assert build("ACAACG$").bwt == "GC$AAAC"


def test_bwt_random_texts():
    random.seed(1)
    for _ in range(20):
        text = "".join(random.choice("ACGT") for _ in range(50)) + "$"
        assert build(text).bwt == naive_bwt(text)


def test_first_occurrence():
    index = build("ACAACG$")
    # 1 '$', 3 'A', 2 'C', 1 'G'
    assert index.first_occurrence == {"$": 0, "A": 1, "C": 4, "G": 6, "T": 7}


def test_occ():
    index = build("ACAACG$")
    # occ[lettre][i] = nombre de lettre dans bwt[0:i], bwt = "GC$AAAC"
    assert index.occ["C"] == [0, 0, 1, 1, 1, 1, 1, 2]
    assert index.occ["A"] == [0, 0, 0, 0, 1, 2, 3, 3]


def test_backward_search_present():
    index = build("ACAACG$")
    top, bottom = index.backward_search("AC")
    assert bottom - top + 1 == 2
    top, bottom = index.backward_search("ACAACG")
    assert bottom - top + 1 == 1


def test_backward_search_absent():
    index = build("ACAACG$")
    assert index.backward_search("GG") is None
    assert index.backward_search("T") is None
    assert index.backward_search("ACAACGA") is None


def test_backward_search_empty_pattern():
    # Le motif vide est présent partout : toutes les lignes de la BWT
    index = build("ACAACG$")
    assert index.backward_search("") == (0, 6)


def test_backward_search_counts_random():
    random.seed(2)
    text = "".join(random.choice("ACGT") for _ in range(300)) + "$"
    index = build(text)
    for _ in range(200):
        length = random.randint(1, 6)
        pattern = "".join(random.choice("ACGT") for _ in range(length))
        expected = naive_count(text, pattern)
        interval = index.backward_search(pattern)
        if expected == 0:
            assert interval is None
        else:
            top, bottom = interval
            assert bottom - top + 1 == expected
