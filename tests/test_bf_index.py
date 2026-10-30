from bf_index import BloomFMIndex


def build(genomes, k=4):
    index = BloomFMIndex(k, 1000, 3)
    for contigs in genomes:
        index.add_genome(contigs)
    return index


def test_one_value_per_genome():
    index = build([["ACGTACGT"], ["TTTTGGGG"], ["CCCCAAAA"]])
    assert len(index.bloom_filters) == 3
    assert len(index.fm_indexes) == 3
    assert len(index.query("ACGT")) == 3


def test_query_present_and_absent():
    index = build([["ACGTACGT"], ["TTTTGGGG"], ["ACGTTTTT"]])
    assert index.query("ACGTA") == [1, 0, 0]
    assert index.query("TTTT") == [0, 1, 1]
    assert index.query("GGGGA") == [0, 0, 0]


def test_order_follows_fof():
    # La i-ème valeur correspond au i-ème génome ajouté
    index = build([["TTTTGGGG"], ["ACGTACGT"]])
    assert index.query("ACGTA") == [0, 1]


def test_query_shorter_than_k():
    # Requête plus courte que k : seul le FM-index est interrogé
    index = build([["ACGTACGT"], ["TTTTGGGG"]], k=4)
    assert index.query("AC") == [1, 0]
    assert index.query("G") == [1, 1]


def test_query_not_across_contigs():
    # "GGTT" n'existe qu'à la jonction des deux contigs : il est absent
    index = build([["AAGG", "TTCC"]])
    assert index.query("GGTT") == [0]
    assert index.query("AAGG") == [1]
    assert index.query("TTCC") == [1]


def test_bloom_filter_false_positive_does_not_matter():
    # Filtre de 1 bit : il répond toujours oui, le FM-index doit corriger
    index = BloomFMIndex(4, 1, 1)
    index.add_genome(["ACGTACGT"])
    index.add_genome(["AAAACCCC"])
    assert index.bloom_filters[1].contains("ACGT")
    assert index.query("ACGTA") == [1, 0]
