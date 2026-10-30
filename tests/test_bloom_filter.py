import random

from bloom_filter import BloomFilter


def random_kmer(k):
    return "".join(random.choice("ACGT") for _ in range(k))


def test_empty_filter_contains_nothing():
    bloom_filter = BloomFilter(1000, 3, 5)
    assert not bloom_filter.contains("ACGTA")
    assert bloom_filter.bits.count(1) == 0


def test_added_kmer_is_contained():
    bloom_filter = BloomFilter(1000, 3, 5)
    bloom_filter.add("ACGTA")
    assert bloom_filter.contains("ACGTA")


def test_no_false_negative():
    # Un filtre de Bloom ne doit jamais oublier un k-mer ajouté
    random.seed(3)
    bloom_filter = BloomFilter(5000, 3, 21)
    kmers = [random_kmer(21) for _ in range(500)]
    for kmer in kmers:
        bloom_filter.add(kmer)
    assert all(bloom_filter.contains(kmer) for kmer in kmers)


def test_add_sequence_adds_all_kmers():
    sequence = "ACGTTGCAAC"
    bloom_filter = BloomFilter(1000, 3, 4)
    bloom_filter.add_sequence(sequence)
    for i in range(len(sequence) - 4 + 1):
        assert bloom_filter.contains(sequence[i:i + 4])


def test_add_sequence_shorter_than_k():
    # Pas de k-mer dans une séquence plus courte que k
    bloom_filter = BloomFilter(1000, 3, 5)
    bloom_filter.add_sequence("ACG")
    assert bloom_filter.bits.count(1) == 0


def test_false_positive_rate_is_low():
    # Avec 100 bits par k-mer et 3 hachages, le taux théorique est ~0.003 %
    random.seed(4)
    bloom_filter = BloomFilter(100000, 3, 21)
    added = {random_kmer(21) for _ in range(1000)}
    for kmer in added:
        bloom_filter.add(kmer)
    tested = [random_kmer(21) for _ in range(10000)]
    false_positives = sum(1 for kmer in tested
                          if kmer not in added and bloom_filter.contains(kmer))
    assert false_positives / len(tested) < 0.01
