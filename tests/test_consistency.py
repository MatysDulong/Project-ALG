"""Vérifie que les deux approches renvoient exactement les mêmes résultats."""
import os
import random
import subprocess
import sys

import pytest

import io_utils
from bf_index import BloomFMIndex
from generalized_index import GeneralizedIndex

FOF_COLI = "data/test/fof_coli_10k.txt"
QUERIES_COLI = "data/test/queries/queries_coli_10k_len100bp.fa"


def build_both(genomes, k=21):
    generalized = GeneralizedIndex(genomes)
    bf = BloomFMIndex(k, 100000, 3)
    for contigs in genomes:
        bf.add_genome(contigs)
    return generalized, bf


def naive_query(genomes, sequence):
    """Réponse attendue, calculée directement sur les contigs."""
    return [int(any(sequence in contig for contig in contigs))
            for contigs in genomes]


def test_random_genomes():
    random.seed(5)
    base = "".join(random.choice("ACGT") for _ in range(2000))
    genomes = [[base[:1500]], [base[500:]], [base[:700], base[1200:]]]
    # Génomes partiellement communs pour avoir des requêtes à 1 ou 2 couleurs
    genomes.append(["".join(random.choice("ACGT") for _ in range(1000))])
    generalized, bf = build_both(genomes, k=11)

    queries = []
    for _ in range(300):
        contigs = random.choice(genomes)
        contig = random.choice(contigs)
        length = random.randint(1, 60)
        start = random.randint(0, len(contig) - length)
        queries.append(contig[start:start + length])
    queries += ["".join(random.choice("ACGT") for _ in range(30))
                for _ in range(50)]

    for sequence in queries:
        expected = naive_query(genomes, sequence)
        assert generalized.query(sequence) == expected
        assert bf.query(sequence) == expected


@pytest.mark.skipif(not os.path.isfile(FOF_COLI),
                    reason="jeux de données de test absents")
def test_coli_10k():
    genomes = [io_utils.read_genome(path)
               for path in io_utils.read_fof(FOF_COLI)]
    generalized, bf = build_both(genomes)
    queries = io_utils.read_queries(QUERIES_COLI)[:200]
    for _, sequence in queries:
        expected = naive_query(genomes, sequence)
        assert generalized.query(sequence) == expected
        assert bf.query(sequence) == expected


@pytest.mark.skipif(not os.path.isfile(FOF_COLI),
                    reason="jeux de données de test absents")
def test_command_line_same_output(tmp_path):
    # Teste les deux programmes imposés, de bout en bout
    outputs = []
    for method in ["generalized", "bf"]:
        index_file = tmp_path / f"{method}.dumped"
        result_file = tmp_path / f"{method}.txt"
        build = subprocess.run(
            [sys.executable, "index_colored_fm_index.py", "-i", FOF_COLI,
             "-o", str(index_file), "-m", method],
            capture_output=True, text=True, check=True)
        assert "OUT TIME_BUILD:" in build.stdout
        assert "OUT TIME_SERIALISATION:" in build.stdout
        query = subprocess.run(
            [sys.executable, "query.py", "-i", str(index_file),
             "-q", QUERIES_COLI, "-o", str(result_file)],
            capture_output=True, text=True, check=True)
        assert "OUT TIME_DESERIALISATION:" in query.stdout
        assert "OUT TIME_QUERY:" in query.stdout
        outputs.append(result_file.read_text())
    assert outputs[0] == outputs[1]
