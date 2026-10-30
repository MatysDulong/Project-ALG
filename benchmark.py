"""Benchmark des deux approches, résultats écrits en CSV.

Deux expériences :
    compare : lance index_colored_fm_index.py et query.py (comme un
        utilisateur) pour les deux méthodes, en faisant varier le jeu de
        données, la taille des génomes, le nombre N de génomes et la longueur
        des requêtes. Mesure aussi la place de l'index en mémoire vive
        (pympler), une seule fois par configuration car c'est lent.
    bloom : fait varier les paramètres des filtres de Bloom (taille, nombre
        de hachages, k) et mesure leur effet sur le temps, la taille et le
        taux de faux positifs.

Usage :
    python benchmark.py [-e compare|bloom|all] [-r REPEATS] [-o OUT_DIR]
"""
import argparse
import csv
import multiprocessing
import os
import pickle
import subprocess
import sys
import tempfile
import time

from pympler import asizeof

import io_utils
from bf_index import BloomFMIndex

DATA_DIR = "data/test"
KINDS = ["coli", "random"]
SIZES = ["10k", "100k"]
NB_GENOMES = [5, 10, 20, 35, 49]
QUERY_LENGTHS = ["100bp", "1000bp"]

# Paramètres de référence de l'expérience bloom (un seul varie à la fois)
BLOOM_KIND, BLOOM_SIZE_DATA = "coli", "100k"
DEFAULT_K, DEFAULT_BLOOM_SIZE, DEFAULT_NUM_HASHES = 21, 2**20, 3
BLOOM_SIZES = [2**12, 2**14, 2**16, 2**18, 2**20, 2**22, 2**24]
NUM_HASHES = [1, 2, 3, 4, 6, 8]
K_VALUES = [11, 15, 21, 31, 51]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare les temps et la taille des deux index.")
    parser.add_argument("-e", "--experiment", default="all",
                        choices=["compare", "bloom", "all"],
                        help="expérience à lancer (défaut : all)")
    parser.add_argument("-r", "--repeats", type=int, default=1,
                        help="nombre de répétitions, on garde la moyenne "
                             "(défaut : 1)")
    parser.add_argument("-o", "--out_dir", default="results",
                        help="dossier des fichiers CSV (défaut : results)")
    args = parser.parse_args()
    if args.repeats <= 0:
        parser.error("--repeats doit être > 0")
    if not os.path.isdir(DATA_DIR):
        parser.error(f"le dossier de données '{DATA_DIR}' n'existe pas")
    return args


def run_program(command):
    """Lance un programme et renvoie ses lignes 'OUT NOM: valeur'."""
    output = subprocess.run([sys.executable] + command, capture_output=True,
                            text=True, check=True).stdout
    values = {}
    for line in output.splitlines():
        if line.startswith("OUT "):
            name, value = line[4:].split(":")
            values[name] = float(value)
    return values


def mean_of_runs(command, repeats):
    """Lance repeats fois le programme et fait la moyenne de chaque valeur."""
    runs = [run_program(command) for _ in range(repeats)]
    return {name: sum(run[name] for run in runs) / repeats
            for name in runs[0]}


def memory_of_index(index_file):
    """Relit l'index et renvoie sa taille en mémoire vive (octets)."""
    with open(index_file, "rb") as f:
        index = pickle.load(f)
    return asizeof.asizeof(index)


def measure_memory(index_file):
    """Lance memory_of_index dans un processus à part.

    Python ne rend pas toujours au système la mémoire libérée : en mesurant
    dans un autre processus, la mémoire ne s'accumule pas d'un index à
    l'autre pendant le benchmark.
    """
    with multiprocessing.Pool(1) as pool:
        return pool.apply(memory_of_index, (index_file,))


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} lignes écrites dans {path}")


def experiment_compare(repeats, out_dir, tmp_dir):
    """Compare generalized et bf avec les programmes en ligne de commande."""
    rows = []
    for kind in KINDS:
        for size in SIZES:
            genome_paths = io_utils.read_fof(
                f"{DATA_DIR}/fof_{kind}_{size}.txt")
            for nb_genomes in NB_GENOMES:
                # FOF réduit aux nb_genomes premiers génomes
                fof = os.path.join(tmp_dir, "fof.txt")
                with open(fof, "w") as f:
                    f.write("\n".join(genome_paths[:nb_genomes]) + "\n")

                for method in ["generalized", "bf"]:
                    index_file = os.path.join(tmp_dir, "index.dumped")
                    build = mean_of_runs(
                        ["index_colored_fm_index.py", "-i", fof,
                         "-o", index_file, "-m", method], repeats)
                    memory_size = measure_memory(index_file)

                    for length in QUERY_LENGTHS:
                        queries = (f"{DATA_DIR}/queries/queries_{kind}_"
                                   f"{size}_len{length}.fa")
                        query = mean_of_runs(
                            ["query.py", "-i", index_file, "-q", queries,
                             "-o", os.path.join(tmp_dir, "result.txt")],
                            repeats)
                        rows.append({
                            "kind": kind,
                            "genome_size": size,
                            "nb_genomes": nb_genomes,
                            "method": method,
                            "query_length": length,
                            "nb_queries": int(query["NB_QUERIES"]),
                            "time_build": build["TIME_BUILD"],
                            "time_serialisation":
                                build["TIME_SERIALISATION"],
                            "index_size": int(build["INDEX_SIZE"]),
                            "memory_size": memory_size,
                            "time_deserialisation":
                                query["TIME_DESERIALISATION"],
                            "time_query": query["TIME_QUERY"],
                        })
                    print(f"compare : {kind} {size} N={nb_genomes} "
                          f"{method} ok")
    write_csv(os.path.join(out_dir, "compare.csv"), rows)


def measure_bloom(genomes, queries, k, bloom_size, num_hashes):
    """Construit un BloomFMIndex et mesure temps, taille et faux positifs."""
    start = time.perf_counter()
    index = BloomFMIndex(k, bloom_size, num_hashes)
    for contigs in genomes:
        index.add_genome(contigs)
    time_build = time.perf_counter() - start

    start = time.perf_counter()
    for _, sequence in queries:
        index.query(sequence)
    time_query = time.perf_counter() - start

    # Faux positif du filtre : il dit oui alors que le premier k-mer est
    # vraiment absent du génome (vérifié directement sur les contigs).
    # fm_queried_rate : part des FM-index interrogés (filtre qui dit oui).
    # useless_rate : FM-index interrogés pour rien (requête finalement
    # absente), faux positifs du filtre compris.
    nb_yes = nb_false_positives = nb_kmer_absent = nb_useless = 0
    for _, sequence in queries:
        first_kmer = sequence[:k]
        presence = index.query(sequence)
        for bloom_filter, contigs, present in zip(index.bloom_filters,
                                                  genomes, presence):
            says_yes = bloom_filter.contains(first_kmer)
            nb_yes += says_yes
            nb_useless += says_yes and not present
            if not any(first_kmer in contig for contig in contigs):
                nb_kmer_absent += 1
                nb_false_positives += says_yes
    nb_pairs = len(queries) * len(genomes)

    return {
        "k": k,
        "bloom_size": bloom_size,
        "num_hashes": num_hashes,
        "time_build": time_build,
        "time_query": time_query,
        "index_size": len(pickle.dumps(index)),
        "false_positive_rate": nb_false_positives / max(nb_kmer_absent, 1),
        "fm_queried_rate": nb_yes / nb_pairs,
        "useless_rate": nb_useless / nb_pairs,
    }


def experiment_bloom(repeats, out_dir):
    """Fait varier un paramètre des filtres de Bloom à la fois."""
    fof = f"{DATA_DIR}/fof_{BLOOM_KIND}_{BLOOM_SIZE_DATA}.txt"
    genomes = [io_utils.read_genome(path) for path in io_utils.read_fof(fof)]
    queries = io_utils.read_queries(
        f"{DATA_DIR}/queries/queries_{BLOOM_KIND}_{BLOOM_SIZE_DATA}"
        f"_len100bp.fa")

    configurations = []
    for bloom_size in BLOOM_SIZES:
        configurations.append(("bloom_size", DEFAULT_K, bloom_size,
                               DEFAULT_NUM_HASHES))
    for num_hashes in NUM_HASHES:
        configurations.append(("num_hashes", DEFAULT_K, DEFAULT_BLOOM_SIZE,
                               num_hashes))
    for k in K_VALUES:
        configurations.append(("k", k, DEFAULT_BLOOM_SIZE,
                               DEFAULT_NUM_HASHES))

    rows = []
    for varied, k, bloom_size, num_hashes in configurations:
        runs = [measure_bloom(genomes, queries, k, bloom_size, num_hashes)
                for _ in range(repeats)]
        row = {"varied": varied}
        for name in runs[0]:
            row[name] = sum(run[name] for run in runs) / repeats
        rows.append(row)
        print(f"bloom : k={k} size={bloom_size} h={num_hashes} ok")
    write_csv(os.path.join(out_dir, "bloom.csv"), rows)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp_dir:
        if args.experiment in ["compare", "all"]:
            experiment_compare(args.repeats, args.out_dir, tmp_dir)
        if args.experiment in ["bloom", "all"]:
            experiment_bloom(args.repeats, args.out_dir)


if __name__ == "__main__":
    main()
