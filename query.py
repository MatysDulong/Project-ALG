import argparse
import os
import pickle
import time

import io_utils
from bf_index import BloomFMIndex


def parse_args():
    """Lit les arguments de la ligne de commande et les vérifie."""
    parser = argparse.ArgumentParser(
        description="Recherche des séquences requêtes dans un index "
                    "sérialisé par index_colored_fm_index.py.")
    parser.add_argument("-i", "--input", required=True,
                        help="fichier de l'index sérialisé")
    parser.add_argument("-q", "--queries", required=True,
                        help="fichier fasta des séquences requêtes")
    parser.add_argument("-o", "--output", required=True,
                        help="fichier texte des résultats")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        parser.error(f"l'index '{args.input}' n'existe pas")
    if not os.path.isfile(args.queries):
        parser.error(f"le fichier de requêtes '{args.queries}' n'existe pas")
    return args


def get_nb_genomes(index):
    """Renvoie le nombre de génomes indexés, quelle que soit la méthode."""
    if isinstance(index, BloomFMIndex):
        return len(index.fm_indexes)
    return index.nb_genomes


def main():
    args = parse_args()

    start = time.perf_counter()
    with open(args.input, "rb") as f:
        index = pickle.load(f)
    time_deserialisation = time.perf_counter() - start

    queries = io_utils.read_queries(args.queries)

    # Les deux index ont la même méthode query() : pas besoin de les distinguer
    start = time.perf_counter()
    results = [(name, index.query(sequence)) for name, sequence in queries]
    time_query = time.perf_counter() - start

    io_utils.write_results(args.output, get_nb_genomes(index), results)

    print(f"OUT TIME_DESERIALISATION: {time_deserialisation}")
    print(f"OUT TIME_QUERY: {time_query}")
    print(f"OUT NB_QUERIES: {len(queries)}")


if __name__ == "__main__":
    main()
