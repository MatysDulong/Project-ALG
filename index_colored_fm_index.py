import argparse
import os
import pickle
import sys
import time

import io_utils
from bf_index import BloomFMIndex
from generalized_index import GeneralizedIndex


def parse_args():
    """Lit les arguments de la ligne de commande et les vérifie."""
    parser = argparse.ArgumentParser(
        description="Construit un FM-index coloré à partir d'un FOF de "
                    "génomes et l'enregistre avec pickle.")
    parser.add_argument("-i", "--input", required=True,
                        help="fichier FOF (un chemin de fasta par ligne)")
    parser.add_argument("-o", "--output", required=True,
                        help="fichier de sortie de l'index sérialisé")
    parser.add_argument("-m", "--method", required=True,
                        choices=["generalized", "bf"],
                        help="generalized : FM-index généralisé + color "
                             "array, bf : N FM-index + N filtres de Bloom")
    # Paramètres utilisés seulement par la méthode bf
    parser.add_argument("-k", type=int, default=21,
                        help="taille des k-mers du filtre de Bloom "
                             "(bf, défaut 21)")
    parser.add_argument("--bloom_size", type=int, default=2**24,
                        help="nombre de bits de chaque filtre de Bloom "
                             "(bf, défaut 2^24)")
    parser.add_argument("--num_hashes", type=int, default=3,
                        help="nombre de fonctions de hachage "
                             "(bf, défaut 3)")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        parser.error(f"le fichier FOF '{args.input}' n'existe pas")
    if args.k <= 0 or args.bloom_size <= 0 or args.num_hashes <= 0:
        parser.error("-k, --bloom_size et --num_hashes doivent être > 0")
    return args


def load_genomes(fof_path):
    """Renvoie la liste des génomes (listes de contigs) dans l'ordre du FOF."""
    genome_paths = io_utils.read_fof(fof_path)
    if len(genome_paths) == 0:
        sys.exit(f"Erreur : le FOF '{fof_path}' est vide")
    for path in genome_paths:
        if not os.path.isfile(path):
            sys.exit(f"Erreur : le génome '{path}' du FOF n'existe pas")
    return [io_utils.read_genome(path) for path in genome_paths]


def build_index(genomes, args):
    """Construit l'index avec la méthode choisie (-m)."""
    if args.method == "generalized":
        return GeneralizedIndex(genomes)
    index = BloomFMIndex(args.k, args.bloom_size, args.num_hashes)
    for genome in genomes:
        index.add_genome(genome)
    return index


def main():
    args = parse_args()

    # La lecture des fasta n'est pas comptée dans TIME_BUILD
    genomes = load_genomes(args.input)

    start = time.perf_counter()
    index = build_index(genomes, args)
    time_build = time.perf_counter() - start

    start = time.perf_counter()
    with open(args.output, "wb") as f:
        pickle.dump(index, f)
    time_serialisation = time.perf_counter() - start

    print(f"OUT TIME_BUILD: {time_build}")
    print(f"OUT TIME_SERIALISATION: {time_serialisation}")
    print(f"OUT NB_GENOMES: {len(genomes)}")
    print(f"OUT INDEX_SIZE: {os.path.getsize(args.output)}")


if __name__ == "__main__":
    main()
