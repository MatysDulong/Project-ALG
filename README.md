# Projet ALG 2026 : FM-index coloré

Matys Dulong et Corentin Perrot, M2 Bioinformatique Rennes

Ce projet indexe un ensemble de N génomes pour savoir, pour une séquence
requête, dans quel(s) génome(s) elle apparaît sans erreur. Deux approches
sont implémentées et comparées :

- **`generalized`** : un seul FM-index (BWT, Count, Rank) de la
  concaténation de tous les génomes, séparés par `$`, et un *color array*
  qui donne pour chaque ligne de la BWT le génome du suffixe ;
- **`bf`** : un FM-index et un filtre de Bloom par génome. Pour une
  requête, seuls les FM-index dont le filtre contient le premier k-mer de
  la requête sont interrogés.

Les deux approches donnent exactement les mêmes résultats (vérifié par les
tests). Le rapport (`report/rapport.pdf`) compare leurs temps et leur
occupation mémoire.


## Installation

Python 3.12 et les bibliothèques de `requirements.txt` (pyfastx,
pydivsufsort, bitarray, mmh3, pympler, pytest). Avec conda :

```
conda create -n alg_projet python=3.12
conda activate alg_projet
pip install -r requirements.txt
```


## Utilisation

Toutes les commandes se lancent depuis la racine du dépôt.

### Construire l'index

```
python index_colored_fm_index.py -i fof_file -o serialized_file -m <generalized|bf> [-h] [options bf]
```

| Option | Rôle |
|---|---|
| `-i` | FOF : un chemin de fichier FASTA par ligne, un génome par fichier |
| `-o` | fichier de sortie de l'index sérialisé (pickle) |
| `-m` | `generalized` (approche 1) ou `bf` (approche 2) |
| `-h` | affiche l'aide |
| `-k` | (bf) longueur des k-mers des filtres de Bloom, défaut 21 |
| `--bloom_size` | (bf) nombre de bits de chaque filtre, défaut 2^24 |
| `--num_hashes` | (bf) nombre de fonctions de hachage, défaut 3 |

Chaque fichier FASTA peut contenir plusieurs séquences (contigs) : elles
sont indexées séparément, une requête à cheval sur deux contigs n'est donc
pas trouvée. L'ordre des génomes dans le FOF donne l'ordre des colonnes du
résultat.

Exemple :

```
python index_colored_fm_index.py -i data/test/fof_coli_100k.txt -o coli.dumped -m bf --bloom_size 1048576
OUT TIME_BUILD: 9.079290577999927
OUT TIME_SERIALISATION: 1.3142513820021122
OUT NB_GENOMES: 49
OUT INDEX_SIZE: 78273780
```

`TIME_BUILD` ne compte pas la lecture des fichiers FASTA, seulement la
construction de l'index. `INDEX_SIZE` est la taille du fichier produit, en
octets.

### Interroger l'index

```
python query.py -i serialized_file -q query.fa -o result.txt [-h]
```

`query.fa` est un fichier FASTA de séquences sur l'alphabet A, C, G, T.
Le fichier `result.txt` contient un en-tête, une ligne vide, puis une ligne
par requête (colonnes séparées par des tabulations) : 1 si la requête est
présente dans le génome i, 0 sinon.

```
python query.py -i coli.dumped -q data/test/queries/queries_coli_100k_len1000bp.fa -o result.txt
OUT TIME_DESERIALISATION: 1.5164086200020392
OUT TIME_QUERY: 0.25634961900141207
OUT NB_QUERIES: 196
```

Début de `result.txt` (7 premières colonnes sur 49) :

```
i_gen	1	2	3	4	5	6	7

random_0_head_042_GCA_000027125.fa	0	1	0	0	0	0	0
random_1_head_042_GCA_000027125.fa	0	1	0	0	0	0	0
```

Les temps sont en secondes.


## Paramètres conseillés

**Quelle approche ?** (détails dans le rapport)

- `generalized` se construit environ 2 fois plus vite et fait une seule
  recherche par requête (contre 4 en moyenne pour `bf` sur les génomes
  E. coli) : c'est le choix conseillé pour le temps ;
- `bf`, avec des filtres bien dimensionnés, produit un fichier plus petit
  (78 Mo contre 118 Mo pour 49 génomes de 100 kb) et permet d'ajouter un
  génome sans tout reconstruire.

**Paramètres des filtres de Bloom (`bf`)**, mesurés sur 49 génomes d'E. coli
de 100 kb :

- `--bloom_size` : la valeur par défaut (2^24 bits = 2 Mo par filtre) est
  **trop grande pour les petits génomes**. Nous conseillons environ 10 bits
  par k-mer du plus gros génome, soit 2^20 (1048576) pour des génomes de
  100 kb (1,5 % de faux positifs) et 2^17 (131072) pour des génomes de
  10 kb. Le taux de faux positifs mesuré suit la formule théorique
  (1 - e^(-hn/m))^h, qui permet de choisir m pour d'autres tailles.
  Jusqu'à 2^16 bits pour 100 kb, les filtres sont saturés et les requêtes
  deviennent 2 à 3 fois plus lentes.
- `--num_hashes` : 3 (ou 4). Au-delà, les faux positifs baissent peu
  (1,5 % pour h = 3, 0,55 % pour h = 8) mais la construction est environ
  1,5 fois plus longue (8,6 s pour h = 3, 12,5 s pour h = 8).
- `-k` : 21 à 31. Un k plus grand rend le premier k-mer plus spécifique et
  évite d'interroger des FM-index pour rien (10 % des FM-index interrogés
  pour rien avec k = 11, 4 % avec k = 31). **k ne doit pas dépasser la
  longueur des requêtes** : pour une requête plus courte que k, le filtre
  n'est pas utilisé et tous les FM-index sont interrogés (le résultat reste
  correct).

Les faux positifs du filtre ne changent jamais le résultat : ils ne font
qu'ajouter des recherches inutiles dans le FM-index, qui les corrige.


## Tests

```
pytest
```

23 tests (`tests/`) : BWT, tables Count et Rank et recherche arrière du
FM-index, filtre de Bloom (pas de faux négatifs, taux de faux positifs),
index `bf` (contigs, requêtes plus courtes que k) et surtout
`test_consistency.py`, qui vérifie que les deux approches et une recherche
naïve donnent le même résultat, sur des génomes aléatoires et sur les
données E. coli, et que les deux programmes en ligne de commande produisent
des fichiers de résultats identiques. Les tests sur les données réelles
sont sautés si `data/` est absent.


## Reproduire les résultats du rapport

1. **Données** : placer les archives et les requêtes fournies sur Moodle
   dans `data/test/` (dossier ignoré par git) :

   ```
   data/test/coli_genomes/10k_coli_genomes.zip
   data/test/coli_genomes/100k_coli_genomes.zip
   data/test/random_genomes/10k_random_genomes.zip
   data/test/random_genomes/100k_random_genomes.zip
   data/test/queries/queries_{coli,random}_{10k,100k}_len{100,1000}bp.fa
   ```

   puis extraire les archives et créer un FOF par collection :

   ```
   for z in data/test/*_genomes/*.zip; do unzip -q "$z" -d "$(dirname "$z")"; done
   for d in data/test/*_genomes/*_genomes; do
       name=$(basename "$d")                     # ex. 10k_coli_genomes
       size=${name%%_*}; kind=$(echo "$name" | cut -d_ -f2)
       ls "$d"/*.fa | sort -V > "data/test/fof_${kind}_${size}.txt"
   done
   ```

2. **Benchmark** (environ 30 min) :

   ```
   python benchmark.py -r 3
   ```

   - `compare` lance les deux programmes pour chaque collection, avec
     N = 5, 10, 20, 35 et 49 génomes (les N premiers du FOF), et les deux
     longueurs de requêtes, avec les paramètres par défaut (k = 21, h = 3,
     m = 2^24). Il mesure les quatre temps, la taille du fichier et la
     taille de l'index en mémoire vive (`pympler.asizeof`, une seule fois
     par configuration) ;
   - `bloom` fait varier un paramètre des filtres à la fois autour de
     k = 21, m = 2^20, h = 3, sur les 49 génomes E. coli de 100 kb et les
     1960 requêtes de 100 pb.

   Les temps sont la moyenne des 3 répétitions (`-r`). Résultats dans
   `results/compare.csv` et `results/bloom.csv`. Une seule expérience :
   `-e compare` ou `-e bloom`.

3. **Rapport** : `report/rapport.tex`, compilé avec
   [tectonic](https://tectonic-typesetting.github.io)
   (`tectonic report/rapport.tex`).

Machine utilisée : Linux, 8 cœurs, 8 Go de RAM, Python 3.12.16. Les temps
absolus dépendent de la machine et de sa charge ; ce sont les comparaisons
entre approches qui comptent.


## Organisation du code

| Fichier | Contenu |
|---|---|
| `index_colored_fm_index.py` | programme de construction de l'index |
| `query.py` | programme de requête |
| `fm_index.py` | classe `FMIndex` : BWT, Count, Rank (`occ`), recherche arrière |
| `generalized_index.py` | approche 1 : `GeneralizedIndex` (FM-index + color array) |
| `bloom_filter.py` | classe `BloomFilter` (bitarray + mmh3) |
| `bf_index.py` | approche 2 : `BloomFMIndex` (un filtre et un FM-index par génome) |
| `io_utils.py` | lecture du FOF, des FASTA, écriture des résultats |
| `benchmark.py` | expériences du rapport |
| `tests/` | tests pytest |
| `results/` | résultats bruts du benchmark (CSV) |
| `report/` | rapport (PDF et source LaTeX) et figures |
