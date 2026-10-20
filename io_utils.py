import pyfastx

def read_fof(path):
    with open(path) as f:
        return [line.strip() for line in f.read().splitlines() if line.strip()]



def read_genome(path):
    """Lit un fichier fasta et renvoie la liste de ses contigs.

    Chaque séquence du fichier est un élément de la liste, en majuscules.
    Un fichier à une seule séquence donne donc une liste d'un élément.
    Les contigs ne sont pas collés : c'est à l'index de les séparer par '$'.
    """
    return [seq.upper() for _, seq in pyfastx.Fastx(path)]

def read_queries(path):
    return  [(name ,seq.upper()) for name, seq in pyfastx.Fastx(path)]

def write_results(path, nb_genom, results):
    with open(path, mode="w") as f :
        header = "i_gen\t" + "\t".join(str(i+1) for i in range(nb_genom))
        f.write(header+"\n \n")
        for result in results :
            f.write(f"{result[0]} {"\t".join(str(i) for i in result[1])} \n")

#print(read_queries("data/q.fasta"))   # [('q1', 'ACGTA'), ('q2', 'TTTTG')]
write_results("data/test_result_bis.txt", 3, [("q1", [1, 0, 1]), ("q2", [0, 0, 1])])