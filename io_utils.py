import pyfastx

def read_fof(path):
    with open(path) as f:
        return [line.strip() for line in f.read().splitlines() if line.strip()]



def read_genome(path):
    sequences = [seq for _, seq in pyfastx.Fastx(path)]
    return "".join(sequences).upper()