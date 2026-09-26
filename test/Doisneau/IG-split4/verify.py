import sys, math, re
from pathlib import Path

def read_density(path):
    """Cell-centred density and n_p of every zone, each concatenated (DATAPACKING=BLOCK layout)."""
    rho, npart = [], []
    zones = re.split(r'^\s*ZONE', Path(path).read_text(), flags=re.IGNORECASE | re.MULTILINE)[1:]
    for z in zones:
        header, body = z.split('\n', 1)
        I, J, K = (int(re.search(rf'\b{k}\s*=\s*(\d+)', header, re.IGNORECASE).group(1)) for k in 'IJK')
        numbers = []
        for token in body.split():
            try:
                numbers.append(float(token))
            except ValueError:
                pass
        n_node = I * J * K
        n_cell = (I - 1) * (J - 1) * max(1, K - 1)
        rho += numbers[3 * n_node : 3 * n_node + n_cell]
        npart += numbers[len(numbers) - n_cell :]
    return rho, npart


out_rho, out_n = read_density('OUTPUT/part-field.tec')
ref_rho, ref_n = read_density('reference/part-field.tec')

n  = len(out_rho)
l2 = math.sqrt(sum((a - b)**2 for a, b in zip(out_rho, ref_rho)) / n) if n == len(ref_rho) else math.inf
# n_p (the last variable) carries the table density through the inlet. It moves with the density along each
# inlet's streamlines, so it is held to the density's tolerance relative to each field's own scale
l2_n = (math.sqrt(sum((a - b)**2 for a, b in zip(out_n, ref_n)) / sum(b * b for b in ref_n))
        if len(out_n) == len(ref_n) else math.inf)

GREEN, RED, RESET = '\033[92m', '\033[91m', '\033[0m'
passed = l2 <= 1e-4 and l2_n <= 1e-4 / math.sqrt(sum(b * b for b in ref_rho) / len(ref_rho))
result = f'{GREEN}PASS{RESET}' if passed else f'{RED}FAIL{RESET}'
print(f'Doisneau IG-split4   -->  {result}')
sys.exit(0 if passed else 1)
