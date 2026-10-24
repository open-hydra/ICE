import sys, math, re
from pathlib import Path

def read_cell_vars(path):
    with open(path) as f:
        lines = f.readlines()

    # Parse I, J, K from ZONE line
    I = J = K = 0
    zone_idx = 0
    for idx, line in enumerate(lines):
        if re.match(r'\s*ZONE', line, re.IGNORECASE):
            zone_idx = idx
            for key, var in [('I', 'I'), ('J', 'J'), ('K', 'K')]:
                m = re.search(rf'\b{key}\s*=\s*(\d+)', line, re.IGNORECASE)
                if m:
                    exec(f'{var} = int(m.group(1))')
            I = int(re.search(r'\bI\s*=\s*(\d+)', line, re.IGNORECASE).group(1))
            J = int(re.search(r'\bJ\s*=\s*(\d+)', line, re.IGNORECASE).group(1))
            K = int(re.search(r'\bK\s*=\s*(\d+)', line, re.IGNORECASE).group(1))
            break

    # Read all numbers after the ZONE line
    numbers = []
    for line in lines[zone_idx + 1:]:
        for token in line.split():
            try:
                numbers.append(float(token))
            except ValueError:
                pass

    # BLOCK layout: x(I*J*K), y(I*J*K), z(I*J*K), then one (I-1)*(J-1)*max(1,K-1) block per variable
    N_node = I * J * K
    N_cell = (I - 1) * (J - 1) * max(1, K - 1)
    start  = 3 * N_node
    nvar   = (len(numbers) - start) // N_cell
    # Whole blocks only: a file cut short or padded holds no variable
    if len(numbers) != start + nvar * N_cell:
        nvar = 0
    return (I, J, K, nvar), [numbers[start + q * N_cell : start + (q + 1) * N_cell] for q in range(nvar)]


out_layout, out = read_cell_vars(Path('OUTPUT/part-field.tec'))
ref_layout, ref = read_cell_vars(Path('reference/part-field.tec'))
ref_rho = ref[0]

# zip() stops at the shorter list: an output on another grid, cut short or with other variables fails uncompared
if out_layout == ref_layout:
    l2 = math.sqrt(sum((a - b)**2 for a, b in zip(out[0], ref_rho)) / len(ref_rho))
    # n_p (the last variable) carries the table density through the inlet. It moves with the density along each
    # inlet's streamlines, so it is held to the density's tolerance relative to each field's own scale
    l2_n = math.sqrt(sum((a - b)**2 for a, b in zip(out[-1], ref[-1])) / sum(b * b for b in ref[-1]))
else:
    print(f'OUTPUT/part-field.tec: I, J, K, variables {out_layout}; reference/part-field.tec {ref_layout}')
    l2 = l2_n = math.inf

GREEN, RED, RESET = '\033[92m', '\033[91m', '\033[0m'
passed = l2 <= 1e-4 and l2_n <= 1e-4 / math.sqrt(sum(b * b for b in ref_rho) / len(ref_rho))
result = f'{GREEN}PASS{RESET}' if passed else f'{RED}FAIL{RESET}'
print(f'Doisneau MK  -->  {result}')
sys.exit(0 if passed else 1)
