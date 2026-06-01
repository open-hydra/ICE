import sys, math, re
from pathlib import Path

def read_density(path):
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

    # BLOCK layout: x(I*J*K), y(I*J*K), z(I*J*K), rho((I-1)*(J-1)*max(1,K-1)), ...
    N_node = I * J * K
    N_cell = (I - 1) * (J - 1) * max(1, K - 1)
    start  = 3 * N_node
    return numbers[start : start + N_cell]


out_rho = read_density(Path('OUTPUT/part-field.tec'))
ref_rho = read_density(Path('reference/part-field.tec'))

n  = len(out_rho)
l2 = math.sqrt(sum((a - b)**2 for a, b in zip(out_rho, ref_rho)) / n)

GREEN, RED, RESET = '\033[92m', '\033[91m', '\033[0m'
passed = l2 <= 1e-10
result = f'{GREEN}PASS{RESET}' if passed else f'{RED}FAIL{RESET}'
print(f'Doisneau AG  -->  {result}')
sys.exit(0 if passed else 1)
