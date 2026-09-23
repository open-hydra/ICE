import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from berthon import verify

# Two rarefactions astride the contact. Measured worst error 8.0e-3.
sys.exit(verify('RCR', tol=1.5e-2))
