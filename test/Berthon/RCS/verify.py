import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from berthon import verify

# Rarefaction, contact and shock. Measured worst error 7.4e-3.
sys.exit(verify('RCS', tol=1.5e-2))
