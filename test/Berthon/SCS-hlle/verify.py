import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from berthon import verify

# The same problem as SCS, solved with the HLLE flux instead of the per-closure
# default. HLLE is less dissipative on the contact, so the density error is lower
# than Rusanov's 4.3e-3; the tolerance is SCS's, which it must not exceed.
sys.exit(verify('SCS', tol=2.5e-2))
