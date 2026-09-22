import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from berthon import verify

# Two shocks and a contact. The largest error is on det(P), whose contact is
# the hardest feature on the mesh; measured 1.6e-2.
sys.exit(verify('SCS', tol=2.5e-2))
