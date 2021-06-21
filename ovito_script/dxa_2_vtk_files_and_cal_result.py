# Boilerplate code generated by OVITO Pro 3.3.0
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *

# Data import:
pipeline = import_file('./dump.atom.gz',multiple_frames=True)

# Dislocation analysis (DXA):
pipeline.modifiers.append(DislocationAnalysisModifier(input_crystal_structure = DislocationAnalysisModifier.Lattice.BCC))

export_file(pipeline, "./dislocations.vtk", "vtk/disloc" ,frame = 40001)

