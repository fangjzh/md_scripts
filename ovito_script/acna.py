# python script for acna analysis by using ovitos
# Syntax: ovitos acna.py

from ovito.io import *
from ovito.data import *
from ovito.modifiers import *
import numpy as np

node = import_file("dump2.atom", multiple_frames=True)

# Perform acna analysis:
acna = CommonNeighborAnalysisModifier(
    mode = CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff,
    )


node.modifiers.append(acna)

select = SelectExpressionModifier(expression="StructureType==3")
node.modifiers.append(select)

# Let OVITO do the computation and export the number of identified 
# defects as a function of simulation time to a text file:
#print('SelectExpression.num_selected')
#export_file(node, "Y-vac2.txt", "txt", 
#    columns = ['Timestep', 'SelectExpression.num_selected'],
#    multiple_frames = True)

# Export the XYZ coordinates of just the antisites by removing all other atoms.
#node.modifiers.append(InvertSelectionModifier())

node.modifiers.append(DeleteSelectedParticlesModifier())
export_file(node, "Pos_acna.atom", "lammps_dump", 
    columns = ['Position.X', 'Position.Y', 'Position.Z'],
    multiple_frames = True)

