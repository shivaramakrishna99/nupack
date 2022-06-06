"""
Predict the folded structure of an RNA sequence
"""

import subprocess
from enum import Enum

from enum import Enum
from pathlib import Path
import subprocess
from typing import Optional

from latch import small_task, workflow
from latch.types import LatchFile


from nupack import *  # Import NUPACK

class Material(Enum):
    dna = "dna"
    rna = "rna"
    rna95 = "rna95"

class Ensemble(Enum):
    stacking = "stacking"
    nostacking = "nostacking"
    
@small_task
def model_spec(
    loop: str = "AAAA",
    structure: str = "....",
    material: Material = Material.rna,
    temperature: float = 37,
    ensemble: Ensemble = Ensemble.stacking,
    sodium: Optional[float] = 1.0,
    magnesium: Optional[float] = 0.0,
    outputFile: Optional[str] = None
) -> LatchFile:

    if not outputFile:
        out = Path("output.txt").resolve()
    else:
        out = outputFile

    nt_model = Model(material=material, ensemble=ensemble, temperature=temperature, sodium=sodium, magnesium=magnesium)

    dGloop = nt_model.loop_energy(loop=loop, structure=structure)
    # stackEnergies = nt_model.stack_energies(loop=loop, structure=structure)

    with open(f"/root/{out}", "w") as f:
        f.write(dGloop)

    return LatchFile("/root/output.txt", "latch:///output.txt")

@workflow
def nupack_loops(
    loop: str = "AAAA",
    structure: str = "....",
    material: Material = Material.rna,
    temperature: float = 37,
    ensemble: Ensemble = Ensemble.stacking,
    sodium: Optional[float] = 1.0,
    magnesium: Optional[float] = 0.0,
    outputFile: Optional[str] = None
) -> LatchFile:
    """Description...

    # NUPACK
    
    ---
    ## About

    <What does NUPACK do?> 

    __metadata__:
        display_name: Predict the folded structure of an RNA sequence
        author: 
            name: Name
            email: 
            github: 
        repository:
        license:
            id: MIT

    Args:
    
        material:
            Specify nucleic acid type as DNA or RNA
            __metadata__:
                display_name: "Nucleic Acid Type"
                appearance:
    """
    return model_spec(
    loop=loop,
    structure=structure,
    material=material,
    temperature=temperature,
    ensemble=ensemble,
    sodium=sodium,
    magnesium=magnesium,
    outputFile=outputFile
    )
