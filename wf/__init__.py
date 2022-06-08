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
    material: str = Material.rna,
    ensemble: str = Ensemble.stacking,
    temperature: float = 37,
    sodium: Optional[float] = 1.0,
    magnesium: Optional[float] = 0.0,
) -> LatchFile:

    nt_model = Model(material=material, ensemble=ensemble, celsius=temperature, sodium=sodium, magnesium=magnesium)

    dGloop = nt_model.loop_energy(loop=loop, structure=structure)
    dGloop = str(dGloop)
    # stackEnergies = nt_model.stack_energies(loop=loop, structure=structure)

    with open("/root/output", "w") as f:
        f.write(dGloop)

    return LatchFile("/root/output", "latch:///loop_energy/output.txt")

@workflow
def nupack_loops(
    loop: str = "AAAA",
    structure: str = "....",
    material: Material = Material.rna,
    temperature: float = 37.0,
    ensemble: Ensemble = Ensemble.stacking,
    sodium: Optional[float] = 1.0,
    magnesium: Optional[float] = 0.0,
) -> LatchFile:
    """Description...

    # NUPACK
    
    ---
    ## About

    <What does NUPACK do?> 

    __metadata__:
        display_name: Create DNA/RNA models and perform analysis for loop structures
        
        author:
            name: NUPACK
            email: support@nupack.org
            github: https://github.com/beliveau-lab/NUPACK
        
        repository: https://github.com/beliveau-lab/NUPACK
        
        license:
            id: BSD-3-Clause

    Args:
    
        material:
            __metadata__:
                display_name: "Nucleic Acid Type"

        temperature:
            __metadata__:
                display_name: "Temperature (in Celsius)"

        ensemble:
            __metadata__:
                display_name: "Ensemble Stacking Type"

        sodium:
            __metadata__:
                display_name: "Sodium concentration (in nM)"

        magnesium:
            __metadata__:
                display_name: "Magnesium concentration (in nM)"

        loop:
            __metadata__:
                display_name: "Loop sequence"

        structure:
            __metadata__:
                display_name: "Loop structure (in dot-bracket notation)"
    """
    return model_spec(
    loop=loop,
    structure=structure,
    material=material,
    temperature=temperature,
    ensemble=ensemble,
    sodium=sodium,
    magnesium=magnesium,
    )
