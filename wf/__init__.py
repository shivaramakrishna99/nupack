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
from natsort import as_ascii

from nupack import *  # Import NUPACK

class Material(Enum):
    dna = "DNA"
    rna = "RNA"
    rna95 = "RNA95"


# Define Model() object and parameters   

@small_task
def loop_stack(
    loop: str = "AAAA",
    structure: str = "....",
    material: str = Material.rna,
    temperature: float = 37,
    sodium: float = 1.0,
    magnesium: float = 0.0,
    outputFile: Optional[str] = "output"
) -> LatchFile:

    nt_model = Model(material=material, celsius=temperature, sodium=sodium, magnesium=magnesium)

    loopE = nt_model.loop_energy(loop=loop, structure=structure)
    loopE = f"{loopE}"

    stackE = nt_model.stack_energies(loop=loop, structure=structure)
    stackE = f"{stackE}"

    content = f"""
        ----------OUTPUT----------
    
        Loop Free Energy: 
        {loopE} kcal/mol
    
        Stacking State Free Energy:
        {stackE} 

        ----------INPUT SPECIFICATIONS----------

        Nucleotide Sequence: {loop}
        Dot Bracket Structure: {structure}

        Energy parameter: {material}
        Temperature: {temperature} °C
        Na+: {sodium} M
        Mg++: {magnesium} M

        ----------END----------
    """

    with open(f"/root/{outputFile}", "w") as f:
        f.write(content)

    return LatchFile(f"/root/{outputFile}", f"latch:///{outputFile}.txt")

@workflow
def nupack_loop_stack(
    loop: str = "AU+AU+AU",
    structure: str = "((+)(+))",
    material: Material = Material.rna,
    temperature: float = 37.0,
    sodium: float = 1.0,
    magnesium: float = 0.0,
    outputFile: Optional[str] = "output"
) -> LatchFile:
    """Analyse loop free energy and stacking state free energies for single and multiloop structures using NUPACK

    # NUPACK - Loop Free Energy and Stacking State Energies
    ---

    ## **About**
    ---
    [NUPACK](https://docs.nupack.org/#about) is a growing software suite for the analysis and design of nucleic acid structures, devices, and systems serving the needs of researchers in the fields of nucleic acid nanotechnology, molecular programming, synthetic biology, and across the life sciences more broadly.

    ## **Citations**
    ---
    ### NUPACK Analysis Algorithms
    **Complex analysis and test tube analysis**
	- M.E. Fornace, N.J. Porubsky, and N.A. Pierce (2020). A unified dynamic programming framework for the analysis of interacting nucleic acid strands: enhanced models, scalability, and speed.  [ACS Synth Biol](https://pubs.acs.org/doi/abs/10.1021/acssynbio.9b00523) , 9:2665-2678, 2020. ( [pdf](http://www.nupack.org/downloads/serve_public_file/fornace20.pdf?type=pdf) ,  [supp info](http://www.nupack.org/downloads/serve_public_file/fornace20_supp.pdf?type=pdf) )
	- R. M. Dirks, J. S. Bois, J. M. Schaeffer, E. Winfree, and N. A. Pierce. Thermodynamic analysis of interacting nucleic acid strands.  [SIAM Rev](http://epubs.siam.org/doi/abs/10.1137/060651100) , 49:65-88, 2007. ( [pdf](http://www.nupack.org/downloads/serve_public_file/sirev07.pdf?type=pdf) )
    **Pseudoknot analysis**
	- R. M. Dirks and N. A. Pierce. An algorithm for computing nucleic acid base-pairing probabilities including pseudoknots.  [J Comput Chem](http://onlinelibrary.wiley.com/doi/10.1002/jcc.10296/abstract) , 25:1295-1304, 2004. ( [pdf](http://www.nupack.org/downloads/serve_public_file/jcc04.pdf?type=pdf) )
	- R. M. Dirks and N. A. Pierce. A partition function algorithm for nucleic acid secondary structure including pseudoknots.  [J Comput Chem](http://onlinelibrary.wiley.com/doi/10.1002/jcc.20057/abstract) , 24:1664-1677, 2003. ( [pdf](http://www.nupack.org/downloads/serve_public_file/jcc03.pdf?type=pdf) ,  [supp info](http://www.nupack.org/downloads/serve_public_file/jcc03_supp.pdf?type=pdf) )

    **Workflow Repository** - https://github.com/shivaramakrishna99/nupack/
    **Acknowledgements** - https://docs.nupack.org/#acknowledgments
    **License** - https://docs.nupack.org/#license
    ---

    __metadata__:
        display_name: NUPACK - Loop Free Energy and Stacking State Energies
        
        author:
            name: The NUPACK Team
            email: support@nupack.org
            github: https://github.com/beliveau-lab/NUPACK
        
        repository: https://github.com/beliveau-lab/NUPACK
        
        license:
            id: BSD-3-Clause

    Args:
    
        material:
            __metadata__:
                display_name: "Nucleic Acid Type"
                _tmp:
                    section_title: Model Specification
                appearance:
                    comment: "Choose between DNA and RNA free energy parameter sets. Default is 'rna', based on Matthews et al., 1999"

        temperature:
            __metadata__:
                display_name: "Temperature (in degree Celsius)"
                appearance:
                    comment: "Temperature of system. Default is 37 °C"

        loop:
            __metadata__:
                display_name: "Loop Sequence(s) as nucleotides"
                _tmp:
                    section_title: Loop Details
                appearance:
                    comment: "Enter the nucleotide sequence of a loop. Separate sequences using the + symbol"

        structure:
            __metadata__:
                display_name: "Loop structure (in dot-bracket notation)"
                appearance:
                    comment: "Enter the dot bracket notation of a loop. Separate sequences using the + symbol"

        sodium:
            __metadata__:
                display_name: "Na+ concentration (in M)"
                _tmp:
                    hidden: true
                appearance:
                    comment: "The sum of the concentrations of (monovalent) sodium, potassium, and ammonium ions, is specified in units of molar. Default: 1.0, Range: [0.05,1.1]"

        magnesium:
            __metadata__:
                display_name: "Mg++ (in nM). Default is 0 nM"
                _tmp:
                    hidden: true
                appearance:
                    comment: "The concentration of (divalent) magnesium ions, is specified in units of molar. Default: 0.0, Range: [0.0,0.2]"

        outputFile:
            __metadata__:
                display_name: "Output Name"
                _tmp:
                    section_title: Output Specification
                appearance:
                    comment: "Specify the name of your output file. Default name given is 'output.txt'"
    """
    return loop_stack(
    loop=loop,
    structure=structure,
    material=material,
    temperature=temperature,
    sodium=sodium,
    magnesium=magnesium,
    outputFile=outputFile
    )
