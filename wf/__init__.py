"""
Specify a NUPACK model for a DNA/RNA loop and calculate two states of energy  
"""

from enum import Enum

from latch import small_task, workflow
from latch.types import LatchFile, LatchMetadata, LatchAuthor, LatchParameter, LatchAppearanceType, LatchRule
from natsort import as_ascii

from nupack import *  # Import NUPACK

class Material(Enum):
    dna = "DNA"
    rna = "RNA"
    rna95 = "RNA95"


# Define Model() object and parameters   

@small_task
def loopStackAnalysis(
    loop: str = "AAAA",
    structure: str = "....",
    material: str = Material.rna,
    temperature: float = 37,
    sodium: float = 1.0,
    magnesium: float = 0.0,
    out: str = "output"
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

metadata = LatchMetadata(
    display_name="NUPACK - Loop and Stack Energies",
    documentation="https://docs.nupack.org/model/#compute-loop-free-energy",
    author=LatchAuthor(
        name="Shivaramakrishna Srinivasan",
        email="shivaramakrishna.srinivasan@gmail.com",
        github="https://github.com/shivaramakrishna99",
    ),
    repository="https://github.com/shivaramakrishna99/nupack-loop-stack",
    license="BSD-3-Clause",
)

metadata.parameters["loop"] = LatchParameter(
    display_name="Loop Sequence(s)",
    description="Enter the nucleotide sequence of a loop. Separate two or more sequences using the '+' symbol",
    section_title= 'Loop Details',
)
metadata.parameters["structure"] = LatchParameter(
    display_name="Loop Structure(s)",
    description="Enter the dot bracket notation of a loop. Separate sequences using the + symbol",
)

metadata.parameters["material"] = LatchParameter(
    display_name="Nucleic Acid Type",
    description="Choose between DNA and RNA free energy parameter sets. Default is 'rna', based on Matthews et al., 1999",
    section_title="Model Specification",
    hidden=True,
)

metadata.parameters["temperature"] = LatchParameter(
    display_name="Temperature (in °C)",
    description="Temperature of system. Default: 37.0",
    hidden=True,
)
metadata.parameters["sodium"] = LatchParameter(
    display_name="Na+ (in M)",
    description="The total concentration of (monovalent) sodium, potassium, and ammonium ions, specified as molarity. Default: 1.0, Range: [0.05,1.1]",
    hidden=True,
    section_title="Additional Model Specification"
)
metadata.parameters["magnesium"] = LatchParameter(
    display_name="Mg++ (in nM)",
    description="The total concentration of (divalent) magnesium ions, specified as molarity. Default: 0.0, Range: [0.0,0.2]",
    hidden=True
)
metadata.parameters["out"] = LatchParameter(
    display_name="Output File Name",
    section_title="Output"
)

@workflow(metadata)
def loopStackAnalysisNUPACK(
    loop: str = "AU+AU+AU",
    structure: str = "((+)(+))",
    material: Material = Material.rna,
    temperature: float = 37.0,
    sodium: float = 1.0,
    magnesium: float = 0.0,
    out: str = "output"
) -> LatchFile:
    """Analyse loop free energy and stacking state free energies, for single and multiloop structures using NUPACK

    # NUPACK - Loop Free Energy and Stacking State Energies

    ## **About**
    ---

    [NUPACK](https://docs.nupack.org/#about) is a growing software suite for the analysis and design of nucleic acid structures, devices, and systems serving the needs of researchers in the fields of nucleic acid nanotechnology, molecular programming, synthetic biology, and across the life sciences more broadly.

    ## **How to use**
    ---

    1. Specify the model by choosing the nucleic acid type. This parameter is based on parameter sets obtained from different research papers. Check them out [here](https://docs.nupack.org/model/#material)
    
    2. Specify the details of the loop's sequence and structural details.

    3. Specify any other changes in the construction of the Model() object using the hidden parameters such as temperature and ion concentrations. 

    4. Run the workflow!

    ## **Citations**
    ---

    ### NUPACK Analysis Algorithms

    **Complex analysis and test tube analysis**
	
    - M.E. Fornace, N.J. Porubsky, and N.A. Pierce (2020). A unified dynamic programming framework for the analysis of interacting nucleic acid strands: enhanced models, scalability, and speed.  [ACS Synth Biol](https://pubs.acs.org/doi/abs/10.1021/acssynbio.9b00523) , 9:2665-2678, 2020. ( [pdf](http://www.nupack.org/downloads/serve_public_file/fornace20.pdf?type=pdf) ,  [supp info](http://www.nupack.org/downloads/serve_public_file/fornace20_supp.pdf?type=pdf) )
	
    - R. M. Dirks, J. S. Bois, J. M. Schaeffer, E. Winfree, and N. A. Pierce. Thermodynamic analysis of interacting nucleic acid strands.  [SIAM Rev](http://epubs.siam.org/doi/abs/10.1137/060651100) , 49:65-88, 2007. ( [pdf](http://www.nupack.org/downloads/serve_public_file/sirev07.pdf?type=pdf) )
    
    **Pseudoknot analysis**
	
    - R. M. Dirks and N. A. Pierce. An algorithm for computing nucleic acid base-pairing probabilities including pseudoknots.  [J Comput Chem](http://onlinelibrary.wiley.com/doi/10.1002/jcc.10296/abstract) , 25:1295-1304, 2004. ( [pdf](http://www.nupack.org/downloads/serve_public_file/jcc04.pdf?type=pdf) )
	
    - R. M. Dirks and N. A. Pierce. A partition function algorithm for nucleic acid secondary structure including pseudoknots.  [J Comput Chem](http://onlinelibrary.wiley.com/doi/10.1002/jcc.20057/abstract) , 24:1664-1677, 2003. ( [pdf](http://www.nupack.org/downloads/serve_public_file/jcc03.pdf?type=pdf) ,  [supp info](http://www.nupack.org/downloads/serve_public_file/jcc03_supp.pdf?type=pdf) )

    **Workflow Repository** - (https://github.com/shivaramakrishna99/nupack-loop-stack/)
    
    **Acknowledgements** - (https://docs.nupack.org/#acknowledgments)

    *Authored by Shivaramakrishna Srinivasan. Feel free to reach out to me at shivaramakrishna.srinivasan@gmail.com*
    ---
    """
    return loopStackAnalysis(
    loop=loop,
    structure=structure,
    material=material,
    temperature=temperature,
    sodium=sodium,
    magnesium=magnesium,
    out=out
    )
