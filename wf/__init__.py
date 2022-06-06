"""
Predict the folded structure of an RNA sequence
"""

import subprocess
from pathlib import Path

from latch import small_task, workflow
from latch.types import LatchFile
from flytekit.core.with_metadata import FlyteMetadata

from nupack import *  # Import NUPACK

@small_task
def nupack_task(
    input_seq: str,
    output_name: str
) -> LatchFile:

    fold = RNA.fold(input_seq)

    out = Path(f"/root/{output_name}.txt")

    return LatchFile(str(out), f"latch://{out}")

@workflow
def nupack( 
    input_seq: str,
    output_name: str,
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

        input_seq:
          RNA sequence for which structure and MFE is to be calculated

          __metadata__:
            display_name: Input RNA sequence
    """
    return nupack_task(
        input_seq=input_seq,
        output_name=output_name
    )
