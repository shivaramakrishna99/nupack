"""
Assemble and sort some COVID reads...
"""

import subprocess
from pathlib import Path

from latch import small_task, workflow
from latch.types import LatchFile
from flytekit.core.with_metadata import FlyteMetadata

@small_task
def RNAfold_task(
    input_seq: str,
    output_name: str
) -> LatchFile:

    fold = RNA.fold(input_seq)

    out = Path(f"/root/{output_name}.txt")

    return LatchFile(str(out), f"latch://{out}")

@workflow
def RNAfold(
    input_seq: str,
    output_name: str,
) -> LatchFile:
    """Description...

    # RNAfold
    
    ---
    ## About

    <What does RNAfold do?> 

    __metadata__:
        display_name: Predict the folded structure of an RNA sequence
        author: 
            name:
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
    return RNAfold_task(
        input_seq=input_seq,
        output_name=output_name
    )
