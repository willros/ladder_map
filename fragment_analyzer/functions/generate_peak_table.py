import pandas as pd
from pathlib import Path
import fragment_analyzer


def generate_peak_table(
    in_files: str | list,
    ladder: str, 
    peak_model: str,
    min_height: int = 100,
    cutoff: int = 175,
    min_ratio: float = 0.3,
    trace_channel: str = "DATA9",
) -> pd.DataFrame:
    """
    Generate a combined dataframe of peaks for all input files.

    Args:
        in_files (Union[str, List[str]]): Path(s) to the input .fsa files, or a directory containing them.
        ladder (str): The ladder used for the fragment analysis.
        peak_model (str): The peak model used for peak area calculations.
        min_height (int, optional): The minimum peak height required for detection. Defaults to 100.
        cutoff (int, optional): The cutoff value for peak area calculations. Defaults to 175.
        min_ratio (float, optional): The minimum ratio required for peak detection. Defaults to 0.3.
        trace_channel (str, optional): The channel used for trace data. Defaults to "DATA9".

    Returns:
        pandas.DataFrame: A combined dataframe containing the peak positions and their corresponding areas.

    Reads all the input .fsa files and generates a combined dataframe of all the peaks detected in each file using the 
    specified ladder and peak model. If a directory is provided instead of specific file paths, all .fsa files within 
    the directory will be used. Peak detection is based on the specified peak model and other optional parameters, 
    such as minimum peak height and cutoff value. 

    Example:
    --------
    peak_df = generate_peak_table(
        in_files="my_folder", ladder="LIZ", peak_model="gauss"
    )
    """
    
    if isinstance(in_files, str):
        in_files = [x for x in Path(in_files).iterdir() if x.suffix == ".fsa"]
        
    peak_dfs = []
    for x in in_files:
        try:
            fsa = fragment_analyzer.FsaFile(
                x,
                ladder,
                min_height=min_height,
                trace_channel=trace_channel
            )
            pla = fragment_analyzer.PeakLadderAssigner(fsa)
            model = fragment_analyzer.FitLadderModel(pla)
            pam = fragment_analyzer.PeakAreaDeMultiplex(
                model,
                cutoff=cutoff, 
                min_ratio=min_ratio
            )
            peak_dfs.append(pam.assays_dataframe(peak_model))
        except:
            print(f"FAILED: {fsa.file_name}")

    return pd.concat(peak_dfs).reset_index(drop=True)
