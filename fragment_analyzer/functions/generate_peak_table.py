import pandas as pd
from pathlib import Path
import fragment_analyzer

def generate_peak_table(folder: str, ladder: str, peak_model: str) -> pd.DataFrame:
    """
    Generates a combined dataframe of all peaks for files in the given folder.

    :param folder: A string representing the path of the folder containing the .fsa files.
    :param ladder: A string representing the name of the ladder used for the fragment analysis.
    :param peak_model: A string representing the peak model used for peak area calculations.
    :return: A Pandas dataframe containing the peak positions and their corresponding areas.

    The function reads all the .fsa files in the folder, and generates a combined dataframe
    of all the peaks detected in each file using the specified ladder and peak model.
    
    ############ Example usage ############
    
    peak_df = generate_peak_table(
        folder="my_folder", ladder="LIZ", peak_model="gauss"
    )
    ############ End Example ##############
    """
    files = [x for x in Path(folder).iterdir() if x.suffix == ".fsa"]
    peak_dfs = []
    for x in files:
        try:
            fsa = fragment_analyzer.FsaFile(x, ladder)
            pla = fragment_analyzer.PeakLadderAssigner(fsa)
            model = fragment_analyzer.FitLadderModel(pla)
            pa = fragment_analyzer.PeakArea(model, peak_model)
            peak_dfs.append(pa.peak_position_area_dataframe)
        except:
            print(f"FAILED: {fsa.file_name}")

    return pd.concat(peak_dfs)
