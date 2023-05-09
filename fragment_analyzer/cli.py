import fire
import pandas as pd
from pathlib import Path
import traceback
import sys

import fragment_analyzer


def report(
    in_path: str,
    out_folder: str,
    ladder: str = "LIZ",
    peak_model: str = "gauss",
    min_ratio: float = 0.3,
    min_height: int = 100,
    cutoff: int = 175,
    trace_channel: str = "DATA9",
    peak_height: int = 200,
) -> None:
    """
    Generate a peak area report for all input files.
    """
    # If in_path is a directory, get a list of all .fsa files in it
    if Path(in_path).is_dir():
        files = [x for x in Path(in_path).iterdir() if x.suffix == ".fsa"]
    else:
        files = [Path(in_path)]

    # Log parameters
    print(f"Runned command:")
    print(f"{" ".join(sys.argv)}")
    print("")
    print("Generating report with the following parameters:")
    print(f"    In path: {in_path}")
    print(f"    Out folder: {out_folder}")
    print(f"    Ladder: {ladder}")
    print(f"    Peak model: {peak_model}")
    print(f"    Min ratio: {min_ratio}")
    print(f"    Min height: {min_height}")
    print(f"    Cutoff: {cutoff}")
    print(f"    Trace channel: {trace_channel}")
    print(f"    Peak Height: {peak_height}")

    # Generate a peak area report for each file
    for file in files:
        try:
            fragment_analyzer.peak_area_report(
                fsa_file=file,
                ladder=ladder,
                folder=out_folder,
                peak_model=peak_model,
                min_ratio=min_ratio,
                min_height=min_height,
                cutoff=cutoff,
                trace_channel=trace_channel,
                peak_height=peak_height,
            )
        except Exception as e:
            print(f"ERROR: {file}")
            print(traceback.format_exc())
            print("")


def peak_table(
    in_path: str,
    out_name: str,
    ladder: str = "LIZ",
    peak_model: str = "gauss",
    min_height: int = 100,
    cutoff: int = 175,
    min_ratio: float = 0.3,
    trace_channel: str = "DATA9",
    peak_height: int = 200,
    excel: bool = False,
) -> pd.DataFrame:
    """
    Generate a combined dataframe of peaks for all input files.
    """
    
    # Logging
    print(f"Runned command:")
    print(f"{" ".join(sys.argv)}")
    print("")
    print("Generating peak table with the following parameters:")
    print(f"    In path: {in_path}")
    print(f"    Out name: {out_name}")
    print(f"    Ladder: {ladder}")
    print(f"    Peak model: {peak_model}")
    print(f"    Min ratio: {min_ratio}")
    print(f"    Min height: {min_height}")
    print(f"    Cutoff: {cutoff}")
    print(f"    Trace channel: {trace_channel}")
    print(f"    Peak Height: {peak_height}")
    
    # If in_path is a directory, get a list of all .fsa files in it
    if Path(in_path).is_dir():
        files = [x for x in Path(in_path).iterdir() if x.suffix == ".fsa"]
    else:
        files = [Path(in_path)]

    peak_dfs = []
    for file in files:
        try:
            fsa = fragment_analyzer.FsaFile(
                file,
                ladder,
                min_height=min_height,
                trace_channel=trace_channel,
            )
            pla = fragment_analyzer.PeakLadderAssigner(fsa)
            model = fragment_analyzer.FitLadderModel(pla)
            pam = fragment_analyzer.PeakAreaDeMultiplex(
                model, cutoff=cutoff, min_ratio=min_ratio, peak_height=peak_height,
            )
            peak_dfs.append(pam.assays_dataframe(peak_model))
        except Exception as e:
            print(f"ERROR: {file}")
            print(traceback.format_exc())
            print("")

    # Combine peak dataframes into a single dataframe
    df = pd.concat(peak_dfs).reset_index(drop=True)

    # Save combined dataframe as a CSV file
    if excel:
        df.to_excel(f"{out_name}.xlsx", index=False)
    else:
        df.to_csv(f"{out_name}.csv", index=False)


def run():
    """
    Run the command-line interface using Fire.
    """
    fire.Fire(
        {
            "report": report,
            "peak_table": peak_table,
        }
    )
