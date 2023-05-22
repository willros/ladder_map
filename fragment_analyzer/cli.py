import fire
import pandas as pd
from pathlib import Path
import sys
import logging

import fragment_analyzer

def setup_logging(outdir: str) -> None:
    """
    Set up the logging object and saves the log file to the same dir as the results files.
    """
    if not (outdir := Path(outdir)).exists():
        outdir.mkdir(parents=True)
        
    LOG_FILE = f"{outdir}/fraggler.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] \n%(message)s",
        datefmt='%Y-%m-%d %I:%M:%S',
        handlers=[
            logging.FileHandler(LOG_FILE),
            logging.StreamHandler()
        ],
    )
    
def get_files(in_path: str) -> list[Path]:
    # If in_path is a directory, get a list of all .fsa files in it
    if Path(in_path).is_dir():
        files = [x for x in Path(in_path).iterdir() if x.suffix == ".fsa"]
    else:
        files = [Path(in_path)]
    return files
    

def report(
    in_path: str,
    out_folder: str,
    ladder: str = "LIZ",
    peak_model: str = "gauss",
    min_ratio: float = 0.2,
    min_height: int = 100,
    cutoff: int = 175,
    trace_channel: str = "DATA9",
    peak_height: int = 200,
    custom_peaks: str = None,
) -> None:
    """
    Generate a peak area report for all input files.
    """

    # Logging 
    setup_logging(out_folder)
    INFO = f"""
    Runned command:
    {' '.join(sys.argv)}

    Generating report with the following parameters:
        In path: {in_path}
        Out folder: {out_folder}
        Ladder: {ladder}
        Peak model: {peak_model}
        Min ratio: {min_ratio}
        Min height: {min_height}
        Cutoff: {cutoff}
        Trace channel: {trace_channel}
        Peak Height: {peak_height}
        Custom Peaks: {custom_peaks}
    """
    logging.info(INFO)
    
    # Files
    files = get_files(in_path)
    
    # Generate a peak area report for each file
    failed_files = []
    for file in files:
        logging.info(f"Processing file: {file}")
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
                custom_peaks=custom_peaks,
            )
        except Exception as e:
            logging.error(f"""
            Not able to process file. Reason:
            {e}
            """)
            failed_files.append(file.stem)
    
    failed_files = "\n".join(failed_files)
    logging.info(f"""Following files were not processed:
    {failed_files}
    """)


def peak_table(
    in_path: str,
    out_name: str,
    ladder: str = "LIZ",
    peak_model: str = "gauss",
    min_height: int = 100,
    cutoff: int = 175,
    min_ratio: float = 0.2,
    trace_channel: str = "DATA9",
    peak_height: int = 200,
    custom_peaks: str = None,
    excel: bool = True,
) -> pd.DataFrame:
    """
    Generate a combined dataframe of peaks for all input files.
    """
    # Logging 
    setup_logging(out_folder)
    INFO = f"""
    Runned command:
    {' '.join(sys.argv)}

    Generating peak table with the following parameters:
        In path: {in_path}
        Out folder: {out_folder}
        Ladder: {ladder}
        Peak model: {peak_model}
        Min ratio: {min_ratio}
        Min height: {min_height}
        Cutoff: {cutoff}
        Trace channel: {trace_channel}
        Peak Height: {peak_height}
        Custom Peaks: {custom_peaks}
        Excel: {excel}
    """
    logging.info(INFO)
    
    # Files
    files = get_files(in_path)
    
    peak_dfs = []
    failed_files = []
    for file in files:
        logging.info(f"Processing file: {file}")
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
                model,
                cutoff=cutoff,
                min_ratio=min_ratio,
                peak_height=peak_height,
                custom_peaks=custom_peaks,
            )
            peak_dfs.append(pam.assays_dataframe(peak_model))
        except Exception as e:
            logging.error(f"Following file did not work: {file}")
            logging.error(f"Reason: {e}")
            failed_files.append(file.stem)
            
    failed_files = "\n".join(failed_files)
    logging.info(f"""Following files were not processed:
    {failed_files}
    """)

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
