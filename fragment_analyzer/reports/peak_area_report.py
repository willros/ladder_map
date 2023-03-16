import panel as pn
from pathlib import Path
import datetime
import pandas as pd

import fragment_analyzer

pn.extension("tabulator")
pn.extension("vega", sizing_mode="stretch_width", template="fast")
pn.widgets.Tabulator.theme = "modern"


def header(
    text: str,
    bg_color: str = "#04c273",
    height: int = 300,
    fontsize: str = "px20",
    textalign: str = "center",
):
    """
    Template for markdown header like block
    """
    return pn.pane.Markdown(
        f"""
        {text}
        """,
        background=bg_color,
        height=height,
        margin=10,
        style={
            "color": "white",
            "padding": "10px",
            "text-align": f"{textalign}",
            "font-size": f"{fontsize}",
        },
    )


def generate_peak_area_no_peaks(name, date, plot_raw):
    """
    Docstring
    """

    head = header(
        text=f"""
        # Fragment Analysis Report
        ## Report of {name}
        ## Date: {date}
        """,
        fontsize="20px",
        bg_color="#03a1fc",
        height=250,
    )

    no_peaks_markdown = header(
        "# No peaks could be generated. Please look at the raw data.", height=100
    )
    raw_plot = pn.pane.Matplotlib(plot_raw.plot_raw_data)
    return pn.Column(
        head,
        no_peaks_markdown,
        raw_plot,
    )


## ------ NEW ------ ###
def peak_area_report(
    fsa_file: str, 
    ladder: str,
    folder: str, 
    peak_model: str,
    min_height: int = 100,
    min_ratio: float = 0.1,
    trace_channel: str = "DATA1"
) -> int:
    """
    Generates an HTML report for the fragment analysis of an FSA file, including peak area data and plots.

    Parameters:
    -----------
    fsa_file : str
        The path to the FSA file to be analyzed.
    ladder : str
        The name of the ladder used in the FSA file.
    folder : str
        The path to the output folder where the report will be saved.
    peak_model : str
        The peak finding model used to identify peaks.

    Returns:
    --------
    int
        An integer representing the status of the report generation:
        - 0 if successful
        - 1 if no peaks were found

    Raises:
    -------
    FileNotFoundError
        If the specified FSA file cannot be found.
    IOError
        If the report file cannot be saved.
    """
    # FSA File and FsaFile Object
    fsa = fragment_analyzer.FsaFile(
        fsa_file,
        ladder, 
        min_height=min_height,
        trace_channel=trace_channel,
    )
    file_name = fsa.file_name
    date = fsa.fsa["ABID1"].decode("utf-8")
    date = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M:%S.%f").strftime("%Y-%m-%d")

    # LadderAssigner, Model and PeakArea
    ladder_assigner = fragment_analyzer.PeakLadderAssigner(fsa)
    model = fragment_analyzer.FitLadderModel(ladder_assigner)
    raw_plots = fragment_analyzer.PlotRawData(model)
    ladder_plots = fragment_analyzer.PlotLadder(model)
    peak_areas = fragment_analyzer.PeakAreaDeMultiplex(
        model, 
        min_ratio=min_ratio,
    )
    peak_plots = fragment_analyzer.PlotPeakArea(peak_areas)

    # create the output folder if it doesn't exist
    outpath = Path(folder)
    if not outpath.exists():
        outpath.mkdir(parents=True)

    # If no peaks could be found
    if not peak_areas.found_peaks:
        outname = outpath / f"FAILED-fragment_analysis-report-{file_name}-{date}.html"
        generate_peak_area_no_peaks(
            file_name, 
            date, 
            raw_plots
        ).save(
            outname,
            title=file_name
        )

        return 1

    # If peaks
    outname = outpath / f"fragment_analysis-report-{file_name}-{date}.html"
    generate_peak_area_report(
        file_name, 
        date, 
        peak_model,
        raw_plots, 
        ladder_plots,
        peak_plots, 
        peak_areas,
    ).save(
        outname,
        title=file_name,
    )

    return 0


def generate_peak_area_report(
    name: str,
    date: str,
    peak_model: str,
    plot_raw, 
    plot_ladder,
    plot_peaks, 
    peak_area,
):
    ### ----- HEADER ----- ###
    head = header(
        text=f"""
        # Fragment Analysis Report
        ## Report of {name} 
        ## Date: {date}
        """,
        fontsize="20px",
        bg_color="#03a1fc",
        height=250,
    )
    ### ----- Raw Data ----- ### 
    
    # Header for this section
    raw_header = header(
        text=f"## Plot of Raw Data",
        bg_color="#04c273",
        height=80,
        textalign="left",
    )
    # PLOT
    raw_data_plot = pn.pane.Matplotlib(
        plot_raw.plot_raw_data, 
        name="Raw Data",
    )
    
    # Section
    raw_tab = pn.Tabs(
        raw_data_plot,
    )
    raw_section = pn.Column(raw_header, raw_tab)
    
    ### ----- Ladder Information ----- ### 
    
    # Header for this section
    ladder_header = header(
        text=f"## Information about the ladder",
        bg_color="#04c273",
        height=80,
        textalign="left",
    )
    # Ladder peak plot
    ladder_peak_plot = pn.pane.Matplotlib(
        plot_ladder.plot_ladder_peaks,
        name="Ladder Peak Plot",
    )
    # Ladder Correlation
    ladder_correlation_plot = pn.pane.Matplotlib(
        plot_ladder.plot_model_fit,
        name="Ladder Correlation Plot",
    )
    
    # Section
    ladder_tab = pn.Tabs(
        ladder_peak_plot, ladder_correlation_plot,
    )
    ladder_section = pn.Column(ladder_header, ladder_tab)
    
    ### ----- Peaks Information ----- ### 
    # Header for this section
    peaks_header = header(
        text=f"## Information about the peaks",
        bg_color="#04c273",
        height=80,
        textalign="left",
    )
    
    # All peaks plot
    all_peaks_plot = pn.pane.Matplotlib(
        plot_peaks.plot_peaks(),
        name="All Peaks",
    )
    
    # Individual peaks and fitting of the model for assays in sample
    # append all_peaks_plot first
    peaks_tab = pn.Tabs(all_peaks_plot)
    for assay in peak_area:
        plot = plot_peaks.plot_areas(peak_model, assay)
        name = f"Assay {assay + 1}"
        plot_pane = pn.pane.Matplotlib(plot, name=name)
        peaks_tab.append(plot_pane)
    
    # Section
    peaks_section = pn.Column(peaks_header, peaks_tab)
    
    
    ### ----- Peaks DataFrame ----- ### 
    # Header for this section
    dataframe_header = header(
        text=f"## Peaks Table",
        bg_color="#04c273",
        height=80,
        textalign="left",
    )
    # Create dataframe
    df = []
    for i in peak_area:
        peak_area.fit_assay_peaks(peak_model, i)
        df.append(peak_area.assay_peak_area_df)

    df = pd.concat(df).reset_index(drop=True)

    # DataFrame Tabulator
    peaks_df_tab = pn.widgets.Tabulator(
    df,
    #editors={"sequence": {"type": "editable", "value": False}},
    layout="fit_columns",
    pagination="local",
    page_size=15,
    show_index=False,
    name="Peaks Table",
)
    
    # Section
    dataframe_tab = pn.Tabs(
        peaks_df_tab
    )
    dataframe_section = pn.Column(dataframe_header, dataframe_tab)
    
    ### CREATE REPORT ###

    all_tabs = pn.Tabs(
        ("Raw Data", raw_section),
        ("Ladder Information", ladder_section),
        ("Peaks", peaks_section),
        ("Table", dataframe_section),
        tabs_location="left",
    )
    report = pn.Column(
        head,
        pn.layout.Divider(),
        all_tabs,
    )
    
    return report


