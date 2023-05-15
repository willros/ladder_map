## ladder map
![logo](examples/logo.png)

Matches ladders to peaks by correlation for fragment analysis. The strategy resembles the one used by [Fragman](https://cran.r-project.org/web/packages/Fragman/index.html) for R.

One difference is that combinations of peaks are generated using [NetworkX](https://networkx.org/) to eliminate impossible combinations. This reduces complexity substantially and allows for an exhaustive search to identify the best match.

## Install

```bash
pip install peak-a-py

conda install -c bioconda peak-a-py
```

## Usage

The library exists of two main classes, namely `LadderMap` and `PeakArea`. 
`LadderMap` matches ladders to peaks by correlation and stores peak and ladder information. The have methods to reassign timeseries data to basepair steps using linear regression. 

`PeakArea` calculates peak area.

### Example Usage

#### Area under curve using naive integrals
```python
from fragment_analyzer.ladder_map import LadderMap, PeakArea

data = "demo/4071_Dx 230113_PRT1_PRT3_rn/PRT3_NA18507_4071_D11_Dx.fsa"

laddermap = LadderMap(data)

peak_area = PeakArea(
    laddermap.adjusted_step_dataframe(),
    start=190, 
    end=200,
    rel_height=.97
)

peak_area.plot_peak_widths()
```

#### Output
![peak_area](examples/peak_area.png)

#### Visualization of best sample ladder peaks
```python
fig = laddermap.plot_best_sample_ladder()
```
#### Output
![sample_ladder](examples/best_sample_ladder.png)

#### Fitting model to the data
##### Voigt Distribution
```python
peak_area.plot_lmfit_model("voigt")
```
#### Output
![voigt_model](examples/voigt_model.png)

##### Gauss Distribution
```python
peak_area.plot_lmfit_model("gauss")
```
#### Output
![gauss_model](examples/gauss_model.png)


#### Looking at more than two peaks
```python
peak_area = PeakArea(
    laddermap.adjusted_step_dataframe(channel="DATA1"),
    start=250, 
    end=300,
    num_peaks=4,
    padding=2,
    model="voigt"
)


peak_area.plot_lmfit_model()
```
The last peak is divided by the mean of the peaks to the left of it:
#### Output
![four_peaks](examples/four_peaks.png)

#### If data needs baseline correction and normalization:
```python
laddermap = LadderMap(data, normalize_peaks=False)
```
Messy output
#### Output
![messy](examples/needs_normalization.png)

#### Normalized data:
```python
laddermap = LadderMap(data, normalize_peaks=True)
```
Normalized peaks:
#### Output
![normalized](examples/normalized.png)


## Usage: peak_area_report

The `peak_area_report` function generates an HTML report for the fragment analysis of an FSA file, including peak area data and plots.

### Parameters:

- `fsa_file` (str): The path to the FSA file to be analyzed.
- `ladder` (str): The name of the ladder used in the FSA file.
- `folder` (str): The path to the output folder where the report will be saved.
- `peak_model` (str): The peak finding model used to identify peaks.
- `min_interpeak_distance` (int, optional, default=30): Minimum distance between peaks.
- `min_height` (int, optional, default=100): Minimum peak height for inclusion in the analysis.
- `min_ratio` (float, optional, default=0.1): Minimum peak area ratio for multiplexing.
- `trace_channel` (str, optional, default="DATA1"): The trace channel in the FSA file.
- `search_peaks_start` (int, optional, default=100): The starting point for peak search.
- `cutoff` (float, optional): Cutoff value for peak area de-multiplexing.

### Returns:

- An integer representing the status of the report generation:
  - 0 if successful
  - 1 if no peaks were found

### Raises:

- FileNotFoundError: If the specified FSA file cannot be found.
- IOError: If the report file cannot be saved.

### Example usage:

```python
result = peak_area_report(
    fsa_file="path/to/fsa_file.fsa",
    ladder="LIZ",
    folder="output_folder",
    peak_model="gauss"
)
```
This example will generate an HTML report for the fragment analysis of the specified FSA file, using the LIZ ladder and a Gaussian peak model. The report will be saved in the output_folder directory.


## Usage: generate_peak_table

The `generate_peak_table` function generates a combined dataframe of all peaks for .fsa files in the given folder, using the specified ladder and peak model.

### Parameters:

- `folder` (str): A string representing the path of the folder containing the .fsa files.
- `ladder` (str): A string representing the name of the ladder used for the fragment analysis.
- `peak_model` (str): A string representing the peak model used for peak area calculations.
- `min_height` (int, optional, default=100): Minimum peak height for inclusion in the analysis.
- `cutoff` (int, optional, default=175): Cutoff value for peak area de-multiplexing.

### Returns:

- A Pandas dataframe containing the peak positions and their corresponding areas.

### Example usage:

```python
peak_df = generate_peak_table(
    folder="my_folder", ladder="LIZ", peak_model="gauss"
)
```
This example will generate a dataframe containing the peak positions and their corresponding areas for all .fsa files in the my_folder directory, using the LIZ ladder and a Gaussian peak model.

-----------------------------------------------------------------------------------------------------------


# Peak Area DeMultiplex

The **PeakAreaDeMultiplex** class is designed for finding peak areas and quotients of peaks in a given data set. It provides methods for peak finding, peak width calculation, peak division by assay, fitting peak models, and calculating peak quotients.

## Dependencies

The project requires the following dependencies:

- `re`
- `pandas`
- `numpy`
- `scipy.signal.find_peaks`
- `scipy.signal.peak_widths`
- `lmfit.models.VoigtModel`
- `lmfit.models.GaussianModel`
- `lmfit.models.LorentzianModel`
- `fragment_analyzer.ladder_fitting.fit_ladder_model.FitLadderModel`

## Custom Errors

The project includes two custom error classes:

- `OverlappingIntervalError`: Raised when there are overlapping intervals in the custom peaks table.
- `WrongColumnsError`: Raised when the input DataFrame does not have the required columns.

## Validation Functions

The project includes two validation functions:

- `is_overlapping(df: pd.DataFrame) -> bool`: Checks if there are overlapping intervals in the DataFrame.
- `has_columns(df) -> bool`: Checks if the DataFrame has the required columns.

## PeakAreaDeMultiplex Class

### Constructor

The `PeakAreaDeMultiplex` class constructor takes the following parameters:

- `model: FitLadderModel`: A `FitLadderModel` object containing the data set to be analyzed.
- `min_ratio: float` (optional, default=0.15): The minimum ratio of peak height to highest peak height required to consider a peak as valid.
- `search_peaks_start: int` (optional, default=110): The starting point in base pairs for the search for peaks.
- `peak_height: int` (optional, default=350): The minimum peak height required to consider a peak as valid.
- `distance_between_assays: int` (optional, default=15): The distance in base pairs between peaks to be considered as belonging to different assays.
- `cutoff: float` (optional, default=None): The cutoff value for dividing peaks into left and right groups based on their base pair positions.
- `custom_peaks: str | pd.DataFrame` (optional, default=None): Custom peaks data either as a path to a CSV file or a DataFrame.

### Attributes

The `PeakAreaDeMultiplex` class has the following attributes:

- `model: FitLadderModel`: A `FitLadderModel` object containing the data set to be analyzed.
- `raw_data: pd.DataFrame`: The raw data from the `FitLadderModel` object.
- `file_name: str`: The name of the file associated with the `FitLadderModel` object.
- `search_peaks_start: int`: The starting point in base pairs for the search for peaks.
- `found_peaks: bool`: A flag indicating whether any peaks were found.
- `peaks_index: np.ndarray`: An array of the indices of the peaks found.
- `peaks_dataframe: pd.DataFrame`: A DataFrame of the peaks found, with base pairs and peak heights.
- `peak_information: pd.DataFrame`: A DataFrame of the peaks found, with base pairs, peak heights, ratios, and peak names.
- `peak_widths: pd.DataFrame`: A DataFrame of the peaks found, with base pairs, peak heights, start and end indices, and peak names.
- `divided_peaks: List[pd.DataFrame]`: A list of DataFrames, each containing a single peak and its associated data.
- `fit_df: List[pd.DataFrame]`: A list of DataFrames, each containing the fitted model parameters for a peak.
- `peak_quotients: pd.DataFrame`: A DataFrame of the peak quotients calculated from the divided peaks.

### Methods

The `PeakAreaDeMultiplex` class provides the following methods:

- `find_peaks() -> None`: Finds peaks in the raw data.
- `calculate_peak_widths() -> None`: Calculates the widths of the peaks.
- `divide_peaks() -> None`: Divides the peaks into left and right groups based on their base pair positions.
- `fit_peaks() -> None`: Fits peak models to the divided peaks.
- `calculate_peak_quotients() -> None`: Calculates peak quotients based on the fitted peak models.

## Usage Example

```python
# Import required libraries
import pandas as pd
from fragment_analyzer.ladder_fitting.fit_ladder_model import FitLadderModel
from peak_area_demultiplex import PeakAreaDeMultiplex

# Create FitLadderModel object
data = pd.read_csv('data.csv')
model = FitLadderModel(data)

# Create PeakAreaDeMultiplex object
demux = PeakAreaDeMultiplex(model, min_ratio=0.1, search_peaks_start=100, peak_height=300)

# Find peaks
demux.find_peaks()

# Calculate peak widths
demux.calculate_peak_widths()

# Divide peaks
demux.divide_peaks()

# Fit peaks
demux.fit_peaks()

# Calculate peak quotients
demux.calculate_peak_quotients()
```



## CLI Tool Example Usage:
### Fraggler report

#### Usage
To generate a peak area report for all input files, use the fraggler report command followed by the required positional arguments and any optional flags.


```console
fraggler report IN_PATH OUT_FOLDER <flags>
```
#### Description
The fraggler report command generates a peak area report for all input files.


#### Positional Arguments
The following positional arguments are required:

- `IN_PATH`: Type `str`. Specifies the input path.
- `OUT_FOLDER`: Type `str`. Specifies the output folder.

#### Flags
The following flags can be used with the `fraggler report` command:

- `-l, --ladder=LADDER`: Type `str`. Specifies the ladder. Default value: 'LIZ'.
- `--peak_model=PEAK_MODEL`: Type `str`. Specifies the peak model. Default value: 'gauss'.
- `--min_ratio=MIN_RATIO`: Type `float`. Specifies the minimum ratio. Default value: 0.3.
- `--min_height=MIN_HEIGHT`: Type `int`. Specifies the minimum height. Default value: 100.
- `--cutoff=CUTOFF`: Type `int`. Specifies the cutoff value. Default value: 175.
- `-t, --trace_channel=TRACE_CHANNEL`: Type `str`. Specifies the trace channel. Default value: 'DATA9'.
- `--peak_height=PEAK_HEIGHT`: Type `int`. Specifies the peak height. Default value: 200.
- `--custom_peaks=CUSTOM_PEAKS`: Type `Optional[str]`. Specifies custom peaks. Default value: None.





### Fraggler peak_table

#### Usage

This project provides a command-line tool, `fraggler peak_table`, that generates a combined dataframe of peaks for all input files.

```console
fraggler peak_table IN_PATH OUT_NAME <flags>
```
#### Description
The fraggler peak_table command generates a combined peak_table for all input files.

- If not specified, fraggler finds peaks agnostic in the `fsa file`. To specifiy custom assays with certain peaks and intervals, the user can add a .csv file to the `--custom_peaks` argument. The csv file MUST have the following shape:
```
| name | start | stop | amount |
|------|-------|------|--------|
| prt1 | 140   | 150  |        |

```
If `amount` if left emtpy, `fraggler` will take all peaks inside the interval. If amount is not empty, fraggler will include the top `N` peaks in the interval, based on height.

#### Positional Arguments
The following positional arguments are required:

- `IN_PATH`: Type: `str`
- `OUT_NAME`: Type: `str`

#### Flags
The following flags can be used with the `fraggler peak_table` command:

- `-l, --ladder=LADDER`: Type: `str`, Default: 'LIZ'
- `--peak_model=PEAK_MODEL`: Type: `str`, Default: 'gauss'
- `--min_height=MIN_HEIGHT`: Type: `int`, Default: 100
- `--cutoff=CUTOFF`: Type: `int`, Default: 175
- `--min_ratio=MIN_RATIO`: Type: `float`, Default: 0.3
- `-t, --trace_channel=TRACE_CHANNEL`: Type: `str`, Default: 'DATA9'
- `--peak_height=PEAK_HEIGHT`: Type: `int`, Default: 200
- `--custom_peaks=CUSTOM_PEAKS`: Type: `Optional[str]`, Default: None
- `-e, --excel=EXCEL`: Type: `bool`, Default: False




# TODO

* Add log class instead of printing inside functions







