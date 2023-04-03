import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_widths
from lmfit.models import VoigtModel, GaussianModel, LorentzianModel
from fragment_analyzer.ladder_fitting.fit_ladder_model import FitLadderModel


class PeakAreaDeMultiplexIterator:
    def __init__(self, number_of_assays):
        self.number_of_assays = number_of_assays
        self.current = 0
        
    def __next__(self):
        if self.current >= self.number_of_assays:
            raise StopIteration
        else:
            result = self.current
            self.current += 1
            return result


class PeakAreaDeMultiplex:
    """
    Class for finding peak areas and quotients of peaks in a given data set.

    Parameters:
    -----------
    model: FitLadderModel
        A FitLadderModel object containing the data set to be analyzed.
    peak_finding_model: str
        The name of the peak-finding model to be used.
    min_ratio: float, optional (default=0.2)
        The minimum ratio of peak height to highest peak height required to consider a peak as valid.
    search_peaks_start: int, optional (default=50)
        The starting point in basepairs for the search for peaks.

    Attributes:
    -----------
    model: FitLadderModel
        A FitLadderModel object containing the data set to be analyzed.
    raw_data: pd.DataFrame
        The raw data from the FitLadderModel object.
    file_name: str
        The name of the file associated with the FitLadderModel object.
    search_peaks_start: int
        The starting point in basepairs for the search for peaks.
    found_peaks: bool
        A flag indicating whether any peaks were found.
    peaks_index: np.ndarray
        An array of the indices of the peaks found.
    peaks_dataframe: pd.DataFrame
        A DataFrame of the peaks found, with basepairs and peak heights.
    peak_information: pd.DataFrame
        A DataFrame of the peaks found, with basepairs, peak heights, ratios, and peak names.
    peak_widths: pd.DataFrame
        A DataFrame of the peaks found, with basepairs, peak heights, start and end indices, and peak names.
    divided_peaks: List[pd.DataFrame]
        A list of DataFrames, each containing a single peak and its associated data.
    fit_df: List[pd.DataFrame]
        A list of DataFrames, each containing the raw data and the best-fit curve for a single peak.
    fit_params: List[dict]
        A list of dictionaries, each containing the parameters of the best-fit curve for a single peak.
    fit_report: List[str]
        A list of strings, each containing the report of the best-fit curve for a single peak.
    quotient: float
        The quotient of the areas of the peaks, calculated as the last peak divided by the mean of the peaks to the left of it.
    """

    def __init__(
        self,
        model: FitLadderModel,
        min_ratio: float = 0.15,
        # TODO
        # Change search_peaks_start to something else
        search_peaks_start: int = 110,
        # TODO
        # Change the peak_height number to something else? 
        peak_height: int = 500,
        distance_between_assays: int = 15,
        cutoff: float = None
    ) -> None:
        self.model = model
        self.raw_data = self.model.adjusted_baisepair_df
        self.file_name = self.model.fsa_file.file_name
        self.search_peaks_start = search_peaks_start
        self.cutoff = cutoff or None
        
        # find peaks
        self.find_peaks_agnostic(
            peak_height=peak_height,
            min_ratio=min_ratio,
            distance_between_assays=distance_between_assays,
        )

        # if no peaks could be found
        self.found_peaks = True
        if self.peak_information.shape[0] == 0:
            self.found_peaks = False
            print(
                f"No peaks could be found in {self.file_name}. Please look at the raw data."
            )

        # if peaks are found
        if self.found_peaks:
            # find peak widths
            self.find_peak_widths()
            # divade peaks based on their assay they belonging
            self.divided_assays = self.divide_peak_assays()
            # how many assays does this sample contain?
            self.number_of_assays = len(self.divided_assays)
            # divide all peaks in each assay into separate dataframes
            self.divided_peaks = [self.divide_peaks(x) for x in self.divided_assays]
            # Print information to the user
            print(f"{self.peak_information.shape[0]} peaks found in {self.file_name}")
            print(f"{self.number_of_assays} assays in {self.file_name}")
            
    def __iter__(self):
        return PeakAreaDeMultiplexIterator(self.number_of_assays)
    

    def find_peaks_agnostic(
        self,
        peak_height: int,
        min_ratio: float,
        distance_between_assays: int,
    ) -> None:
        peaks_dataframe = self.raw_data.loc[
            lambda x: x.basepairs > self.search_peaks_start
        ]
        peaks_index, _ = find_peaks(peaks_dataframe.peaks, height=peak_height)

        peak_information = (
            peaks_dataframe
            .iloc[peaks_index]
            .assign(peaks_index=peaks_index)
            .assign(peak_name=lambda x: range(1, x.shape[0] + 1))
            # separate the peaks into different assay groups depending on the distance
            # between the peaks
            .assign(difference=lambda x: x.basepairs.diff())
            .fillna(100)
            .assign(
                assay=lambda x: np.select(
                    [x.difference > distance_between_assays],
                    [x.peak_name * 10],
                    default=pd.NA,
                )
            )
            .fillna(method="ffill")
            .assign(
                max_peak=lambda x: x.groupby("assay")["peaks"]
                .transform(np.max)
            )
            .assign(ratio=lambda x: x.peaks / x.max_peak)
            .loc[lambda x: x.ratio > min_ratio]
            .assign(peak_name=lambda x: range(1, x.shape[0] + 1))
        )
        # TODO
        # Remove old
        #peak_information = (
        #    peaks_dataframe.iloc[peaks_index]
        #    .assign(peaks_index=peaks_index)
        #    .assign(ratio=lambda x: x.peaks / x.peaks.max())
        #    .loc[lambda x: x.ratio > min_ratio]
        #    .assign(peak_name=lambda x: range(1, x.shape[0] + 1))
        #    # separate the peaks into different assay groups depending on the distance
        #    # between the peaks
        #    .assign(difference=lambda x: x.basepairs.diff())
        #    .fillna(100)
        #    .assign(
        #        assay=lambda x: np.select(
        #            [x.difference > distance_between_assays],
        #           [x.peak_name * 10],
        #            default=pd.NA,
        #        )
        #    )
        #    .fillna(method="ffill")
        #)

        # update peaks_index based on the above filtering
        peaks_index = peak_information.peaks_index.to_numpy()

        # update class attributes
        self.peaks_index = peaks_index
        self.peaks_dataframe = peaks_dataframe
        self.peak_information = peak_information

    def find_peak_widths(self, rel_height: float = 0.95):
        widths = peak_widths(
            self.peaks_dataframe.peaks,
            self.peaks_index,
            rel_height=rel_height,
        )

        df = pd.DataFrame(widths).T
        df.columns = ["x", "peak_height", "peak_start", "peak_end"]

        self.peak_widths = (
            df.assign(peak_start=lambda x: np.floor(x.peak_start).astype(int))
            .assign(peak_end=lambda x: np.ceil(x.peak_end).astype(int))
            .assign(peak_name=lambda x: range(1, x.shape[0] + 1))
            .merge(self.peak_information, on="peak_name")
        )

    def divide_peak_assays(self) -> list[pd.DataFrame]:
        """
        Divide the peaks belonging to different assays based on their assay number
        """
        df = self.peak_widths
        return [df.loc[df.assay == x] for x in df.assay.unique()]

    def divide_peaks(self, assay: pd.DataFrame, padding: int = 4) -> list[pd.DataFrame]:
        # add some padding to the left and right to be sure to include everything in the peak
        return [
            self.peaks_dataframe.iloc[x.peak_start - padding : x.peak_end + padding]
            for x in assay.itertuples()
        ]

    def fit_lmfit_model(self, peak_finding_model: str, assay_number: int):
        if assay_number >= self.number_of_assays:
            raise IndexError(
                f"""
                The sample only contains {self.number_of_assays} assays. 
                Use a number inside of the range.
                Indexing starts at 0.
                """
            )

        if peak_finding_model == "gauss":
            model = GaussianModel()
        elif peak_finding_model == "voigt":
            model = VoigtModel()
        elif peak_finding_model == "lorentzian":
            model = LorentzianModel()
        else:
            raise NotImplementedError(
                f"""
                {peak_finding_model} is not implemented! 
                Options: [gauss, voigt, lorentzian]
                """
            )

        fit_df = []
        fit_params = []
        fit_report = []
        for df in self.divided_peaks[assay_number]:
            df = df.copy()
            y = df.peaks.to_numpy()
            x = df.basepairs.to_numpy()

            params = model.guess(y, x)
            out = model.fit(y, params, x=x)

            fit_df.append(df.assign(fitted=out.best_fit, model=peak_finding_model))
            fit_params.append(out.values)
            fit_report.append(out.fit_report())

        # Update the instances of the model fit
        self.fit_df = fit_df
        self.fit_params = fit_params
        self.fit_report = fit_report

    def calculate_quotient(
        self, 
    ) -> None:
        
        """
        :Params:
        """
        areas = np.array([x["amplitude"] for x in self.fit_params])
        
        right_by_left = True
        if self.cutoff is not None:
            if pd.concat(self.fit_df).basepairs.mean() < self.cutoff:
                right_by_left = False

        # if there only is 1 peak, return 0
        if len(areas) == 1:
            self.quotient = 0
            return

        # if there only are 2 peaks, return the quotient
        if len(areas) == 2:
            # left peak divided by right peak
            if not right_by_left:
                self.quotient = areas[0] / areas[1]
                return
            # right peak divided by left peak
            self.quotient = areas[1] / areas[0]
            return

        # TODO change this to the proper assay
        # return the last peak divided by the mean of the peaks to the left of it
        self.quotient = areas[-1] / areas[:-1].mean()

    def peak_position_area_dataframe(self, assay_number: int) -> pd.DataFrame:
        """
        Returns a DataFrame of each peak and its properties
        """
        dataframes = []
        for i, _ in enumerate(self.fit_df):
            df = (
                self.fit_df[i]
                .loc[lambda x: x.peaks == x.peaks.max()]
                .assign(area=self.fit_params[i]["amplitude"])
                .assign(peak_name=f"Peak {i + 1}")
                .drop(columns="time")
                .reset_index(drop=True)
                .rename(
                    columns={
                        "peaks": "peak_height",
                        "fitted": "fitted_peak_height",
                    }
                )
                .drop_duplicates("peak_name")
                .assign(file_name=self.file_name)
            )
            dataframes.append(df)

        self.assay_peak_area_df = pd.concat(dataframes).assign(
            quotient=self.quotient,
            peak_number=lambda x: x.shape[0],
            assay_number=assay_number + 1,
        )

    def fit_assay_peaks(
        self, 
        peak_finding_model: str,
        assay_number: int,
    ) -> None:
        """
        Runs fit_lmfit_model, calculate_quotient and peak_position_area_dataframe
        """
        self.fit_lmfit_model(peak_finding_model, assay_number)
        self.calculate_quotient()
        self.peak_position_area_dataframe(assay_number)
