import pandas as pd
import numpy as np
import networkx as nx
from scipy import signal
from sklearn.metrics import r2_score
from sklearn.preprocessing import minmax_scale
from scipy.interpolate import UnivariateSpline

from fragment_analyzer.utils.fsa_file import FsaFile


class PeakLadderAssigner:
    """
    A class that assigns peak sizes to a ladder for use in DNA fragment size analysis.
    The peak sizes are generated based on a given size standard and the input parameters.

    Attributes:
    fsa_file (FsaFile): an object of the FsaFile class that contains parameters for analysis

    Methods:
    assign_ladder_peak_sizes(): Assigns peak sizes to a ladder based on the given parameters and returns the best combination.
    get_peaks(size_standard): Finds peaks in the size standard based on given parameters and returns the peaks.
    generate_graph(peaks): Generates a graph based on the peaks and returns the graph.
    generate_combinations(graph): Generates combinations of peaks based on the graph and returns the combinations.
    get_best_fit(combinations, method): Finds the best fit for the given combinations and returns the best combination.
    _polynomial_model_inv_r2_score(ladder, comb): Returns a score based on the polynomial model inverse r2 score.
    _max_fractional_deviation_score(ladder, comb): Returns a score based on the maximum fractional deviation.
    _max_first_derivative_score(combination): Returns a score based on the maximum first derivative.
    _max_second_derivative_score(combination): Returns a score based on the maximum second derivative.
    _max_spline_second_derivative_score(combination): Returns a score based on the maximum spline second derivative.
    """

    def __init__(self, fsa_file: FsaFile) -> None:
        """
        Initializes the class with an object of the FsaFile class.

        Args:
        fsa_file (FsaFile): an object of the FsaFile class that contains parameters for analysis
        """

        self.fsa_file = fsa_file

    def assign_ladder_peak_sizes(self):
        """
        Assigns peak sizes to a ladder based on the given parameters and returns the best combination.

        Returns:
        np.array: the best combination of peak sizes
        """
        peaks = self.get_peaks(self.fsa_file.size_standard)
        graph = self.generate_graph(peaks)
        combinations = self.generate_combinations(graph)
        best_combination = self.get_best_fit(combinations)

        return best_combination

    def get_peaks(self, size_standard) -> np.array:
        """
        Finds peaks in the size standard based on given parameters and returns the peaks.

        Args:
        size_standard (np.array): the size standard array

        Returns:
        np.array: the peaks array
        """

        peaks_obj = signal.find_peaks(
            size_standard,
            distance=self.fsa_file.min_interpeak_distance,
            height=self.fsa_file.min_height,
        )

        peaks = peaks_obj[0]
        heights = peaks_obj[1]["peak_heights"]

        df = pd.DataFrame({"peaks": peaks, "heights": heights})

        idxmax = df["heights"].idxmax()
        df = df.drop(idxmax)

        peaks_adj = df.nlargest(self.fsa_file.max_peak_count, ["heights"])

        return peaks_adj["peaks"].sort_values().to_numpy()

    def generate_graph(self, peaks) -> nx.DiGraph:
        """
        Generates a graph based on the peaks and returns the graph.

        Args:
        peaks (np.array): the peaks array

        Returns:
        nx.DiGraph: the graph object
        """
        G = nx.DiGraph()

        for p in peaks:
            G.add_node(p)

        i = 0
        while i < peaks.size:
            j = i + 1
            while j < peaks.size:
                diff = peaks[j] - peaks[i]
                if diff <= self.fsa_file.max_ladder_trace_distance:
                    G.add_edge(peaks[i], peaks[j], length=diff)
                j += 1
            i += 1

        return G

    def generate_combinations(self, graph):
        """
        Generates combinations of nodes from a given directed acyclic graph (DAG)
        that satisfy certain conditions.

        Parameters:
        graph (networkx.DiGraph): a DAG representing a state machine

        Returns:
        A numpy array representing a sequence of nodes that satisfies the
        conditions of the DAG.
        """

        # Get start nodes that have zero in-degree
        start_nodes = [node for node in graph.nodes if graph.in_degree(node) == 0]

        # Get end nodes that have zero out-degree
        end_nodes = [node for node in graph.nodes if graph.out_degree(node) == 0]

        # Get all simple paths from start node to end node
        all_paths = list(nx.all_simple_paths(graph, start_nodes[0], end_nodes[0]))

        # Generate combinations of nodes that satisfy certain conditions
        for p_arr in all_paths:
            for i in range(0, len(p_arr) - self.fsa_file.ref_count + 1):
                yield np.array(p_arr[i : i + self.fsa_file.ref_count])

    def get_best_fit(self, combinations, method="2nd_derivative"):
        """
        Calculates the best fit for a given set of combinations using a specified method.

        Parameters:
        combinations (iterable): an iterable containing combinations of nodes
        method (str): a string indicating the method to use. The default method is
            "2nd_derivative".

        Returns:
        A numpy array representing the combination of nodes with the best fit.
        """

        # Create an empty dataframe to store combinations and scores
        df = pd.DataFrame()

        # Store the combinations in the dataframe
        df["combination"] = list(combinations)

        # Calculate the score for each combination using the specified method
        if method == "2nd_derivative":
            df["score"] = np.vectorize(self._max_spline_second_derivative_score)(
                df["combination"]
            )

            # Sort the dataframe by score in ascending order
            df = df.sort_values(by="score", ascending=True)

            # Get the best combination (i.e., the one with the lowest score)
            best = df.head(1)

            # Return the best combination as a numpy array
            return best.combination.squeeze()

    @staticmethod
    def _polynomial_model_inv_r2_score(ladder: np.array, comb: np.array) -> float:
        """
        Calculates the inverse of the R^2 (coefficient of determination) score for a polynomial regression model.

        Parameters:
        - ladder (np.array): array containing the ladder positions of a dataset
        - comb (np.array): array containing the corresponding comb heights of a dataset

        Returns:
        - float: the inverse of the R^2 score for the polynomial model

        Notes:
        - This method fits a polynomial regression model of degree 3 to the dataset using NumPy's polyfit function.
        - It then calculates the predicted values of the comb heights using NumPy's poly1d function.
        - Finally, it calculates the R^2 score between the predicted comb heights and the actual comb heights using the r2_score function
          from the scikit-learn library, and returns the inverse of that score.

        """
        # fit polynomial regression model
        fit = np.polyfit(ladder, comb, 3)

        # create predicted values using fitted model
        predict = np.poly1d(fit)

        # calculate R^2 score between predicted and actual values, and return its inverse
        return 1 - r2_score(comb, predict(ladder))

    @staticmethod
    def _max_fractional_deviation_score(ladder: np.ndarray, comb: np.ndarray):
        """
        Calculates the maximum fractional deviation score for a dataset.

        Parameters:
        - ladder (np.ndarray): array containing the ladder positions of a dataset
        - comb (np.ndarray): array containing the corresponding comb heights of a dataset

        Returns:
        - float: the maximum fractional deviation score

        Notes:
        - This method calculates the interval lengths between adjacent ladder positions and corresponding comb heights using NumPy's
          diff function and min-max scaling with a feature range equal to the range of ladder positions.
        - It then calculates the fractional deviation between the interval lengths and the scaled comb heights, and returns the
          maximum deviation.

        """
        # calculate interval lengths and scaled comb height intervals
        l_intervals = np.diff(ladder)
        c_intervals = np.diff(minmax_scale(comb, feature_range=(ladder[0], ladder[-1])))

        # calculate fractional deviation between interval lengths and scaled comb height intervals
        frac_deviation = np.abs(l_intervals - c_intervals) / l_intervals

        # find index of maximum deviation and return the corresponding value
        max_frac_deviation_idx = np.argmax(frac_deviation)
        return frac_deviation[max_frac_deviation_idx]

    def _max_first_derivative_score(self, combination: np.ndarray) -> float:
        """
        Calculates the maximum first derivative score for the given combination of features.

        Args:
            combination (np.ndarray): A 1-dimensional numpy array containing a combination of features.

        Returns:
            float: The maximum first derivative score for the given combination of features.
        """
        # Scale the combination to the range of reference sizes
        comb_scaled = minmax_scale(
            combination,
            feature_range=(self.fsa_file.ref_sizes[0], self.fsa_file.ref_sizes[-1]),
        )

        # Calculate the differences between consecutive values of the scaled combination and the reference sizes
        diff_intervals = np.diff(comb_scaled) - np.diff(self.fsa_file.ref_sizes)

        # Calculate the absolute gradient of the differences between consecutive values
        abs_diff_intervals_gradient = np.abs(np.gradient(diff_intervals))

        # Find the index of the maximum absolute gradient
        max_abs_diff_intervals_gradient_idx = np.argmax(abs_diff_intervals_gradient)

        # Return the maximum absolute gradient
        return abs_diff_intervals_gradient[max_abs_diff_intervals_gradient_idx]

    def _max_second_derivative_score(self, combination: np.ndarray) -> float:
        """
        Calculates the maximum absolute second derivative score of the given combination.

        Args:
            combination (np.ndarray): An array of feature values for a given combination.

        Returns:
            float: The maximum absolute second derivative score of the given combination.

        Raises:
            None.

        """
        # Scale the combination to match the range of reference sizes.
        comb_scaled = minmax_scale(
            combination,
            feature_range=(self.fsa_file.ref_sizes[0], self.fsa_file.ref_sizes[-1]),
        )

        # Calculate the difference between the scaled combination and reference sizes.
        diff_intervals = np.diff(comb_scaled) - np.diff(self.fsa_file.ref_sizes)

        # Calculate the absolute second derivative of the difference intervals.
        abs_second_derivative = np.abs(np.gradient(np.gradient(diff_intervals)))

        # Find the index of the maximum absolute second derivative score.
        max_second_derivative_idx = np.argmax(abs_second_derivative)

        # Return the maximum absolute second derivative score.
        return abs_second_derivative[max_second_derivative_idx]

    def _max_spline_second_derivative_score(self, combination: np.ndarray) -> float:
        """
        Calculates the maximum absolute value of the second derivative of a given combination of values using a UnivariateSpline.

        Args:
        - self (object): an instance of the class containing the function
        - combination (np.ndarray): an array of values representing a combination

        Returns:
        - max_value (float): the maximum absolute value of the second derivative of the given combination

        """
        # create a UnivariateSpline object using the reference sizes and the given combination
        spl = UnivariateSpline(self.fsa_file.ref_sizes, combination, s=0)

        # calculate the second derivative of the spline
        der2 = spl.derivative(n=2)

        # evaluate the second derivative at each reference size and return the maximum absolute value
        max_value = max(abs(der2(self.fsa_file.ref_sizes)))

        return max_value
