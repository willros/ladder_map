"""
Class to fit model to ladder and size standard
"""

from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import SplineTransformer
from sklearn.metrics import mean_squared_error, r2_score
import pandas as pd

from fragment_analyzer.ladder_fitting.peak_ladder_assigner import PeakLadderAssigner
from fragment_analyzer.utils.fsa_file import FsaFile


class FitLadderModel:
    def __init__(self, ladder_assigner: PeakLadderAssigner) -> None:
        """
        Initialize FitLadderModel object with ladder_assigner object and its attributes.

        Args:
            ladder_assigner (PeakLadderAssigner): Object of PeakLadderAssigner class.

        Returns:
            None.
        """
        self.ladder_assigner = ladder_assigner
        self.fsa_file = self.ladder_assigner.fsa_file
        self.best_combination = ladder_assigner.assign_ladder_peak_sizes().reshape(
            -1, 1
        )

        self.fit_model()
        self.mse, self.r2 = self.model_score()
        self.adjusted_baisepair_df = self.generate_adjusted_trace_df()

    def fit_model(self) -> None:
        """
        Fit model based on ladder type.

        Returns:
            None.
        """
        match self.fsa_file.ladder:
            case "ROX":
                self._fit_ROX_ladder()
            case "LIZ":
                self._fit_LIZ_ladder()
            case _:
                print("Ladder not found")

    def model_score(self) -> tuple[float, float]:
        """
        Calculate mean squared error and R-squared score of the fitted model.

        Returns:
            Tuple (mse, r2).
        """
        true_Y = self.fsa_file.ref_sizes
        predicted_Y = self.model.predict(self.best_combination)
        mse = mean_squared_error(true_Y, predicted_Y)
        r2 = r2_score(true_Y, predicted_Y)

        return mse, r2

    def generate_adjusted_trace_df(self) -> pd.DataFrame:
        """
        Generate a dataframe with adjusted basepairs and peaks.

        Returns:
            Dataframe with columns time, peaks and basepairs.
        """
        df = (
            pd.DataFrame({"peaks": self.fsa_file.trace})
            .reset_index()
            .rename(columns={"index": "time"})
            .assign(
                basepairs=lambda x: self.model.predict(x.time.to_numpy().reshape(-1, 1))
            )
            .loc[lambda x: x.basepairs >= 0]
        )

        return df

    def _fit_ROX_ladder(self) -> None:
        """
        Fit model with ROX ladder.

        Returns:
            None.
        """
        self.model = make_pipeline(
            SplineTransformer(degree=4, n_knots=6, extrapolation="continue"),
            LinearRegression(fit_intercept=True),
        )

        X = self.best_combination
        y = self.fsa_file.ref_sizes

        self.model.fit(X, y)

    def _fit_LIZ_ladder(self) -> None:
        """
        Fit model with LIZ ladder.

        Returns:
            None.
        """
        self.model = make_pipeline(
            SplineTransformer(degree=3, n_knots=3, extrapolation="continue"),
            LinearRegression(fit_intercept=True),
        )

        X = self.best_combination
        y = self.fsa_file.ref_sizes

        self.model.fit(X, y)
