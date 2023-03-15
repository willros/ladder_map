import matplotlib.pyplot as plt
from fragment_analyzer.ladder_fitting.fit_ladder_model import FitLadderModel


class PlotRawData:
    def __init__(self, model: FitLadderModel):
        self.model = model

    @property
    def plot_raw_data(self):
        data = self.model.adjusted_baisepair_df
        fig = plt.figure(figsize=(20, 10))

        plt.plot(data.basepairs, data.peaks)
        plt.xlabel("Basepairs")
        plt.ylabel("Intensity")

        plt.close()
        return fig
