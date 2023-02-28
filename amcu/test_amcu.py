from amcu.automcu import AutoMCU
import pytest


# @pytest.fixture
def test_amcu():
    em_csvs = ["data/EM1.csv", "data/EM2.csv"]
    em_counts = [2]
    em_names = ["EM1", "EM2"]
    band_ranges = [[0, 1], [1, 2]]
    amcu = AutoMCU()

    output = amcu.apply_image("../data/input_img", "../data/output_img.tif")
    assert output.shape == (10, 1)
