from amcu.automcu import AutoMCU
import pytest


def test_automcu():
    """
    test_automcu is a test to run the code without bug and get sure that there's no bug in the main body of the code.
    """
    amcu = AutoMCU()
    input_img = "\\gdcs-stor.rc.asu.edu\gdcs-stor\CarbonMapper\AlgDevel\Agriculture\Florida_Apr2022\east_agriculture\GAO20220411t125900p0000"
    output_path = "data\GAO20220411t125900p0000"
    scale = 0.1

    amcu.apply_image(
        refl_path=input_img,
        input_scale=scale,
        output_path=output_path,
        hdrpath=input_img,
    )
