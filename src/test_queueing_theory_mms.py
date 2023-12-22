import math as ma
import queueing_theory as qt
import pytest as pyt


class TestMMs:
    """Class to test the M/M/s queuing theory model"""

    # Test for valid inputs and stable system
    def test_mms_model_compute_p_zero_valid_input_stable_system(self):
        lam, miu, s = 10, 10, 2
        expected_probability = 1 / (
            1 + (lam / miu) + ((lam / miu) ** 2 / (2 * (1 - (lam / (s * miu)))))
        )
        assert qt.mms_model_compute_p_zero(lam, miu, s) == pyt.approx(
            expected_probability
        )

    # Test for unstable system (rho >= 1)
    def test_mms_model_compute_p_zero_unstable_system(self):
        lam, miu, s = 20, 10, 1
        with pyt.raises(NameError) as exception_info:
            qt.mms_model_compute_p_zero(lam, miu, s)
        assert str(exception_info.value) == qt.UNSTABLE_MESSAGE

    # Test for negative or zero input values
    def test_mms_model_compute_p_zero_negative_zero_input(self):
        test_cases = [(-10, 10, 1), (10, -10, 1), (10, 10, 0)]
        for lam, miu, s in test_cases:
            with pyt.raises(NameError) as exception_info:
                qt.mms_model_compute_p_zero(lam, miu, s)
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR
