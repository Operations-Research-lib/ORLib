import queueing_theory as qt
import math as m
from pytest import approx, raises


class TestMMs:
    """Class to test the M/M/s queuing theory model"""

    # Test for valid inputs and stable system
    def test_mms_model_compute_p_zero_valid_input_stable_system_returns_expected(self):
        # arrange
        lam, miu, s = 10, 10, 2
        expected_probability = 1 / (
            1 + (lam / miu) + ((lam / miu) ** 2 / (2 * (1 - (lam / (s * miu)))))
        )
        # act => assert
        assert qt.mms_model_compute_p_zero(lam, miu, s) == approx(expected_probability)

    # Test for unstable system (rho >= 1)
    def test_mms_model_compute_p_zero_unstable_system_raises_exception(self):
        # arrange
        lam, miu, s = 20, 10, 1
        with raises(NameError) as exception_info:
            # act
            qt.mms_model_compute_p_zero(lam, miu, s)
            # assert
        assert str(exception_info.value) == qt.UNSTABLE_MESSAGE

    # Test for negative or zero input values
    def test_mms_model_compute_p_zero_negative_zero_input_raises_exception(self):
        # arrange
        test_cases = [(-10, 10, 1), (10, -10, 1), (10, 10, 0)]
        for lam, miu, s in test_cases:
            with raises(NameError) as exception_info:
                # act
                qt.mms_model_compute_p_zero(lam, miu, s)
                # assert
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR

    # Test for unstable system (rho >= 1)
    def test_mms_model_compute_lq_unstable_system_raises_exception(self):
        # arrange
        lam, miu, s = 20, 10, 1
        # act
        with raises(NameError) as exception_info:
            qt.mms_model_compute_lq(lam, miu, s)
        # assert
        assert str(exception_info.value) == qt.UNSTABLE_MESSAGE

    def test_mms_model_compute_lq_negative_zero_input_raises_exception(self):
        test_cases = [(-10, 10, 1), (10, -10, 1), (10, 10, 0)]
        for lam, miu, s in test_cases:
            with raises(NameError) as exception_info:
                # act
                qt.mms_model_compute_lq(lam, miu, s)
                # assert
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR

    # test stable system valid input
    def test_mms_model_compute_lq_valid_input_stable_system_returns_expected(self):
        # arrange
        lam = 25
        miu = 20
        s = 2
        lq_expected = 0.801282051
        # act => assert
        assert qt.mms_model_compute_lq(lam, miu, s) == approx(lq_expected)
