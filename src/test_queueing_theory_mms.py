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
        assert qt.mms_model_compute_p_zero(lam, miu, s) == approx(
            expected_probability
        ), "Did not return expected result in compute p_zero."

    # Test for unstable system (rho >= 1)
    def test_mms_model_compute_p_zero_unstable_system_raises_exception(self):
        # arrange
        lam, miu, s = 20, 10, 1
        with raises(NameError) as exception_info:
            # act
            qt.mms_model_compute_p_zero(lam, miu, s)
            # assert
        assert (
            str(exception_info.value) == qt.UNSTABLE_MESSAGE
        ), "Did not raise error when system is unstable in compute p_zero."

    # Test for negative or zero input values
    def test_mms_model_compute_p_zero_negative_zero_input_raises_exception(self):
        # arrange
        test_cases = [(-10, 10, 1), (10, -10, 1), (10, 10, 0)]
        for lam, miu, s in test_cases:
            with raises(NameError) as exception_info:
                # act
                qt.mms_model_compute_p_zero(lam, miu, s)
                # assert
            assert (
                str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR
            ), "Did not raise error when invalid inputs in compute p_zero."

    # Test for unstable system (rho >= 1)
    def test_mms_model_compute_lq_unstable_system_raises_exception(self):
        # arrange
        lam, miu, s = 20, 10, 1
        # act
        with raises(NameError) as exception_info:
            qt.mms_model_compute_lq(lam, miu, s)
        # assert
        assert (
            str(exception_info.value) == qt.UNSTABLE_MESSAGE
        ), "Did not raise error when system is unstable in compute lq."

    def test_mms_model_compute_lq_negative_zero_input_raises_exception(self):
        # arrange
        test_cases = [(-10, 10, 1), (10, -10, 1), (10, 10, 0)]
        for lam, miu, s in test_cases:
            with raises(NameError) as exception_info:
                # act
                qt.mms_model_compute_lq(lam, miu, s)
                # assert
            assert (
                str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR
            ), "Did not raise error when invalid input in compute lq."

    # test stable system valid input
    def test_mms_model_compute_lq_valid_input_stable_system_returns_expected(self):
        # arrange
        lam = 25
        miu = 20
        s = 2
        lq_expected = 0.801282051
        # act => assert
        assert qt.mms_model_compute_lq(lam, miu, s) == approx(
            lq_expected
        ), "Did not give expected result in compute wq."

    # test compute wq
    # testing for unstable systems
    def test_mms_model_compute_wq_unstable_system_raises_exception(self):
        # arrange
        lam, miu, s = 20, 10, 1
        # act
        with raises(NameError) as exception_info:
            qt.mms_model_compute_wq(lam, miu, s)
        # assert
        assert (
            str(exception_info.value) == qt.UNSTABLE_MESSAGE
        ), "Did not raise error when system is unstable in compute wq."

    # testing for invalid inputs
    def test_mms_model_compute_wq_negative_zero_input_raises_exception(self):
        # arrange
        test_cases = [(-10, 10, 1), (10, -10, 1), (10, 10, 0)]
        for lam, miu, s in test_cases:
            with raises(NameError) as exception_info:
                # act
                qt.mms_model_compute_wq(lam, miu, s)
                # assert
            assert (
                str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR
            ), "Did not raise error with invalid inputs in compute wq."

    # test compute w
    # testing with valid inputs
    def test_mms_model_wq_valid_input_stable_system_expected_result(self):
        # arrange
        lam = 25
        miu = 20
        s = 2
        wq_expected = 0.032051282
        # act => assert
        assert qt.mms_model_compute_wq(lam, miu, s) == approx(
            wq_expected
        ), "Did not return expected result in compute wq."

    # testing with invalid inputs
    def test_mms_model_compute_w_negative_zero_input_raises_exception(self):
        # arrange
        test_cases = [(-10, 10, 1), (10, -10, 1), (10, 10, 0)]
        for lam, miu, s in test_cases:
            with raises(NameError) as exception_info:
                # act
                qt.mms_model_compute_w(lam, miu, s)
                # assert
            assert (
                str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR
            ), "Did not raise error with invalid inputs in compute w."

    # testing for unstable systems
    def test_mms_model_compute_w_valid_input_stable_system_expected_result(self):
        # arrange
        lam = 25
        miu = 20
        s = 2
        w_expected = 0.082051282
        # act => assert
        assert qt.mms_model_compute_w(lam, miu, s) == approx(w_expected)
