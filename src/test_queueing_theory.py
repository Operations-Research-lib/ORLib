"""Module providing unit test for queueing_theory"""
import math as ma
from pytest import raises
import queueing_theory as qt


class TestMm1:
    """Class to test the M/M/1 queueing theory model"""

    def test_mm1_model_compute_rho_expected_valid_output(self):
        """
        Test case for mm1_model_compute_rho function.

        This function tests the calculation of rho (traffic intensity) in the M/M/1
        queueing model with valid inputs.
        """
        # arrange
        lam_input = 20
        miu_input = 30
        expected_output = 2/3

        # act
        rho_calc = qt.mm1_model_compute_rho(lam = lam_input, miu = miu_input)

        # assert
        assert rho_calc != 0, "rho is equal to 0"
        assert ma.isclose(rho_calc, expected_output,
                          rel_tol=1e-9), "rho is not approximately equal to 0.8"


    def test_mm1_model_compute_rho_negative_zero_input_raises_exception(self):
        """
        Test case for mm1_model_compute_rho function.

        This function tests if some input have negative or zero values. Expected
        output is an exception.
        """
        # arrange
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_rho(lam = lam_input, miu = miu_input)

            # assert
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR


    def test_mm1_model_compute_wq_expected_valid_output(self):
        """
        Test case for mm1_model_compute_wq function.

        This function tests the calculation of Wq (expected waiting time in the
        queue) in the M/M/1 queueing model.
        """
        # arrange
        lam_input = 20
        miu_input = 30
        expected_output = 1/15

        # act
        wq_calc = qt.mm1_model_compute_wq(lam = lam_input, miu = miu_input)

        # assert
        assert wq_calc != 0, "Wq is equal to 0"
        assert ma.isclose(wq_calc, expected_output,
                          rel_tol=1e-9), "Wq is not approximately equal to 0.033"


    def test_mm1_model_compute_wq_unstable_system_input_raises_exception(self):
        """
        Test case for mm1_model_compute_wq function.

        This function tests the if there is input

        """
        # arrange
        input_cases = [(10,10), (10,5)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_wq(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.UNSTABLE_MESSAGE


    def test_mm1_model_compute_wq_negative_zero_input_raises_exception(self):
        """
        Test case for mm1_model_compute_wq function.

        This function tests if some input have negative or zero values.
        """

        # arrange
        input_cases = [(-10,10), (10,-10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_wq(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR


    def test_mm1_model_compute_w_expected_valid_output(self):
        """
        Test case for mm1_model_compute_w function.

        This function tests the calculation of Wq (expected waiting time in the
        queue) in the M/M/1 queueing model.
        """
        # arrange
        lam_input = 20
        miu_input = 30
        expected_output = 1/10

        # act
        wq_calc = qt.mm1_model_compute_w(lam = lam_input, miu = miu_input)

        # assert
        assert wq_calc != 0, "Wq is equal to 0"
        assert ma.isclose(wq_calc, expected_output,
                          rel_tol=1e-9), "Wq is not approximately equal to 0.033"


    def test_mm1_model_compute_w_unstable_system_input_raises_exception(self):
        """
        Test case for mm1_model_compute_w function.

        This function tests the if there is input

        """
        # arrange
        input_cases = [(10,10), (10,5)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_w(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.UNSTABLE_MESSAGE


    def test_mm1_model_compute_w_negative_zero_input_raises_exception(self):
        """
        Test case for mm1_model_compute_w function.

        This function tests if some input have negative or zero values.
        """

        # arrange
        input_cases = [(-10,10), (10,-10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_w(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR
