"""Module providing unit test for queueing_theory"""
import math as ma

import queueing_theory as qt


class TestMm1:
    """Class to test the M/M/1 queueing theory model"""

    def test_mm1_model_compute_rho_expected_valid_output(self):
        """
        Test case for mm1_model_compute_rho function.

        This function tests the calculation of rho (traffic intensity) in the M/M/1
        queueing model.
        It verifies that the calculated value of rho is not equal to 0 and is
        approximately equal to 0.8.

        """
        lam_input = 8
        miu_input = 10
        expected_output = 0.8

        rho_calc = qt.mm1_model_compute_rho(lam = lam_input, miu = miu_input)

        assert rho_calc != 0, "rho is equal to 0"
        assert ma.isclose(rho_calc, expected_output,
                          rel_tol=1e-9), "rho is not approximately equal to 0.8"
