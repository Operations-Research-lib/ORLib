"""Module providing unit test for queueing_theory"""
import math as ma
from pytest import approx, raises
import queueing_theory as qt


class TestMm1:
    """Class to test the M/M/1 queueing theory model"""


    def test_mm1_model_info_expected_valid_output(self):
        """
        Test case for mm1_model_info function.

        This function tests the calculation of the M/M/1 queueing model with valid
        inputs.
        """
        # arrange
        lam_input = 20
        miu_input = 30
        expected_outputs = [2/3, 4/3, 2, 1/15, 1/10]

        # act
        mm1_expected_info = qt.mm1_model_info(lam = lam_input, miu = miu_input)

        # assert
        assert len(mm1_expected_info) == 5
        for i, expected_info_output in enumerate(mm1_expected_info):
            assert expected_info_output != 0, "One of the outputs is equal to 0"
            assert expected_info_output == approx(expected_outputs[i], rel=1e-9)



    def test_mm1_model_info_unstable_system_input_raises_exception(self):
        """
        Test case for mm1_model_compute_lq function.

        This function test if the input raise an exception when the system is
        unstable.
        """
        # arrange
        input_cases = [(10,10), (10,5)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_info(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.UNSTABLE_MESSAGE


    def test_mm1_model_info_negative_zero_input_raises_exception(self):
        """
        Test case for mm1_model_compute_lq function.

        This function tests if some input have negative or zero values.
        """

        # arrange
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_info(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR


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



    def test_mm1_model_compute_lq_expected_valid_output(self):
        """
        Test case for mm1_model_compute_lq function.

        This function tests the calculation of Lq in the M/M/1 queueing model.
        """
        lam_input = 20
        miu_input = 30
        expected_output = 4/3

        lq_calc = qt.mm1_model_compute_lq(lam = lam_input, miu = miu_input)

        assert lq_calc != 0, "Lq is equal to 0"
        assert ma.isclose(lq_calc, expected_output,
                          rel_tol=1e-9), "Lq is not approximately equal to 1,333..."


    def test_mm1_model_compute_lq_unstable_system_input_raises_exception(self):
        """
        Test case for mm1_model_compute_lq function.

        This function test if the input raise an exception when the system is
        unstable.
        """
        # arrange
        input_cases = [(10,10), (10,5)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_lq(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.UNSTABLE_MESSAGE


    def test_mm1_model_compute_lq_negative_zero_input_raises_exception(self):
        """
        Test case for mm1_model_compute_lq function.

        This function tests if some input have negative or zero values.
        """

        # arrange
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_lq(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR


    def test_mm1_model_compute_l_expected_valid_output(self):
        """
        Test case for mm1_model_compute_l function.

        This function tests the calculation of L in the M/M/1 queueing model.
        """
        lam_input = 20
        miu_input = 30
        expected_output = 2

        lq_calc = qt.mm1_model_compute_l(lam = lam_input, miu = miu_input)

        assert lq_calc != 0, "L is equal to 0"
        assert ma.isclose(lq_calc, expected_output,
                          rel_tol=1e-9), "L is not approximately equal to 2"


    def test_mm1_model_compute_l_unstable_system_input_raises_exception(self):
        """
        Test case for mm1_model_compute_l function.

        This function test if the input raise an exception when the system is
        unstable.
        """
        # arrange
        input_cases = [(10,10), (10,5)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_l(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.UNSTABLE_MESSAGE


    def test_mm1_model_compute_l_negative_zero_input_raises_exception(self):
        """
        Test case for mm1_model_compute_l function.

        This function tests if some input have negative or zero values.
        """

        # arrange
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_l(lam = lam_input, miu = miu_input)
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
                          rel_tol=1e-9), "Wq is not approximately equal to 0.666..."


    def test_mm1_model_compute_wq_unstable_system_input_raises_exception(self):
        """
        Test case for mm1_model_compute_wq function.

        This function test if the input raise an exception when the system is
        unstable.
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
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_wq(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR


    def test_mm1_model_compute_w_expected_valid_output(self):
        """
        Test case for mm1_model_compute_w function.

        This function tests the calculation of w (expected waiting time in the
        queue) in the M/M/1 queueing model.
        """
        # arrange
        lam_input = 20
        miu_input = 30
        expected_output = 1/10

        # act
        wq_calc = qt.mm1_model_compute_w(lam = lam_input, miu = miu_input)

        # assert
        assert wq_calc != 0, "w is equal to 0"
        assert ma.isclose(wq_calc, expected_output,
                          rel_tol=1e-9), "w is not approximately equal to 0.1"


    def test_mm1_model_compute_w_unstable_system_input_raises_exception(self):
        """
        Test case for mm1_model_compute_w function.

        This function tests if the input raise an exception when the system is
        unstable.
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
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_w(lam = lam_input, miu = miu_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR


    def test_mm1_model_compute_pn_expected_valid_output(self):
        """
        Test case for mm1_model_compute_pn function.

        This function tests the calculation of pn, n = 0,1,2.
        """
        # arrange
        input_output_cases = [(20,30,0,1/3), (20,30,1,2/9), (20,30,2,4/27)]

        for lam_input, miu_input, n_input, expected_output in input_output_cases:
            # act
            wq_calc = qt.mm1_model_compute_pn(lam = lam_input, miu = miu_input, n = n_input)
            # assert
            assert wq_calc != 0, "pn is equal to 0"
            assert ma.isclose(wq_calc, expected_output,
                            rel_tol=1e-9), "pn is not approximately equal"


    def test_mm1_model_compute_pn_unstable_system_input_raises_exception(self):
        """
        Test case for mm1_model_compute_pn function.

        This function tests if the input raise an exception when the system is
        unstable.
        """
        # arrange
        input_cases = [(10,10,1), (10,5,1)]

        for lam_input, miu_input, n_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_pn(lam = lam_input, miu = miu_input, n = n_input)
            # assert
            assert str(exception_info.value) == qt.UNSTABLE_MESSAGE


    def test_mm1_model_compute_pn_negative_zero_input_raises_exception(self):
        """
        Test case for mm1_model_compute_pn function.

        This function tests if some input have negative or zero values.
        """

        # arrange
        input_cases = [(-10,10,1), (10,-10,1), (10,0,1), (0,10,1), (10,10,-1)]

        for lam_input, miu_input, n_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1_model_compute_pn(lam = lam_input, miu = miu_input, n = n_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_INPUT_ERROR