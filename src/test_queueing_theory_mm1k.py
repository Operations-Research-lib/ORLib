"""Module providing unit test for queueing_theory"""
from pytest import approx, raises
import queueing_theory as qt


class TestMm1k:
    """Class to test the M/M/1/K queueing theory model"""

    def test_mm1k_model_compute_p_zero_expected_valid_output(self):
        """
        Test case for mm1k_model_compute_p_zero function.

        This function tests the calculation of P0 in the M/M/1/K queueing model.
        """
        # arrange
        lam_input = 20
        miu_input = 30
        k_input = 2
        expected_output = 9/19

        # act
        p_zero_calc = qt.mm1k_model_compute_p_zero(lam = lam_input, miu = miu_input, k = k_input)

        # assert
        assert p_zero_calc != 0, "P0 is equal to 0"
        assert p_zero_calc == approx(
            expected_output, rel=1e-2
        ), "Did not return expected result in compute p0."


    def test_mm1k_model_compute_p_zero_negative_zero_lam_miu_inputs_raises_exception(self):
        """
        Test case for mm1k_model_compute_p_zero function.

        This function tests if lambda or miu have negative or zero values.
        """
        # arrange
        k_input = 2
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_p_zero(lam = lam_input, miu = miu_input, k = k_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_LAM_MIU_ERROR


    def test_mm1k_model_compute_p_zero_negative_zero_k_input_raises_expection(self):
        """
        Test case for mm1k_model_compute_p_zero function.

        This function tests if the k input raise an exception when is zero or negative.
        """
        k_cases = [0,-1]
        n_input = 1
        lam_input = 20
        miu_input = 30

        for k_input in k_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_pn(lam = lam_input, miu = miu_input, k = k_input, n = n_input)

            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_K_ERROR


    def test_mm1k_model_compute_lq_expected_valid_output(self):
        """
        Test case for mm1k_model_compute_lq function.

        This function tests the calculation of Lq in the M/M/1/K queueing model.
        """
        # arrange
        lam_input = 20
        miu_input = 30
        k_input = 2
        expected_output = 0.21053

        # act
        lq_calc = qt.mm1k_model_compute_lq(lam = lam_input, miu = miu_input, k = k_input)

        # assert
        assert lq_calc != 0, "Lq is equal to 0"
        assert lq_calc == approx(
            expected_output, rel=1e-2
        ), "Did not return expected result in compute lq."


    def test_mm1k_model_compute_lq_negative_zero_lambda_miu_inputs_raises_exception(self):
        """
        Test case for mm1k_model_compute_lq function.

        This function tests if lambfa or miu input have negative or zero values.
        """
        # arrange
        k_input = 2
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_lq(lam = lam_input, miu = miu_input, k = k_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_LAM_MIU_ERROR


    def test_mm1k_model_compute_lq_negative_zero_k_input_raises_expection(self):
        """
        Test case for mm1k_model_compute_lq function.

        This function tests if the k input raise an exception when is zero or negative.
        """
        k_cases = [0,-1]
        lam_input = 20
        miu_input = 30

        for k_input in k_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_lq(lam = lam_input, miu = miu_input, k = k_input)

            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_K_ERROR


    def test_mm1k_model_compute_l_expected_valid_output(self):
        """
        Test case for mm1k_model_compute_l function.

        This function tests the calculation of L in the M/M/1/K queueing model.
        """
        k_input = 2
        lam_input = 20
        miu_input = 30
        expected_output = 0.73684

        l_calc = qt.mm1k_model_compute_l(lam = lam_input, miu = miu_input, k = k_input)

        assert l_calc != 0, "L is equal to 0"
        assert l_calc == approx(
            expected_output, rel=1e-2
        ), "Did not return expected result in compute l."


    def test_mm1k_model_compute_l_negative_zero_lambda_miu_inputs_raises_exception(self):
        """
        Test case for mm1k_model_compute_l function.

        This function tests if lamda or miu inputs have negative or zero values.
        """
        # arrange
        k_input = 2
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_l(lam = lam_input, miu = miu_input, k = k_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_LAM_MIU_ERROR


    def test_mm1k_model_compute_l_negative_zero_k_input_raises_expection(self):
        """
        Test case for mm1k_model_compute_l function.

        This function tests if the k input raise an exception when is zero or negative.
        """
        k_cases = [0,-1]
        lam_input = 20
        miu_input = 30

        for k_input in k_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_l(lam = lam_input, miu = miu_input, k = k_input)

            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_K_ERROR


    def test_mm1k_model_compute_wq_expected_valid_output(self):
        """
        Test case for mm1k_model_compute_wq function.

        This function tests the calculation of Wq (expected waiting time in the
        queue) in the M/M/1/K queueing model.
        """
        # arrange
        lam_input = 20
        miu_input = 30
        k_input = 2
        expected_output = 0.01333

        # act
        wq_calc = qt.mm1k_model_compute_wq(lam = lam_input, miu = miu_input, k = k_input)

        # assert
        assert wq_calc != 0, "Wq is equal to 0"
        assert wq_calc == approx(
            expected_output, rel=1e-2
        ), "Did not return expected result in compute lq."



    def test_mm1k_model_compute_wq_negative_zero_lambda_miu_input_raises_exception(self):
        """
        Test case for mm1k_model_compute_wq function.

        This function tests if lambda or miu inputs have negative or zero values.
        """
        # arrange
        k_input = 2
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_wq(lam = lam_input, miu = miu_input, k = k_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_LAM_MIU_ERROR


    def test_mm1k_model_compute_wq_negative_zero_k_input_raises_expection(self):
        """
        Test case for mm1k_model_compute_wq function.

        This function tests if the k input raise an exception when is zero or negative.
        """
        k_cases = [0,-1]
        lam_input = 20
        miu_input = 30

        for k_input in k_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_wq(lam = lam_input, miu = miu_input, k = k_input)

            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_K_ERROR


    def test_mm1k_model_compute_w_expected_valid_output(self):
        """
        Test case for mm1k_model_compute_w function.

        This function tests the calculation of w (expected waiting time in the
        queue) in the M/M/1/K queueing model.
        """
        # arrange
        lam_input = 20
        miu_input = 30
        k_input = 2
        expected_output = 0.0467

        # act
        w_calc = qt.mm1k_model_compute_w(lam = lam_input, miu = miu_input, k = k_input)

        # assert
        assert w_calc != 0, "w is equal to 0"
        assert w_calc == approx(
            expected_output, rel=1e-2
        ), "Did not return expected result in compute w."


    def test_mm1k_model_compute_w_negative_zero_lambda_miu_inputs_raises_exception(self):
        """
        Test case for mm1k_model_compute_w function.

        This function tests if lambda or miu inputs have negative or zero values.
        """
        # arrange
        k_input = 2
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_w(lam = lam_input, miu = miu_input, k = k_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_LAM_MIU_ERROR


    def test_mm1k_model_compute_w_negative_zero_k_input_raises_expection(self):
        """
        Test case for mm1k_model_compute_p_zero function.

        This function tests if the k input raise an exception when is zero or negative.
        """
        k_cases = [0,-1]
        lam_input = 20
        miu_input = 30

        for k_input in k_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_w(lam = lam_input, miu = miu_input, k = k_input)

            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_K_ERROR


    def test_mm1k_model_compute_pn_expected_valid_output(self):
        """
        Test case for mm1k_model_compute_pn function.

        This function tests the calculation of pn, n = 0,1,2.
        """
        # arrange
        k_input = 2
        lam_input = 20
        miu_input = 30
        n_input = [0,1,2,3]
        expected_output_p0 = 9/19
        expected_output_p1 = 6/19
        expected_output_p2 = 4/19
        expected_output_p3 = 0.0

        input_output_cases = [(lam_input, miu_input, n_input[0], expected_output_p0),
                              (lam_input, miu_input, n_input[1], expected_output_p1),
                              (lam_input, miu_input, n_input[2], expected_output_p2),
                              (lam_input, miu_input, n_input[3], expected_output_p3)]

        for lam_input, miu_input, n_input, expected_output in input_output_cases:
            # act
            pn_calc = qt.mm1k_model_compute_pn(lam = lam_input, miu = miu_input, k = k_input,
                                               n = n_input)
            # assert
            assert pn_calc != 0, "pn is equal to 0"
            assert pn_calc == approx(
                expected_output, rel=1e-2
            ), "Did not return expected result in compute pn."



    def test_mm1k_model_compute_pn_negative_zero_lambda_miu_inputs_raises_exception(self):
        """
        Test case for mm1k_model_compute_pn function.

        This function tests if some input have negative or zero values.
        """
        # arrange
        k_input = 2
        n_input = 1
        input_cases = [(-10,10), (10,-10), (10,0), (0,10)]

        for lam_input, miu_input in input_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_pn(lam = lam_input, miu = miu_input, k = k_input, n = n_input)
            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_LAM_MIU_ERROR


    def test_mm1k_model_compute_pn_negative_zero_k_input_raises_expection(self):
        """
        Test case for mm1k_model_compute_pn function.

        This function tests if the input raise an exception when the system is
        unstable.
        """
        # arrange
        k_cases = [0,-1]
        n_input = 1
        lam_input = 20
        miu_input = 30

        for k_input in k_cases:
            with raises(Exception) as exception_info:
                # act
                qt.mm1k_model_compute_pn(lam = lam_input, miu = miu_input, k = k_input, n = n_input)

            # assert
            assert str(exception_info.value) == qt.NEGATIVE_ZERO_K_ERROR


    def test_mm1k_model_compute_pn_negative_n_raises_exception(self):
        """
        Test case for mm1k_model_compute_pn function.

        This function tests if n input have negative value.
        """
        # arrange
        k_input = 2
        n_input = -1
        lam_input = 20
        miu_input = 30

        # act
        with raises(Exception) as exception_info:
            qt.mm1k_model_compute_pn(lam = lam_input, miu = miu_input, k = k_input, n = n_input)

        # assert
        assert str(exception_info.value) == qt.NEGATIVE_N_ERROR
