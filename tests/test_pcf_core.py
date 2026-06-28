import numpy as np
import pytest
from code.pcf_core import (
    compute_first_difference,
    compute_second_difference,
    lag_1_autocorrelation,
    compute_pcf,
    test_against_white_noise,
    detrend_series,
    windowed_pcf
)


class TestDifferences:
    """Test difference operators"""
    
    def test_first_difference_length(self):
        """First difference should reduce length by 1"""
        X = np.array([1, 2, 4, 7, 11])
        d1 = compute_first_difference(X)
        assert len(d1) == len(X) - 1
        assert np.allclose(d1, [1, 2, 3, 4])
    
    def test_second_difference_length(self):
        """Second difference should reduce length by 2"""
        X = np.array([1, 2, 4, 7, 11])
        d2 = compute_second_difference(X)
        assert len(d2) == len(X) - 2
        assert np.allclose(d2, [1, 1, 1])


class TestAutocorrelation:
    """Test lag-1 autocorrelation computation"""
    
    def test_white_noise_uncorrelated(self):
        """White noise should have rho1 ≈ 0"""
        np.random.seed(42)
        wn = np.random.randn(10000)
        rho1 = lag_1_autocorrelation(wn)
        assert abs(rho1) < 0.05


class TestPCF:
    """Test PCF metric computation"""
    
    def test_white_noise_floor(self):
        """White noise should have δ ≈ 0"""
        np.random.seed(42)
        wn = np.random.randn(5000)
        result = compute_pcf(wn)
        assert abs(result['delta']) < 0.05
        assert result['white_noise_floor'] == -2/3
    
    def test_structured_data_above_floor(self):
        """AR(1) process should have δ > 0"""
        np.random.seed(42)
        ar1 = np.zeros(5000)
        ar1[0] = np.random.randn()
        for t in range(1, len(ar1)):
            ar1[t] = 0.7 * ar1[t-1] + np.random.randn()
        result = compute_pcf(ar1)
        assert result['delta'] > 0
        assert result['above_floor'] is True
    
    def test_pcf_output_keys(self):
        """PCF output should contain all required keys"""
        X = np.random.randn(100)
        result = compute_pcf(X)
        required_keys = ['delta', 'rho1', 'ci_lower', 'ci_upper', 
                        'white_noise_floor', 'above_floor', 'second_differences']
        for key in required_keys:
            assert key in result


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
