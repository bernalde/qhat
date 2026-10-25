"""Altered-target scaling preserves supplied optional metadata without invention."""
from copy import deepcopy
from fractions import Fraction
from study import scaled_preflight_model


def model():
    term=[[[0,'Z']],'2']
    return dict(n=1,electrons=1,case='metadata_fixture',denominator='5',coefficients=[term],blocks=[[term]])


def test_absent_scalar_metadata_stays_absent():
    supplied=model();before=deepcopy(supplied);altered=scaled_preflight_model(supplied)
    assert supplied==before and 'removed_scalar_exact' not in altered
    assert Fraction(int(altered['coefficients'][0][1]),int(altered['denominator']))==Fraction(2,5)*Fraction(7,5)
    assert altered['blocks'][0][0]==altered['coefficients'][0]


def test_supplied_scalar_metadata_is_scaled_exactly():
    supplied=model();supplied['removed_scalar_exact']='3/7';altered=scaled_preflight_model(supplied)
    assert altered['removed_scalar_exact']=='3/5' and supplied['removed_scalar_exact']=='3/7'
