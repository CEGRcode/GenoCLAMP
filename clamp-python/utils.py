import numpy as np
from typing import Union

def boltzmann(arr: np.ndarray, alpha: float, axis: Union[None, int, tuple] = None):
    return np.sum(arr * np.exp(alpha * arr), axis=axis) / np.sum(np.exp(alpha * arr), axis=axis)
    