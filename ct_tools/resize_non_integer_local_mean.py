# Copyright 2018 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
from typing import Tuple, Iterable


def _reflect_breaks(size: int) -> np.ndarray:
    """Calculate cell boundaries with reflecting boundary conditions."""
    result = np.concatenate([[0], 0.5 + np.arange(size - 1), [size - 1]])
    assert len(result) == size + 1
    return result


def _interval_overlap(first_breaks: np.ndarray, second_breaks: np.ndarray) -> np.ndarray:
    """Return the overlap distance between all pairs of intervals.

    Args:
    first_breaks: breaks between entries in the first set of intervals, with
      shape (N+1,). Must be a non-decreasing sequence.
    second_breaks: breaks between entries in the second set of intervals, with
      shape (M+1,). Must be a non-decreasing sequence.

    Returns:
    Array with shape (N, M) giving the size of the overlapping region between
    each pair of intervals.
    """
    first_upper = first_breaks[1:]
    second_upper = second_breaks[1:]
    upper = np.minimum(first_upper[:, np.newaxis], second_upper[np.newaxis, :])

    first_lower = first_breaks[:-1]
    second_lower = second_breaks[:-1]
    lower = np.maximum(first_lower[:, np.newaxis], second_lower[np.newaxis, :])

    return np.maximum(upper - lower, 0)


def _resize_weights(old_size: int, new_size: int, reflect: bool = False) -> np.ndarray:
    """Create a weight matrix for resizing with the local mean along an axis.

    Args:
    old_size: old size.
    new_size: new size.
    reflect: whether or not there are reflecting boundary conditions.

    Returns:
    NumPy array with shape (new_size, old_size). Rows sum to 1.
    """
    if not reflect:
        old_breaks = np.linspace(0, old_size, num=old_size + 1)
        new_breaks = np.linspace(0, old_size, num=new_size + 1)
    else:
        old_breaks = _reflect_breaks(old_size)
        new_breaks = (old_size - 1) / (new_size - 1) * _reflect_breaks(new_size)

    weights = _interval_overlap(new_breaks, old_breaks)
    weights /= np.sum(weights, axis=1, keepdims=True)
    assert weights.shape == (new_size, old_size)
    return weights


def resize(array: np.ndarray, shape: Tuple[int, ...], reflect_axes: Iterable[int] = ()) -> np.ndarray:
    """Resize an array with the local mean / bilinear scaling.

    Works for both upsampling and downsampling in a fashion equivalent to
    block_mean and zoom, but allows for resizing by non-integer multiples. Prefer
    block_mean and zoom when possible, as this implementation is probably slower.

    Args:
    array: array to resize.
    shape: shape of the resized array.
    reflect_axes: iterable of axis numbers with reflecting boundary conditions,
      mirrored over the center of the first and last cell.

    Returns:
    Array resized to shape.

    Raises:
    ValueError: if any values in reflect_axes fall outside the interval
      [-array.ndim, array.ndim).
    """
    reflect_axes_set = set()
    for axis in reflect_axes:
        if not -array.ndim <= axis < array.ndim:
            raise ValueError('invalid axis: {}'.format(axis))
        reflect_axes_set.add(axis % array.ndim)

    output = array
    for axis, (old_size, new_size) in enumerate(zip(array.shape, shape)):
        reflect = axis in reflect_axes_set
        weights = _resize_weights(old_size, new_size, reflect=reflect)
        product = np.tensordot(output, weights, [[axis], [-1]])
        output = np.moveaxis(product, -1, axis)
    return output
