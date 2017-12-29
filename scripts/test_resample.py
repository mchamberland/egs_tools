#!/home/mchamber/venvs/py-3.6/bin/python
import ct_tools.resample as ctr
import numpy as np

weights = ctr.fast_downscale(5, 2)
test = [1, 2, 3, 4, 5]
new = np.matmul(weights, test)
print('Done!')
