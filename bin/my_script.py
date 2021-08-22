#!/usr/bin/env python
import numpy as np
for data in ["sc", "vol", "len", "sc_edge_normalized", "sc_vol_normalized"]:
curr_data = np.load(data + ".npy")
np.savetxt(data + ".csv", curr_data, delimiter=",")
