import cProfile
import galacticops as go
import random
ds, zs = go.dist_lookup(cosmology=True, H_0=69.6, W_m=0.286, W_v=0.714, z_max=8.0)

rd = [random.random() for d in range(100000)]

def test():
    for d in rd:
        go.interpolate_z(d, ds, zs)

test()
