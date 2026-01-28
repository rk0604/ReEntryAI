'''
{
  "t": 183.2,
  "triangles": [
    {"id": 0, "qdot": 615000, "Q": 2.41e8},
    {"id": 1, "qdot": 602000, "Q": 2.35e8},
    {"id": 2, "qdot": 540000, "Q": 2.00e8}
    // ...
  ]
}  

Better because it's a fixed-length, control-relevant summary of the full triangle mesh (and stochastic atmosphere), 
so model can learn stable policies without drowning in hundreds of weakly-informative per-cell values
RL observation (recommended): [h, v, gamma, s, delta_rho,
                               qdot_max, Q_max, qdot_mean,
                               frac_over_qdot_limit, frac_over_Q_limit,
                               center_ring_qdot_mean, shoulder_ring_qdot_mean, nonuniformity]

'''