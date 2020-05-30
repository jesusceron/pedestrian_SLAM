import numpy as np
import scipy.stats
from filterpy.kalman import KalmanFilter
from filterpy.monte_carlo import stratified_resample
from filterpy.monte_carlo import systematic_resample


def tracker1(x_initial, R_std, Q_std):
    tracker = KalmanFilter(dim_x=1, dim_z=1)

    tracker.F = np.array([1])
    tracker.u = 0.
    tracker.H = np.array([1])

    tracker.R = R_std
    tracker.Q = Q_std
    tracker.x = np.array([x_initial]).T
    tracker.P = 50

    return tracker


def rssi_filter(data, R_std, Q_std):

    x_initial = data[0]

    # run filter
    robot_tracker = tracker1(x_initial, R_std, Q_std)
    mu, cov, _, _ = robot_tracker.batch_filter(data)

    res_nan = np.isnan(mu)
    mu = mu[~res_nan]
    mu = mu.tolist()
    
    return mu

def calculate_weight(particle_weight, dist, covV, distance_rssi):

    particle_weight *= (scipy.stats.norm(dist, covV).pdf(distance_rssi) + 1.e-300)

    return particle_weight
