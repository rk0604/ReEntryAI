# main entry to the sim
import numpy as np
import re
import constants
import math
import StochasticEntrySim.math_3d as math_3d
import matplotlib.pyplot as plt


def parse_input_file(input_path: str):
    """
    Parse a simple key = value input file into a dictionary.
    Lines starting with # or blank lines are ignored.
    Numeric values (including negatives and scientific notation)
    are converted to float; everything else is kept as string.
    """
    state_dict = {}

    with open(input_path, "r") as file:
        input_list = file.readlines()

    extract_pattern = r"^\s*([A-Za-z_]\w*)\s*=\s*([^#\n]+)"
    num_pattern = re.compile(r"^[+-]?\d+(\.\d+)?([eE][+-]?\d+)?$")

    for line in input_list:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        m = re.match(extract_pattern, line)
        if not m:
            continue

        key = m.group(1)
        raw_value = m.group(2).strip()

        if num_pattern.fullmatch(raw_value):
            try:
                value = float(raw_value)
                state_dict[key] = value
            except ValueError:
                state_dict[key] = raw_value
        else:
            state_dict[key] = raw_value

    return state_dict


def constant_bank_sigma_fn_builder(state_dict):
    sigma_const_deg = state_dict.get("sigma_const_deg", 0.0)
    sigma_const_rad = float(sigma_const_deg) * math.pi / 180.0

    def sigma_fn(t, x):
        return sigma_const_rad

    return sigma_fn


def main():
    # 1. Load input
    input_path = "input.txt"
    state_dict = parse_input_file(input_path)

    # 2. Vehicle and aero params
    params = {
        "m": float(state_dict.get("mass_kg", 2000.0)),
        "S": float(state_dict.get("ref_area_m2", 15.0)),
        "CD": float(state_dict.get("CD_const", 1.0)),
        "CL": float(state_dict.get("CL_const", 0.3)),
    }

    # 3. Initial state (x = [r, phi, lam, V, gamma, chi])
    h0 = float(state_dict.get("init_altitude_m", 120000.0))
    r0 = constants.RADIUS_EARTH + h0

    phi0_deg = float(state_dict.get("init_lat_deg", 0.0))
    lam0_deg = float(state_dict.get("init_lon_deg", 0.0))
    V0 = float(state_dict.get("init_speed_mps", 7700.0))
    gamma0_deg = float(state_dict.get("init_fpa_deg", -5.0))
    chi0_deg = float(state_dict.get("init_heading_deg", 90.0))

    phi0 = phi0_deg * math.pi / 180.0
    lam0 = lam0_deg * math.pi / 180.0
    gamma0 = gamma0_deg * math.pi / 180.0
    chi0 = chi0_deg * math.pi / 180.0

    x0 = [r0, phi0, lam0, V0, gamma0, chi0]

    # 4. Time settings
    t0 = float(state_dict.get("t0_s", 0.0))
    t_final = float(state_dict.get("t_final_s", 2000.0))
    dt = float(state_dict.get("dt_s", 0.25))

    sigma_fn = constant_bank_sigma_fn_builder(state_dict)

    # stopping conditions
    min_altitude = float(state_dict.get("min_altitude_m", 0.0))
    min_speed = float(state_dict.get("min_speed_mps", 10.0))

    # 5. Integrate with Euler
    t = t0
    x = x0[:]
    history = []

    while t <= t_final:
        h = math_3d.altitude(x[0])
        V = x[3]

        if h <= min_altitude or V <= min_speed:
            break

        history.append((t, *x))
        dx = math_3d.eom_3d(t, x, params, sigma_fn)

        # force all derivatives to be floats, and update
        dx = [float(d) for d in dx]
        x = [xi + dt * dxi for xi, dxi in zip(x, dx)]
        t += dt

    if len(history) == 0:
        print("No trajectory points recorded. Check initial conditions or dt.")
        return

    # 6. Plot
    history = np.array(history)
    t_arr = history[:, 0]
    r_arr = history[:, 1]
    V_arr = history[:, 4]

    h_arr = (r_arr - constants.RADIUS_EARTH) / 1000.0

    plt.figure()
    plt.plot(t_arr, h_arr)
    plt.xlabel("Time (s)")
    plt.ylabel("Altitude (km)")
    plt.title("Reentry Trajectory: Altitude vs Time")
    plt.grid(True)

    plt.figure()
    plt.plot(t_arr, V_arr)
    plt.xlabel("Time (s)")
    plt.ylabel("Speed (m/s)")
    plt.title("Reentry Trajectory: Speed vs Time")
    plt.grid(True)

    plt.show()


if __name__ == "__main__":
    main()
