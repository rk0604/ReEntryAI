import AtmosphereModel

def verify_density():
    test_altitudes_km = [0, 50, 80, 120, 400, 900, 1000]

    print(f"{'Alt (km)':<10} | {'PyMSIS (kg/m3)':<20}")
    print("-" * 33)

    for h_km in test_altitudes_km:
        h_m = h_km * 1000.0
        try:
            rho_py = AtmosphereModel.get_atmosphere_density(h_m)
            rho_py_str = f"{rho_py:.5e}" if rho_py is not None else "None"
        except Exception as e:
            rho_py_str = f"Error: {e}"

        print(f"{h_km:<10} | {rho_py_str:<20}")

if __name__ == "__main__":
    verify_density()
