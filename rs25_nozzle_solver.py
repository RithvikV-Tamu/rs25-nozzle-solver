from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

PSI_TO_PA = 6894.757
LBF_TO_N = 4.4482216152605
INCH_TO_M = 0.0254
G0 = 9.80665
UNIVERSAL_GAS_CONSTANT = 8314.46

# Model inputs carried over from the original AERO 351 script.
GAMMA = 1.20
GAS_CONSTANT = 355.0
CHAMBER_PRESSURE = 2994.0 * PSI_TO_PA
CHAMBER_TEMPERATURE = 3600.0
EXIT_AREA_RATIO = 69.0
EXIT_DIAMETER = 2.30
NOZZLE_LENGTH = 3.07
NUM_POINTS = 400

# Public RS-25 values from the Wikipedia RS-25 page (see README for sources).
REFERENCE_VACUUM_THRUST = 512300.0 * LBF_TO_N
REFERENCE_VACUUM_ISP = 452.3
REFERENCE_MASS_FLOW = 514.49
REFERENCE_THROAT_DIAMETER = 10.3 * INCH_TO_M
REFERENCE_EXIT_DIAMETER = 90.7 * INCH_TO_M
REFERENCE_MIXTURE_RATIO = 6.03

FIGURE_FOLDER = Path(__file__).parent / "figures"


def area_mach_ratio(mach, gamma):
    term = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * mach**2)
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return term**exponent / mach


def solve_supersonic_mach(area_ratio, gamma, mach_low=1.000001, mach_high=10.0):
    if area_ratio < 1.0:
        raise ValueError("Area ratio A/A* must be at least 1.")
    if abs(area_ratio - 1.0) < 1.0e-10:
        return 1.0

    def residual(mach):
        return area_mach_ratio(mach, gamma) - area_ratio

    return brentq(residual, mach_low, mach_high)


def nozzle_geometry(exit_diameter, exit_area_ratio, nozzle_length, num_points):
    exit_area = np.pi * exit_diameter**2 / 4.0
    throat_area = exit_area / exit_area_ratio
    x = np.linspace(0.0, nozzle_length, num_points)
    xi = x / nozzle_length
    smooth_shape = 3.0 * xi**2 - 2.0 * xi**3
    area = throat_area * (1.0 + (exit_area_ratio - 1.0) * smooth_shape)
    return x, area, throat_area, exit_area


def choked_mass_flow(throat_area, chamber_pressure, chamber_temperature, gamma, gas_constant):
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return (
        throat_area * chamber_pressure / np.sqrt(chamber_temperature)
        * np.sqrt(gamma / gas_constant)
        * (2.0 / (gamma + 1.0))**exponent
    )


def solve_nozzle(
    gamma=GAMMA,
    gas_constant=GAS_CONSTANT,
    chamber_pressure=CHAMBER_PRESSURE,
    chamber_temperature=CHAMBER_TEMPERATURE,
    exit_area_ratio=EXIT_AREA_RATIO,
    exit_diameter=EXIT_DIAMETER,
    nozzle_length=NOZZLE_LENGTH,
    num_points=NUM_POINTS,
):
    x, area, throat_area, exit_area = nozzle_geometry(
        exit_diameter, exit_area_ratio, nozzle_length, num_points
    )
    area_ratio = area / throat_area
    mach = np.array([solve_supersonic_mach(ratio, gamma) for ratio in area_ratio])

    temperature_ratio = 1.0 + 0.5 * (gamma - 1.0) * mach**2
    temperature = chamber_temperature / temperature_ratio
    pressure = chamber_pressure / temperature_ratio**(gamma / (gamma - 1.0))
    density = pressure / (gas_constant * temperature)
    velocity = mach * np.sqrt(gamma * gas_constant * temperature)

    mass_flow = choked_mass_flow(
        throat_area, chamber_pressure, chamber_temperature, gamma, gas_constant
    )
    vacuum_thrust = mass_flow * velocity[-1] + pressure[-1] * exit_area
    vacuum_isp = vacuum_thrust / (mass_flow * G0)

    return {
        "x": x,
        "area": area,
        "radius": np.sqrt(area / np.pi),
        "area_ratio": area_ratio,
        "mach": mach,
        "temperature": temperature,
        "pressure": pressure,
        "density": density,
        "velocity": velocity,
        "throat_area": throat_area,
        "throat_diameter": np.sqrt(4.0 * throat_area / np.pi),
        "exit_area": exit_area,
        "mass_flow": mass_flow,
        "vacuum_thrust": vacuum_thrust,
        "vacuum_isp": vacuum_isp,
    }


def percent_difference(model_value, reference_value):
    return 100.0 * (model_value - reference_value) / reference_value


def complete_combustion_gas_constant(mixture_ratio):
    # Assumes 2 H2 + O2 -> 2 H2O goes to completion with the extra H2 left over.
    # No dissociation, so this is only a rough estimate for comparison.
    molar_mass_h2 = 2.016
    molar_mass_o2 = 31.998
    molar_mass_h2o = 18.015
    mass_o2 = mixture_ratio
    mass_h2 = 1.0
    moles_o2 = mass_o2 / molar_mass_o2
    moles_h2o = 2.0 * moles_o2
    moles_h2_left = mass_h2 / molar_mass_h2 - 2.0 * moles_o2
    total_mass = mass_o2 + mass_h2
    mean_molar_mass = total_mass / (moles_h2o + moles_h2_left)
    return UNIVERSAL_GAS_CONSTANT / mean_molar_mass, mean_molar_mass


def print_inputs():
    print("RS-25 simplified 1-D isentropic nozzle model")
    print("--------------------------------------------")
    print(f"gamma                    = {GAMMA:.3f}")
    print(f"Gas constant R           = {GAS_CONSTANT:.1f} J/(kg K)")
    print(f"Chamber pressure P0      = {CHAMBER_PRESSURE / 1e6:.3f} MPa")
    print(f"Chamber temperature T0   = {CHAMBER_TEMPERATURE:.1f} K")
    print(f"Exit area ratio Ae/At    = {EXIT_AREA_RATIO:.1f}")
    print(f"Exit diameter De         = {EXIT_DIAMETER:.3f} m")
    print(f"Nozzle length L          = {NOZZLE_LENGTH:.3f} m")
    print()


def print_exit_values(result):
    print("Computed exit-plane values (no public reference found for these)")
    print(f"  Exit Mach number       = {result['mach'][-1]:.4f}")
    print(f"  Exit pressure          = {result['pressure'][-1] / 1000:.3f} kPa")
    print(f"  Exit temperature       = {result['temperature'][-1]:.2f} K")
    print(f"  Exit density           = {result['density'][-1]:.5f} kg/m^3")
    print(f"  Exit velocity          = {result['velocity'][-1]:.2f} m/s")
    print()


def comparison_rows(result):
    return [
        ("Vacuum thrust (kN)", result["vacuum_thrust"] / 1000, REFERENCE_VACUUM_THRUST / 1000),
        ("Vacuum Isp (s)", result["vacuum_isp"], REFERENCE_VACUUM_ISP),
        ("Mass flow (kg/s)", result["mass_flow"], REFERENCE_MASS_FLOW),
        ("Throat diameter (m)", result["throat_diameter"], REFERENCE_THROAT_DIAMETER),
    ]


def print_comparison(result):
    print("Comparison with public RS-25 reference values")
    print(f"  {'Quantity':<22}{'Model':>12}{'Reference':>12}{'Diff (%)':>10}")
    for name, model_value, reference_value in comparison_rows(result):
        difference = percent_difference(model_value, reference_value)
        print(f"  {name:<22}{model_value:>12.4g}{reference_value:>12.4g}{difference:>+10.1f}")
    print()


def print_sensitivity():
    combustion_gas_constant, mean_molar_mass = complete_combustion_gas_constant(
        REFERENCE_MIXTURE_RATIO
    )
    published_area_ratio = (REFERENCE_EXIT_DIAMETER / REFERENCE_THROAT_DIAMETER)**2
    cases = [
        ("Baseline (original inputs)", {}),
        (f"R = {combustion_gas_constant:.1f} (complete combustion, O/F 6.03)",
         {"gas_constant": combustion_gas_constant}),
        (f"Ae/At = {published_area_ratio:.1f} (published diameters)",
         {"exit_area_ratio": published_area_ratio}),
    ]
    print("Sensitivity check: change one input at a time")
    print(f"  Complete-combustion molar mass = {mean_molar_mass:.2f} kg/kmol")
    print(f"  Original R = {GAS_CONSTANT:.1f} implies molar mass "
          f"{UNIVERSAL_GAS_CONSTANT / GAS_CONSTANT:.2f} kg/kmol")
    print(f"  {'Case':<46}{'Thrust %':>10}{'Isp %':>8}{'mdot %':>8}")
    for label, changes in cases:
        result = solve_nozzle(**changes)
        thrust_difference = percent_difference(result["vacuum_thrust"], REFERENCE_VACUUM_THRUST)
        isp_difference = percent_difference(result["vacuum_isp"], REFERENCE_VACUUM_ISP)
        mass_flow_difference = percent_difference(result["mass_flow"], REFERENCE_MASS_FLOW)
        print(f"  {label:<46}{thrust_difference:>+10.1f}{isp_difference:>+8.1f}"
              f"{mass_flow_difference:>+8.1f}")
    print()


def plot_nozzle_shape(result):
    x = result["x"]
    radius = result["radius"]
    plt.figure(figsize=(8, 4.5))
    plt.plot(x, radius, color="tab:blue", linewidth=2)
    plt.plot(x, -radius, color="tab:blue", linewidth=2)
    plt.fill_between(x, -radius, radius, color="tab:blue", alpha=0.15)
    plt.xlabel("Distance from throat, x (m)")
    plt.ylabel("Radius, r (m)")
    plt.title("Approximate RS-25 nozzle shape (throat to exit)")
    plt.axis("equal")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(FIGURE_FOLDER / "nozzle_shape.png", dpi=200)
    plt.close()


def plot_bracketing(result):
    sample_index = np.argmin(np.abs(result["x"] - 0.60 * NOZZLE_LENGTH))
    sample_area_ratio = result["area_ratio"][sample_index]
    sample_mach = result["mach"][sample_index]
    mach_values = np.linspace(1.0001, 8.0, 600)
    residual = area_mach_ratio(mach_values, GAMMA) - sample_area_ratio

    plt.figure(figsize=(7, 5))
    plt.plot(mach_values, residual, label=r"$f(M) = A/A^*(M) - A/A_t$")
    plt.axhline(0.0, color="black", linewidth=1)
    plt.axvline(1.0001, linestyle="--", color="gray", label="Lower bracket")
    plt.axvline(8.0, linestyle="--", color="gray", label="Upper bracket")
    plt.scatter([sample_mach], [0.0], color="tab:red", zorder=5,
                label=f"Root, M = {sample_mach:.3f}")
    plt.xlabel("Mach number, M")
    plt.ylabel("f(M)")
    plt.title(f"Bracketed root at x = 0.6 L (A/At = {sample_area_ratio:.1f})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURE_FOLDER / "bracketing_method.png", dpi=200)
    plt.close()


def plot_flow_properties(result):
    x = result["x"]
    figure, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    axes[0].plot(x, result["mach"], color="tab:blue")
    axes[0].set_ylabel("Mach number")
    axes[1].semilogy(x, result["pressure"] / 1000, color="tab:green")
    axes[1].set_ylabel("Pressure (kPa, log scale)")
    axes[2].plot(x, result["temperature"], color="tab:red")
    axes[2].set_ylabel("Temperature (K)")
    axes[2].set_xlabel("Distance from throat, x (m)")
    for axis in axes:
        axis.grid(True)
    axes[0].set_title("Mach number, pressure, and temperature along the nozzle")
    figure.tight_layout()
    figure.savefig(FIGURE_FOLDER / "flow_properties.png", dpi=200)
    plt.close(figure)


def plot_density_velocity(result):
    x = result["x"]
    figure, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    axes[0].plot(x, result["density"], color="tab:purple")
    axes[0].set_ylabel("Density (kg/m^3)")
    axes[1].plot(x, result["velocity"], color="tab:orange")
    axes[1].set_ylabel("Velocity (m/s)")
    axes[1].set_xlabel("Distance from throat, x (m)")
    for axis in axes:
        axis.grid(True)
    axes[0].set_title("Density and velocity along the nozzle")
    figure.tight_layout()
    figure.savefig(FIGURE_FOLDER / "density_velocity.png", dpi=200)
    plt.close(figure)


def main():
    result = solve_nozzle()
    print_inputs()
    print_exit_values(result)
    print_comparison(result)
    print_sensitivity()

    FIGURE_FOLDER.mkdir(exist_ok=True)
    plot_nozzle_shape(result)
    plot_bracketing(result)
    plot_flow_properties(result)
    plot_density_velocity(result)
    print(f"Figures saved to {FIGURE_FOLDER}")


if __name__ == "__main__":
    main()
