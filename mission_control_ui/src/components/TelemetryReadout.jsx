import Panel from "./Panel.jsx";
import { fmt } from "../format.js";

/**
 * Snapshot of the FINAL logged trajectory row — the touchdown state.
 * Shows every channel a flight controller would scan in last seconds.
 */
export default function TelemetryReadout({ trajectory }) {
  const last = trajectory[trajectory.length - 1] || {};
  const altKm = (last.alt_m ?? 0) / 1000;
  const dynP = last.q_pa ?? 0;
  const heat = last.nominal_heat_qdot_max_hi ?? 0;

  return (
    <Panel title="Final Telemetry Snapshot" corner="T-FINAL">
      <div className="kvgrid">
        <span className="k">T (s)</span>
        <span className="v">{fmt.num(last.t_s, 2)}</span>

        <span className="k">Altitude</span>
        <span className="v amber">
          {altKm < 0.001 ? fmt.m(last.alt_m, 2) : fmt.alt(last.alt_m)}
        </span>

        <span className="k">Velocity</span>
        <span className="v green">{fmt.vel(last.V_mps)}</span>

        <span className="k">Gamma</span>
        <span className="v">{fmt.deg((last.gamma_rad ?? 0) * (180 / Math.PI), 2)}</span>

        <span className="k">Heading (chi)</span>
        <span className="v">{fmt.deg((last.chi_rad ?? 0) * (180 / Math.PI), 2)}</span>

        <span className="k">Latitude</span>
        <span className="v cyan">{fmt.deg((last.phi_rad ?? 0) * (180 / Math.PI), 4)}</span>

        <span className="k">Longitude</span>
        <span className="v cyan">{fmt.deg((last.lam_rad ?? 0) * (180 / Math.PI), 4)}</span>

        <span className="k">Density</span>
        <span className="v">{fmt.num(last.rho_kgm3, 4)} kg/m³</span>

        <span className="k">Dyn. pressure</span>
        <span className="v">{fmt.num(dynP, 1)} Pa</span>

        <span className="k">CD / CL / L:D</span>
        <span className="v">
          {fmt.num(last.CD, 3)} / {fmt.num(last.CL, 3)} / {fmt.num(last.LD, 3)}
        </span>

        <span className="k">Heat rate (peak)</span>
        <span className={`v ${heat / 1e6 > 12 ? "amber" : ""}`}>
          {fmt.num(heat / 1e6, 2)} MW/m²
        </span>

        <span className="k">Range to go</span>
        <span className="v amber">
          {fmt.km((last.range_to_go_m ?? 0) / 1000)}
        </span>

        <span className="k">Sigma actual</span>
        <span className="v">{fmt.deg((last.sigma_actual_rad ?? 0) * (180 / Math.PI), 2)}</span>

        <span className="k">Mass</span>
        <span className="v">{fmt.num(last.cpas_mass_kg_now, 0)} kg</span>

        <span className="k">CPAS phase</span>
        <span className="v green">{last.cpas_phase}</span>
      </div>
    </Panel>
  );
}
