import Panel, { Pill } from "./Panel.jsx";
import { fmt } from "../format.js";

export default function MissionStatusPanel({ runSummary }) {
  const ok = runSummary.termination_reason === "nominal_ground_reached";
  return (
    <Panel
      title="Mission Status"
      corner="SYS-01"
      status={
        <Pill kind={ok ? "nominal" : "alert"}>
          {ok ? "GROUND" : "ABORT"}
        </Pill>
      }
    >
      <div className="readouts">
        <div className="readout">
          <div className="label">Termination</div>
          <div className="value cyan" style={{ fontSize: 13 }}>
            {runSummary.termination_reason}
          </div>
        </div>
        <div className="readout">
          <div className="label">Final MET</div>
          <div className="value">
            {fmt.met(runSummary.final_t_s)}
          </div>
        </div>
        <div className="readout">
          <div className="label">Final Altitude</div>
          <div className="value amber">
            {fmt.alt(runSummary.final_alt_m)}
          </div>
        </div>
        <div className="readout">
          <div className="label">Final Velocity</div>
          <div className={`value ${runSummary.final_V_mps <= 12 ? "" : "red"}`}>
            {fmt.vel(runSummary.final_V_mps)}
          </div>
        </div>
      </div>

      <div className="kvgrid">
        <span className="k">Steps logged</span>
        <span className="v">{fmt.int(runSummary.num_logged_steps)}</span>
        <span className="k">Guidance updates</span>
        <span className="v">{fmt.int(runSummary.guidance_updates)}</span>
        <span className="k">Thruster fires</span>
        <span className="v">{fmt.int(runSummary.thruster_fire_events)}</span>
        <span className="k">Interval recenters</span>
        <span className="v amber">
          {fmt.int(runSummary.interval_recenter_count)}
        </span>
        <span className="k">Interval splits</span>
        <span className="v">{fmt.int(runSummary.interval_split_count)}</span>
        <span className="k">Final gamma</span>
        <span className="v">{fmt.deg(runSummary.final_gamma_deg, 2)}</span>
      </div>
    </Panel>
  );
}
