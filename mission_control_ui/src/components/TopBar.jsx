import { Pill } from "./Panel.jsx";
import { fmt } from "../format.js";

export default function TopBar({ runSummary, landingSummary }) {
  const ok = runSummary.termination_reason === "nominal_ground_reached";
  const intervalFailed = !!runSummary.interval_failed_detected;

  let statusKind, statusText;
  if (!ok) {
    statusKind = "alert";
    statusText = "ABORT";
  } else if (intervalFailed) {
    statusKind = "warn";
    statusText = "DEGRADED";
  } else {
    statusKind = "nominal";
    statusText = "NOMINAL";
  }

  return (
    <div className="topbar">
      <div className="brand">
        REENTRY AI
        <span className="sub">
          MISSION CONTROL // {runSummary.trajectory_id}
        </span>
      </div>
      <div className="spacer" />

      <Pill kind={statusKind}>{statusText}</Pill>

      <div className="stat">
        <span className="label">MET</span>
        <span className="value">{fmt.met(runSummary.final_t_s)}</span>
      </div>
      <div className="stat">
        <span className="label">AERO MODEL</span>
        <span className="value" style={{ fontSize: 12 }}>
          {runSummary.aero_model}
        </span>
      </div>
      <div className="stat">
        <span className="label">MISS</span>
        <span className="value" style={{
          color: landingSummary.great_circle_error_km > 100
            ? "var(--red)"
            : "var(--green)",
        }}>
          {fmt.km(landingSummary.great_circle_error_km)}
        </span>
      </div>
    </div>
  );
}
