import Panel, { Pill } from "./Panel.jsx";
import { fmt } from "../format.js";

export default function LandingPanel({ landingSummary, runSummary }) {
  const miss = landingSummary.great_circle_error_km;
  const onTarget = miss <= 5;
  const closeMiss = miss <= 50;

  let kind, txt;
  if (onTarget) { kind = "nominal"; txt = "ON TARGET"; }
  else if (closeMiss) { kind = "warn"; txt = "OFF TARGET"; }
  else { kind = "alert"; txt = "DIVERGED"; }

  return (
    <Panel title="Landing Solution" corner="SYS-02"
           status={<Pill kind={kind}>{txt}</Pill>}>
      <div className="readouts">
        <div className="readout">
          <div className="label">Great-Circle Miss</div>
          <div className={`value ${onTarget ? "" : closeMiss ? "amber" : "red"}`}>
            {fmt.km(miss)}
          </div>
        </div>
        <div className="readout">
          <div className="label">Flat-Ground Miss</div>
          <div className="value cyan">
            {fmt.km(landingSummary.flat_ground_error_km)}
          </div>
        </div>
        <div className="readout">
          <div className="label">North Error</div>
          <div className="value amber">
            {fmt.km((landingSummary.north_error_km ?? 0))}
          </div>
        </div>
        <div className="readout">
          <div className="label">East Error</div>
          <div className="value amber">
            {fmt.km(landingSummary.east_error_km)}
          </div>
        </div>
      </div>

      <div className="kvgrid">
        <span className="k">Target latitude</span>
        <span className="v">{fmt.deg(runSummary.target_phi_deg, 3)}</span>
        <span className="k">Target longitude</span>
        <span className="v">{fmt.deg(runSummary.target_lam_deg, 3)}</span>
        <span className="k">Final latitude</span>
        <span className="v">{fmt.deg(landingSummary.final_phi_deg, 3)}</span>
        <span className="k">Final longitude</span>
        <span className="v">{fmt.deg(landingSummary.final_lam_deg, 3)}</span>
        <span className="k">Final heading</span>
        <span className="v">{fmt.deg(landingSummary.final_chi_deg, 2)}</span>
        <span className="k">RCS total fires</span>
        <span className="v amber">{fmt.int(landingSummary.rcs_total_fires)}</span>
      </div>
    </Panel>
  );
}
