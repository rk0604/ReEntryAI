import Panel, { Pill } from "./Panel.jsx";
import { fmt } from "../format.js";

export default function CpasStatusPanel({ runSummary }) {
  const c = runSummary.cpas_summary || {};
  const operational = c.mains_operational_mask || [false, false, false];
  const nMainsOk = operational.filter(Boolean).length;

  let kind, txt;
  if (c.phase === "landed") { kind = "nominal"; txt = "DESCENT COMPLETE"; }
  else if (c.phase === "main") { kind = "info"; txt = "UNDER MAINS"; }
  else if (c.phase === "drogue" || c.phase === "pilot") { kind = "info"; txt = c.phase.toUpperCase(); }
  else if (c.phase === "stowed") { kind = "warn"; txt = "STOWED"; }
  else { kind = "warn"; txt = (c.phase || "UNKNOWN").toUpperCase(); }

  return (
    <Panel title="CPAS Parachute System" corner="SYS-03"
           status={<Pill kind={kind}>{txt}</Pill>}>
      <div className="readouts">
        <div className="readout">
          <div className="label">FBC Jettison</div>
          <div className={`value ${c.fbc_jettisoned ? "" : "amber"}`} style={{ fontSize: 16 }}>
            {c.fbc_jettisoned ? `T+${fmt.t(c.fbc_jettison_t_s, 1)}` : "ARMED"}
          </div>
        </div>
        <div className="readout">
          <div className="label">Mains Operational</div>
          <div className={`value ${nMainsOk === 3 ? "" : "red"}`}>
            {nMainsOk} / 3
          </div>
        </div>
        <div className="readout">
          <div className="label">Mass Shed</div>
          <div className="value amber">
            {fmt.num(c.mass_shed_kg_total, 0)} kg
          </div>
        </div>
        <div className="readout">
          <div className="label">Pendulum Azimuth</div>
          <div className="value cyan">
            {fmt.deg(c.pendulum_azimuth_deg, 1)}
          </div>
        </div>
      </div>

      <div className="kvgrid">
        <span className="k">Drogue deploy</span>
        <span className="v">{fmt.t(c.drogue_deploy_t_s)}</span>
        <span className="k">Pilot deploy</span>
        <span className="v">{fmt.t(c.pilot_deploy_t_s)}</span>
        <span className="k">Main deploy</span>
        <span className="v">{fmt.t(c.main_deploy_t_s)}</span>
        <span className="k">Landed</span>
        <span className="v green">{fmt.t(c.landed_t_s)}</span>
      </div>
    </Panel>
  );
}
