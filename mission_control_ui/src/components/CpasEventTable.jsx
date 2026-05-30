import Panel from "./Panel.jsx";
import { classifyCpasEvent, fmt } from "../format.js";

export default function CpasEventTable({ cpasEvents }) {
  return (
    <Panel title="CPAS Event Log" corner={`${cpasEvents.length} events`}>
      <div className="scroll">
        <table className="table">
          <thead>
            <tr>
              <th>T+ (s)</th>
              <th>Step</th>
              <th>Phase</th>
              <th>Event</th>
              <th className="right">Alt</th>
              <th className="right">V</th>
              <th className="right">Reef</th>
              <th className="right">Mains</th>
            </tr>
          </thead>
          <tbody>
            {cpasEvents.map((e, i) => (
              <tr key={i}>
                <td>{fmt.num(e.t_s, 1)}</td>
                <td className="dim">{e.step}</td>
                <td className="dim">{e.phase}</td>
                <td className={classifyCpasEvent(e.event)}>{e.event}</td>
                <td className="right">{fmt.alt(e.alt_m)}</td>
                <td className="right">{fmt.vel(e.V_mps)}</td>
                <td className="right dim">{e.reefing_stage}</td>
                <td className="right dim">{e.num_mains_active}/3</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </Panel>
  );
}
