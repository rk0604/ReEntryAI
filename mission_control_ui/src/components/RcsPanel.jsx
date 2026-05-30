import { useMemo } from "react";
import Panel, { Pill } from "./Panel.jsx";
import { fmt } from "../format.js";

/**
 * Aggregate the thruster_fires.csv into per-thruster counts and total firings.
 */
export default function RcsPanel({ thrusterFires }) {
  const stats = useMemo(() => {
    const counts = {};
    for (const r of thrusterFires) {
      counts[r.thruster] = (counts[r.thruster] || 0) + 1;
    }
    const entries = Object.entries(counts).sort((a, b) => b[1] - a[1]);
    return {
      total: thrusterFires.length,
      entries,
      maxCount: entries.length ? entries[0][1] : 1,
    };
  }, [thrusterFires]);

  return (
    <Panel
      title="RCS Thruster Activity"
      corner={`${fmt.int(stats.total)} fires`}
      status={<Pill kind="info">8 PODS</Pill>}
    >
      <table className="table">
        <thead>
          <tr>
            <th>Thruster</th>
            <th>Fires</th>
            <th>Utilization</th>
          </tr>
        </thead>
        <tbody>
          {stats.entries.map(([name, count]) => {
            const pct = (count / stats.maxCount) * 100;
            const isPos = name.endsWith("RPOS");
            return (
              <tr key={name}>
                <td>{name}</td>
                <td className="right">{fmt.int(count)}</td>
                <td style={{ width: "55%" }}>
                  <div
                    style={{
                      height: 6,
                      background: isPos ? "var(--cyan)" : "var(--amber)",
                      width: `${pct}%`,
                      boxShadow: `0 0 6px ${
                        isPos ? "var(--cyan-soft)" : "var(--amber-soft)"
                      }`,
                      borderRadius: 1,
                    }}
                  />
                </td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </Panel>
  );
}
