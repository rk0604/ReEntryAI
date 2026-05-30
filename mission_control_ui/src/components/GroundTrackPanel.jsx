import { useMemo } from "react";
import Panel from "./Panel.jsx";
import { fmt } from "../format.js";

/**
 * 2-D ground track in lat/lon space. Plots:
 *   - Start (cyan square)
 *   - Trajectory (green line, downsampled)
 *   - End (amber circle)
 *   - Target (red X)
 */
export default function GroundTrackPanel({ trajectory, runSummary, landingSummary }) {
  const path = useMemo(() => {
    if (!trajectory || trajectory.length === 0) return null;

    const stride = Math.max(1, Math.floor(trajectory.length / 600));
    const pts = [];
    for (let i = 0; i < trajectory.length; i += stride) {
      const r = trajectory[i];
      pts.push({
        lon: (r.lam_rad ?? 0) * (180 / Math.PI),
        lat: (r.phi_rad ?? 0) * (180 / Math.PI),
      });
    }

    const targetLat = runSummary.target_phi_deg;
    const targetLon = runSummary.target_lam_deg;
    const startLat = pts[0].lat, startLon = pts[0].lon;
    const endLat = landingSummary.final_phi_deg;
    const endLon = landingSummary.final_lam_deg;

    const allLats = [...pts.map((p) => p.lat), targetLat, startLat, endLat];
    const allLons = [...pts.map((p) => p.lon), targetLon, startLon, endLon];

    return {
      pts,
      latMin: Math.min(...allLats),
      latMax: Math.max(...allLats),
      lonMin: Math.min(...allLons),
      lonMax: Math.max(...allLons),
      startLat, startLon,
      endLat, endLon,
      targetLat, targetLon,
    };
  }, [trajectory, runSummary, landingSummary]);

  if (!path) {
    return (
      <Panel title="Ground Track" corner="MAP">
        <div className="plot" />
      </Panel>
    );
  }

  const W = 600, H = 360;
  const padL = 40, padR = 12, padT = 12, padB = 30;
  const lonSpan = Math.max(1e-6, path.lonMax - path.lonMin);
  const latSpan = Math.max(1e-6, path.latMax - path.latMin);
  // pad bounds by 10%
  const lonMin = path.lonMin - lonSpan * 0.08;
  const lonMax = path.lonMax + lonSpan * 0.08;
  const latMin = path.latMin - latSpan * 0.08;
  const latMax = path.latMax + latSpan * 0.08;

  const sx = (lon) => padL + ((lon - lonMin) / (lonMax - lonMin)) * (W - padL - padR);
  const sy = (lat) => padT + (1 - (lat - latMin) / (latMax - latMin)) * (H - padT - padB);

  const d = path.pts
    .map((p, i) => `${i === 0 ? "M" : "L"}${sx(p.lon).toFixed(1)},${sy(p.lat).toFixed(1)}`)
    .join(" ");

  const missKm = landingSummary.great_circle_error_km;
  const onTgt = missKm <= 5;

  return (
    <Panel
      title="Ground Track"
      corner={`Miss ${fmt.km(missKm)}`}
    >
      <div className="plot" style={{ minHeight: 360 }}>
        <svg viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="xMidYMid meet">
          {/* axes */}
          <line className="axes" x1={padL} y1={padT} x2={padL} y2={H - padB} />
          <line className="axes" x1={padL} y1={H - padB} x2={W - padR} y2={H - padB} />
          {/* gridlines (5x5) */}
          {Array.from({ length: 5 }, (_, i) => {
            const lat = latMin + (i / 4) * (latMax - latMin);
            return (
              <g key={`gy${i}`}>
                <line
                  className="grid-line"
                  x1={padL} x2={W - padR}
                  y1={sy(lat)} y2={sy(lat)}
                />
                <text className="label" x={padL - 4} y={sy(lat) + 3} textAnchor="end">
                  {lat.toFixed(1)}°
                </text>
              </g>
            );
          })}
          {Array.from({ length: 5 }, (_, i) => {
            const lon = lonMin + (i / 4) * (lonMax - lonMin);
            return (
              <g key={`gx${i}`}>
                <line
                  className="grid-line"
                  x1={sx(lon)} x2={sx(lon)}
                  y1={padT} y2={H - padB}
                />
                <text className="label" x={sx(lon)} y={H - padB + 14} textAnchor="middle">
                  {lon.toFixed(1)}°
                </text>
              </g>
            );
          })}

          {/* trajectory */}
          <path className="trace green" d={d} />

          {/* target */}
          <g transform={`translate(${sx(path.targetLon)}, ${sy(path.targetLat)})`}>
            <line x1={-8} y1={-8} x2={8} y2={8} stroke="var(--red)" strokeWidth="2" />
            <line x1={-8} y1={8} x2={8} y2={-8} stroke="var(--red)" strokeWidth="2" />
            <text x={12} y={4} className="label" fill="var(--red)">TARGET</text>
          </g>

          {/* start */}
          <g transform={`translate(${sx(path.startLon)}, ${sy(path.startLat)})`}>
            <rect x={-5} y={-5} width="10" height="10" fill="var(--cyan)" />
            <text x={9} y={4} className="label" fill="var(--cyan)">EI</text>
          </g>

          {/* end */}
          <g transform={`translate(${sx(path.endLon)}, ${sy(path.endLat)})`}>
            <circle r={6} fill={onTgt ? "var(--green)" : "var(--amber)"} />
            <text x={10} y={4} className="label" fill={onTgt ? "var(--green)" : "var(--amber)"}>
              SPLASH
            </text>
          </g>

          {/* axis labels */}
          <text x={W / 2} y={H - 4} className="axis-label" textAnchor="middle">
            LONGITUDE
          </text>
          <text
            x={14} y={padT + (H - padT - padB) / 2}
            className="axis-label"
            textAnchor="middle"
            transform={`rotate(-90, 14, ${padT + (H - padT - padB) / 2})`}
          >
            LATITUDE
          </text>
        </svg>
      </div>
    </Panel>
  );
}
