import { useMemo } from "react";

/**
 * Lightweight SVG sparkline for a single channel from the trajectory.
 *
 * Props:
 *   data      -- array of {x, y} numeric pairs (already downsampled)
 *   events    -- optional array of {t, label} markers
 *   color     -- "green" | "cyan" | "amber" | "red"
 *   yLabel    -- string
 *   yFormat   -- function(y) -> string for y-axis tick labels
 */
export default function TelemetryPlot({
  data,
  events = [],
  color = "green",
  yLabel = "",
  yFormat = (y) => y.toFixed(0),
}) {
  const layout = useMemo(() => {
    if (!data || data.length === 0) return null;
    const xs = data.map((d) => d.x);
    const ys = data.map((d) => d.y);
    const xMin = Math.min(...xs);
    const xMax = Math.max(...xs);
    const yMin = Math.min(...ys);
    const yMax = Math.max(...ys);
    return { xMin, xMax, yMin, yMax };
  }, [data]);

  if (!layout) return <div className="plot" />;

  const W = 600;
  const H = 220;
  const padL = 46, padR = 12, padT = 14, padB = 28;

  const { xMin, xMax, yMin, yMax } = layout;
  const xSpan = Math.max(1e-9, xMax - xMin);
  const ySpan = Math.max(1e-9, yMax - yMin);

  const sx = (x) => padL + ((x - xMin) / xSpan) * (W - padL - padR);
  const sy = (y) => padT + (1 - (y - yMin) / ySpan) * (H - padT - padB);

  const path = data
    .map((d, i) => `${i === 0 ? "M" : "L"}${sx(d.x).toFixed(1)},${sy(d.y).toFixed(1)}`)
    .join(" ");

  // 5 y-axis tick marks
  const yTicks = Array.from({ length: 5 }, (_, i) =>
    yMin + (i / 4) * (yMax - yMin)
  );
  const xTicks = Array.from({ length: 5 }, (_, i) =>
    xMin + (i / 4) * (xMax - xMin)
  );

  return (
    <div className="plot">
      <svg viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="none">
        {/* grid */}
        {yTicks.map((t, i) => (
          <line
            key={`gy-${i}`}
            className="grid-line"
            x1={padL} x2={W - padR}
            y1={sy(t)} y2={sy(t)}
          />
        ))}
        {xTicks.map((t, i) => (
          <line
            key={`gx-${i}`}
            className="grid-line"
            x1={sx(t)} x2={sx(t)}
            y1={padT} y2={H - padB}
          />
        ))}

        {/* axes */}
        <line className="axes" x1={padL} y1={padT} x2={padL} y2={H - padB} />
        <line className="axes" x1={padL} y1={H - padB} x2={W - padR} y2={H - padB} />

        {/* y tick labels */}
        {yTicks.map((t, i) => (
          <text
            key={`yt-${i}`}
            className="label"
            x={padL - 6} y={sy(t) + 3}
            textAnchor="end"
          >
            {yFormat(t)}
          </text>
        ))}

        {/* x tick labels */}
        {xTicks.map((t, i) => (
          <text
            key={`xt-${i}`}
            className="label"
            x={sx(t)} y={H - padB + 14}
            textAnchor="middle"
          >
            {t.toFixed(0)}s
          </text>
        ))}

        {/* events */}
        {events.map((e, i) => (
          <g key={`ev-${i}`}>
            <line
              className="event-line"
              x1={sx(e.t)} x2={sx(e.t)}
              y1={padT} y2={H - padB}
            />
            <text
              className="event-label"
              x={sx(e.t) + 2}
              y={padT + 8 + (i % 3) * 9}
            >
              {e.label}
            </text>
          </g>
        ))}

        {/* y axis label */}
        <text
          className="axis-label"
          x={10} y={padT + (H - padT - padB) / 2}
          transform={`rotate(-90, 10, ${padT + (H - padT - padB) / 2})`}
          textAnchor="middle"
        >
          {yLabel}
        </text>

        {/* trace */}
        <path className={`trace ${color}`} d={path} />
      </svg>
    </div>
  );
}
