import { useMemo } from "react";
import "./App.css";
import { useMissionData } from "./useMissionData.js";

import TopBar from "./components/TopBar.jsx";
import MissionStatusPanel from "./components/MissionStatusPanel.jsx";
import LandingPanel from "./components/LandingPanel.jsx";
import CpasStatusPanel from "./components/CpasStatusPanel.jsx";
import GroundTrackPanel from "./components/GroundTrackPanel.jsx";
import TelemetryPlot from "./components/TelemetryPlot.jsx";
import Panel from "./components/Panel.jsx";
import CpasEventTable from "./components/CpasEventTable.jsx";
import TelemetryReadout from "./components/TelemetryReadout.jsx";
import RcsPanel from "./components/RcsPanel.jsx";
import Trajectory3D from "./components/Trajectory3D.jsx";

const PI = Math.PI;

/**
 * Downsample to a fixed max number of points by simple stride.
 * Keeps the first and last point. Good enough for visual sparklines.
 */
function downsample(trajectory, channelExtractor, maxPoints = 500) {
  if (!trajectory || trajectory.length === 0) return [];
  const stride = Math.max(1, Math.floor(trajectory.length / maxPoints));
  const out = [];
  for (let i = 0; i < trajectory.length; i += stride) {
    const r = trajectory[i];
    const y = channelExtractor(r);
    if (Number.isFinite(y)) out.push({ x: r.t_s, y });
  }
  // ensure the last point is included
  const lastR = trajectory[trajectory.length - 1];
  const lastY = channelExtractor(lastR);
  if (Number.isFinite(lastY) && out[out.length - 1]?.x !== lastR.t_s) {
    out.push({ x: lastR.t_s, y: lastY });
  }
  return out;
}

function App() {
  const data = useMissionData();

  // Pre-compute event markers shared by the plots. Only the main parachute
  // milestones -- reefing/disreef sub-stages cluster together and clutter the
  // chart. Map raw event names to short, readable labels.
  const eventMarkers = useMemo(() => {
    if (data.status !== "ready") return [];
    const majorLabels = {
      drogue_deploy: "Drogue",
      pilot_deploy: "Pilot",
      main_deploy: "Main",
      landed: "Landed",
    };
    return data.cpasEvents
      .filter((e) => e.event in majorLabels)
      .map((e) => ({ t: e.t_s, label: majorLabels[e.event] }));
  }, [data]);

  const altSeries  = useMemo(() => data.trajectory ? downsample(data.trajectory, r => (r.alt_m ?? 0) / 1000) : [], [data]);
  const velSeries  = useMemo(() => data.trajectory ? downsample(data.trajectory, r => r.V_mps) : [], [data]);
  const gammaSeries= useMemo(() => data.trajectory ? downsample(data.trajectory, r => (r.gamma_rad ?? 0) * 180 / PI) : [], [data]);
  const heatSeries = useMemo(() => data.trajectory ? downsample(data.trajectory, r => (r.nominal_heat_qdot_max_hi ?? 0) / 1e6) : [], [data]);

  if (data.status === "loading") {
    return (
      <div className="center">
        <div className="blink">LOADING TELEMETRY...</div>
        <div style={{ color: "var(--text-dim)", fontSize: 11 }}>
          /data/run_summary.json · trajectory.csv · cpas_events.csv · thruster_fires.csv
        </div>
      </div>
    );
  }

  if (data.status === "error") {
    return (
      <div className="center error">
        <div>TELEMETRY LINK FAILED</div>
        <div style={{ color: "var(--text-dim)", fontSize: 11, maxWidth: 600 }}>
          {data.error}
        </div>
      </div>
    );
  }

  const { runSummary, landingSummary, cpasEvents, trajectory, thrusterFires } = data;

  return (
    <div className="app">
      <TopBar runSummary={runSummary} landingSummary={landingSummary} />

      {/* Hero: live 3D trajectory playback */}
      <Trajectory3D
        trajectory={trajectory}
        cpasEvents={cpasEvents}
        runSummary={runSummary}
        landingSummary={landingSummary}
      />

      {/* Row 1: 3 status panels */}
      <div className="grid cols-3">
        <MissionStatusPanel runSummary={runSummary} />
        <LandingPanel landingSummary={landingSummary} runSummary={runSummary} />
        <CpasStatusPanel runSummary={runSummary} />
      </div>

      {/* Row 2: ground track + altitude plot */}
      <div className="grid cols-2">
        <GroundTrackPanel
          trajectory={trajectory}
          runSummary={runSummary}
          landingSummary={landingSummary}
        />
        <Panel
          title="Altitude vs Time"
          corner={`${trajectory.length} samples`}
        >
          <TelemetryPlot
            data={altSeries}
            events={eventMarkers}
            color="green"
            yLabel="ALT (km)"
            yFormat={(y) => y.toFixed(0)}
          />
        </Panel>
      </div>

      {/* Row 3: velocity, gamma, heat-rate plots */}
      <div className="grid cols-3">
        <Panel title="Velocity vs Time" corner="V_mps">
          <TelemetryPlot
            data={velSeries}
            events={eventMarkers}
            color="cyan"
            yLabel="V (m/s)"
            yFormat={(y) => y.toFixed(0)}
          />
        </Panel>
        <Panel title="Flight-Path Angle" corner="GAMMA">
          <TelemetryPlot
            data={gammaSeries}
            events={eventMarkers}
            color="amber"
            yLabel="γ (deg)"
            yFormat={(y) => y.toFixed(0)}
          />
        </Panel>
        <Panel title="Peak Heat Rate" corner="MW/m²">
          <TelemetryPlot
            data={heatSeries}
            events={eventMarkers}
            color="red"
            yLabel="q̇ (MW/m²)"
            yFormat={(y) => y.toFixed(1)}
          />
        </Panel>
      </div>

      {/* Row 4: CPAS event log + RCS table + final telemetry snapshot */}
      <div className="grid cols-3">
        <CpasEventTable cpasEvents={cpasEvents} />
        <RcsPanel thrusterFires={thrusterFires} />
        <TelemetryReadout trajectory={trajectory} />
      </div>
    </div>
  );
}

export default App;
