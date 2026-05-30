import { useMemo, useRef, useState, useEffect } from "react";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls, Stars, Line } from "@react-three/drei";
import * as THREE from "three";

import Panel, { Pill } from "./Panel.jsx";
import { fmt } from "../format.js";

/**
 * Convert geodetic (lat, lon, alt_m) to a 3D point in km, with the
 * Earth centred at the origin and Z = north pole. Same convention the
 * sim itself uses.
 *
 * altScale lets the visualization exaggerate altitude only (Earth's
 * radius is unchanged), so the entry corridor is visible at to-scale
 * camera distances. altScale = 1 is true to scale.
 */
function geodeticToVec3(phi_rad, lam_rad, alt_m, R_earth_km = 6378.1, altScale = 1) {
  const r = R_earth_km + (alt_m / 1000.0) * altScale;
  const x = r * Math.cos(phi_rad) * Math.cos(lam_rad);
  const y = r * Math.cos(phi_rad) * Math.sin(lam_rad);
  const z = r * Math.sin(phi_rad);
  return new THREE.Vector3(x, y, z);
}

// -----------------------------------------------------------------------------
// Earth: dark blue sphere with a wireframe overlay
// -----------------------------------------------------------------------------
function Earth({ radius = 6378.1 }) {
  return (
    <group>
      {/* solid sphere */}
      <mesh>
        <sphereGeometry args={[radius, 64, 64]} />
        <meshStandardMaterial
          color="#0a1e3a"
          emissive="#031024"
          roughness={0.9}
          metalness={0.1}
        />
      </mesh>
      {/* wireframe overlay slightly above the sphere */}
      <mesh>
        <sphereGeometry args={[radius * 1.001, 48, 32]} />
        <meshBasicMaterial
          color="#1e4a7a"
          wireframe={true}
          transparent={true}
          opacity={0.35}
        />
      </mesh>
      {/* equator highlight ring */}
      <mesh rotation={[Math.PI / 2, 0, 0]}>
        <torusGeometry args={[radius * 1.002, 4, 8, 96]} />
        <meshBasicMaterial color="#5ddbe0" transparent opacity={0.4} />
      </mesh>
    </group>
  );
}

// -----------------------------------------------------------------------------
// Trajectory: the full path as a glowing line
// -----------------------------------------------------------------------------
function TrajectoryPath({ points, currentIdx }) {
  // Past segment (already flown) — bright green
  const flown = useMemo(
    () => points.slice(0, Math.max(2, currentIdx + 1)),
    [points, currentIdx]
  );
  // Future segment — dimmer
  const remaining = useMemo(
    () => points.slice(Math.max(0, currentIdx)),
    [points, currentIdx]
  );

  return (
    <>
      <Line points={flown} color="#46d160" lineWidth={2.5} transparent opacity={0.95} />
      {remaining.length > 1 && (
        <Line points={remaining} color="#46d160" lineWidth={1.0} transparent opacity={0.25} />
      )}
    </>
  );
}

// -----------------------------------------------------------------------------
// CPAS event markers along the trajectory
// -----------------------------------------------------------------------------
function EventMarkers({ events, points, times }) {
  return events.map((e, i) => {
    // Find nearest trajectory index to this event time
    let idx = 0;
    let best = Infinity;
    for (let k = 0; k < times.length; k++) {
      const d = Math.abs(times[k] - e.t_s);
      if (d < best) {
        best = d;
        idx = k;
      }
    }
    const p = points[idx];
    if (!p) return null;
    return (
      <mesh key={i} position={p}>
        <sphereGeometry args={[40, 8, 8]} />
        <meshBasicMaterial color="#f6c043" />
      </mesh>
    );
  });
}

// -----------------------------------------------------------------------------
// Capsule marker, animated along the trajectory and oriented heat-shield
// forward (apex trailing along the anti-velocity direction).
// -----------------------------------------------------------------------------
function Capsule({ position, velocityDir, scale = 1 }) {
  const groupRef = useRef();

  // coneGeometry's natural axis is +Y (apex at +Y, base at -Y).
  // For heat-shield-forward flight, the BASE leads (faces +velocity)
  // and the APEX trails. So +Y of the cone should point opposite to V.
  useFrame(() => {
    if (!groupRef.current || !velocityDir) return;
    const apex = velocityDir.clone().negate().normalize();
    const q = new THREE.Quaternion().setFromUnitVectors(
      new THREE.Vector3(0, 1, 0),
      apex
    );
    groupRef.current.quaternion.copy(q);
  });

  // Squat Orion-ish proportions: ~5 m wide base, ~3 m tall
  // (cone radius 80, height 60 in scene units works visually)
  return (
    <group ref={groupRef} position={position}>
      <mesh>
        <coneGeometry args={[80 * scale, 60 * scale, 16]} />
        <meshStandardMaterial
          color="#ff4d4d"
          emissive="#ff4d4d"
          emissiveIntensity={0.55}
          roughness={0.4}
        />
      </mesh>
      {/* glowing heat-shield disk at the base (anti-apex end) */}
      <mesh position={[0, -30 * scale, 0]} rotation={[Math.PI / 2, 0, 0]}>
        <ringGeometry args={[60 * scale, 90 * scale, 24]} />
        <meshBasicMaterial color="#f6c043" transparent opacity={0.7} side={THREE.DoubleSide} />
      </mesh>
      <mesh>
        <sphereGeometry args={[120 * scale, 16, 16]} />
        <meshBasicMaterial color="#ff4d4d" transparent opacity={0.10} />
      </mesh>
    </group>
  );
}

// -----------------------------------------------------------------------------
// Scene root + animation loop
// -----------------------------------------------------------------------------
function Scene({
  points, times, events, playState, setIdx, currentIdx,
  targetVec, startVec, cameraMode, controlsRef, camera,
}) {
  // advance the playback index over time
  useFrame((_, dt) => {
    if (!playState.playing) return;
    setIdx((i) => {
      const next = i + playState.speed * dt * 4;
      if (next >= points.length - 1) {
        playState.playing = false;
        return points.length - 1;
      }
      return next;
    });
  });

  const idxFloor = Math.floor(currentIdx);
  const capsulePos = points[idxFloor] || points[0];

  // Velocity direction at the current point -- finite difference along the path
  let velocityDir = null;
  if (idxFloor + 1 < points.length) {
    velocityDir = points[idxFloor + 1].clone().sub(points[idxFloor]);
  } else if (idxFloor > 0) {
    velocityDir = points[idxFloor].clone().sub(points[idxFloor - 1]);
  }

  // In follow mode, keep the OrbitControls target glued to the capsule
  useFrame(() => {
    if (cameraMode === "follow" && controlsRef.current && capsulePos) {
      controlsRef.current.target.copy(capsulePos);
    }
  });

  return (
    <>
      <ambientLight intensity={0.4} />
      <directionalLight position={[15000, 15000, 8000]} intensity={1.2} color="#ffffff" />
      <directionalLight position={[-12000, -8000, -6000]} intensity={0.3} color="#5ddbe0" />

      <Stars radius={400000} depth={80000} count={4000} factor={400} fade speed={0.5} />

      <Earth />

      <TrajectoryPath points={points} currentIdx={Math.floor(currentIdx)} />

      <EventMarkers events={events} points={points} times={times} />

      {/* Start marker (cyan) */}
      {startVec && (
        <mesh position={startVec}>
          <octahedronGeometry args={[80, 0]} />
          <meshBasicMaterial color="#5ddbe0" />
        </mesh>
      )}
      {/* Target marker (red cross via two boxes) */}
      {targetVec && (
        <group position={targetVec}>
          <mesh>
            <boxGeometry args={[200, 30, 30]} />
            <meshBasicMaterial color="#ff4d4d" />
          </mesh>
          <mesh>
            <boxGeometry args={[30, 200, 30]} />
            <meshBasicMaterial color="#ff4d4d" />
          </mesh>
          <mesh>
            <boxGeometry args={[30, 30, 200]} />
            <meshBasicMaterial color="#ff4d4d" />
          </mesh>
        </group>
      )}

      <Capsule position={capsulePos} velocityDir={velocityDir} scale={1.5} />

      <OrbitControls
        ref={controlsRef}
        enablePan={true}
        enableZoom={true}
        enableRotate={true}
        zoomSpeed={0.8}
        rotateSpeed={0.5}
        minDistance={50}
        maxDistance={200000}
      />
    </>
  );
}

// -----------------------------------------------------------------------------
// Public component
// -----------------------------------------------------------------------------
export default function Trajectory3D({ trajectory, cpasEvents, runSummary, landingSummary }) {
  // Visual altitude exaggeration multiplier. Real Earth radius is 6378 km;
  // the entry corridor of 120 km is only ~2% of that, so at to-scale the
  // capsule looks like it's hugging the surface. Multiplying altitude
  // (only above-ground portion) makes the descent profile visible.
  const [altScale, setAltScale] = useState(30);

  // Downsample the trajectory to keep frame rate up — ~600 points is plenty.
  const { points, times } = useMemo(() => {
    const stride = Math.max(1, Math.floor(trajectory.length / 600));
    const pts = [];
    const ts = [];
    for (let i = 0; i < trajectory.length; i += stride) {
      const r = trajectory[i];
      pts.push(geodeticToVec3(r.phi_rad, r.lam_rad, r.alt_m, 6378.1, altScale));
      ts.push(r.t_s);
    }
    const lastR = trajectory[trajectory.length - 1];
    pts.push(geodeticToVec3(lastR.phi_rad, lastR.lam_rad, lastR.alt_m, 6378.1, altScale));
    ts.push(lastR.t_s);
    return { points: pts, times: ts };
  }, [trajectory, altScale]);

  const startVec = points[0];
  const targetVec = useMemo(() => {
    if (!runSummary) return null;
    // Target is at altitude 0, so altScale doesn't affect it.
    return geodeticToVec3(
      (runSummary.target_phi_deg * Math.PI) / 180,
      (runSummary.target_lam_deg * Math.PI) / 180,
      0,
      6378.1,
      altScale,
    );
  }, [runSummary, altScale]);

  // Filter to major events for markers
  const majorEvents = useMemo(() => {
    const majors = new Set([
      "fbc_jettison", "drogue_deploy", "drogue_disreef",
      "pilot_deploy", "drogue_cut", "main_deploy",
      "main_disreef_1", "main_disreef_2", "landed",
    ]);
    return cpasEvents.filter((e) => majors.has(e.event));
  }, [cpasEvents]);

  // Animation state — playing, speed, current index
  const [currentIdx, setCurrentIdx] = useState(0);
  const playState = useRef({ playing: false, speed: 25 });
  const [, force] = useState(0);

  // Camera focus mode + OrbitControls ref
  const [cameraMode, setCameraMode] = useState("trajectory");
  const controlsRef = useRef();

  // Compute trajectory centroid + bounding-box radius so we can frame it.
  const trajFrame = useMemo(() => {
    if (!points.length) return null;
    const c = new THREE.Vector3();
    points.forEach((p) => c.add(p));
    c.multiplyScalar(1 / points.length);
    let maxDist = 0;
    points.forEach((p) => { maxDist = Math.max(maxDist, p.distanceTo(c)); });
    return { center: c, radius: maxDist };
  }, [points]);

  // When the user switches focus mode, snap the OrbitControls target and
  // camera distance to the right place for that mode.
  useEffect(() => {
    if (!controlsRef.current) return;
    const controls = controlsRef.current;
    const cam = controls.object;

    if (cameraMode === "earth") {
      controls.target.set(0, 0, 0);
      cam.position.set(16000, 12000, 16000);
    } else if (cameraMode === "trajectory" && trajFrame) {
      controls.target.copy(trajFrame.center);
      // Sit the camera at ~3x the trajectory's bounding radius from its
      // centroid, looking back toward Earth's centre.
      const outward = trajFrame.center.clone().normalize();
      const dist = Math.max(800, trajFrame.radius * 3.0);
      cam.position.copy(
        trajFrame.center.clone().add(outward.multiplyScalar(dist * 0.7))
          .add(new THREE.Vector3(0, dist * 0.5, dist * 0.5))
      );
    } else if (cameraMode === "follow") {
      const p = points[Math.floor(currentIdx)] || points[0];
      controls.target.copy(p);
      const outward = p.clone().normalize();
      cam.position.copy(p.clone().add(outward.multiplyScalar(400))
        .add(new THREE.Vector3(0, 300, 300)));
    }
    controls.update();
  }, [cameraMode, trajFrame, points, altScale]);  // re-frame when scale changes

  const togglePlay = () => {
    if (currentIdx >= points.length - 1) setCurrentIdx(0);
    playState.current.playing = !playState.current.playing;
    force((x) => x + 1);
  };
  const setSpeed = (s) => {
    playState.current.speed = s;
    force((x) => x + 1);
  };
  const onScrub = (e) => {
    setCurrentIdx(Number(e.target.value));
    playState.current.playing = false;
    force((x) => x + 1);
  };
  const reset = () => {
    setCurrentIdx(0);
    playState.current.playing = false;
    force((x) => x + 1);
  };

  const currentRow = trajectory[Math.floor((currentIdx / (points.length - 1)) * (trajectory.length - 1))];
  const currentT = currentRow?.t_s ?? 0;
  const currentAlt = currentRow?.alt_m ?? 0;
  const currentV = currentRow?.V_mps ?? 0;
  const currentPhase = currentRow?.cpas_phase ?? "—";

  return (
    <Panel
      title="3D Trajectory — Live Playback"
      corner={`${points.length} samples`}
      status={
        <Pill kind={playState.current.playing ? "info" : "nominal"}>
          {playState.current.playing ? "PLAYING" : "PAUSED"}
        </Pill>
      }
    >
      <div style={{ position: "relative", height: 500, background: "var(--bg-0)", borderRadius: 3, border: "1px solid var(--border)" }}>
        <Canvas
          camera={{ position: [16000, 12000, 16000], fov: 45, near: 10, far: 1_000_000 }}
          style={{ background: "transparent" }}
        >
          <Scene
            points={points}
            times={times}
            events={majorEvents}
            playState={playState.current}
            setIdx={setCurrentIdx}
            currentIdx={currentIdx}
            startVec={startVec}
            targetVec={targetVec}
            cameraMode={cameraMode}
            controlsRef={controlsRef}
          />
        </Canvas>

        {/* HUD overlay */}
        <div style={{
          position: "absolute", top: 10, left: 10,
          background: "rgba(5,7,9,0.78)", border: "1px solid var(--border)",
          padding: "8px 12px", borderRadius: 3, fontSize: 11, lineHeight: 1.7,
          fontFamily: "var(--mono)",
        }}>
          <div style={{ color: "var(--text-dim)" }}>T+ <span style={{ color: "var(--cyan)" }}>{fmt.num(currentT, 1)}s</span></div>
          <div style={{ color: "var(--text-dim)" }}>ALT <span style={{ color: "var(--amber)" }}>{fmt.alt(currentAlt)}</span></div>
          <div style={{ color: "var(--text-dim)" }}>V <span style={{ color: "var(--green)" }}>{fmt.vel(currentV)}</span></div>
          <div style={{ color: "var(--text-dim)" }}>PHASE <span style={{ color: "var(--magenta)" }}>{currentPhase}</span></div>
        </div>

        {/* Legend */}
        <div style={{
          position: "absolute", top: 10, right: 10,
          background: "rgba(5,7,9,0.78)", border: "1px solid var(--border)",
          padding: "8px 12px", borderRadius: 3, fontSize: 10, lineHeight: 1.6,
          fontFamily: "var(--mono)", color: "var(--text-dim)",
        }}>
          <div><span style={{color:"#5ddbe0"}}>◆</span> entry interface</div>
          <div><span style={{color:"#46d160"}}>━</span> trajectory</div>
          <div><span style={{color:"#f6c043"}}>●</span> CPAS event</div>
          <div><span style={{color:"#ff4d4d"}}>▼</span> capsule</div>
          <div><span style={{color:"#ff4d4d"}}>✕</span> target</div>
        </div>
      </div>

      {/* Playback controls */}
      <div style={{ display: "flex", gap: 12, alignItems: "center", flexWrap: "wrap", marginTop: 6 }}>
        <button onClick={togglePlay} style={btnStyle}>
          {playState.current.playing ? "⏸ PAUSE" : "▶ PLAY"}
        </button>
        <button onClick={reset} style={btnStyle}>⏮ RESET</button>
        <div style={{ display: "flex", gap: 4 }}>
          {[1, 5, 25, 100, 500].map((s) => (
            <button
              key={s}
              onClick={() => setSpeed(s)}
              style={{
                ...btnStyle,
                color: playState.current.speed === s ? "var(--cyan)" : "var(--text-dim)",
                borderColor: playState.current.speed === s ? "var(--cyan)" : "var(--border)",
              }}
            >
              {s}×
            </button>
          ))}
        </div>
        <input
          type="range"
          min={0}
          max={points.length - 1}
          value={currentIdx}
          onChange={onScrub}
          style={{ flex: 1, accentColor: "var(--green)" }}
        />
        <div style={{ color: "var(--text-dim)", fontSize: 11, minWidth: 90, textAlign: "right" }}>
          {Math.floor(currentIdx)} / {points.length - 1}
        </div>
      </div>

      {/* Camera focus toggle */}
      <div style={{ display: "flex", gap: 12, alignItems: "center", marginTop: 4 }}>
        <span style={{ color: "var(--text-dim)", fontSize: 10, letterSpacing: "0.16em" }}>
          CAMERA
        </span>
        <div style={{ display: "flex", gap: 4 }}>
          {[
            { v: "earth",      label: "EARTH" },
            { v: "trajectory", label: "TRAJECTORY" },
            { v: "follow",     label: "FOLLOW CAPSULE" },
          ].map(({ v, label }) => (
            <button
              key={v}
              onClick={() => setCameraMode(v)}
              style={{
                ...btnStyle,
                color: cameraMode === v ? "var(--cyan)" : "var(--text-dim)",
                borderColor: cameraMode === v ? "var(--cyan)" : "var(--border)",
              }}
            >
              {label}
            </button>
          ))}
        </div>
        <span style={{ color: "var(--text-dim)", fontSize: 10, marginLeft: "auto" }}>
          {cameraMode === "earth"      && "orbiting Earth's centre"}
          {cameraMode === "trajectory" && "framed on the trajectory centroid"}
          {cameraMode === "follow"     && "locked to capsule — orbit/zoom still work"}
        </span>
      </div>

      {/* Altitude exaggeration toggle */}
      <div style={{ display: "flex", gap: 12, alignItems: "center", marginTop: 4 }}>
        <span style={{ color: "var(--text-dim)", fontSize: 10, letterSpacing: "0.16em" }}>
          ALT SCALE
        </span>
        <div style={{ display: "flex", gap: 4 }}>
          {[
            { v: 1,   label: "1× (true)" },
            { v: 10,  label: "10×" },
            { v: 30,  label: "30×" },
            { v: 100, label: "100×" },
            { v: 300, label: "300×" },
          ].map(({ v, label }) => (
            <button
              key={v}
              onClick={() => setAltScale(v)}
              style={{
                ...btnStyle,
                color: altScale === v ? "var(--amber)" : "var(--text-dim)",
                borderColor: altScale === v ? "var(--amber)" : "var(--border)",
              }}
            >
              {label}
            </button>
          ))}
        </div>
        <span style={{ color: "var(--text-dim)", fontSize: 10, marginLeft: "auto" }}>
          {altScale === 1
            ? "true scale — entry corridor invisible at this zoom"
            : `altitude above surface visually multiplied ×${altScale}`}
        </span>
      </div>
    </Panel>
  );
}

const btnStyle = {
  background: "var(--bg-2)",
  border: "1px solid var(--border)",
  color: "var(--text)",
  fontFamily: "var(--mono)",
  fontSize: 11,
  letterSpacing: "0.12em",
  padding: "4px 12px",
  cursor: "pointer",
  borderRadius: 2,
};
