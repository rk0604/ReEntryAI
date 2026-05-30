import { useEffect, useState } from "react";
import Papa from "papaparse";

/**
 * Fetch a CSV from /data/ and parse it with PapaParse.
 */
async function loadCSV(path) {
  const res = await fetch(path);
  if (!res.ok) throw new Error(`Failed to load ${path}: ${res.status}`);
  const text = await res.text();
  return new Promise((resolve, reject) => {
    Papa.parse(text, {
      header: true,
      dynamicTyping: true,
      skipEmptyLines: true,
      complete: (r) => resolve(r.data),
      error: reject,
    });
  });
}

async function loadJSON(path) {
  const res = await fetch(path);
  if (!res.ok) throw new Error(`Failed to load ${path}: ${res.status}`);
  return res.json();
}

/**
 * Load everything in revision_v1 in parallel.
 *
 * Returns:
 *   { status: 'loading' | 'ready' | 'error',
 *     runSummary, landingSummary, cpasEvents, trajectory, thrusterFires,
 *     error }
 */
export function useMissionData() {
  const [state, setState] = useState({ status: "loading" });

  useEffect(() => {
    let cancelled = false;
    (async () => {
      try {
        const [
          runSummary,
          landingSummary,
          cpasEvents,
          trajectory,
          thrusterFires,
        ] = await Promise.all([
          loadJSON("/data/run_summary.json"),
          loadJSON("/data/landing_summary.json"),
          loadCSV("/data/cpas_events.csv"),
          loadCSV("/data/trajectory.csv"),
          loadCSV("/data/thruster_fires.csv"),
        ]);
        if (cancelled) return;
        setState({
          status: "ready",
          runSummary,
          landingSummary,
          cpasEvents,
          trajectory,
          thrusterFires,
        });
      } catch (err) {
        if (!cancelled) {
          setState({ status: "error", error: err.message || String(err) });
        }
      }
    })();
    return () => {
      cancelled = true;
    };
  }, []);

  return state;
}
