/**
 * Display formatters for mission readouts.
 */

const isNum = (x) => typeof x === "number" && Number.isFinite(x);

export const fmt = {
  // Mission Elapsed Time as HH:MM:SS
  met(sec) {
    if (!isNum(sec)) return "--:--:--";
    sec = Math.max(0, Math.round(sec));
    const h = Math.floor(sec / 3600);
    const m = Math.floor((sec % 3600) / 60);
    const s = sec % 60;
    return [h, m, s].map((n) => String(n).padStart(2, "0")).join(":");
  },

  t(sec, digits = 2) {
    return isNum(sec) ? `${sec.toFixed(digits)}s` : "—";
  },

  alt(m) {
    if (!isNum(m)) return "—";
    if (Math.abs(m) >= 1000) return `${(m / 1000).toFixed(2)} km`;
    return `${m.toFixed(1)} m`;
  },

  vel(mps) {
    return isNum(mps) ? `${mps.toFixed(2)} m/s` : "—";
  },

  km(km) {
    return isNum(km) ? `${km.toFixed(2)} km` : "—";
  },

  m(m, digits = 1) {
    return isNum(m) ? `${m.toFixed(digits)} m` : "—";
  },

  deg(d, digits = 3) {
    return isNum(d) ? `${d.toFixed(digits)}°` : "—";
  },

  int(n) {
    return isNum(n) ? n.toLocaleString() : "—";
  },

  num(n, digits = 3) {
    return isNum(n) ? n.toFixed(digits) : "—";
  },

  bool(b) {
    return b ? "TRUE" : "FALSE";
  },

  truthClass(b, good = true) {
    if (b === good) return "green";
    return "amber";
  },
};

/**
 * Classify a CPAS event name for table colour-coding.
 */
export function classifyCpasEvent(evt) {
  if (!evt) return "";
  if (evt === "landed") return "evt-landed";
  if (evt.includes("squid")) return "evt-squid";
  if (evt.includes("deploy") || evt === "fbc_jettison") return "evt-deploy";
  if (evt.includes("disreef")) return "evt-major";
  return "";
}
