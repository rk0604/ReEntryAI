export default function Panel({ title, corner, children, status }) {
  return (
    <div className="panel">
      <div className="head">
        <span className="title">{title}</span>
        {status}
        <span className="corner">{corner}</span>
      </div>
      <div className="body">{children}</div>
    </div>
  );
}

export function Pill({ kind = "nominal", children }) {
  return (
    <span className={`pill ${kind}`}>
      <span className="dot" />
      {children}
    </span>
  );
}
