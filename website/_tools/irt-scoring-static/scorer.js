const NORMAL_LOG_CONST = -0.5 * Math.log(2 * Math.PI);

export function normalizeBinaryResponse(value) {
  if (value == null) return null;
  if (typeof value === "boolean") return value ? 1 : 0;
  if (typeof value === "number") return value === 0 || value === 1 ? value : null;

  const raw = String(value).trim();
  if (!raw) return null;
  const low = raw.toLowerCase();

  if (["na", "n/a", "missing"].includes(low)) return null;
  if (["yes", "y", "true", "1"].includes(low)) return 1;
  if (["no", "n", "false", "0"].includes(low)) return 0;

  const numeric = Number(raw);
  return numeric === 0 || numeric === 1 ? numeric : null;
}

export function standardizeSex(value) {
  if (value == null) return null;
  const raw = String(value).trim();
  if (!raw) return null;
  const low = raw.toLowerCase();
  if (["male", "m", "boy"].includes(low)) return "Male";
  if (["female", "f", "girl"].includes(low)) return "Female";
  return raw;
}

export function standardizeAgeGroup(value) {
  if (value == null) return null;
  const raw = String(value).trim();
  if (!raw) return null;
  const low = raw.toLowerCase();
  if (["age_17_19", "17-19", "17_19", "17 to 19"].includes(low)) return "Age_17_19";
  if (["age_20_24", "20-24", "20_24", "20 to 24"].includes(low)) return "Age_20_24";
  return raw;
}

export function prettyAgeGroup(value) {
  return ({
    Age_17_19: "17-19",
    Age_20_24: "20-24"
  })[value] || value;
}

export function prettyGroupLevel(value) {
  if (value == null || value === "") return "-";
  if (value === "Overall") return "Overall";
  if (value.includes("__")) {
    const [sex, ageGroup] = value.split("__");
    return [sex, prettyAgeGroup(ageGroup)].filter(Boolean).join(" ");
  }
  return prettyAgeGroup(value);
}

export function ensureGroups(record) {
  const out = { ...record };
  if ("SEX" in out) out.SEX = standardizeSex(out.SEX);
  if ("AGEGRP" in out) out.AGEGRP = standardizeAgeGroup(out.AGEGRP);
  if (!("AGEGRP" in out) && out.AGE != null && out.AGE !== "") {
    const age = Number(out.AGE);
    out.AGEGRP = Number.isFinite(age) ? (age < 20 ? "Age_17_19" : "Age_20_24") : null;
  }
  if (!("SEX_AGEGRP" in out)) {
    out.SEX_AGEGRP = out.SEX && out.AGEGRP ? `${out.SEX}__${out.AGEGRP}` : null;
  }
  return out;
}

function logistic(x) {
  if (x >= 0) {
    const z = Math.exp(-x);
    return 1 / (1 + z);
  }
  const z = Math.exp(x);
  return z / (1 + z);
}

function normalLogDensity(theta) {
  return NORMAL_LOG_CONST - 0.5 * theta * theta;
}

function linearInterpolate(x, points, xKey, yKey) {
  const n = points.length;
  if (n === 0) return null;
  if (x <= points[0][xKey]) return Number(points[0][yKey]);
  if (x >= points[n - 1][xKey]) return Number(points[n - 1][yKey]);
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const mid = Math.floor((lo + hi) / 2);
    if (points[mid][xKey] <= x) lo = mid;
    else hi = mid;
  }
  const x0 = Number(points[lo][xKey]);
  const x1 = Number(points[hi][xKey]);
  const y0 = Number(points[lo][yKey]);
  const y1 = Number(points[hi][yKey]);
  const t = (x - x0) / (x1 - x0);
  return y0 + t * (y1 - y0);
}

function itemsAnswered(record, items) {
  return items.reduce((acc, item) => {
    const v = record[item];
    return acc + (v === 0 || v === 1 ? 1 : 0);
  }, 0);
}

function totalScore(record, items) {
  let answered = 0;
  let total = 0;
  for (const item of items) {
    const v = record[item];
    if (v === 0 || v === 1) {
      answered += 1;
      total += v;
    }
  }
  return answered === 0 ? null : total;
}

function eapTheta(record, engine) {
  const answered = engine.items.filter((item) => record[item] === 0 || record[item] === 1);
  if (answered.length === 0) return null;

  const grid = engine.theta_map.map((row) => Number(row.Theta));
  const logPost = grid.map((theta) => {
    let ll = normalLogDensity(theta);
    for (const item of answered) {
      const params = engine.parameters[item];
      const p = logistic(params.a * theta + params.d);
      ll += record[item] === 1 ? Math.log(Math.max(p, 1e-12)) : Math.log(Math.max(1 - p, 1e-12));
    }
    return ll;
  });

  const maxLog = Math.max(...logPost);
  const weights = logPost.map((v) => Math.exp(v - maxLog));
  const weightSum = weights.reduce((a, b) => a + b, 0);
  if (!Number.isFinite(weightSum) || weightSum <= 0) return null;

  let mean = 0;
  for (let i = 0; i < grid.length; i += 1) {
    mean += grid[i] * weights[i];
  }
  return mean / weightSum;
}

function resolveCutInfo(record, theta, engine) {
  if (theta == null) {
    return {
      cutoff: null,
      grouping: null,
      level: null,
      depression: null
    };
  }

  const table = engine.cutoff_table || [];
  const overall = table.find((row) => row.Grouping === "none" && row.Group_Level === "Overall") || null;
  const overallCutoff = overall ? Number(overall.Cutoff_Theta) : null;
  const grouped = ensureGroups(record);
  const mode = engine.cutoff_apply_mode || "none";

  const matchCut = (grouping, level) => {
    const row = table.find((entry) => entry.Grouping === grouping && entry.Group_Level === level);
    if (!row) return null;
    return {
      cutoff: Number(row.Cutoff_Theta),
      grouping,
      level
    };
  };

  let cutInfo = null;
  if (mode === "none") {
    cutInfo = overall ? { cutoff: overallCutoff, grouping: "none", level: "Overall" } : null;
  } else {
    const level = grouped[mode];
    if (level != null && level !== "") {
      cutInfo = matchCut(mode, String(level));
    }
    if (!cutInfo && overall) {
      cutInfo = { cutoff: overallCutoff, grouping: "none", level: "Overall" };
    }
  }

  if (!cutInfo || !Number.isFinite(cutInfo.cutoff)) {
    return {
      cutoff: null,
      grouping: null,
      level: null,
      depression: null
    };
  }

  return {
    ...cutInfo,
    depression: theta >= cutInfo.cutoff ? 1 : 0
  };
}

export function scoreRecord(rawRecord, engine) {
  const record = { ...rawRecord };
  for (const item of engine.items) {
    record[item] = normalizeBinaryResponse(record[item]);
  }

  const answered = itemsAnswered(record, engine.items);
  const ssqTotal = totalScore(record, engine.items);
  const theta = eapTheta(record, engine);
  const phqExpected = theta == null ? null : linearInterpolate(theta, engine.theta_map, "Theta", "PHQ_Expected_Sum");
  const phqProb = theta == null ? null : linearInterpolate(theta, engine.theta_map, "Theta", "PHQ_Prob_GE10");
  const cut = resolveCutInfo(record, theta, engine);

  const grouped = ensureGroups(record);

  return {
    ...grouped,
    SSQ10_Total_Score: ssqTotal,
    Theta_Harmonized: theta,
    Depression_Binary: cut.depression,
    PHQ_Expected_Sum: phqExpected,
    PHQ_Prob_GE10: phqProb,
    Cutoff_Theta_Applied: cut.cutoff,
    Cutoff_Grouping_Used: cut.grouping,
    Cutoff_Group_Level_Used: cut.level,
    Cutoff_Apply_Mode: engine.cutoff_apply_mode,
    Cutoff_Method: engine.cutoff_method,
    Engine_ID: engine.engine_id,
    Engine_Label: engine.label,
    Items_Answered: answered
  };
}

export function scoreBatch(rows, engine) {
  return rows.map((row) => scoreRecord(row, engine));
}

export function parseCsv(text) {
  const rows = [];
  let row = [];
  let cell = "";
  let inQuotes = false;

  const pushCell = () => {
    row.push(cell);
    cell = "";
  };

  const pushRow = () => {
    rows.push(row);
    row = [];
  };

  for (let i = 0; i < text.length; i += 1) {
    const ch = text[i];
    const next = text[i + 1];

    if (ch === "\"") {
      if (inQuotes && next === "\"") {
        cell += "\"";
        i += 1;
      } else {
        inQuotes = !inQuotes;
      }
    } else if (ch === "," && !inQuotes) {
      pushCell();
    } else if ((ch === "\n" || ch === "\r") && !inQuotes) {
      if (ch === "\r" && next === "\n") i += 1;
      pushCell();
      pushRow();
    } else {
      cell += ch;
    }
  }

  if (cell.length > 0 || row.length > 0) {
    pushCell();
    pushRow();
  }

  if (rows.length === 0) return [];
  const headers = rows[0].map((x) => x.trim());
  return rows.slice(1).filter((r) => r.some((x) => x !== "")).map((r) => {
    const out = {};
    headers.forEach((header, idx) => {
      out[header] = r[idx] ?? "";
    });
    return out;
  });
}

function escapeCsv(value) {
  if (value == null) return "";
  const text = String(value);
  if (/[",\n\r]/.test(text)) {
    return `"${text.replace(/"/g, "\"\"")}"`;
  }
  return text;
}

export function toCsv(rows) {
  if (!rows || rows.length === 0) return "";
  const headers = Array.from(new Set(rows.flatMap((row) => Object.keys(row))));
  const lines = [
    headers.join(","),
    ...rows.map((row) => headers.map((header) => escapeCsv(row[header])).join(","))
  ];
  return lines.join("\n");
}

export function prettyMode(value) {
  return ({
    none: "Overall",
    SEX: "Sex",
    AGEGRP: "Age group",
    SEX_AGEGRP: "Sex + age group"
  })[value] || value;
}

export function prettyMethod(value) {
  return ({
    sens_at_spec: "Sensitivity at target specificity",
    spec_at_sens: "Specificity at target sensitivity",
    youden: "Youden index threshold",
    model_tcc: "Model-based TCC threshold"
  })[value] || value;
}
