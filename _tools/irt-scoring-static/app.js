import {
  normalizeBinaryResponse,
  parseCsv,
  prettyAgeGroup,
  prettyGroupLevel,
  prettyMethod,
  prettyMode,
  scoreBatch,
  scoreRecord,
  toCsv
} from "./scorer.js";

const bundlePath = "./engines/production_engines.json";
const exampleCsvPath = "./engines/example_batch.csv";
const DASH = "-";

const tabs = [
  { id: "overview", label: "Overview" },
  { id: "manual", label: "Manual Scoring" },
  { id: "batch", label: "Batch Scoring" },
  { id: "results", label: "Results" },
  { id: "privacy", label: "Privacy" },
  { id: "documentation", label: "Documentation" }
];

const state = {
  bundle: null,
  engine: null,
  activeTab: "overview",
  latestBatch: null,
  latestBatchSourceRows: null,
  latestManualInput: null,
  latestRun: null
};

const dom = {
  bundleStatus: document.getElementById("bannerEngine"),
  bannerTab: document.getElementById("bannerTab"),
  bannerLatestSummary: document.getElementById("bannerLatestSummary"),
  openBatchBtn: document.getElementById("openBatchBtn"),
  manualEngineSelect: document.getElementById("manualEngineSelect"),
  batchEngineSelect: document.getElementById("batchEngineSelect"),
  manualItems: document.getElementById("manualItems"),
  manualSex: document.getElementById("manualSex"),
  manualAge: document.getElementById("manualAge"),
  scoreManualBtn: document.getElementById("scoreManualBtn"),
  resetManualBtn: document.getElementById("resetManualBtn"),
  manualHighlights: document.getElementById("manualHighlights"),
  manualStatus: document.getElementById("manualStatus"),
  manualResultTable: document.getElementById("manualResultTable"),
  manualCompareNote: document.getElementById("manualCompareNote"),
  manualEngineCompareTable: document.getElementById("manualEngineCompareTable"),
  batchFile: document.getElementById("batchFile"),
  scoreBatchBtn: document.getElementById("scoreBatchBtn"),
  downloadExampleBtn: document.getElementById("downloadExampleBtn"),
  downloadBatchBtn: document.getElementById("downloadBatchBtn"),
  batchHighlights: document.getElementById("batchHighlights"),
  batchStatus: document.getElementById("batchStatus"),
  batchPreviewTable: document.getElementById("batchPreviewTable"),
  latestSummaryCards: document.getElementById("latestSummaryCards"),
  latestSummaryText: document.getElementById("latestSummaryText"),
  latestResultPreviewTable: document.getElementById("latestResultPreviewTable"),
  overviewStats: document.getElementById("overviewStats"),
  overviewCurrent: document.getElementById("overviewCurrent"),
  overviewSentence: document.getElementById("overviewSentence"),
  engineTable: document.getElementById("engineTable"),
  schemaTableOverview: document.getElementById("schemaTableOverview"),
  schemaTableDoc: document.getElementById("schemaTableDoc"),
  tabButtons: Array.from(document.querySelectorAll("[data-tab-target]")),
  tabPanels: Array.from(document.querySelectorAll("[data-tab-panel]")),
  gotoButtons: Array.from(document.querySelectorAll("[data-goto]"))
};

const toolHelpConfig = {
  title: "IRT Scoring Workspace",
  intro: "Use this help panel for orientation while scoring SSQ-10 records. Switch Guide on to reveal short hover notes on fields, buttons, tables, and output labels.",
  sections: [
    {
      heading: "Recommended flow",
      items: [
        "Start in Overview to confirm the active engine, threshold scope, and expected input schema.",
        "Use Manual Scoring for one participant or test profile before running a larger file.",
        "Use Batch Scoring for a CSV cohort file, then review validation notes before downloading outputs.",
        "Use Results and Documentation to confirm what the returned fields mean before sharing scored data."
      ]
    },
    {
      heading: "Before exporting",
      items: [
        "Confirm that the selected engine and threshold strategy match the intended analysis.",
        "Check validation messages and the preview table before downloading a scored CSV.",
        "Remember that the depression indicator is for harmonised research scoring, not standalone clinical diagnosis."
      ]
    },
    {
      heading: "Guide mode",
      items: [
        "Guide Off keeps the interface clean.",
        "Guide On highlights supported labels and controls. Hover or focus them to see a short explanation.",
        "Guide notes are interpretive aids; they do not change the scoring engine or exported data."
      ]
    }
  ]
};

const guideState = {
  enabled: false,
  tooltip: null,
  queued: false
};

const toolGuideDefinitions = new Map([
  ["Browser-local scoring", "The scoring run happens in this browser session rather than on a remote server."],
  ["No raw data retention", "Uploaded raw data are not written back to project storage by this browser tool."],
  ["Open Batch Scoring", "Jump directly to the CSV upload workflow."],
  ["Selected engine", "The packaged IRT engine currently used for scoring and threshold application."],
  ["Current workspace", "The active tab or workspace section you are viewing."],
  ["Latest run", "A short summary of the most recent successful manual or batch scoring run."],
  ["Privacy mode", "A reminder that uploaded data remain in the active browser session unless you choose to download an output."],
  ["Overview", "Orientation to the active engine, workflow, expected schema, and limitations."],
  ["Manual Scoring", "Use this tab to score one participant profile at a time."],
  ["Batch Scoring", "Use this tab to upload, validate, score, preview, and download a CSV file."],
  ["Results", "Shows the latest scoring summary and a compact preview from the current browser session."],
  ["Privacy", "Explains what the browser tool does and does not retain."],
  ["Documentation", "Defines expected inputs, returned fields, abbreviations, and statistical terms."],
  ["Threshold scope", "The grouping level used when applying the stored decision threshold."],
  ["Threshold approach", "The rule used to select the threshold, such as sensitivity-oriented or specificity-oriented scoring."],
  ["Overall sensitivity", "The proportion of reference-positive participants correctly identified in the packaging sample."],
  ["Overall specificity", "The proportion of reference-negative participants correctly identified in the packaging sample."],
  ["Overall accuracy", "The overall agreement between engine classification and reference classification in the packaging sample."],
  ["Current selection", "The engine or result currently driving the visible output summary."],
  ["Engine version", "Choose which packaged production-safe engine to apply."],
  ["SEX (optional)", "Optional subgroup field used when the selected engine has sex-specific thresholds."],
  ["AGE (optional)", "Optional age field used when the selected engine has age-specific thresholds."],
  ["SSQ-10 items", "The ten binary SSQ item responses used by the scoring engine."],
  ["Blank", "Missing response. Blank values are not counted as endorsed symptoms."],
  ["0", "No or not endorsed."],
  ["1", "Yes or endorsed."],
  ["Score participant", "Run the selected engine on the current manual profile."],
  ["View Latest Results", "Move to the Results tab for the current session summary."],
  ["Reset", "Clear the manual scoring inputs and result preview."],
  ["Upload CSV", "Choose a cohort CSV file with the required SSQ-10 item columns."],
  ["Score uploaded file", "Validate and score the selected CSV in the browser."],
  ["Download scored CSV", "Save the scored output from the current browser session."],
  ["Download example CSV", "Download a small template file for testing the batch workflow."],
  ["Rows scored", "Number of records successfully scored in the current run."],
  ["Mean theta", "Average harmonised depression theta across the scored records."],
  ["Depression prevalence", "Proportion of scored records classified as depression-positive by the selected threshold."],
  ["Harmonized Theta", "The harmonised latent depression score returned by the IRT engine."],
  ["Depression", "Binary research-use classification based on the applied threshold."],
  ["Expected PHQ-9 Sum", "A PHQ-oriented score implied by the harmonised theta."],
  ["Probability PHQ-9 >= 10", "Estimated probability of meeting or exceeding a PHQ-9 score of 10."],
  ["Cutoff Target", "The subgroup or overall threshold actually applied to the record."],
  ["Items Answered", "Number of SSQ-10 item responses available for scoring."],
  ["Total SSQ-10 Score", "Observed raw sum of endorsed SSQ-10 items."]
]);

const toolGuideRules = [
  { test: (label) => /sensitivity/i.test(label), text: "Sensitivity summarizes how well reference-positive cases are detected." },
  { test: (label) => /specificity/i.test(label), text: "Specificity summarizes how well reference-negative cases are excluded." },
  { test: (label) => /accuracy/i.test(label), text: "Accuracy is the overall proportion of correct classifications in the packaging sample." },
  { test: (label) => /engine/i.test(label), text: "Engine labels identify the packaged scoring and threshold scenario being applied." },
  { test: (label) => /threshold|cutoff/i.test(label), text: "Threshold fields describe the decision rule used to convert theta into a binary indicator." },
  { test: (label) => /schema|required|expected values/i.test(label), text: "Schema fields describe what the uploaded data should contain before scoring." },
  { test: (label) => /download/i.test(label), text: "Download buttons create local files from the current browser session." }
];

function makeToolText(tag, className, text) {
  const element = document.createElement(tag);
  if (className) element.className = className;
  element.textContent = text;
  return element;
}

function normaliseGuideLabel(value) {
  return String(value || "").replace(/\s+/g, " ").trim().replace(/:$/, "");
}

function guideForLabel(label) {
  const cleanLabel = normaliseGuideLabel(label);
  if (!cleanLabel) return "";
  if (toolGuideDefinitions.has(cleanLabel)) return toolGuideDefinitions.get(cleanLabel);
  for (const [key, text] of toolGuideDefinitions) {
    if (cleanLabel.includes(key) || key.includes(cleanLabel)) return text;
  }
  const rule = toolGuideRules.find((item) => item.test(cleanLabel));
  return rule ? rule.text : "";
}

function annotateToolGuides() {
  const selector = [
    "h1",
    "h2",
    "h3",
    ".mini-label",
    ".field-label",
    ".tab-btn",
    ".primary-btn",
    ".secondary-btn",
    ".ghost-btn",
    ".status-pill",
    ".workflow-step strong",
    ".guide-card h3",
    ".stat-card span",
    ".choice-value",
    ".choice-label",
    "th"
  ].join(",");

  document.querySelectorAll(selector).forEach((element) => {
    if (element.closest(".tool-help-root") || element.closest(".tool-help-panel")) return;
    const guide = guideForLabel(element.textContent);
    if (!guide) return;
    element.classList.add("tool-guide-target");
    element.dataset.toolGuideText = guide;
    if (!element.hasAttribute("tabindex") && !/^(BUTTON|A|INPUT|SELECT|TEXTAREA)$/.test(element.tagName)) {
      element.setAttribute("tabindex", "0");
    }
  });
}

function queueToolGuideAnnotation() {
  if (guideState.queued) return;
  guideState.queued = true;
  window.requestAnimationFrame(() => {
    guideState.queued = false;
    annotateToolGuides();
    updateToolGuide();
  });
}

function updateToolGuide() {
  document.documentElement.classList.toggle("tool-guide-on", guideState.enabled);
  document.querySelectorAll("[data-tool-guide-switch]").forEach((control) => {
    control.classList.toggle("is-active", guideState.enabled);
    control.setAttribute("aria-pressed", guideState.enabled ? "true" : "false");
  });
  document.querySelectorAll("[data-tool-guide-label]").forEach((label) => {
    label.textContent = guideState.enabled ? "On" : "Off";
  });
  if (!guideState.enabled) hideToolGuideTooltip();
}

function moveToolGuideTooltip(event) {
  if (!guideState.tooltip) return;
  const width = guideState.tooltip.offsetWidth || 300;
  const height = guideState.tooltip.offsetHeight || 90;
  const left = Math.min(window.innerWidth - width - 12, event.clientX + 14);
  const top = Math.min(window.innerHeight - height - 12, event.clientY + 14);
  guideState.tooltip.style.left = `${Math.max(12, left)}px`;
  guideState.tooltip.style.top = `${Math.max(12, top)}px`;
}

function showToolGuideTooltip(target, event = null) {
  if (!guideState.enabled || !target || !target.dataset.toolGuideText || !guideState.tooltip) return;
  guideState.tooltip.textContent = target.dataset.toolGuideText;
  guideState.tooltip.classList.add("is-visible");
  if (event) {
    moveToolGuideTooltip(event);
  } else {
    const rect = target.getBoundingClientRect();
    moveToolGuideTooltip({ clientX: rect.left, clientY: rect.bottom });
  }
}

function hideToolGuideTooltip() {
  if (guideState.tooltip) guideState.tooltip.classList.remove("is-visible");
}

function initToolHelp() {
  if (document.querySelector(".tool-help-root")) return;
  const actions = document.querySelector(".topbar-actions");
  if (!actions) return;

  const root = document.createElement("div");
  root.className = "tool-help-root";

  const helpButton = document.createElement("button");
  helpButton.className = "tool-help-button";
  helpButton.type = "button";
  helpButton.setAttribute("aria-expanded", "false");
  helpButton.setAttribute("aria-label", `Open help for ${toolHelpConfig.title}`);
  helpButton.append(makeToolText("span", "tool-help-button__icon", "?"), makeToolText("span", "", "Help"));

  const guideButton = document.createElement("button");
  guideButton.className = "tool-guide-button";
  guideButton.type = "button";
  guideButton.setAttribute("aria-pressed", "false");
  guideButton.setAttribute("aria-label", "Toggle tool guide");
  guideButton.setAttribute("data-tool-guide-switch", "");
  guideButton.appendChild(makeToolText("span", "", "Guide"));
  const guideLabel = makeToolText("b", "", "Off");
  guideLabel.setAttribute("data-tool-guide-label", "");
  guideButton.appendChild(guideLabel);

  const panel = document.createElement("aside");
  panel.className = "tool-help-panel";
  panel.setAttribute("aria-label", `${toolHelpConfig.title} help`);

  const handle = document.createElement("div");
  handle.className = "tool-help-panel__handle";
  const titleBlock = document.createElement("div");
  titleBlock.append(makeToolText("p", "", "Tool help"), makeToolText("h3", "", toolHelpConfig.title));
  const close = document.createElement("button");
  close.className = "tool-help-close";
  close.type = "button";
  close.setAttribute("aria-label", "Close help");
  close.textContent = "x";
  handle.append(titleBlock, close);

  const body = document.createElement("div");
  body.className = "tool-help-panel__body";
  body.appendChild(makeToolText("p", "tool-help-panel__intro", toolHelpConfig.intro));
  toolHelpConfig.sections.forEach((section) => {
    const block = document.createElement("section");
    block.className = "tool-help-section";
    block.appendChild(makeToolText("h4", "", section.heading));
    const list = document.createElement("ul");
    section.items.forEach((item) => {
      list.appendChild(makeToolText("li", "", item));
    });
    block.appendChild(list);
    body.appendChild(block);
  });
  panel.append(handle, body);

  root.append(helpButton, guideButton);
  actions.prepend(root);
  document.body.appendChild(panel);

  guideState.tooltip = document.createElement("div");
  guideState.tooltip.className = "tool-guide-tooltip";
  document.body.appendChild(guideState.tooltip);

  helpButton.addEventListener("click", () => {
    const open = !panel.classList.contains("is-open");
    panel.classList.toggle("is-open", open);
    helpButton.setAttribute("aria-expanded", open ? "true" : "false");
  });

  close.addEventListener("click", () => {
    panel.classList.remove("is-open");
    helpButton.setAttribute("aria-expanded", "false");
  });

  guideButton.addEventListener("click", () => {
    guideState.enabled = !guideState.enabled;
    annotateToolGuides();
    updateToolGuide();
  });

  document.addEventListener("mouseover", (event) => {
    const target = event.target.closest("[data-tool-guide-text]");
    if (target) showToolGuideTooltip(target, event);
  });
  document.addEventListener("mousemove", (event) => {
    if (guideState.tooltip?.classList.contains("is-visible")) moveToolGuideTooltip(event);
  });
  document.addEventListener("mouseout", (event) => {
    if (event.target.closest("[data-tool-guide-text]")) hideToolGuideTooltip();
  });
  document.addEventListener("focusin", (event) => {
    const target = event.target.closest("[data-tool-guide-text]");
    if (target) showToolGuideTooltip(target);
  });
  document.addEventListener("focusout", (event) => {
    if (event.target.closest("[data-tool-guide-text]")) hideToolGuideTooltip();
  });
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape") {
      panel.classList.remove("is-open");
      helpButton.setAttribute("aria-expanded", "false");
      hideToolGuideTooltip();
    }
  });

  const appShell = document.querySelector(".app-shell");
  if (appShell && "MutationObserver" in window) {
    const observer = new MutationObserver(queueToolGuideAnnotation);
    observer.observe(appShell, { childList: true, subtree: true });
  }

  annotateToolGuides();
  updateToolGuide();
}

function fmtNumber(value, digits = 3) {
  if (value == null || Number.isNaN(value)) return DASH;
  return Number(value).toFixed(digits);
}

function fmtPercent(value, digits = 1) {
  if (value == null || Number.isNaN(value)) return DASH;
  return `${(100 * Number(value)).toFixed(digits)}%`;
}

function metricCard(label, value, note = "") {
  const div = document.createElement("div");
  div.className = "stat-card";
  div.innerHTML = `<span>${label}</span><strong>${value}</strong>${note ? `<small>${note}</small>` : ""}`;
  return div;
}

function setBlockMessage(el, message, type = "empty") {
  el.className = type === "status" ? "data-status" : "empty-note";
  el.textContent = message;
}

function clearTable(tableEl) {
  tableEl.querySelector("thead").innerHTML = "";
  tableEl.querySelector("tbody").innerHTML = "";
}

function setTable(tableEl, rows, columns) {
  const thead = tableEl.querySelector("thead");
  const tbody = tableEl.querySelector("tbody");
  thead.innerHTML = "";
  tbody.innerHTML = "";

  if (!rows || rows.length === 0) return;

  const headRow = document.createElement("tr");
  columns.forEach(({ key, label }) => {
    const th = document.createElement("th");
    th.dataset.key = key;
    th.textContent = label;
    headRow.appendChild(th);
  });
  thead.appendChild(headRow);

  rows.forEach((row) => {
    const tr = document.createElement("tr");
    columns.forEach(({ key, format }) => {
      const td = document.createElement("td");
      const value = row[key];
      if (typeof format === "function") {
        td.textContent = format(value, row);
      } else if (value == null || value === "") {
        td.textContent = DASH;
      } else if (typeof value === "number") {
        td.textContent = Number.isInteger(value) ? String(value) : Number(value).toFixed(3);
      } else {
        td.textContent = String(value);
      }
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });
}

function mean(values) {
  const ok = values.filter((value) => value != null && !Number.isNaN(value));
  if (ok.length === 0) return null;
  return ok.reduce((acc, value) => acc + Number(value), 0) / ok.length;
}

function prevalence(values) {
  const ok = values.filter((value) => value === 0 || value === 1);
  if (ok.length === 0) return null;
  return ok.reduce((acc, value) => acc + value, 0) / ok.length;
}

function currentOverallMetrics(engine) {
  return (engine.metrics_applied || []).find(
    (row) => row.Grouping === "none" && row.Group_Level === "Overall"
  ) || null;
}

function currentTabLabel() {
  const match = tabs.find((tab) => tab.id === state.activeTab);
  return match ? match.label : "Overview";
}

function selectTab(tabId) {
  state.activeTab = tabId;
  dom.tabButtons.forEach((btn) => {
    const active = btn.dataset.tabTarget === tabId;
    btn.classList.toggle("active", active);
    btn.setAttribute("aria-selected", active ? "true" : "false");
  });
  dom.tabPanels.forEach((panel) => {
    panel.classList.toggle("is-active", panel.dataset.tabPanel === tabId);
  });
  dom.bannerTab.textContent = currentTabLabel();
}

function idColumns(rows) {
  if (!rows || rows.length === 0) return [];
  const headers = Object.keys(rows[0]);
  return headers
    .filter((name) => /(^id$|_id$|participant|record|study|case)/i.test(name))
    .filter((name) => !/^engine(_|$)/i.test(name))
    .slice(0, 2);
}

function inputSchemaRows(engine) {
  return [
    ...engine.items.map((item) => ({
      field: item,
      required: "Yes",
      expected_values: "0/1 or Yes/No"
    })),
    { field: "SEX", required: "No", expected_values: "Male/Female (optional)" },
    { field: "AGE", required: "No", expected_values: "Numeric age in years (optional)" },
    { field: "AGEGRP", required: "No", expected_values: "17-19 or 20-24 (optional)" }
  ];
}

function compactPreview(rows, maxRows = 10) {
  if (!rows || rows.length === 0) return [];
  const keep = rows.slice(0, maxRows).map((row) => ({ ...row }));
  return keep;
}

function previewColumns(rows) {
  const ids = idColumns(rows);
  return [
    ...ids.map((key) => ({ key, label: key.replace(/_/g, " ") })),
    { key: "SEX", label: "Sex" },
    { key: "AGE", label: "Age" },
    { key: "AGEGRP", label: "Age Group", format: (value) => (value ? prettyAgeGroup(value) : DASH) },
    { key: "Items_Answered", label: "Items Answered" },
    { key: "SSQ10_Total_Score", label: "Total SSQ-10 Score" },
    { key: "Theta_Harmonized", label: "Harmonized Theta", format: (value) => fmtNumber(value) },
    {
      key: "Depression_Binary",
      label: "Depression",
      format: (value) => (value === 1 ? "Yes" : value === 0 ? "No" : DASH)
    },
    { key: "PHQ_Expected_Sum", label: "Expected PHQ-9 Sum", format: (value) => fmtNumber(value) },
    { key: "PHQ_Prob_GE10", label: "Probability PHQ-9 >= 10", format: (value) => fmtPercent(value) },
    {
      key: "Cutoff_Group_Level_Used",
      label: "Cutoff Target",
      format: (value) => prettyGroupLevel(value)
    }
  ];
}

function latestSummaryCards(scoredRows, engine) {
  return [
    metricCard("Rows scored", String(scoredRows.length)),
    metricCard("Mean theta", fmtNumber(mean(scoredRows.map((row) => row.Theta_Harmonized)))),
    metricCard("Depression prevalence", fmtPercent(prevalence(scoredRows.map((row) => row.Depression_Binary)))),
    metricCard("Engine", engine.label, prettyMode(engine.cutoff_apply_mode)),
    metricCard("Threshold scope", prettyMode(engine.cutoff_apply_mode)),
    metricCard("Threshold approach", prettyMethod(engine.cutoff_method))
  ];
}

function refreshLatestRunView() {
  dom.latestSummaryCards.innerHTML = "";
  clearTable(dom.latestResultPreviewTable);

  if (!state.latestRun) {
    setBlockMessage(
      dom.latestSummaryText,
      "Only the latest result from the current browser session is shown here.",
      "status"
    );
    dom.bannerLatestSummary.textContent = "No scoring run yet.";
    return;
  }

  state.latestRun.cards.forEach((card) => dom.latestSummaryCards.appendChild(card));
  setBlockMessage(dom.latestSummaryText, state.latestRun.summaryText, "status");
  dom.bannerLatestSummary.textContent = state.latestRun.summaryText;
  setTable(
    dom.latestResultPreviewTable,
    compactPreview(state.latestRun.scoredRows),
    previewColumns(state.latestRun.scoredRows)
  );
}

function setLatestRun(runType, scoredRows, engine) {
  const summaryText = runType === "manual"
    ? `Manual scoring completed for 1 row using ${engine.label}.`
    : `Batch scoring completed for ${scoredRows.length} row(s) using ${engine.label}.`;

  state.latestRun = {
    type: runType,
    scoredRows,
    engine,
    summaryText,
    cards: latestSummaryCards(scoredRows, engine)
  };
  refreshLatestRunView();
}

function renderOverview(engine) {
  const overall = currentOverallMetrics(engine);
  dom.bundleStatus.textContent = engine.label;
  dom.overviewCurrent.textContent = engine.label;
  dom.overviewSentence.textContent = engine.notes || engine.description || DASH;
  dom.overviewStats.innerHTML = "";

  [
    metricCard("Threshold scope", prettyMode(engine.cutoff_apply_mode)),
    metricCard("Threshold approach", prettyMethod(engine.cutoff_method)),
    metricCard("Overall sensitivity", overall ? fmtPercent(overall.Sensitivity) : DASH),
    metricCard("Overall specificity", overall ? fmtPercent(overall.Specificity) : DASH),
    metricCard("Overall accuracy", overall ? fmtPercent(overall.Accuracy) : DASH)
  ].forEach((card) => dom.overviewStats.appendChild(card));
}

function renderEngineTable() {
  const rows = state.bundle.engines.map((engine) => {
    const overall = currentOverallMetrics(engine);
    return {
      label: engine.label,
      scope: prettyMode(engine.cutoff_apply_mode),
      approach: prettyMethod(engine.cutoff_method),
      sensitivity: overall ? overall.Sensitivity : null,
      specificity: overall ? overall.Specificity : null,
      accuracy: overall ? overall.Accuracy : null
    };
  });

  setTable(
    dom.engineTable,
    rows,
    [
      { key: "label", label: "Engine" },
      { key: "scope", label: "Threshold Scope" },
      { key: "approach", label: "Threshold Approach" },
      { key: "sensitivity", label: "Sensitivity", format: (value) => fmtPercent(value) },
      { key: "specificity", label: "Specificity", format: (value) => fmtPercent(value) },
      { key: "accuracy", label: "Accuracy", format: (value) => fmtPercent(value) }
    ]
  );
}

function renderSchemaTables(engine) {
  const rows = inputSchemaRows(engine);
  const columns = [
    { key: "field", label: "Field" },
    { key: "required", label: "Required" },
    { key: "expected_values", label: "Expected values" }
  ];
  setTable(dom.schemaTableOverview, rows, columns);
  setTable(dom.schemaTableDoc, rows, columns);
}

function renderManualItems(engine, values = {}) {
  dom.manualItems.innerHTML = "";
  engine.items.forEach((item) => {
    const current = values[item] ?? null;
    const card = document.createElement("div");
    card.className = "item-card";
    card.innerHTML = `
      <h4>${item}</h4>
      <div class="option-row">
        <label class="choice-chip">
          <input type="radio" name="${item}" value="" ${current == null ? "checked" : ""}>
          <span class="choice-text"><span class="choice-value">Blank</span><span class="choice-label">Missing</span></span>
        </label>
        <label class="choice-chip">
          <input type="radio" name="${item}" value="0" ${current === 0 ? "checked" : ""}>
          <span class="choice-text"><span class="choice-value">0</span><span class="choice-label">No</span></span>
        </label>
        <label class="choice-chip">
          <input type="radio" name="${item}" value="1" ${current === 1 ? "checked" : ""}>
          <span class="choice-text"><span class="choice-value">1</span><span class="choice-label">Yes</span></span>
        </label>
      </div>
    `;
    dom.manualItems.appendChild(card);
  });
}

function getManualRecord() {
  const record = {
    SEX: dom.manualSex.value || null,
    AGE: dom.manualAge.value ? Number(dom.manualAge.value) : null
  };
  state.engine.items.forEach((item) => {
    const selected = document.querySelector(`input[name="${item}"]:checked`);
    record[item] = selected ? normalizeBinaryResponse(selected.value) : null;
  });
  return record;
}

function renderManualResult(result) {
  dom.manualHighlights.innerHTML = "";
  [
    metricCard("Total SSQ-10 Score", result.SSQ10_Total_Score ?? DASH),
    metricCard("Harmonized Theta", fmtNumber(result.Theta_Harmonized)),
    metricCard("Depression", result.Depression_Binary === 1 ? "Yes" : result.Depression_Binary === 0 ? "No" : DASH),
    metricCard("Expected PHQ-9 Sum", fmtNumber(result.PHQ_Expected_Sum)),
    metricCard("Probability PHQ-9 >= 10", fmtPercent(result.PHQ_Prob_GE10)),
    metricCard("Cutoff Target", prettyGroupLevel(result.Cutoff_Group_Level_Used))
  ].forEach((card) => dom.manualHighlights.appendChild(card));

  setBlockMessage(dom.manualStatus, "Manual scoring completed in the browser for the current session.", "status");
  setTable(dom.manualResultTable, [result], previewColumns([result]));
}

function renderManualCompareTable(record) {
  const scoredRows = state.bundle.engines.map((engine) => scoreRecord(record, engine));
  setTable(
    dom.manualEngineCompareTable,
    scoredRows,
    [
      { key: "Engine_Label", label: "Engine" },
      { key: "Theta_Harmonized", label: "Harmonized Theta", format: (value) => fmtNumber(value) },
      {
        key: "Depression_Binary",
        label: "Depression",
        format: (value) => (value === 1 ? "Yes" : value === 0 ? "No" : DASH)
      },
      { key: "PHQ_Expected_Sum", label: "Expected PHQ-9 Sum", format: (value) => fmtNumber(value) },
      { key: "PHQ_Prob_GE10", label: "Probability PHQ-9 >= 10", format: (value) => fmtPercent(value) },
      { key: "Cutoff_Theta_Applied", label: "Applied Cutoff Theta", format: (value) => fmtNumber(value) },
      {
        key: "Cutoff_Group_Level_Used",
        label: "Cutoff Target",
        format: (value) => prettyGroupLevel(value)
      }
    ]
  );
  setBlockMessage(
    dom.manualCompareNote,
    "The current profile has been rescored across all production-safe browser engines.",
    "status"
  );
}

function validateBatchRows(rows, engine) {
  const issues = [];
  const warnings = [];

  if (!rows || rows.length === 0) {
    issues.push("The uploaded file does not contain any data rows.");
    return { valid: false, issues, warnings };
  }

  const headers = Object.keys(rows[0]);
  const missing = engine.items.filter((item) => !headers.includes(item));
  if (missing.length > 0) {
    issues.push(`Missing required SSQ-10 item columns: ${missing.join(", ")}`);
  }

  if (issues.length === 0) {
    const badCells = [];
    rows.forEach((row, rowIndex) => {
      engine.items.forEach((item) => {
        const raw = row[item];
        if (raw == null || String(raw).trim() === "") return;
        const normalized = normalizeBinaryResponse(raw);
        if (!(normalized === 0 || normalized === 1)) {
          badCells.push(`${item} row ${rowIndex + 1}`);
        }
      });
    });

    if (badCells.length > 0) {
      issues.push(`Values outside the supported 0/1 or Yes/No coding were found in: ${badCells.slice(0, 10).join(", ")}`);
    }
  }

  const allMissingRows = rows
    .map((row, index) => ({
      index: index + 1,
      answered: engine.items.reduce((acc, item) => acc + ((normalizeBinaryResponse(row[item]) === 0 || normalizeBinaryResponse(row[item]) === 1) ? 1 : 0), 0)
    }))
    .filter((entry) => entry.answered === 0)
    .map((entry) => entry.index);

  if (allMissingRows.length > 0) {
    warnings.push(`Rows with all SSQ-10 items missing: ${allMissingRows.slice(0, 10).join(", ")}`);
  }

  return {
    valid: issues.length === 0,
    issues,
    warnings
  };
}

function renderBatchResult(scored, validationMessage) {
  dom.batchHighlights.innerHTML = "";
  [
    metricCard("Rows scored", String(scored.length)),
    metricCard("Mean theta", fmtNumber(mean(scored.map((row) => row.Theta_Harmonized)))),
    metricCard("Depression prevalence", fmtPercent(prevalence(scored.map((row) => row.Depression_Binary)))),
    metricCard("Engine", state.engine.label, prettyMode(state.engine.cutoff_apply_mode)),
    metricCard("Threshold scope", prettyMode(state.engine.cutoff_apply_mode)),
    metricCard("Threshold approach", prettyMethod(state.engine.cutoff_method))
  ].forEach((card) => dom.batchHighlights.appendChild(card));

  setBlockMessage(dom.batchStatus, validationMessage, "status");
  setTable(dom.batchPreviewTable, compactPreview(scored), previewColumns(scored));
}

function updateEngineSelectors(engineId) {
  dom.manualEngineSelect.value = engineId;
  dom.batchEngineSelect.value = engineId;
}

function selectEngine(engineId, { reRenderItems = true, rescore = true } = {}) {
  const engine = state.bundle.engines.find((entry) => entry.engine_id === engineId);
  if (!engine) return;

  state.engine = engine;
  updateEngineSelectors(engineId);
  renderOverview(engine);
  renderSchemaTables(engine);

  if (reRenderItems) {
    renderManualItems(engine, state.latestManualInput || {});
  }

  if (rescore && state.latestManualInput) {
    const rescored = scoreRecord(state.latestManualInput, engine);
    renderManualResult(rescored);
    renderManualCompareTable(state.latestManualInput);
    setLatestRun("manual", [rescored], engine);
  }

  if (rescore && state.latestBatchSourceRows) {
    const rescoredBatch = scoreBatch(state.latestBatchSourceRows, engine);
    state.latestBatch = rescoredBatch;
    renderBatchResult(
      rescoredBatch,
      `Validation passed for ${rescoredBatch.length} row(s). Data were scored in the browser only. Changing the engine rescored the active batch for this session.`
    );
    setLatestRun("batch", rescoredBatch, engine);
  }
}

async function scoreBatchFile() {
  const file = dom.batchFile.files?.[0];
  if (!file) {
    setBlockMessage(dom.batchStatus, "Choose a CSV file first.", "status");
    return;
  }

  const text = await file.text();
  const parsed = parseCsv(text);
  const validation = validateBatchRows(parsed, state.engine);

  dom.batchHighlights.innerHTML = "";
  clearTable(dom.batchPreviewTable);
  dom.downloadBatchBtn.disabled = true;

  if (!validation.valid) {
    setBlockMessage(dom.batchStatus, validation.issues.join(" "), "status");
    return;
  }

  state.latestBatchSourceRows = parsed;
  state.latestBatch = scoreBatch(parsed, state.engine);

  const pieces = [`Validation passed for ${state.latestBatch.length} row(s). Data were scored in the browser only.`];
  if (validation.warnings.length > 0) {
    pieces.push(validation.warnings.join(" "));
  }
  renderBatchResult(state.latestBatch, pieces.join(" "));
  dom.downloadBatchBtn.disabled = state.latestBatch.length === 0;
  setLatestRun("batch", state.latestBatch, state.engine);
}

function downloadText(filename, text) {
  const blob = new Blob([text], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

function wireEvents() {
  initToolHelp();

  dom.tabButtons.forEach((btn) => {
    btn.addEventListener("click", () => selectTab(btn.dataset.tabTarget));
  });

  dom.gotoButtons.forEach((btn) => {
    btn.addEventListener("click", () => selectTab(btn.dataset.goto));
  });

  dom.openBatchBtn.addEventListener("click", () => selectTab("batch"));

  dom.manualEngineSelect.addEventListener("change", (event) => {
    selectEngine(event.target.value);
  });

  dom.batchEngineSelect.addEventListener("change", (event) => {
    selectEngine(event.target.value);
  });

  dom.scoreManualBtn.addEventListener("click", () => {
    const record = getManualRecord();
    state.latestManualInput = { ...record };
    const result = scoreRecord(record, state.engine);

    if (result.Items_Answered === 0) {
      dom.manualHighlights.innerHTML = "";
      clearTable(dom.manualResultTable);
      clearTable(dom.manualEngineCompareTable);
      setBlockMessage(dom.manualStatus, "At least one SSQ-10 item is required for scoring.", "status");
      setBlockMessage(dom.manualCompareNote, "Score a participant to compare the current profile across engines.", "empty");
      return;
    }

    renderManualResult(result);
    renderManualCompareTable(record);
    setLatestRun("manual", [result], state.engine);
  });

  dom.resetManualBtn.addEventListener("click", () => {
    state.latestManualInput = null;
    dom.manualSex.value = "";
    dom.manualAge.value = "";
    renderManualItems(state.engine);
    dom.manualHighlights.innerHTML = "";
    clearTable(dom.manualResultTable);
    clearTable(dom.manualEngineCompareTable);
    setBlockMessage(dom.manualStatus, "Score a participant to see the harmonised outputs.", "empty");
    setBlockMessage(dom.manualCompareNote, "Score a participant to compare the current profile across engines.", "empty");
  });

  dom.scoreBatchBtn.addEventListener("click", () => {
    scoreBatchFile().catch((error) => {
      setBlockMessage(dom.batchStatus, `Batch scoring failed: ${error.message}`, "status");
    });
  });

  dom.downloadBatchBtn.addEventListener("click", () => {
    if (!state.latestBatch) return;
    downloadText("portable_irt_browser_scores.csv", toCsv(state.latestBatch));
  });

  dom.downloadExampleBtn.addEventListener("click", async () => {
    const response = await fetch(exampleCsvPath);
    const text = await response.text();
    downloadText("portable_irt_example_batch.csv", text);
  });
}

function populateEngineSelectors(bundle) {
  [dom.manualEngineSelect, dom.batchEngineSelect].forEach((selectEl) => {
    selectEl.innerHTML = "";
    bundle.engines.forEach((engine) => {
      const option = document.createElement("option");
      option.value = engine.engine_id;
      option.textContent = engine.label;
      selectEl.appendChild(option);
    });
  });
}

async function init() {
  const response = await fetch(bundlePath);
  if (!response.ok) {
    throw new Error(`Unable to load ${bundlePath}`);
  }

  state.bundle = await response.json();
  populateEngineSelectors(state.bundle);
  renderEngineTable();
  wireEvents();

  const defaultEngine = state.bundle.engines.find((entry) => entry.engine_id === state.bundle.default_engine_id)
    || state.bundle.engines[0];

  selectEngine(defaultEngine.engine_id, { reRenderItems: true, rescore: false });
  setBlockMessage(dom.manualStatus, "Score a participant to see the harmonised outputs.", "empty");
  setBlockMessage(dom.manualCompareNote, "Score a participant to compare the current profile across engines.", "empty");
  setBlockMessage(dom.batchStatus, "No file scored yet.", "empty");
  refreshLatestRunView();
  selectTab("overview");
}

init().catch((error) => {
  dom.bundleStatus.textContent = `Failed to load engine bundle: ${error.message}`;
});
