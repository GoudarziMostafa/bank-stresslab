"use strict";

const el = (id) => document.getElementById(id);
const format = (n, digits = 1) => new Intl.NumberFormat("en-US", {
  minimumFractionDigits: digits, maximumFractionDigits: digits,
}).format(n).replace("-", "−");
const bp = (n) => `${format(n)} bp`;
const pct = (n, digits = 2) => `${format(n * 100, digits)}%`;
const signed = (n, digits = 1) => (n > 0 ? "+" : "") + format(n, digits);
const names = { synthetic: "Synthetic baseline", basel: "Basel benchmark", gems: "GEMs private lending" };
const inputFields = [
  ["pd", "Annual PD", pct],
  ["lgd", "Fixed LGD", pct],
  ["rho_gaussian", "Gaussian ρ", (n) => format(n, 3)],
  ["df", "Student-t ν", (n) => format(n, 0)],
];
const metric = (label, value) => `<div><dt>${label}</dt><dd>${value}</dd></div>`;
const badge = (feasible, text) => `<span class="status${feasible ? "" : " blocked"}">${text}</span>`;
const interval = (lo, hi, digits) => `${format(lo, digits)} to ${format(hi, digits)}`;
const pointLabel = (df) => `ν = ${format(df, 1)}`;

function segments(id, labels, active, onChange) {
  const group = el(id);
  const buttons = Object.entries(labels).map(([key, label]) => {
    const button = document.createElement("button");
    button.type = "button";
    button.textContent = label;
    button.setAttribute("aria-pressed", String(key === active));
    button.addEventListener("click", () => {
      for (const sibling of group.children) sibling.setAttribute("aria-pressed", String(sibling === button));
      onChange(key);
    });
    return button;
  });
  group.replaceChildren(...buttons);
}

function plot(canvas, limits, ticks, band = null) {
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  if (!width || !height) return null;
  const ratio = window.devicePixelRatio || 1;
  canvas.width = Math.round(width * ratio);
  canvas.height = Math.round(height * ratio);
  const ctx = canvas.getContext("2d");
  ctx.scale(ratio, ratio);

  const [left, top, right, bottom] = [58, 12, width - 16, height - 34];
  const x = (n) => left + (right - left) * (n - limits.x[0]) / (limits.x[1] - limits.x[0]);
  const y = (n) => bottom - (bottom - top) * (n - limits.y[0]) / (limits.y[1] - limits.y[0]);
  const line = (points, color, weight = 1, dash = []) => {
    ctx.strokeStyle = color;
    ctx.lineWidth = weight;
    ctx.setLineDash(dash);
    ctx.beginPath();
    ctx.moveTo(x(points[0][0]), y(points[0][1]));
    for (const [a, b] of points.slice(1)) ctx.lineTo(x(a), y(b));
    ctx.stroke();
    ctx.setLineDash([]);
  };
  const dot = (a, b, radius, color) => {
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(x(a), y(b), radius, 0, 2 * Math.PI);
    ctx.fill();
  };
  if (band) {
    ctx.fillStyle = "#e4edfe";
    ctx.fillRect(x(band[0]), top, x(band[1]) - x(band[0]), bottom - top);
  }
  ctx.font = "12px system-ui";
  ctx.textBaseline = "middle";
  ctx.textAlign = "right";
  for (const [value, label] of ticks.y) {
    line([[limits.x[0], value], [limits.x[1], value]], "#dbe5ed");
    ctx.fillStyle = "#536679";
    ctx.fillText(label, left - 7, y(value));
  }
  ctx.textAlign = "center";
  for (const [value, label] of ticks.x) {
    ctx.fillStyle = "#536679";
    ctx.fillText(label, x(value), bottom + 15);
  }
  line([[limits.x[0], limits.y[0]], [limits.x[1], limits.y[0]]], "#9bafbf");
  return { line, dot };
}

function floorChart(canvas, scenario, selected, crossing) {
  const log = Math.log10;
  const maximum = scenario.curve.at(-1)[0];
  const upper = Math.ceil(scenario.curve[0][1] * 100 / .05) * .0005;
  const xTicks = maximum > 200 ? [4, 10, 25, 60, 150, maximum] : [4, 10, 20, 50, 100, 200];
  const band = crossing && [log(crossing.lower_95 ?? maximum), log(crossing.upper_95 ?? maximum)];
  const chart = plot(canvas, { x: [log(4), log(maximum)], y: [0, upper] }, {
    x: xTicks.map((v) => [log(v), v]),
    y: Array.from({ length: 5 }, (_, i) => [i * upper / 4, pct(i * upper / 4, 3)]),
  }, band);
  if (!chart) return;
  chart.line([[log(4), scenario.target], [log(maximum), scenario.target]], "#315ac0", 2, [5, 4]);
  chart.line(scenario.curve.map(([df, floor]) => [log(df), floor]), "#087b72", 2.5);
  chart.dot(log(scenario.crossing), scenario.target, 4, "#315ac0");
  chart.dot(log(selected.df), selected.floor, 5, "#a44b16");
}

function bootstrapChart(bootstrap) {
  const bounds = [0, .0045];
  const ticks = [0, .001, .002, .003, .004].map((n) => [n, pct(n, 2)]);
  const chart = plot(el("bootstrap-chart"), { x: bounds, y: bounds }, { x: ticks, y: ticks });
  if (!chart) return;
  chart.line([[0, 0], [.0045, .0045]], "#667c91", 1.2, [4, 4]);
  for (const [pair, floor] of bootstrap.points) chart.dot(pair, floor, 2, "rgba(0, 107, 101, .28)");
}

function renderOverview(data) {
  const { gems } = data.anchors;
  el("headline-note").textContent =
    `GEMs calibration uses ${gems.last_year - gems.first_year + 1} annual observations. Both books use one credit factor and fixed loss given default (LGD). Student-t adds a shared scale. Losses are measured against total exposure.`;

  el("overview").innerHTML = Object.entries(data.anchors).map(([key, row]) => {
    const relative = pct(Math.abs(row.gap_bp) / row.gaussian_es_bp, 1);
    const examples = data.gems_range.map((example) => `<div>
      <small>ν = ${example.df}</small><strong>${bp(Math.abs(example.gap_bp))}</strong></div>`).join("");
    const result = key === "gems"
      ? `<p>Illustrative Student-t ES gaps below Gaussian</p>
         <div class="examples">${examples}</div>
         <p class="note">At ν = ${row.df}, the gap is ${relative} of Gaussian ES.</p>`
      : `<strong class="big-number">${format(Math.abs(row.gap_bp), 0)} bp</strong>
         <p>Student-t 99.9% ES is ${format(Math.abs(row.gap_bp), 0)} bp below Gaussian.</p>
         <p class="note">Matched at ν = ${row.df} · ${relative} of Gaussian ES</p>`;
    return `<article class="metric-card ${key}">
      <h2>${names[key]}</h2>
      ${result}
      <dl class="metrics">${inputFields.map(([field, label, display]) =>
        metric(key === "gems" && field === "df" ? "Example ν" : label, display(row[field]))).join("")}</dl>
    </article>`;
  }).join("");
}

function crossingInterval(crossing) {
  if (crossing.lower_95 === null) return "no finite crossing";
  const upper = crossing.upper_95 === null ? "no finite crossing" : format(crossing.upper_95, 0);
  return `${format(crossing.lower_95, 0)} to ${upper}`;
}

function setupFloor(data) {
  let activeCase = "gems";
  const selector = el("floor-df");
  const presets = data.floors.gems.presets;
  const dfs = presets.map(({ df }) => df).join(",");
  if (Object.values(data.floors).some(({ presets }) => presets.map(({ df }) => df).join(",") !== dfs)) {
    throw new Error("Floor scenarios must share the same degrees-of-freedom presets");
  }
  const crossing = data.bootstrap.crossing;
  selector.replaceChildren(...presets.map(
    ({ df }) => new Option(`ν = ${df}`, df),
  ));

  function selected() {
    const scenario = data.floors[activeCase];
    return [scenario, scenario.presets.find(({ df }) => df === Number(selector.value))];
  }
  function draw() {
    floorChart(el("floor-chart"), ...selected(), activeCase === "gems" ? crossing : null);
  }
  function render() {
    const [scenario, choice] = selected();
    const feasible = choice.floor <= scenario.target;
    const showInterval = activeCase === "gems";
    const uncertainty = showInterval ? `<p class="note crossing-note">Approximate 95% interval for the crossing: <strong>${crossingInterval(crossing)}</strong>. The band shows how the crossing moves as PD and pair risk are resampled together.</p>` : "";
    el("crossing-legend").hidden = !showInterval;
    el("floor-readout").innerHTML = `${badge(feasible, feasible ? "Feasible match" : "No feasible match")}
      <dl class="readout">
        ${metric("Annual PD", pct(scenario.pd))}
        ${metric("Gaussian ρ", format(scenario.rho, 3))}
        ${metric("Pair target", pct(scenario.target, 4))}
        ${metric("Student-t floor", pct(choice.floor, 4))}
        ${metric("Point crossing at ν", format(scenario.crossing))}
      </dl>
      <p class="note">${feasible ? "The floor is below the target." : "The floor exceeds the target."}</p>${uncertainty}`;
    el("floor-chart").setAttribute("aria-label", `${names[activeCase]}. The minimum pair-default probability crosses the target at ν = ${format(scenario.crossing)}. The selected ν is ${choice.df}.${showInterval ? ` Approximate 95% crossing interval ${crossingInterval(crossing)}.` : ""}`);
    draw();
  }
  segments("floor-cases", names, activeCase, (value) => {
    activeCase = value;
    render();
  });
  selector.addEventListener("change", render);
  const output = data.gems_range.map((row) => `<tr><td>${row.df}</td><td>${bp(row.gaussian_es_bp)}</td>
    <td>${bp(row.student_es_bp)}</td><td>${signed(row.gap_bp)} bp</td></tr>`).join("");
  el("gems-range").innerHTML = `<p class="note">Illustrative choices at the fitted GEMs PD, LGD and Gaussian correlation.
    At the point crossing (${pointLabel(data.floors.gems.crossing)}), the ES gap is ${signed(data.gems_boundary_gap_bp)} bp.</p>
    <div class="table-box"><table><thead><tr><th scope="col">ν</th><th scope="col">Gaussian ES</th>
    <th scope="col">Student-t ES</th><th scope="col">Difference</th></tr></thead><tbody>${output}</tbody></table></div>`;
  render();
  new ResizeObserver(draw).observe(el("floor-chart").parentElement);
}

function setupSwaps(data) {
  const controls = el("swap-controls");
  controls.innerHTML = inputFields.map(([field, label, display]) => `<fieldset>
    <legend>${label}</legend>${data.choices[field].map((value, i) => `<label class="choice">
      <input type="radio" name="${field}" value="${i}" ${i === 0 ? "checked" : ""}>
      <span>${i === 0 ? "Synthetic" : "GEMs"}<strong>${display(value)}</strong></span></label>`).join("")}
    </fieldset>`).join("");
  const mode = el("swap-mode");
  function render() {
    const key = [...controls.querySelectorAll("input:checked")].map((input) => input.value).join("");
    const row = data.swaps[key];
    const matched = mode.value === "matched";
    const feasible = row.status === "matched";
    const student = matched ? row.student_es_bp : row.fixed_rho_es_bp;
    const rho = matched ? row.rho_student : row.rho_gaussian;
    const delta = student === null ? "—" : `${signed(student - row.gaussian_es_bp)} bp`;
    const diagnostics = matched
      ? `Gaussian pair target ${pct(row.pair_target, 4)} · Student-t floor ${pct(row.pair_floor, 4)}`
      : `Gaussian pair probability ${pct(row.pair_target, 4)} · Student-t pair probability ${pct(row.fixed_rho_pair, 4)}`;
    el("swap-results").innerHTML = `${badge(!matched || feasible,
      matched ? (feasible ? "Pair probabilities matched" : "No feasible match") : "Same latent correlation")}
      <dl class="swap-metrics">
        ${metric("Gaussian 99.9% ES", bp(row.gaussian_es_bp))}
        ${metric("Student-t 99.9% ES", student === null ? "—" : bp(student))}
        ${metric("Student-t minus Gaussian", delta)}
      </dl>
      <p class="note">${diagnostics} ${rho === null ? "" : `· Student-t ρ ${format(rho, 4)}`}</p>`;
  }
  controls.addEventListener("change", render);
  mode.addEventListener("change", render);
  el("swap-presets").addEventListener("click", (event) => {
    const key = event.target.dataset.key;
    if (!key) return;
    inputFields.forEach(([field], index) => {
      controls.querySelector(`input[name="${field}"][value="${key[index]}"]`).checked = true;
    });
    render();
  });
  const rows = [["Synthetic baseline", data.swaps["0000"]], ...inputFields.map(([, label], i) =>
    [`GEMs ${label}`, data.swaps[inputFields.map((_, j) => i === j ? "1" : "0").join("")]])];
  el("one-input-results").innerHTML = `<table><thead><tr><th scope="col">Input changed</th>
    <th scope="col">Student-t minus Gaussian ES</th></tr></thead>
    <tbody>${rows.map(([label, row]) => `<tr><td>${label}</td>
      <td>${row.gap_bp === null ? "No feasible match" : `${signed(row.gap_bp)} bp`}</td></tr>`).join("")}</tbody></table>`;
  const [syntheticLgd, gemsLgd] = data.choices.lgd;
  el("lgd-note").textContent = `Swapping fixed LGD alone multiplies both loss distributions by ${format(gemsLgd / syntheticLgd, 3)}. The pair probabilities and the relative ES gap stay the same.`;
  render();
}

function setupConcentration(data) {
  const labels = { var_bp: "Value at Risk", es_bp: "Expected shortfall" };
  function render(field) {
    el("concentration-results").innerHTML = ["synthetic", "gems"].map((key) => {
      const rows = data.concentration.filter((row) => row.case === key);
      const values = ["Gaussian", "Student-t"].map((model) => {
        const [equal, concentrated] = rows.filter((row) => row.model === model);
        return `<tr><td>${model}</td><td>${bp(equal[field])}</td><td>${bp(concentrated[field])}</td>
          <td class="up">${signed(concentrated[field] - equal[field])} bp</td></tr>`;
      });
      return `<article><h3>${names[key]}</h3>
        <p class="note">99.9% ${labels[field]} · ν = ${data.anchors[key].df}${key === "gems" ? " (illustrative)" : ""}</p>
        <div class="table-box"><table><thead><tr><th scope="col">Model</th><th scope="col">Equal</th>
          <th scope="col">10% name</th><th scope="col">Change</th></tr></thead>
          <tbody>${values.join("")}</tbody></table></div></article>`;
    }).join("");
  }
  const [gaussian, student] = data.concentration.filter((row) =>
    row.case === "gems" && row.allocation === "concentrated");
  el("concentration-tail").textContent =
    `At GEMs inputs, the large name defaults in ${pct(gaussian.large_name_tail_share, 3)} of the worst 0.1% of Gaussian losses and ${pct(student.large_name_tail_share, 3)} of the matched Student-t tail. Both VaRs land on the same attainable loss, ${bp(gaussian.var_bp)}. Above the cutoff, loss probabilities differ and ES still separates the models.`;
  segments("concentration-measures", labels, "var_bp", render);
  render("var_bp");
}

function setupEvidence(data) {
  const { bootstrap, simulation } = data;
  const difference = bootstrap.intervals.student_floor_minus_pair;
  const pair = bootstrap.intervals.pair_probability;
  const excess = bootstrap.crossing.excess_pair;
  const range = (row) => interval(row.lower_95 * 100, row.upper_95 * 100, 3);
  el("bootstrap-summary").innerHTML = `<p><strong>GEMs pair-default probability</strong><br>Approximate 95% interval ${range(pair)}%.</p>
    <p><strong>Student-t floor minus pair estimate</strong><br>Approximate 95% interval ${range(difference)} percentage points at ν = 10.</p>
    <p><strong>GEMs pair probability minus PD²</strong><br>Approximate 95% interval ${interval(excess.lower_95 * 100, excess.upper_95 * 100, 4)} percentage points.</p>
    <p>The floor and the pair estimate move together across resamples (correlation ${format(bootstrap.correlation, 3)}). The chart displays a subset of the ${format(difference.replicates, 0)} draws and the line of equality.</p>
    <p><strong>Synthetic ES difference</strong><br>95% Monte Carlo interval ${interval(simulation.low95_bp, simulation.high95_bp, 1)} bp.</p>
    <p><a href="../reports/sensitivity.md#gems-private-lending">Bootstrap design</a> · <a href="../reports/matched_dependence.md#numerical-evidence">Simulation design</a></p>`;
  new ResizeObserver(() => bootstrapChart(bootstrap)).observe(el("bootstrap-chart").parentElement);
}

async function start() {
  try {
    const response = await fetch("data.json");
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    const data = await response.json();
    el("dashboard").hidden = false;
    renderOverview(data);
    setupFloor(data);
    setupSwaps(data);
    setupConcentration(data);
    setupEvidence(data);
    el("load-status").hidden = true;
  } catch (error) {
    console.error(error);
    el("dashboard").hidden = true;
    el("load-status").textContent = "The dashboard data could not load. Serve the repository with a local HTTP server and open /dashboard/.";
  }
}

start();
