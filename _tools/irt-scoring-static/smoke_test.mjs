import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { parseCsv, scoreBatch } from "./scorer.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

async function main() {
  const bundlePath = process.argv[2]
    ? path.resolve(process.argv[2])
    : path.join(__dirname, "engines", "production_engines.json");
  const csvPath = process.argv[3]
    ? path.resolve(process.argv[3])
    : path.join(__dirname, "engines", "example_batch.csv");
  const engineId = process.argv[4] || null;

  const bundle = JSON.parse(await fs.readFile(bundlePath, "utf8"));
  const csvText = await fs.readFile(csvPath, "utf8");
  const rows = parseCsv(csvText);

  const engine = bundle.engines.find((entry) => entry.engine_id === (engineId || bundle.default_engine_id));
  if (!engine) {
    throw new Error(`Engine not found: ${engineId || bundle.default_engine_id}`);
  }

  const scored = scoreBatch(rows, engine);
  process.stdout.write(JSON.stringify({
    engine_id: engine.engine_id,
    n_rows: scored.length,
    scored
  }));
}

main().catch((error) => {
  console.error(error.message);
  process.exit(1);
});
