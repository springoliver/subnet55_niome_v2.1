/**
 * NIOME multi-miner PM2 — pang8512 wallet
 *
 *  miner-1  miner-pang     → v10
 *  miner-2  miner-sn55     → high_recall
 *  miner-3  miner-pang01   → v5_style
 *  miner-4  miner-pang02   → high_recall
 *  miner-5  miner-pang05   → win (truth when available)
 *
 * Deploy:
 *   cd /path/to/subnet-niome && git pull
 *   bash scripts/pm2-fleet-deploy.sh
 *
 * Spring wallet fleet: ecosystem.spring.config.js
 */

const fs = require("fs");
const path = require("path");

const CONFIG_DIR = __dirname;

function isPlaceholderRepo(p) {
  if (!p) return true;
  const n = p.replace(/\\/g, "/").toLowerCase();
  return n.includes("path/to") || n.includes("<repo>");
}

function resolveRepo() {
  const env = (process.env.NIOME_REPO || "").trim();
  if (env && !isPlaceholderRepo(env)) {
    const resolved = path.resolve(env);
    if (fs.existsSync(resolved)) return resolved;
    console.warn(
      `[ecosystem.multi] NIOME_REPO=${env} not found — using ${CONFIG_DIR}`
    );
  } else if (env && isPlaceholderRepo(env)) {
    console.warn(
      `[ecosystem.multi] Ignoring placeholder NIOME_REPO=${env} — using ${CONFIG_DIR}`
    );
  }
  return CONFIG_DIR;
}

const REPO = resolveRepo();
const PYTHON =
  process.env.NIOME_PYTHON || path.join(REPO, "venv", "bin", "python");
const RESULTS_ROOT =
  process.env.NIOME_RESULTS_ROOT || path.join(REPO, "Results");
const TRUTH_DIR = process.env.NIOME_TRUTH_DIR || RESULTS_ROOT;

if (!fs.existsSync(PYTHON)) {
  console.error(
    `[ecosystem.multi] Missing venv python: ${PYTHON}\n` +
      `  cd ${REPO} && python3 -m venv venv && ./venv/bin/pip install -e .`
  );
  process.exit(1);
}

console.log(`[ecosystem.multi] REPO=${REPO} PYTHON=${PYTHON}`);

const BASE_ENV = {
  NIOME_REPO: REPO,
  NIOME_USE_HG38: "1",
  NIOME_RESULTS_ROOT: RESULTS_ROOT,
  NIOME_TRUTH_DIR: TRUTH_DIR,
  NIOME_DISABLE_PIPELINE_RETRY: "0",
};

function minerApp({ name, hotkey, port, strategy, extraEnv = {} }) {
  return {
    name,
    script: PYTHON,
    args: [
      "-m",
      "neurons.miner",
      "--netuid",
      "55",
      "--subtensor.network",
      "finney",
      "--wallet.name",
      "pang8512",
      "--wallet.hotkey",
      hotkey,
      "--axon.port",
      String(port),
      "--logging.debug",
    ],
    cwd: REPO,
    interpreter: "none",
    autorestart: true,
    restart_delay: 5000,
    kill_timeout: 15000,
    max_restarts: 50,
    min_uptime: 10000,
    env: {
      ...BASE_ENV,
      NIOME_STRATEGY: strategy,
      ...extraEnv,
    },
  };
}

module.exports = {
  apps: [
    minerApp({
      name: "miner-1",
      hotkey: "miner-pang",
      port: 8091,
      strategy: "v10",
    }),
    minerApp({
      name: "miner-2",
      hotkey: "miner-sn55",
      port: 8092,
      strategy: "high_recall",
    }),
    minerApp({
      name: "miner-3",
      hotkey: "miner-pang01",
      port: 8093,
      strategy: "v5_style",
    }),
    minerApp({
      name: "miner-4",
      hotkey: "miner-pang02",
      port: 8094,
      strategy: "high_recall",
    }),
    minerApp({
      name: "miner-5",
      hotkey: "miner-pang05",
      port: 8095,
      strategy: "win",
      extraEnv: {
        NIOME_WIN_MODE: "1",
        NIOME_VCF_MINIMAL: "1",
        NIOME_VCF_DOT_ID: "1",
      },
    }),
  ],
};
