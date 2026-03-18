/**
 * NIOME spring wallet fleet (Vultr) — use: pm2 start ecosystem.spring.config.js
 *
 *  springhot → 141 (auto→v10)  spring01 → 50 (high_recall)  spring02 → 99 (v10)
 *  spring03  → 209 (high_recall)  spring04 → 235 (win)  spring05 → 226 (win)
 *  spring06  → 124 (win)  spring07 → 9 (v10)  spring08 → 97 (win)
 *  spring09  → 217 (v10)
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
  } else if (env && isPlaceholderRepo(env)) {
    console.warn(
      `[ecosystem.spring] Ignoring placeholder NIOME_REPO=${env} — using ${CONFIG_DIR}`
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
  console.error(`[ecosystem.spring] Missing venv: ${PYTHON}`);
  process.exit(1);
}

console.log(`[ecosystem.spring] REPO=${REPO} PYTHON=${PYTHON}`);

const BASE_ENV = {
  NIOME_REPO: REPO,
  NIOME_USE_HG38: "1",
  NIOME_RESULTS_ROOT: RESULTS_ROOT,
  NIOME_TRUTH_DIR: TRUTH_DIR,
  NIOME_DISABLE_PIPELINE_RETRY: "0",
};

function minerApp({
  name,
  hotkey,
  uid,
  port,
  strategy,
  extraEnv = {},
  forceValidatorPermit = false,
}) {
  const args = [
    "-m",
    "neurons.miner",
    "--netuid",
    "55",
    "--subtensor.network",
    "finney",
    "--wallet.name",
    "spring",
    "--wallet.hotkey",
    hotkey,
    "--axon.port",
    String(port),
    "--logging.debug",
  ];
  if (forceValidatorPermit) {
    args.push("--blacklist.force_validator_permit");
  }

  const uidStrategy = {};
  uidStrategy[`NIOME_UID_STRATEGY_${uid}`] = strategy;

  return {
    name,
    script: PYTHON,
    args,
    cwd: REPO,
    interpreter: "none",
    autorestart: true,
    restart_delay: 5000,
    kill_timeout: 15000,
    env: {
      ...BASE_ENV,
      NIOME_STRATEGY: strategy,
      ...uidStrategy,
      ...extraEnv,
    },
  };
}

module.exports = {
  apps: [
    minerApp({
      name: "springhot",
      hotkey: "springhot",
      uid: 141,
      port: 8105,
      strategy: "auto",
      forceValidatorPermit: true,
    }),
    // minerApp({
    //   name: "spring01",
    //   hotkey: "spring01",
    //   uid: 50,
    //   port: 8101,
    //   strategy: "high_recall",
    // }),
    minerApp({
      name: "spring02",
      hotkey: "spring02",
      uid: 99,
      port: 8112,
      strategy: "v10",
      forceValidatorPermit: true,
    }),
    // minerApp({
    //   name: "spring03",
    //   hotkey: "spring03",
    //   uid: 209,
    //   port: 8103,
    //   strategy: "high_recall",
    //   forceValidatorPermit: true,
    // }),
    // minerApp({
    //   name: "spring04",
    //   hotkey: "spring04",
    //   uid: 235,
    //   port: 8104,
    //   strategy: "win",
    //   forceValidatorPermit: true,
    // }),
    minerApp({
      name: "spring05",
      hotkey: "spring05",
      uid: 226,
      port: 8107,
      strategy: "win",
      extraEnv: {
        NIOME_WIN_MODE: "1",
        NIOME_VCF_MINIMAL: "1",
        NIOME_VCF_DOT_ID: "1",
      },
    }),
    minerApp({
      name: "spring06",
      hotkey: "spring06",
      uid: 124,
      port: 8110,
      strategy: "win",
      forceValidatorPermit: true,
      extraEnv: {
        NIOME_WIN_MODE: "1",
        NIOME_VCF_MINIMAL: "1",
        NIOME_VCF_DOT_ID: "1",
      },
    }),
    minerApp({
      name: "spring07",
      hotkey: "spring07",
      uid: 9,
      port: 8111,
      strategy: "v10",
    }),
    minerApp({
      name: "spring08",
      hotkey: "spring08",
      uid: 97,
      port: 8108,
      strategy: "win",
      forceValidatorPermit: true,
      extraEnv: {
        NIOME_WIN_MODE: "1",
        NIOME_VCF_MINIMAL: "1",
        NIOME_VCF_DOT_ID: "1",
      },
    }),
    minerApp({
      name: "spring09",
      hotkey: "spring09",
      uid: 217,
      port: 8109,
      strategy: "high_recall",
    }),
  ],
};
