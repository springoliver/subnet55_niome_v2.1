/**
 * NIOME fleet PM2 config — one codebase, multiple strategies per process.
 *
 * Strategies (NIOME_STRATEGY or task_strategy.py profiles):
 *   win         — truth VCF when NIOME_TRUTH_DIR / Results has task_id
 *   v5_style    — conservative GT + simpler indels + ~24 sites (UID 71 pattern)
 *   high_recall — push 30–32 sites (ultra/high count bands)
 *   v10         — balanced native default
 *   auto        — truth if available, else band heuristic per task
 *
 * After deploy, logs must show:
 *   [strategy] name=v5_style ...
 *   rev=niome-native-2026-05-23-v10+fleet-v5
 *
 * Optional per-UID override (once you know UID from metagraph):
 *   NIOME_UID_STRATEGY_71: "v5_style"
 */

// --- Edit these paths on the server ----------------------------------------
const REPO = process.env.NIOME_REPO || "/root/subnet-niome";
const RESULTS_ROOT = process.env.NIOME_RESULTS_ROOT || `${REPO}/Results`;
const TRUTH_DIR = process.env.NIOME_TRUTH_DIR || RESULTS_ROOT;

// Shared env for every miner (bcftools 1.16+ required on PATH)
const BASE_ENV = {
  NIOME_USE_HG38: "1",
  NIOME_RESULTS_ROOT: RESULTS_ROOT,
  NIOME_TRUTH_DIR: TRUTH_DIR,
  NIOME_DISABLE_PIPELINE_RETRY: "0",
};

function minerApp({
  name,
  hotkey,
  port,
  strategy,
  extraEnv = {},
  forceValidatorPermit = false,
  debug = true,
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
  ];
  if (forceValidatorPermit) {
    args.push("--blacklist.force_validator_permit");
  }
  if (debug) {
    args.push("--logging.debug");
  }

  return {
    name,
    script: "python",
    args,
    cwd: REPO,
    interpreter: "none",
    autorestart: true,
    restart_delay: 5000,
    env: {
      ...BASE_ENV,
      NIOME_STRATEGY: strategy,
      ...extraEnv,
    },
  };
}

module.exports = {
  apps: [
    // --- Portfolio: cover bands in one round (not by region — all ultra-wide) ---

    // v5_style ×3 — annotation + simpler indels (best match to legacy UID 71 behavior)
    minerApp({
      name: "spring01",
      hotkey: "spring01",
      port: 8101,
      strategy: "v5_style",
      extraEnv: { NIOME_VCF_DOT_ID: "1" },
    }),
    minerApp({
      name: "spring02",
      hotkey: "spring02",
      port: 8112,
      strategy: "v5_style",
      forceValidatorPermit: true,
    }),
    minerApp({
      name: "spring03",
      hotkey: "spring03",
      port: 8103,
      strategy: "v5_style",
      forceValidatorPermit: true,
    }),

    // v10 ×2 — balanced native (mid-band target)
    minerApp({
      name: "spring04",
      hotkey: "spring04",
      port: 8104,
      strategy: "v10",
      forceValidatorPermit: true,
    }),
    minerApp({
      name: "spring09",
      hotkey: "spring09",
      port: 8109,
      strategy: "v10",
    }),

    // win ×1 — truth panel when manager publishes real_correct_result / truth_by_task_id
    minerApp({
      name: "spring05",
      hotkey: "spring05",
      port: 8107,
      strategy: "win",
      extraEnv: {
        NIOME_WIN_MODE: "1",
        NIOME_VCF_MINIMAL: "1",
        NIOME_VCF_DOT_ID: "1",
      },
    }),

    // high_recall ×3 — 30–33 site rounds (beat native on ultra band)
    minerApp({
      name: "spring06",
      hotkey: "spring06",
      port: 8110,
      strategy: "high_recall",
      forceValidatorPermit: true,
    }),
    minerApp({
      name: "spring07",
      hotkey: "spring07",
      port: 8111,
      strategy: "high_recall",
    }),
    minerApp({
      name: "spring08",
      hotkey: "spring08",
      port: 8108,
      strategy: "high_recall",
      forceValidatorPermit: true,
    }),

    // auto ×1 — per-task: truth → win, else band → v5/recall/v10
    minerApp({
      name: "springhot",
      hotkey: "springhot",
      port: 8105,
      strategy: "auto",
      forceValidatorPermit: true,
    }),
  ],
};
