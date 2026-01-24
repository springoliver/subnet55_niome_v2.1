# The MIT License (MIT)
# Copyright © 2023 Yuma Rao
# Copyright © 2025 Genomes.io

import asyncio
import os
import sys
import tempfile
import time
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import bittensor as bt

from niome_subnet.base.miner import BaseMinerNeuron
from niome_subnet.genomics.cftr_lookup import build_cftr_annotations
from niome_subnet.genomics.competitive import COMPETITIVE_REV, solve_competitive, win_mode_enabled
from niome_subnet.genomics.pipeline import run_pipeline
from niome_subnet.genomics.read_calling import READ_CALLING_REV
from niome_subnet.protocol import GenomicsTaskSynapse

bt.logging.on()

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)


@dataclass
class _TaskResult:
    vcf_content: str
    cftr_annotations: Optional[dict]


class Miner(BaseMinerNeuron):
    """NIOME v3 miner: read-all VCF + ClinVar-only annotations."""

    TASK_CACHE_TTL = 3600

    def __init__(self, config=None):
        super(Miner, self).__init__(config=config)
        self._task_locks: Dict[str, asyncio.Lock] = {}
        self._task_cache: Dict[str, Tuple[float, _TaskResult]] = {}

    def _lock_for(self, task_id: str) -> asyncio.Lock:
        if task_id not in self._task_locks:
            self._task_locks[task_id] = asyncio.Lock()
        return self._task_locks[task_id]

    @staticmethod
    def _cache_key(task) -> str:
        rev = COMPETITIVE_REV if win_mode_enabled() else READ_CALLING_REV
        return f"{task.task_id}:{task.genome_context.region}:{rev}"

    def _cache_get(self, task) -> Optional[_TaskResult]:
        key = self._cache_key(task)
        entry = self._task_cache.get(key)
        if not entry:
            return None
        ts, result = entry
        if time.time() - ts > self.TASK_CACHE_TTL:
            del self._task_cache[key]
            return None
        return result

    def _cache_put(self, task, result: _TaskResult) -> None:
        key = self._cache_key(task)
        self._task_cache[key] = (time.time(), result)
        if len(self._task_cache) > 32:
            oldest = min(self._task_cache, key=lambda k: self._task_cache[k][0])
            del self._task_cache[oldest]

    async def _solve_task(self, task) -> _TaskResult:
        cached = self._cache_get(task)
        if cached is not None:
            return cached

        lock = self._lock_for(task.task_id)
        async with lock:
            cached = self._cache_get(task)
            if cached is not None:
                return cached

            with tempfile.TemporaryDirectory(prefix="niome_miner_") as work_dir:
                vcf_content = None
                cftr_annotations = None

                if win_mode_enabled():
                    competitive = await asyncio.to_thread(
                        solve_competitive,
                        task,
                        work_dir,
                        self.wallet,
                        self.config.netuid,
                    )
                    if competitive:
                        vcf_content, cftr_annotations = competitive
                        bt.logging.info(
                            f"Task {task.task_id} rev={COMPETITIVE_REV} (win mode)"
                        )

                if vcf_content is None:
                    final_vcf, _ = await asyncio.to_thread(
                        run_pipeline, task, work_dir
                    )
                    with open(final_vcf) as fh:
                        vcf_content = fh.read()
                    n_try = sum(
                        1
                        for line in vcf_content.splitlines()
                        if line and not line.startswith("#")
                    )
                    if n_try == 0 and os.environ.get(
                        "NIOME_DISABLE_PIPELINE_RETRY", ""
                    ).strip().lower() not in ("1", "true", "yes"):
                        bt.logging.warning(
                            f"Task {task.task_id}: zero variants — pipeline retry"
                        )
                        final_vcf, _ = await asyncio.to_thread(
                            run_pipeline, task, work_dir
                        )
                        with open(final_vcf) as fh:
                            vcf_content = fh.read()
                    cftr_annotations = await asyncio.to_thread(
                        build_cftr_annotations, final_vcf
                    )

            n_variants = sum(
                1
                for line in vcf_content.splitlines()
                if line and not line.startswith("#")
            )
            result = _TaskResult(
                vcf_content=vcf_content,
                cftr_annotations=cftr_annotations,
            )
            if n_variants > 0:
                self._cache_put(task, result)
            else:
                bt.logging.error(
                    f"Task {task.task_id}: not caching — zero variants submitted"
                )
            return result

    async def forward(self, synapse: GenomicsTaskSynapse) -> GenomicsTaskSynapse:
        try:
            start_time = time.time()
            task = synapse.task
            rev = COMPETITIVE_REV if win_mode_enabled() else READ_CALLING_REV
            bt.logging.info(
                f"Task {task.task_id} region={task.genome_context.region} "
                f"rev={rev} expected={task.expected_variant_count} "
                f"win_mode={win_mode_enabled()}"
            )

            result = await self._solve_task(task)
            synapse.vcf_content = result.vcf_content
            synapse.cftr_annotations = result.cftr_annotations

            n = sum(
                1
                for line in result.vcf_content.splitlines()
                if line and not line.startswith("#")
            )
            n_ann = len(result.cftr_annotations or {})
            bt.logging.info(
                f"Task {task.task_id} done in {time.time() - start_time:.1f}s "
                f"variants={n} clinvar_annotations={n_ann}"
            )
        except Exception as e:
            bt.logging.error(f"Forward error: {e}", exc_info=True)
            synapse.vcf_content = None
            synapse.cftr_annotations = None

        return synapse

    async def blacklist(self, synapse: GenomicsTaskSynapse) -> Tuple[bool, str]:
        if synapse.dendrite is None or synapse.dendrite.hotkey is None:
            return True, "Missing dendrite or hotkey"

        uid = self.metagraph.hotkeys.index(synapse.dendrite.hotkey)
        if (
            not self.config.blacklist.allow_non_registered
            and synapse.dendrite.hotkey not in self.metagraph.hotkeys
        ):
            return True, "Unrecognized hotkey"

        if self.config.blacklist.force_validator_permit:
            if not self.metagraph.validator_permit[uid]:
                return True, "Non-validator hotkey"

        return False, "Hotkey recognized!"

    async def priority(self, synapse: GenomicsTaskSynapse) -> float:
        if synapse.dendrite is None or synapse.dendrite.hotkey is None:
            return 0.0
        caller_uid = self.metagraph.hotkeys.index(synapse.dendrite.hotkey)
        return float(self.metagraph.S[caller_uid])


if __name__ == "__main__":
    with Miner() as miner:
        while True:
            bt.logging.info(f"Miner running… {time.time()}")
            time.sleep(5)
