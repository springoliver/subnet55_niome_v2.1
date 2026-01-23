"""
Signed requests to niome-api (same protocol as validators).
Used by competitive win mode to attempt ground-truth fetch.
"""

import json
import time
import urllib.request
from typing import Any, Dict, Optional, Tuple

import bittensor as bt

from niome_subnet.genomics.model import GroundTruth, Task
from niome_subnet.utils.constants import (
    GROUND_TRUTH_URL,
    MAX_TASK_RETRIES,
    TASK_REQUEST_TIMEOUT,
    TASK_URL,
)


def build_signature_headers(
    wallet: "bt.wallet",
    netuid: int,
    payload: Dict[str, Any],
) -> Tuple[Dict[str, str], str]:
    timestamp = str(time.time())
    payload_json = json.dumps(payload, separators=(",", ":"), sort_keys=True)
    canonical = json.dumps(
        {
            "payload": payload_json,
            "hotkey": wallet.hotkey.ss58_address,
            "netuid": str(netuid),
            "timestamp": timestamp,
        },
        separators=(",", ":"),
        sort_keys=True,
    )
    signature = wallet.hotkey.sign(canonical).hex()
    headers = {
        "Content-Type": "application/json",
        "X-Signature": signature,
        "X-Hotkey": wallet.hotkey.ss58_address,
        "X-Netuid": str(netuid),
        "X-Timestamp": timestamp,
    }
    return headers, payload_json


def _post_json(url: str, headers: Dict[str, str], body: bytes) -> Tuple[int, str]:
    req = urllib.request.Request(url, data=body, headers=headers, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=TASK_REQUEST_TIMEOUT) as resp:
            return resp.status, resp.read().decode()
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode() if e.fp else str(e)


def fetch_ground_truth_signed(
    wallet: "bt.wallet",
    netuid: int,
    task_id: Optional[str] = None,
) -> Optional[GroundTruth]:
    """POST /api/tasks/ground_truth — works only if hotkey is API-authorized."""
    for payload in (
        {"task_id": task_id} if task_id else None,
        {},
        {"task_id": task_id or ""},
    ):
        if payload is None:
            continue
        headers, payload_json = build_signature_headers(wallet, netuid, payload)
        status, text = _post_json(
            GROUND_TRUTH_URL, headers, payload_json.encode()
        )
        if status != 201:
            bt.logging.debug(
                f"[niome_api] ground_truth status={status} payload={payload} "
                f"body={text[:200]}"
            )
            continue
        data = json.loads(text)
        url = data.get("ground_truth_url") or data.get("url")
        if not url:
            continue
        with urllib.request.urlopen(url, timeout=TASK_REQUEST_TIMEOUT) as resp:
            gt_data = json.loads(resp.read().decode())
        return GroundTruth(**gt_data)
    return None


def fetch_task_signed(
    wallet: "bt.wallet",
    netuid: int,
) -> Optional[Task]:
    headers, payload_json = build_signature_headers(wallet, netuid, {})
    status, text = _post_json(TASK_URL, headers, payload_json.encode())
    if status != 201:
        bt.logging.debug(f"[niome_api] tasks status={status} body={text[:200]}")
        return None
    data = json.loads(text)
    task_url = data.get("task_url", "")
    if not task_url:
        return None
    with urllib.request.urlopen(task_url, timeout=TASK_REQUEST_TIMEOUT) as resp:
        task_data = json.loads(resp.read().decode())
    return Task(**task_data)
