"""Thread-safe in-memory cache for rendered GUI scalar fields."""

from __future__ import annotations

import secrets
import sys
import threading
import time
from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np


@dataclass
class _CacheEntry:
    value: Any
    size_bytes: int
    touched_at: float


def _object_size(value: Any, seen: set[int] | None = None) -> int:
    """Estimate retained bytes without counting shared objects twice."""

    if seen is None:
        seen = set()
    object_id = id(value)
    if object_id in seen:
        return 0
    seen.add(object_id)

    if isinstance(value, np.ndarray):
        return int(value.nbytes)
    if isinstance(value, dict):
        return sys.getsizeof(value) + sum(
            _object_size(key, seen) + _object_size(item, seen)
            for key, item in value.items()
        )
    if isinstance(value, (list, tuple, set, frozenset)):
        return sys.getsizeof(value) + sum(_object_size(item, seen) for item in value)
    return sys.getsizeof(value)


class SurfaceCache:
    """Bound GUI surfaces by entry count, memory, and idle lifetime.

    A single surface larger than the byte budget is retained as the sole entry.
    This keeps the active analysis usable while still evicting older windows.
    """

    def __init__(
        self,
        *,
        max_entries: int = 32,
        max_bytes: int = 512 * 1024 * 1024,
        ttl_seconds: float = 2 * 60 * 60,
        clock: Callable[[], float] = time.monotonic,
    ):
        if max_entries < 1:
            raise ValueError("max_entries should be at least one.")
        if max_bytes < 1:
            raise ValueError("max_bytes should be positive.")
        if ttl_seconds <= 0:
            raise ValueError("ttl_seconds should be positive.")
        self.max_entries = int(max_entries)
        self.max_bytes = int(max_bytes)
        self.ttl_seconds = float(ttl_seconds)
        self._clock = clock
        self._entries: OrderedDict[str, _CacheEntry] = OrderedDict()
        self._total_bytes = 0
        self._lock = threading.RLock()

    @property
    def entry_count(self) -> int:
        with self._lock:
            return len(self._entries)

    @property
    def total_bytes(self) -> int:
        with self._lock:
            return self._total_bytes

    def put(self, value: Any) -> str:
        key = secrets.token_hex(16)
        now = self._clock()
        entry = _CacheEntry(
            value=value,
            size_bytes=_object_size(value),
            touched_at=now,
        )
        with self._lock:
            self._purge_expired(now)
            self._entries[key] = entry
            self._total_bytes += entry.size_bytes
            self._enforce_limits()
        return key

    def get(self, key: str | None) -> Any | None:
        if not key:
            return None
        cache_key = str(key)
        now = self._clock()
        with self._lock:
            self._purge_expired(now)
            entry = self._entries.get(cache_key)
            if entry is None:
                return None
            entry.touched_at = now
            self._entries.move_to_end(cache_key)
            return entry.value

    def clear(self) -> None:
        with self._lock:
            self._entries.clear()
            self._total_bytes = 0

    def _purge_expired(self, now: float) -> None:
        expired = [
            key
            for key, entry in self._entries.items()
            if now - entry.touched_at > self.ttl_seconds
        ]
        for key in expired:
            self._discard(key)

    def _enforce_limits(self) -> None:
        while len(self._entries) > self.max_entries:
            self._discard_oldest()
        while self._total_bytes > self.max_bytes and len(self._entries) > 1:
            self._discard_oldest()

    def _discard_oldest(self) -> None:
        key = next(iter(self._entries))
        self._discard(key)

    def _discard(self, key: str) -> None:
        entry = self._entries.pop(key, None)
        if entry is not None:
            self._total_bytes -= entry.size_bytes


__all__ = ["SurfaceCache"]
