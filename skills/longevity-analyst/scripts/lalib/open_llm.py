"""Run one LongevityBench prompt through a named Longevity-LLM.

Weights are loaded once per process. The caller must pass the Hugging Face
repository id from bench_route control metadata. Training splits are never used here.
Called only via `la.py bench --control-model` (对照); never the default Skill answer path.

On Apple Silicon, prefer the MLX 4-bit conversion of L-Qwen3.5-9B when present
(`…/longevity-llm-mlx-4bit`); full bf16 weights OOM on typical MPS budgets.
"""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Dict, Optional, Tuple

_CACHE: Dict[str, object] = {}

_MLX_DEFAULT = (
    Path(os.environ["LONGEVITY_LLM_MLX"])
    if os.environ.get("LONGEVITY_LLM_MLX")
    else None
)


def _mlx_path_for(repo: str) -> Optional[Path]:
    p = Path(repo)
    if p.name.endswith("-mlx-4bit") and p.is_dir():
        return p
    cand = p.parent / f"{p.name}-mlx-4bit"
    if cand.is_dir():
        return cand
    if _MLX_DEFAULT is not None and _MLX_DEFAULT.is_dir() and "longevity-llm" in str(p):
        return _MLX_DEFAULT
    return None


def _load_mlx(path: Path):
    key = f"mlx:{path}"
    if key in _CACHE:
        return _CACHE[key]
    from mlx_lm import load

    model, tok = load(str(path))
    print(f"[open_llm] loaded MLX {path.name}", flush=True)
    _CACHE[key] = ("mlx", tok, model)
    return _CACHE[key]


def _load_hf(repo: str):
    key = f"hf:{repo}"
    if key in _CACHE:
        return _CACHE[key]
    import torch
    from transformers import AutoModelForCausalLM, AutoTokenizer, PreTrainedTokenizerFast

    try:
        tok = AutoTokenizer.from_pretrained(repo, trust_remote_code=True)
    except ValueError:
        tok = PreTrainedTokenizerFast.from_pretrained(repo)
    dtype = torch.float16
    if torch.cuda.is_available():
        device = "cuda"
        device_map = "auto"
    elif getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        device = "mps"
        device_map = {"": "mps"}
        try:
            torch.mps.empty_cache()
        except Exception:
            pass
    else:
        device = "cpu"
        device_map = {"": "cpu"}
        dtype = torch.float32
    kwargs = {
        "trust_remote_code": True,
        "dtype": dtype,
        "device_map": device_map,
        "low_cpu_mem_usage": True,
    }
    model = None
    err = None
    try:
        from transformers import Qwen3_5ForCausalLM

        model = Qwen3_5ForCausalLM.from_pretrained(repo, **kwargs)
    except Exception as e:
        err = e
        try:
            model = AutoModelForCausalLM.from_pretrained(repo, **kwargs)
            err = None
        except Exception as e2:
            err = e2
    if model is None:
        raise RuntimeError(f"failed to load Longevity-LLM on {device}: {err}")
    param_dev = str(next(model.parameters()).device)
    if device in ("mps", "cuda") and device not in param_dev:
        raise RuntimeError(f"wanted device={device} but weights landed on {param_dev}")
    print(f"[open_llm] loaded {type(model).__name__} on {param_dev}", flush=True)
    model.eval()
    _CACHE[key] = ("hf", tok, model)
    return _CACHE[key]


def _load(repo: str):
    mlx = _mlx_path_for(repo)
    if mlx is not None:
        try:
            return _load_mlx(mlx)
        except Exception as e:
            print(f"[open_llm] MLX load failed ({e}); trying HF", flush=True)
    return _load_hf(repo)


def _chat_prompt(tok, prompt: str, suffix: str = "") -> str:
    messages = [
        {
            "role": "system",
            "content": (
                "You are a helpful assistant with expertise in aging biology. "
                "For multiple-choice questions, respond with only the option letter."
            ),
        },
        {"role": "user", "content": prompt},
    ]
    if hasattr(tok, "apply_chat_template"):
        text = tok.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)
    else:
        text = prompt
    return text + suffix


def _generate(repo: str, prompt: str, max_new_tokens: int = 32, suffix: str = "") -> str:
    backend, tok, model = _load(repo)
    text = _chat_prompt(tok, prompt, suffix=suffix)
    if backend == "mlx":
        from mlx_lm import generate

        return generate(model, tok, prompt=text, max_tokens=max_new_tokens, verbose=False).strip()
    import torch

    inputs = tok(text, return_tensors="pt")
    device = next(model.parameters()).device
    inputs = {k: v.to(device) for k, v in inputs.items()}
    with torch.no_grad():
        out = model.generate(
            **inputs,
            max_new_tokens=max_new_tokens,
            do_sample=False,
            pad_token_id=getattr(tok, "eos_token_id", None),
        )
    new_tokens = out[0][inputs["input_ids"].shape[-1] :]
    return tok.decode(new_tokens, skip_special_tokens=True).strip()


def _first_letter(text: str) -> Optional[str]:
    # Prefer letter after think / answer markers (reasoning models).
    for pat in (
        r"(?:final answer|answer)\s*[:=]\s*([A-E])\b",
        r"</think>\s*([A-E])\b",
        r"^[\s\W]*([A-E])\b",
    ):
        m = re.search(pat, text, re.I | re.M)
        if m:
            return m.group(1).upper()
    m = re.search(r"\b([A-E])\b", text)
    return m.group(1) if m else None


def _first_number(text: str) -> Optional[float]:
    m = re.search(r"([+-]?\d+(?:\.\d+)?)\s*%?", text)
    return float(m.group(1)) if m else None


def answer(prompt: str, fmt: str, repo: str) -> Optional[dict]:
    """One exam field from an open Longevity-LLM. None when generation fails."""
    if fmt in ("pairwise", "binary", "multiclass"):
        raw = _generate(repo, prompt, max_new_tokens=32, suffix="Answer: ")
        letter = _first_letter(raw)
        if not letter:
            return None
        return {"answer": letter, "raw": raw}
    if fmt == "regression":
        raw = _generate(repo, prompt, max_new_tokens=64, suffix="The predicted age in years is ")
        value = _first_number(raw)
        if value is None:
            return None
        return {"answer": value, "raw": raw}
    return None
