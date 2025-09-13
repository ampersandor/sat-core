import os
import shlex
import subprocess
import logging
from enum import Enum
from pathlib import Path

import celery
import requests
from dotenv import load_dotenv
from celery import Celery, signals
from wonderwords import RandomWord
from dataclasses import dataclass
from .bluebase import BlueBase, Statistic

load_dotenv()
logger = logging.getLogger(__name__)
broker_url = os.getenv("BROKER_URL")
webhook_url = os.getenv("WEBHOOK_URL")
r = RandomWord()

app = Celery("tasks", broker=broker_url)
r = RandomWord()
PREFIX = Path("/data")


class Tool(Enum):
    MAFFT = "mafft"
    VSEARCH = "vsearch"
    UCLUST = "uclust"


@dataclass
class Response:
    status: str
    task_id: str
    align_file: str
    stat_file: str
    statistic: Statistic
    output_dir: str
    message: str



def _extract_args(request):
    """태스크 원본 인자 안전 추출"""
    k = getattr(request, "kwargs", {}) or {}
    a = tuple(getattr(request, "args", ()) or ())
    dir_name = k.get("dir_name") if k else (a[0] if len(a) > 0 else None)
    return dir_name


def _post_webhook(payload: dict):
    if not webhook_url:
        logger.warning("webhook_url is not set; skip send.")
        return
    try:
        resp = requests.post(webhook_url, json=payload, timeout=5)
        logger.info("Webhook %s %s", resp.status_code, resp.text[:500])
    except Exception as e:
        logger.exception("[Webhook error] %s", e)


def _notify(request, status: str, align_file=None, stat_file=None, statistic=None, message=None):
    task_id = request.id
    output_dir = _extract_args(request)
    response = Response(
        status=status,
        task_id=task_id,
        align_file=align_file,
        stat_file=stat_file,
        statistic=statistic,
        output_dir=output_dir,
        message=message,
    )
    _post_webhook(response.__dict__)


@app.task(bind=True, name="run_tool")
def run_tool(self: celery.Task, dir_name, base_name, tool, options):
    logger.info(f"Starting {tool}: {self.request.id}")
    output_dir = PREFIX / dir_name
    input_path = output_dir / base_name
    random_word = r.word()

    align_file_name = f"{random_word}.aln"
    align_file = output_dir / align_file_name
    log_file = output_dir / f"{random_word}.log"

    cmd = [Tool(tool.lower()).value, *shlex.split(options or ""), str(input_path)]
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(align_file, "w") as f_out, open(log_file, "w", encoding="utf-8") as f_log:
        process = subprocess.run(cmd, stdout=f_out, stderr=f_log, text=True)
        if process.returncode != 0:
            with open(log_file, "r", encoding="utf-8") as f:
                error_msg = f.read()
            logger.error("%s failed: %s", cmd[0], error_msg[:1000])
            raise RuntimeError(f"{cmd[0]} failed with rc={process.returncode}: {error_msg}")

    logger.info("%s completed. Output: %s", cmd[0], align_file)
    stat_file_name, statistic = BlueBase(str(input_path), str(output_dir)).main()

    return align_file_name, stat_file_name, statistic

@signals.task_success.connect
def on_success(sender=None, result=None, **kwargs):
    align_file, stat_file, statistic = result
    _notify(sender.request, "SUCCESS", align_file=align_file, stat_file=stat_file, statistic=statistic)


@signals.task_failure.connect
def on_failure(sender=None, exception=None, **kwargs):
    _notify(sender.request, "ERROR", message=str(exception))