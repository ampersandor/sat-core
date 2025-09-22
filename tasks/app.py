import os
import re
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
from Bio import SeqIO

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
    statistic: str
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
        statistic=statistic.to_dict() if statistic else None,
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

    cmd = create_cmd(Tool(tool.lower()), input_path, align_file, shlex.split(options or ""))
    logger.info(f"Running command:\n{cmd}")

    output_dir.mkdir(parents=True, exist_ok=True)

    with open(align_file, "w") as f_out, open(log_file, "w", encoding="utf-8") as f_log:
        process = subprocess.run(cmd, stdout=f_out, stderr=f_log, text=True, shell=True)

        if process.returncode != 0:
            with open(log_file, "r", encoding="utf-8") as f:
                error_msg = f.read()
            logger.error("%s failed: %s", cmd, error_msg[:1000])
            raise RuntimeError(f"{cmd} failed with rc={process.returncode}: {error_msg}")

    logger.info("%s completed. Output: %s", cmd, align_file)

    if tool in [Tool.VSEARCH, Tool.UCLUST]:
        align_file = clean_align_file(align_file)

    stat_file_name, statistic = BlueBase(str(align_file), str(output_dir)).main()

    return align_file_name, stat_file_name, statistic


def create_cmd(tool: Tool, input_path, output_path, options: list):
    match tool:
        case Tool.MAFFT:
            return create_mafft_cmd(input_path, options)
        case Tool.VSEARCH:
            return create_vsearch_cmd(input_path, output_path, options)
        case Tool.UCLUST:
            return create_uclust_cmd(input_path, output_path, options)
        case _:
            raise ValueError(f"Invalid tool: {tool}") 


def create_mafft_cmd(input_path, options: list):
    cmd = f"{Tool.MAFFT.value} {" ".join(options)} {input_path}"
    return cmd


# vsearch --cluster_smallmen {input_file} --id {identity} --usersort --msaout {output_path}
def create_vsearch_cmd(input_path, output_path, options: list):
    cmd = f"{Tool.VSEARCH.value} --cluster_smallmem {input_path} {" ".join(options)} --msaout {output_path}"
    return cmd

def create_uclust_cmd(input_path, output_path, options: list):
    base_name = os.path.basename(output_path).split(".")[0]
    base_path = os.path.dirname(output_path)
    uc_file = os.path.join(base_path, f"{base_name}_pctid_0.uc")
    temp_fa = os.path.join(base_path, f"{base_name}_pctid_0.fa")

    cmd1 = f"{Tool.UCLUST.value} --input {input_path} --uc {uc_file} {" ".join(options)}"
    cmd2 = f"{Tool.UCLUST.value} --uc2fasta {uc_file} --input {input_path} --output {temp_fa}"
    cmd3 = f"{Tool.UCLUST.value} --staralign {temp_fa} --output {output_path}"
    return f"{cmd1} && {cmd2} && {cmd3}"

def clean_align_file(align_file):
    clean_align_file = align_file.replace(".aln", "_clean.aln")
    max_length = 0
    for record in SeqIO.parse(align_file, "fasta"):
        if record.id == "consensus":
            continue
        max_length = max(max_length, len(record.seq))

    with open(clean_align_file, "w") as f:
        for record in SeqIO.parse(align_file, "fasta"):
            if record.id == "consensus":
                continue
            if record.id.startswith("*"):
                record.id = record.id[1:]
            
            record.seq = re.sub(r'[.~]', '-', record.seq).ljust(max_length, "-").upper()
            SeqIO.write(record, f, "fasta")
    
    os.remove(align_file)

    return clean_align_file


@signals.task_success.connect
def on_success(sender=None, result=None, **kwargs):
    align_file, stat_file, statistic = result
    _notify(sender.request, "SUCCESS", align_file=align_file, stat_file=stat_file, statistic=statistic)


@signals.task_failure.connect
def on_failure(sender=None, exception=None, **kwargs):
    _notify(sender.request, "ERROR", message=str(exception))