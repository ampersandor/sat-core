import os
import re
import shlex
import subprocess
import logging
import time
import threading
import psutil
from enum import Enum
from pathlib import Path
from typing import Optional

import celery
import requests
from dotenv import load_dotenv
from celery import Celery, signals
from wonderwords import RandomWord
from dataclasses import dataclass
from Bio import SeqIO
from Bio.Seq import Seq

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
class ResourceUsage:
    """리소스 사용량 정보"""
    max_memory_mb: float
    avg_memory_mb: float
    max_cpu_percent: float
    avg_cpu_percent: float
    execution_time_seconds: float


@dataclass
class Response:
    status: str
    task_id: str
    align_file: str
    stat_file: str
    statistic: str
    output_dir: str
    message: str
    resource_usage: Optional[dict] = None


class ResourceMonitor:
    """프로세스 리소스 모니터링 클래스"""
    
    def __init__(self, process: subprocess.Popen, interval: float = 0.1):
        self.process = process
        self.interval = interval
        self.memory_samples = []
        self.cpu_samples = []
        self.start_time = time.time()
        self.monitoring = False
        self.monitor_thread = None
        
    def start_monitoring(self):
        """모니터링 시작"""
        self.monitoring = True
        self.start_time = time.time()
        self.monitor_thread = threading.Thread(target=self._monitor_loop)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()
        
    def stop_monitoring(self) -> ResourceUsage:
        """모니터링 중지 및 결과 반환"""
        self.monitoring = False
        if self.monitor_thread:
            self.monitor_thread.join(timeout=1.0)
            
        execution_time = time.time() - self.start_time
        
        # 샘플이 없는 경우 기본값 반환
        if not self.memory_samples or not self.cpu_samples:
            return ResourceUsage(
                max_memory_mb=0.0,
                avg_memory_mb=0.0,
                max_cpu_percent=0.0,
                avg_cpu_percent=0.0,
                execution_time_seconds=execution_time
            )
            
        return ResourceUsage(
            max_memory_mb=max(self.memory_samples),
            avg_memory_mb=sum(self.memory_samples) / len(self.memory_samples),
            max_cpu_percent=max(self.cpu_samples),
            avg_cpu_percent=sum(self.cpu_samples) / len(self.cpu_samples),
            execution_time_seconds=execution_time
        )
        
    def _monitor_loop(self):
        """모니터링 루프"""
        try:
            # psutil Process 객체 생성
            ps_process = psutil.Process(self.process.pid)
            
            while self.monitoring and self.process.poll() is None:
                try:
                    # 메모리 사용량 (MB) - 부모 프로세스 + 모든 자식 프로세스 포함
                    memory_mb = 0.0
                    cpu_percent = 0.0
                    
                    try:
                        # 부모 프로세스 메모리
                        memory_mb = ps_process.memory_info().rss / (1024 * 1024)
                        # 부모 프로세스 CPU (멀티스레드 포함)
                        cpu_percent = ps_process.cpu_percent()
                        
                        # 자식 프로세스들의 리소스 추가 (fork된 경우 대비)
                        for child in ps_process.children(recursive=True):
                            try:
                                print(f"child: {child}, memory_mb: {child.memory_info().rss / (1024 * 1024)}, cpu_percent: {child.cpu_percent()}")
                                memory_mb += child.memory_info().rss / (1024 * 1024)
                                cpu_percent += child.cpu_percent()
                            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                                # 자식 프로세스가 이미 종료된 경우 무시
                                pass
                                
                    except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                        # 부모 프로세스가 종료된 경우
                        break
                    
                    self.memory_samples.append(memory_mb)
                    self.cpu_samples.append(cpu_percent)
                    
                    time.sleep(self.interval)
                    
                except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                    # 프로세스가 종료되었거나 접근할 수 없는 경우
                    break
                except Exception as e:
                    logger.warning(f"Resource monitoring error: {e}")
                    break
                    
        except Exception as e:
            logger.error(f"Failed to start resource monitoring: {e}")


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


def _notify(request, status: str, align_file=None, stat_file=None, statistic=None, message=None, resource_usage=None):
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
        resource_usage=resource_usage.__dict__ if resource_usage else None,
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
    with open(log_file, "w", encoding="utf-8") as f_log:
        f_log.write(f"Running command:\n{cmd}\n")

    resource_usage = None
    
    with open(align_file, "w", encoding="utf-8") as f_out, open(log_file, "a", encoding="utf-8") as f_log:
        # 프로세스 시작
        process = subprocess.Popen(
            cmd, 
            stdout=f_out, 
            stderr=f_log, 
            text=True, 
            shell=True
        )
        
        # 리소스 모니터링 시작
        monitor = ResourceMonitor(process)
        monitor.start_monitoring()
        
        try:
            # 프로세스 완료 대기
            return_code = process.wait()
            
            # 리소스 모니터링 중지 및 결과 수집
            resource_usage = monitor.stop_monitoring()
            
            logger.info(f"Resource usage - Memory: {resource_usage.max_memory_mb:.2f}MB (max), "
                       f"CPU: {resource_usage.max_cpu_percent:.2f}% (max), "
                       f"Time: {resource_usage.execution_time_seconds:.2f}s")
            
            if return_code != 0:
                with open(log_file, "r", encoding="utf-8") as f:
                    error_msg = f.read()
                logger.error("%s failed: %s", cmd, error_msg[:1000])
                raise RuntimeError(f"{cmd} failed with rc={return_code}: {error_msg}")
                
        except Exception as e:
            # 예외 발생 시에도 리소스 사용량 수집
            resource_usage = monitor.stop_monitoring()
            raise

    logger.info("%s completed. Output: %s", cmd, align_file)

    if Tool(tool.lower()) in [Tool.VSEARCH, Tool.UCLUST]:
        logger.info("Cleaning output from %s", Tool(tool.lower()))
        max_length = clean_align_file(align_file)
        if max_length == 0:
            logger.error("max_length: %d", max_length)
            raise ValueError("Alignment file is empty")

    stat_file_name, statistic = BlueBase(str(align_file), str(output_dir)).main()

    return align_file_name, stat_file_name, statistic, resource_usage


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
    cmd = f"{Tool.MAFFT.value} {" ".join(options)} --thread 8 {input_path}"
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
    clean_align_file = align_file.with_name(align_file.stem + "_clean" + align_file.suffix)
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
                record.description = record.id
            
            record.seq = Seq(re.sub(r'[.~]', '-', str(record.seq)).ljust(max_length, "-").upper())
            SeqIO.write(record, f, "fasta")
    
    os.remove(align_file)
    os.rename(clean_align_file, align_file)
    return max_length


@signals.task_success.connect
def on_success(sender=None, result=None, **kwargs):
    align_file, stat_file, statistic, resource_usage = result
    _notify(sender.request, "SUCCESS", align_file=align_file, stat_file=stat_file, 
            statistic=statistic, resource_usage=resource_usage)


@signals.task_failure.connect
def on_failure(sender=None, exception=None, **kwargs):
    _notify(sender.request, "ERROR", message=str(exception))


