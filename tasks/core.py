import os
from dotenv import load_dotenv
import subprocess
from enum import Enum
import logging
import requests
import hashlib
from celery import Celery, signals

load_dotenv()
logger = logging.getLogger(__name__)
broker_url = os.getenv("BROKER_URL", "pyamqp://guest@localhost//")
webhook_url = os.getenv("WEBHOOK_URL", "http://localhost:8080/api/v1/analyze/update")

app = Celery("tasks", broker=broker_url)
data_dir = os.getenv("DATA_DIR", "/data")

class Tool(Enum):
    MAFFT = "mafft"
    VSEARCH = "vsearch"
    UCLUST = "uclust"

@app.task(bind=True)
def run_tool(self, input_path, tool, options):
    file_name = os.path.basename(input_path).split(".")[0]
    input_path = os.path.join(data_dir, input_path)
    hash_output_file = hashlib.sha256((file_name + tool + "".join(options)).encode()).hexdigest()
    res = hash_output_file + ".aln"
    output_file = os.path.join(data_dir, res)
    log_file = os.path.join(data_dir, hash_output_file + ".log")
    cmd = [tool, *options.split(), input_path]
    with open(output_file, "w") as f_out, open(log_file, "w", encoding="utf-8") as f_log:
        process = subprocess.run(cmd, stdout=f_out, stderr=f_log, text=True)
        if process.returncode != 0:
            with open(log_file, "r", encoding="utf-8") as f:
                error_msg = f.read()
            logger.error(f"{tool} alignment failed with error: {error_msg}")
            raise Exception(f"{tool} alignment failed with return code {process.returncode}: {error_msg}")
    logger.info(f"{tool} alignment completed successfully. Output saved to: {output_file}")

    return res


@signals.task_success.connect
def on_success(sender=None, result=None, **kwargs):
    task_id = sender.request.id
    output_file = result
    print(f"on_success: {task_id}, {output_file}")
    print(f"webhook_url: {webhook_url}")
    if webhook_url:
        try:
            payload = {
                "task_id": task_id,
                "status": "SUCCESS",
                "output_file": output_file,
                "message": None,
            }
            print(f"Sending webhook payload: {payload}")
            response = requests.post(webhook_url, json=payload)
            print(f"Webhook response status: {response.status_code}")
            print(f"Webhook response content: {response.text}")
        except Exception as e:
            print(f"[Webhook error] {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"webhook_url is not set")

@signals.task_failure.connect
def on_failure(sender=None, exception=None, **kwargs):
    task_id = sender.request.id
    print(f"on_failure: {task_id}, {exception}")
    print(f"webhook_url: {webhook_url}")
    if webhook_url:
        try:
            payload = {
                "task_id": task_id,
                "status": "ERROR",
                "output_file": None,
                "message": str(exception),
            }
            print(f"Sending webhook payload: {payload}")
            response = requests.post(webhook_url, json=payload)
            print(f"Webhook response status: {response.status_code}")
            print(f"Webhook response content: {response.text}")
        except Exception as e:
            print(f"[Webhook error] {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"webhook_url is not set")