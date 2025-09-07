import os
import subprocess
import logging
from enum import Enum

import requests
from dotenv import load_dotenv
from celery import Celery, signals
from wonderwords import RandomWord

load_dotenv()
logger = logging.getLogger(__name__)
broker_url = os.getenv("BROKER_URL", "pyamqp://guest@localhost//")
webhook_url = os.getenv("WEBHOOK_URL", "http://localhost:8080/api/v1/analyze/update")
r = RandomWord()

app = Celery("tasks", broker=broker_url)
PREFIX = "/data"

class Tool(Enum):
    MAFFT = "mafft"
    VSEARCH = "vsearch"
    UCLUST = "uclust"

@app.task(bind=True)
def run_tool(self, dir_name, base_name, tool, options):
    input_path = os.path.join(PREFIX, dir_name, base_name)
    random_word = r.word()

    output_file = os.path.join(PREFIX, dir_name, random_word + ".aln")
    log_file = os.path.join(PREFIX, dir_name, random_word + ".log")

    cmd = [tool, *options.split(), input_path]
    with open(output_file, "w") as f_out, open(log_file, "w", encoding="utf-8") as f_log:
        process = subprocess.run(cmd, stdout=f_out, stderr=f_log, text=True)
        if process.returncode != 0:
            with open(log_file, "r", encoding="utf-8") as f:
                error_msg = f.read()
            logger.error(f"{tool} alignment failed with error: {error_msg}")
            raise Exception(f"{tool} alignment failed with return code {process.returncode}: {error_msg}")
    logger.info(f"{tool} alignment completed successfully. Output saved to: {output_file}")

    return random_word + ".aln"


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