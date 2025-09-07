FROM python:3.13-slim

RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*


WORKDIR /app
COPY ./main.py .

RUN pip install celery fastapi uvicorn requests

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]