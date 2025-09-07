FROM python:3.13-slim

RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*

RUN pip install poetry
RUN poetry install

WORKDIR /app
COPY ./main.py .


CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]