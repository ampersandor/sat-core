FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*
RUN curl -LsSf https://astral.sh/uv/install.sh | sh

RUN apt-get update && apt-get install -y tzdata \
    && ln -fs /usr/share/zoneinfo/Asia/Seoul /etc/localtime \
    && echo "Asia/Seoul" > /etc/timezone \
    && dpkg-reconfigure -f noninteractive tzdata \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/root/.local/bin:$PATH"
WORKDIR /app
COPY pyproject.toml uv.lock ./

RUN uv sync


FROM python:3.12-slim AS runtime

# 비루트 사용자 생성
RUN addgroup --system --gid 1001 ampersandor && \
    adduser --system --uid 1001 --gid 1001 --shell /bin/bash --home /home/ampersandor ampersandor && \
    mkdir -p /home/ampersandor/.local && \
    chown -R ampersandor:ampersandor /home/ampersandor

ENV PYTHONPATH=/app
ENV PYTHONUNBUFFERED=1

COPY --from=builder /app/.venv/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /root/.local/bin /usr/local/bin

COPY main.py ./app/

# 비루트 사용자로 전환
USER ampersandor
WORKDIR /app

CMD ["uv", "run", "python", "main.py"]
