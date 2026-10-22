FROM python:3.13-slim@sha256:57b29a7c050c4e679e4bea572e14b3965cf4b5839de2da49d55e49ffc4024f39

# terminus 系 Agent 在容器里要起 tmux 会话、用 asciinema 录屏，提前装好
#（否则运行时在线安装可能撞上 harbor 的 120s 工具安装预算直接拖垮 trial，2026-10-22 实测）
RUN apt-get update && apt-get install -y --no-install-recommends tmux && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir asciinema==2.4.0

WORKDIR /app

COPY src/ /app/src/
COPY data/ /app/data/

CMD ["bash"]
