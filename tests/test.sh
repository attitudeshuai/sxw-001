#!/bin/bash
# 评分脚本：reward 只写 /logs/verifier/reward.txt；开头先覆盖为 0；
# EXIT/INT/TERM/HUP 都兜底归零（只有正常跑到打分那一步才保留分数）；写完 chmod 0444。
set -u

mkdir -p /logs/verifier
echo 0 > /logs/verifier/reward.txt

finish() {
  if [ ! -f /tmp/.harbor_scored ]; then
    echo 0 > /logs/verifier/reward.txt
  fi
  chmod 0444 /logs/verifier/reward.txt 2>/dev/null || true
}
trap finish EXIT INT TERM HUP

# 判分工具由本脚本在判分阶段安装（离线 wheels 随 tests/ 挂载，不预装进镜像）
python3 -m pip install --no-index --find-links /tests/wheels pytest >/dev/null 2>&1

timeout 240 python3 -m pytest /tests/test_outputs.py -q --junitxml=/tmp/results.xml >/tmp/pytest.log 2>&1

SCORE=$(python3 /tests/score.py /tmp/results.xml) || SCORE=0
case "$SCORE" in
  ''|*[!0-9.]*) SCORE=0 ;;
esac
echo "$SCORE" > /logs/verifier/reward.txt && touch /tmp/.harbor_scored
