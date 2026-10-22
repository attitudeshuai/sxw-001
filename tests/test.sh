#!/bin/bash
# 评分脚本：reward 只写 /logs/verifier/reward.txt；开头先覆盖为 0；
# EXIT/INT/TERM/HUP 都兜底归零（只有正常跑到打分那一步才保留分数）；写完 chmod 0444。
# 判分依赖结构化 pytest 完成报告（junitxml + score.py 的 scoring-report.json）；
# pytest 没正常完成或完成报告不存在/不可解析时记 INCOMPLETE 并写 0。
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

RESULTS=/tmp/results.xml
timeout 240 python3 -m pytest /tests/test_outputs.py -q --junitxml="$RESULTS" >/tmp/pytest.log 2>&1
PYTEST_RC=$?

# 完成报告存在且可解析才采信（非空检查 + ET.parse 试解析）
REPORT_OK=1
if [ ! -s "$RESULTS" ]; then
  REPORT_OK=0
else
  python3 -c "import xml.etree.ElementTree as ET; ET.parse('$RESULTS')" 2>/dev/null || REPORT_OK=0
fi

if [ $PYTEST_RC -ne 0 ] || [ $REPORT_OK -eq 0 ]; then
  echo "INCOMPLETE: pytest 未正常完成或完成报告不可用（rc=$PYTEST_RC, report_ok=$REPORT_OK）" >&2
  echo 0 > /logs/verifier/reward.txt
  touch /tmp/.harbor_scored
  exit 0
fi

python3 /tests/score.py --junit "$RESULTS" \
  --report-out /tmp/scoring-report.json \
  --reward-out /logs/verifier/reward.txt || echo 0 > /logs/verifier/reward.txt
touch /tmp/.harbor_scored
