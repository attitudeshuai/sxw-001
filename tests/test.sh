#!/bin/bash
# 评分脚本：reward 只写 /logs/verifier/reward.txt；开头先覆盖为 0，崩溃/超时/信号都归零。
set -u

mkdir -p /logs/verifier
echo 0 > /logs/verifier/reward.txt
trap 'echo 0 > /logs/verifier/reward.txt' INT TERM HUP

SCORE=$(python3 /tests/test_outputs.py) || SCORE=0
case "$SCORE" in
  ''|*[!0-9.]*) SCORE=0 ;;
esac
echo "$SCORE" > /logs/verifier/reward.txt
