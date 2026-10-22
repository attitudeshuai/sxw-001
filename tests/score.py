#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""把 pytest 的 junitxml 结果按 WEIGHTS 折算成总分（0~1），打印到 stdout。

用法：python3 /tests/score.py /tmp/results.xml
WEIGHTS 从 /tests/test_outputs.py 里静态读出（ast，不执行被测代码）。
"""
import ast
import os
import sys
import xml.etree.ElementTree as ET


def load_weights(path):
    with open(path, encoding="utf-8") as f:
        tree = ast.parse(f.read())
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and t.id == "WEIGHTS":
                    return ast.literal_eval(node.value)
    return {}


def main():
    xml_path = sys.argv[1]
    weights = load_weights(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "test_outputs.py"))
    if not os.path.exists(xml_path) or not weights:
        print("0")
        return
    try:
        root = ET.parse(xml_path).getroot()
    except ET.ParseError:
        print("0")
        return
    failed = set()
    for tc in root.iter("testcase"):
        name = tc.get("name", "")
        if list(tc):  # 有子节点（failure/error/skipped）算不通过
            failed.add(name)
    total = 0.0
    for name, w in weights.items():
        if name not in failed:
            total += w
        status = "PASS" if name not in failed else "FAIL"
        print("%s: %s (权重 %.3f)" % (name, status, w), file=sys.stderr)
    print("%.6f" % total)


if __name__ == "__main__":
    main()
