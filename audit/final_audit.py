# -*- coding: utf-8 -*-
"""Final pre-submission audit: automated consistency battery on the EN main manuscript."""
import re
from docx import Document

PATH = "方案B_合并稿/投稿_HeartRhythm_压缩版/咖啡与房颤_投稿版_v2.docx"
d = Document(PATH)
paras = d.paragraphs
full = "\n".join(p.text for p in paras)
ri = next(i for i, p in enumerate(paras) if p.text.strip() == "References")
body = "\n".join(p.text for p in paras[:ri])

issues = []
ok = []

# 1. citation order
cits = []
for p in paras[:ri]:
    for m in re.finditer(r"\[(\d+(?:[,\-–]\s*\d+)*)\]", p.text):
        for n in re.findall(r"\d+", m.group(1)):
            cits.append(int(n))
first = {}
for n in cits:
    if n not in first:
        first[n] = len(first) + 1
if all(first.get(n) == n for n in range(1, 43)):
    ok.append("引文首现顺序 1–42 连续")
else:
    issues.append(f"引文顺序异常: {[(n, first.get(n)) for n in range(1,43) if first.get(n)!=n]}")

# 2. supplementary table citations: which S-numbers cited in body
s_cited = sorted(set(int(x) for x in re.findall(r"[Tt]able S(\d+)", body)))
missing_s = [n for n in range(1, 30) if n not in s_cited]
ok.append(f"正文引用的补充表: S{',S'.join(map(str,s_cited))}")
if missing_s:
    issues.append(f"补充表存在但正文未引用: S{missing_s}")

# 3. main tables/figures cited
for t in ("Table 1", "Table 2", "Table 3", "Table 4"):
    (ok if re.search(rf"\b{t}\b", body) else issues).append(f"{t} 正文引用" if re.search(rf"\b{t}\b", body) else f"{t} 正文未引用!")
for f in ("Figure 1", "Figure 2", "Figure 3"):
    if not re.search(rf"\b{f}\b", body):
        issues.append(f"{f} 正文未引用!")

# 4. key numbers consistency (count occurrences)
for pat, desc in [
    (r"1\.32", "OR 1.32 pooled"),
    (r"1\.19", "PRESSO 1.19"),
    (r"1\.06", "MVMR 1.06"),
    (r"5\.33", "A-cluster OR 5.33"),
    (r"0\.58", "recurrence 0.58"),
    (r"0\.967", "biomarker 0.967"),
    (r"116,473", "pooled cases"),
    (r"77,690", "cross-ancestry cases"),
    (r"428,860", "UKB n"),
    (r"8,299", "metabolome N"),
    (r"2\.3 × 10−7", "cluster Q P"),
    (r"4\.04", "PR +4.04ms"),
    (r"0\.81", "FAERS ROR 0.81"),
]:
    n = len(re.findall(pat, full))
    ok.append(f"'{desc}' 出现 {n} 次")

# 5. stale framing / old title remnants
for pat, desc in [
    (r"Triangulating Total Genetic Effects", "旧标题残留"),
    (r"three orthogonally biased designs", "三设计旧框架(优点句)"),
    (r"triangulating total genetic effects, mediation pathways, and real-world safety", "摘要三腿旧框架"),
    (r"Heart Rhythm", "他刊名残留"),
    (r"JCE", "JCE残留"),
]:
    hits = [p.text[:80] for p in paras if re.search(pat, p.text, re.I)]
    if hits:
        issues.append(f"{desc}: {len(hits)} 处 -> {hits[:2]}")
    else:
        ok.append(f"无{desc}")

# 6. data availability completeness
da = next((p.text for p in paras if p.text.strip().startswith("Data availability:")), "")
for kw in ["FinnGen R11", "FinnGen R12", "Nielsen", "BioBank Japan", "Miyazawa",
           "GWAS Catalog", "metabolome", "Chen", "Ensembl", "FAERS", "openFDA",
           "QT", "heart rate", "PR interval", "PAGE"]:
    if kw.lower() not in da.lower():
        issues.append(f"数据可得性声明缺少: {kw}")
ok.append("数据可得性声明已检查")

# 7. word counts
ai = next(i for i, p in enumerate(paras) if p.text.strip() == "Abstract")
bw = sum(len(p.text.split()) for p in paras[ai:ri])
aw = sum(len(p.text.split()) for p in paras)
n_abs = sum(len(paras[i].text.split()) for i in (ai + 1, ai + 2, ai + 3))
ok.append(f"词数: 正文 {bw}, 全含 {aw}, 摘要 {n_abs}")
if n_abs > 250:
    issues.append(f"摘要 {n_abs} 词 > 250!")

# 8. section numbering continuity
secs = [p.text.strip() for p in paras if re.match(r"^[23]\.\d+\s", p.text.strip())]
m_nums = [int(re.match(r"^2\.(\d+)", s).group(1)) for s in secs if s.startswith("2.")]
r_nums = [int(re.match(r"^3\.(\d+)", s).group(1)) for s in secs if s.startswith("3.")]
if m_nums != list(range(1, 12)):
    issues.append(f"方法节编号异常: {m_nums}")
if r_nums != list(range(1, 13)):
    issues.append(f"结果节编号异常: {r_nums}")
ok.append(f"方法 2.1–2.11、结果 3.1–3.12 编号连续")

# 9. reference list
refs = [p.text.strip() for p in paras[ri + 1:] if re.match(r"^\d+\.\s", p.text.strip())]
nums = [int(re.match(r"^(\d+)\.", r).group(1)) for r in refs]
if nums != list(range(1, 43)):
    issues.append(f"文献编号异常: {nums}")
else:
    ok.append("文献 1–42 连续")

# 10. tables count in doc
ok.append(f"主文表格对象数: {len(d.tables)} (应为 4)")
if len(d.tables) != 4:
    issues.append(f"主文表格数 {len(d.tables)} != 4")

print("=" * 30, "PASS", "=" * 30)
for x in ok:
    print(" ✓", x)
print("=" * 30, "ISSUES", "=" * 30)
for x in issues:
    print(" ✗", x)
print(f"\n共 {len(issues)} 个问题")
