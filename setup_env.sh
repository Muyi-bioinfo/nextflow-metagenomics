#!/usr/bin/env bash
# ============================================================================
# setup_env.sh —— nf-meta 环境的创建与安装后修复
#
# 用法:
#     bash setup_env.sh              # 创建 nf-meta 并修复
#     bash setup_env.sh --repair     # 仅对已存在的环境执行修复
#     bash setup_env.sh --name foo   # 用别的环境名
#
# 为什么需要这个脚本 —— conda 解出的环境并不等于可用的环境。
# 有两类问题 environment.yml 表达不了, 必须在安装后修:
#
#   1. humann 的 conda 包内含 bin/bowtie2* 与 bin/diamond 的**自带副本**,
#      安装时会覆盖同名的独立包。因 humann 在依赖图中较晚链接, 覆盖结果是:
#          bowtie2  2.5.5  ->  2.2.3   (2014 年版)
#          diamond  2.1.11 ->  2.0.15
#      全程仅有一条 libmamba warning, 不报错。bowtie2 2.2.3 会被 HOST_REMOVAL
#      直接用上; diamond 2.0.15 则违反 checkm2 硬钉的 2.1.11。
#      修法: 对这两个包 --force-reinstall, 把正确的二进制写回去。
#
#   2. eggnog-mapper 2.1.12 build _0 到自己的包内目录找可执行文件
#      (site-packages/eggnogmapper/bin/), 而该目录是空的 —— 于是 emapper
#      报 "Diamond was not found"。(修好此问题的 _1/_2 build 要求
#      python <3.12, 与本环境的 python 3.12 冲突, 用不了。)
#      修法: 把 PATH 上的对应二进制软链进去。
#
# 两处修复都是幂等的, 可重复执行。
# ============================================================================

set -euo pipefail

ENV_NAME="nf-meta"
REPAIR_ONLY=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repair) REPAIR_ONLY=1; shift ;;
        --name)   ENV_NAME="$2";  shift 2 ;;
        -h|--help) sed -n '2,20p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "未知参数: $1" >&2; exit 1 ;;
    esac
done

# 优先用 mamba, 回退到 conda
if command -v mamba >/dev/null 2>&1; then
    CONDA_BIN=mamba
elif command -v conda >/dev/null 2>&1; then
    CONDA_BIN=conda
else
    echo "ERROR: 未找到 mamba 或 conda。" >&2
    exit 1
fi

# ----------------------------------------------------------------------------
# 1. 创建环境
# ----------------------------------------------------------------------------
if [[ "$REPAIR_ONLY" -eq 0 ]]; then
    if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        echo "ERROR: 环境 '$ENV_NAME' 已存在。"
        echo "       重建请先: conda env remove -n $ENV_NAME"
        echo "       仅修复请: bash setup_env.sh --repair"
        exit 1
    fi
    echo ">>> 创建环境 '$ENV_NAME' (约 8 GB, 含 R 全家桶与 tensorflow, 耗时较长)"
    "$CONDA_BIN" env create -f "${SCRIPT_DIR}/environment.yml" -n "$ENV_NAME" -y
fi

ENV_PREFIX="$(conda env list | awk -v n="$ENV_NAME" '$1==n {print $NF}')"
if [[ -z "$ENV_PREFIX" || ! -d "$ENV_PREFIX" ]]; then
    echo "ERROR: 找不到环境 '$ENV_NAME' 的路径。" >&2
    exit 1
fi
echo ">>> 环境路径: $ENV_PREFIX"

# ----------------------------------------------------------------------------
# 2. 修复 humann 对 bowtie2 / diamond 的二进制覆盖
# ----------------------------------------------------------------------------
echo ">>> 修复 humann 覆盖的 bowtie2 / diamond"
"$CONDA_BIN" install -n "$ENV_NAME" -y --force-reinstall \
    --override-channels -c conda-forge -c bioconda \
    bowtie2=2.5.5 diamond=2.1.11

# ----------------------------------------------------------------------------
# 3. 把可执行文件软链进 eggnog-mapper 的包内 bin 目录
# ----------------------------------------------------------------------------
EGGNOG_BIN="$(find "${ENV_PREFIX}/lib" -maxdepth 4 -type d -path '*/site-packages/eggnogmapper/bin' 2>/dev/null | head -n 1)"
if [[ -n "$EGGNOG_BIN" ]]; then
    echo ">>> 修复 eggnog-mapper 包内 bin 目录: $EGGNOG_BIN"
    for b in diamond mmseqs prodigal hmmscan hmmsearch hmmpress phmmer; do
        if [[ -x "${ENV_PREFIX}/bin/${b}" ]]; then
            ln -sf "${ENV_PREFIX}/bin/${b}" "${EGGNOG_BIN}/${b}"
        fi
    done
else
    echo "!!! 警告: 未找到 eggnogmapper/bin 目录, 跳过该修复。"
fi

# ----------------------------------------------------------------------------
# 4. 验证 —— 每个工具都必须真的能启动
# ----------------------------------------------------------------------------
echo
echo ">>> 验证工具可启动 (求解成功 != 能跑)"
B="${ENV_PREFIX}/bin"
FAILED=0

check() {
    local name="$1"; shift
    if "$@" >/dev/null 2>&1; then
        printf "  %-16s OK\n" "$name"
    else
        printf "  %-16s !! FAILED: %s\n" "$name" "$*"
        FAILED=$((FAILED + 1))
    fi
}

check python      "$B/python"     --version
check fastqc      "$B/fastqc"     --version
check fastp       "$B/fastp"      --version
check bowtie2     "$B/bowtie2"    --version
check samtools    "$B/samtools"   --version
check multiqc     "$B/multiqc"    --version
check kraken2     "$B/kraken2"    --version
check bracken     "$B/bracken"    -v
check metaphlan   "$B/metaphlan"  --version
check humann      "$B/humann"     --version
check megahit     "$B/megahit"    --version
check spades      "$B/spades.py"  --version
check quast       "$B/quast.py"   --version
check metabat2    "$B/metabat2"   --help
check checkm2     "$B/checkm2"    --version
check gtdbtk      "$B/gtdbtk"     --version
check drep        "$B/dRep"       --help
check prodigal    "$B/prodigal"   -v
check diamond     "$B/diamond"    --version
check emapper     "$B/emapper.py" --version
check coverm      "$B/coverm"     --version
check nextflow    "$B/nextflow"   -version
# rgi 无 --version, 需给子命令
check rgi         "$B/rgi"        main --help

# 版本回归检查 —— 这两个曾被 humann 静默覆盖
echo
echo ">>> 覆盖回归检查"
BT_VER="$("$B/bowtie2" --version 2>/dev/null | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
DM_VER="$("$B/diamond" --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
printf "  bowtie2 = %s  (须为 2.5.5, 若为 2.2.3 则 humann 覆盖未修复)\n" "${BT_VER:-未知}"
printf "  diamond = %s  (须为 2.1.11, 若为 2.0.15 则 humann 覆盖未修复)\n" "${DM_VER:-未知}"
[[ "$BT_VER" == "2.5.5"  ]] || { echo "  !! bowtie2 版本错误"; FAILED=$((FAILED + 1)); }
[[ "$DM_VER" == "2.1.11" ]] || { echo "  !! diamond 版本错误"; FAILED=$((FAILED + 1)); }

echo
if [[ "$FAILED" -eq 0 ]]; then
    echo "=== 全部通过。启用: conda activate ${ENV_NAME} ==="
else
    echo "=== 有 ${FAILED} 项失败, 见上。 ==="
    exit 1
fi
