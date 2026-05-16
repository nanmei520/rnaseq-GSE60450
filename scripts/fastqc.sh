#!/bin/bash
# ================================================
# fastqc.sh
# 功能：对原始数据和清洗后数据进行质量控制
# 工具：FastQC MultiQC
# 输入：raw_data/*.fastq.gz 
# 输出：results/fastqc_raw/
# ================================================

cd ~/rnaseq_practice

# 创建输出目录
mkdir -p results/fastqc_raw
mkdir -p results/fastqc_clean

echo "=== 对原始数据进行 FastQC 质控 ==="
fastqc raw_data/*.fastq.gz \           #质控文件路径
    -o results/fastqc_raw/ \              #输出文件路径
    -t 4

echo "=== 原始数据 FastQC 完成 ==="
echo "=== 汇总原始数据质控报告 ==="

multiqc ~/rnaseq_practice/fastqc_raw/ -o ~/rnaseq_practice/results/multiqc_report/

echo "=== 清洗后数据 FastQC 汇总完成 ==="

