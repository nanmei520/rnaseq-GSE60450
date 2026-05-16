#!/bin/bash
# ================================================
# featurecounts.sh
# 功能：统计每个基因的 read 计数，生成表达矩阵
# 工具：featureCounts（subread 包）
# 输入：bam/*_Aligned.sortedByCoord.out.bam
# 输出：counts/all_samples_counts.txt
#       counts/counts_matrix.txt（DESeq2 直接输入）
# ================================================

cd ~/rnaseq_practice

echo "=== 开始 featureCounts 基因计数 ==="

# 整理所有 BAM 文件路径
BAM_FILES=$(ls bam/*_Aligned.sortedByCoord.out.bam | tr '\n' ' ')

featureCounts \
    -a ref/gencode.vM36.chr19.gtf \         #注释GTF文件
    -o counts/all_samples_counts.txt \        #输出文件
    -T 4 \                                            #线程数
    -s 0 \                                             #0是非链特异性（我用的数据GSE60450是非链特异性）                 1是正链特异性，2是反链特异性
    -Q 10 \                                            #设置最低MAPQ质量分数，过滤掉低质量比对
    --minOverlap 1 \                         #设置read与exon至少有1个碱基对重叠才计数
    --largestOverlap \                      #设置多基因重叠时归属于重叠最多的那个基因
    ${BAM_FILES}                             #所有的BAM文件

echo "=== featureCounts 完成 ==="
echo "=== 整理计数矩阵（提取 GeneID + 计数列）==="

# 去掉注释行和注释列，只保留 GeneID + 各样本计数
grep -v '^#' counts/all_samples_counts.txt | \
    cut -f1,7- > counts/counts_matrix.txt

echo "=== 计数矩阵整理完成 ==="
wc -l counts/counts_matrix.txt       查看基因数量

