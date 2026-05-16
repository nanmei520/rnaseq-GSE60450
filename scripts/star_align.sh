#!/bin/bash
# ================================================
# star_align.sh
# 功能：建立 STAR 索引并进行比对
# 工具：STAR
# 参考：chr19 单染色体（GRCm39 / GENCODE M36）
# 输出：star_index/（索引）、bam/（BAM文件）
# ================================================

cd ~/rnaseq_practice

# ─── 第一步：建立 STAR 索引 ─────────────────────
# 注意：--genomeSAindexNbases 11 适配 chr19 小基因组
# 完整基因组去掉此参数即可，STAR 会自动计算

echo "=== 开始建立 STAR 索引 ==="

STAR \
    --runMode genomeGenerate \              #建索引模式
    --genomeDir star_index/ \                     #索引文件的输出目录
    --genomeFastaFiles ref/chr19.fa \                     #参考的基因组文件
    --sjdbGTFfile ref/gencode.vM36.chr19.gtf \        #参考的基因组对应的注释文件
    --sjdbOverhang 99 \                                         #读取序列长度减一，我的是100，故填写99
    --genomeSAindexNbases 11 \                       #构建后缀数组时的基数，我用的小基因组，计算出对应的基数
    --runThreadN 4                                           #线程数

echo "=== 索引建立完成 ==="

# ─── 第二步：比对 ────────────────────────────────
# 循环处理全部 4 个样本

echo "=== 开始 STAR 比对 ==="

for SRR in SRR1552444 SRR1552445 SRR1552446 SRR1552447; do
    echo "--- 比对 ${SRR} ---"
    STAR \
        --runMode alignReads \               #比对模式
        --genomeDir star_index/ \               #索引目录位置
        --readFilesIn clean_data/${SRR}_clean.fastq.gz \   #输入质控后的FASTQ
        --readFilesCommand zcat \                             #告诉STAR输入文件是gzip压缩的
        --outSAMtype BAM SortedByCoordinate \          #输出坐标排序的BAM文件
        --outSAMattributes NH HI AS NM \                  #BAM中包含的额外标签字段
        --outFileNamePrefix bam/${SRR}_ \                  #设置输出文件前缀
        --runThreadN 4 \                                                   #线程数
        --outBAMsortingThreadN 2 \                           #BAM排序线程数
        --quantMode GeneCounts                                  #输出基因计数
    echo "${SRR} 比对完成"
done




# ─── 第三步：再次进行FastQC质控 ──────────────────────────
fastqc ~/rnaseq_practice/clean_data/\
          -o results/fastqc_clean/\
          -t 4


# ─── 第四步：再次对FastQC进行MultiQC汇总 ──────────────────────────
multiqc ~/rnaseq_practice/results/fastqc_clean/\
            -o results/fastqc_clean/
            


# ─── 第五步：BAM 建索引 ──────────────────────────
echo "=== 对 BAM 文件建立索引 ==="

for BAM in bam/*_Aligned.sortedByCoord.out.bam; do
    samtools index ${BAM}
done

echo "=== 全部比对完成 ==="

# ─── 查看比对统计 ────────────────────────────────
echo "=== 比对统计摘要 ==="
for SRR in SRR1552444 SRR1552445 SRR1552446 SRR1552447; do
    echo "--- ${SRR} ---"
    grep "Uniquely mapped reads %" bam/${SRR}_Log.final.out
done
