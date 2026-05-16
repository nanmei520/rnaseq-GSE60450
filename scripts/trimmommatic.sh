#!/bin/bash
# ================================================
# trimmomatic.sh
# 功能：去除接头序列和低质量碱基
# 工具：Trimmomatic（单端 SE 模式）
# 输入：raw_data/SRR155244*.fastq.gz
# 输出：clean_data/SRR155244*_clean.fastq.gz
# ================================================

cd ~/rnaseq_practice

# 定义接头文件路径（随 Trimmomatic 安装自带）
ADAPTERS=$(dirname $(which trimmomatic))/../share/trimmomatic/adapters/TruSeq3-SE.fa

echo "=== 开始 Trimmomatic 质量过滤 ==="

for SRR in SRR1552444 SRR1552445 SRR1552446 SRR1552447; do
    echo "--- 处理 ${SRR} ---"
    trimmomatic SE \                         
        raw_data/${SRR}.fastq.gz \             #输入文件
        clean_data/${SRR}_clean.fastq.gz \     #输出文件
        ILLUMINACLIP:${ADAPTERS}:2:30:10 \    #去接头命令 分别对应（错配数：palindrome阈值：简单clip阈值）
        LEADING:3 \                            #5’端Q<3的碱基删除掉
        TRAILING:3 \                            #3’端Q<3的碱基删除掉
        SLIDINGWINDOW:4:15 \        #滑动窗口4bp（4个碱基对），均值<15的截断
        MINLEN:36 \                          #过滤后长度<36bp的读段丢掉
        -threads 4                     #线程数
    echo "${SRR} 完成"
done
echo "=== Trimmomatic 全部完成 ==="

