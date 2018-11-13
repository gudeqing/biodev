1. raw fastq
2. mapping to get bam
3. introduce mutation to raw fastq based on bam info
4. analysis with pipeline

simulate.py是模拟脚本的串联
simulate_mutation是模拟脚本，一次只能模拟一个样，当前使用v2版

read—base-stat2 是与模拟无关的脚本，用于统计read的不同位置的碱基组成比例
complex—sim是复杂模拟的脚本，CNV部分并未完整测试通过