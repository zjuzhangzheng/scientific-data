#amplicon data
wd=/c/amplicon
    db=/c/public
    PATH=$PATH:${db}/win
    cd ${wd}
mkdir -p temp
mkdir -p seq
gunzip seq/*.gz
ls -sh seq/
zless seq/KO1_1.fq.gz|head -n4
zless seq/KO1_1.fq | head -n20 | cut -c 1-60
seqkit stat seq/KO1_1.fq.gz
# 解压缩并写入新文件，不删除原压缩包
gunzip -c ${db}/usearch/rdp_16s_v18.fa.gz > ${db}/usearch/rdp_16s_v18.fa
seqkit stat ${db}/usearch/rdp_16s_v18.fa
head -n2 ${db}/usearch/rdp_16s_v18.fa
 #依照实验设计批处理并合并双端序列
 time for i in `tail -n+2 result/metadata.txt | cut -f 1`;do
      vsearch --fastq_mergepairs seq/${i}_1.fq.gz --reverse seq/${i}_2.fq.gz \
      --fastqout temp/${i}.merged.fq --relabel ${i}.
    done
#合并所有样品至同一文件
    cat temp/*.merged.fq > temp/all.fq
 #查看文件大小
 ls -lsh temp/all.fq
 # Cut primers and quality filter
    time vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 29 --fastq_stripright 18 \
      --fastq_maxee_rate 0.01 \
      --fastaout temp/filtered.fa
 # Dereplicate and cluster/denoise
vsearch --derep_fulllength temp/filtered.fa \
      --minuniquesize 8 --sizeout --relabel Uni_ \
      --output temp/uniques.fa 
#聚类OTU/去噪ASV
  usearch -cluster_otus temp/uniques.fa \
      -otus temp/otus.fa \
      -relabel OTU_
 # Reference-based chimera detect
 mkdir -p result/raw
 vsearch --uchime_ref temp/otus.fa \
      -db ${db}/usearch/rdp_16s_v18.fa \
      --nonchimeras result/raw/otus.fa
#vsearch生成特征表
time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 0.97 --threads 4 \
    	--otutabout result/raw/otutab.txt 
#去除质体和非细菌
vsearch --sintax result/raw/otus.fa \
      --db ${db}/usearch/rdp_16s_v18.fa \
      --sintax_cutoff 0.1 \
      --tabbedout result/raw/otus.sintax 
# Normlize by subsample
 mkdir -p result/alpha
    Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
      --depth 0 --seed 1 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt
    usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat
#计算多样性指数
usearch -alpha_div result/otutab_rare.txt \
      -output result/alpha/alpha.txt
 #稀释曲线：取1%-100%的序列中OTUs数量，2s
    #Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast / with_replacement / without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
    usearch -alpha_div_rare result/otutab_rare.txt \
      -output result/alpha/alpha_rare.txt \
      -method without_replacement
#Beta多样性

    #结果有多个文件，需要目录
    mkdir -p result/beta/
    #基于OTU构建进化树 Make OTU tree, 4s
    usearch -cluster_agg result/otus.fa -treeout result/otus.tree
    #生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
    usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
      -filename_prefix result/beta/
    #物种注释分类汇总
    #OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

    #OTU对应物种8列格式：注意注释是非整齐
    #生成物种表格OTU/ASV中空白补齐为Unassigned
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt

    #统计门纲目科属，使用 rank参数 p c o f g，为phylum, class, order, family, genus缩写
    mkdir -p result/tax
    for i in p c o f g;do
      usearch -sintax_summary result/otus.sintax \
      -otutabin result/otutab_rare.txt -rank ${i} \
      -output result/tax/sum_${i}.txt
    done
    sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
    # 列出所有文件
    wc -l result/tax/sum_*.txt
    head -n3 result/tax/sum_g.txt
    #有参定量特征表

    # 比对Greengenes，用于PICRUSt/Bugbase功能预测
    mkdir -p result/gg/
    #与GG所有97% OTUs比对，用于功能预测

 

    #  vsearch比对，更准更慢，但并行24-96线程更强
    vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fasta \
     --otutabout result/gg/otutab.txt --id 0.97 --threads 12
    比对率81.04%, 1核30m, 12核7m

    #统计
    usearch -otutab_stats result/gg/otutab.txt  \
      -output result/gg/otutab.stat
    cat result/gg/otutab.stat

#R语言多样性和物种分析
#Alpha多样性箱线图
Rscript ${db}/script/alpha_boxplot.R --alpha_index richness \
      --input result/alpha/vegan.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59
    # 使用循环绘制6种常用指数
    for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
      Rscript ${db}/script/alpha_boxplot.R --alpha_index ${i} \
        --input result/alpha/vegan.txt --design result/metadata.txt \
        --group Group --output result/alpha/ \
        --width 89 --height 59
    done

    # Alpha多样性柱状图+标准差
    Rscript ${db}/script/alpha_barplot.R --alpha_index richness \
      --input result/alpha/vegan.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59
     #稀释曲线

    Rscript ${db}/script/alpha_rare_curve.R \
      --input result/alpha/alpha_rare.txt --design result/metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59

#多样性维恩图  
# 四组比较，图和代码见输入文件目录，运行目录为当前项目根目录
    bash ${db}/script/sp_vennDiagram.sh \
      -f result/alpha/otu_group_exist.txt \
      -a WT -b KO -c OE -d All \
      -w 3 -u 3 \
      -p WT_KO_OE_All
#物种组成Taxonomy
 Rscript ${db}/script/tax_stackplot.R \
      --input result/tax/sum_p.txt --design result/metadata.txt \
      --group Group --output result/tax/sum_p.stackplot \
      --legend 5 --width 89 --height 59
    # 批量绘制输入包括p/c/o/f/g共5级
    for i in p c o f g; do
    Rscript ${db}/script/tax_stackplot.R \
      --input result/tax/sum_${i}.txt --design result/metadata.txt \
      --group Group --output result/tax/sum_${i}.stackplot \
      --legend 8 --width 89 --height 59; done
#构建进化树
 # 起始文件为 result/tree目录中 otus.fa(序列)、annotation.txt(物种和相对丰度)文件
    # Muscle软件进行序列对齐，3s
    muscle -in otus.fa -out otus_aligned.fas
# FastTree快速建树(Linux)
    # 注意FastTree软件输入文件为fasta格式的文件，而不是通常用的Phylip格式。输出文件是Newick格式。
    # 该方法适合于大数据，例如几百个OTUs的系统发育树！
    # Ubuntu上安装fasttree可以使用`apt install fasttree`
    fasttree -gtr -nt otus_aligned.fa > otus.nwk
