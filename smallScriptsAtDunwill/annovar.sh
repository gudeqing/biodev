annovar="/nfs2/software/annovar/ANNOVAR_2019.10.24/table_annovar.pl"
annovardb="/nfs2/database/annovar/humandb/"
 ${annovar} $1 $annovardb -buildver hg19 -out $1 -remove \
-protocol refGene,cytoBand,avsnp150,snp138NonFlagged,AVENIO,EpioneHS,cosmic90,Oncomine,CGI,CBMDB,CIVIC,DOCM,CHASMplus,IntOGen,TCGA_PCDM,TCGA,icgc21,CancerHotspots,nci60,clinvar_20200316,esp6500siv2_all,1000g2015aug_all,exac03,gnomad_exome,gme,cg69,hrcr1,intervar_20180118,dbnsfp33a \
-operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
--argument '-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs' \
-nastring . -vcfinput --dot2underline --thread 12 --maxgenethread 20
