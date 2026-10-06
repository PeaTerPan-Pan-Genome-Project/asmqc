# Patched copy of src/runLR.sh of CRAQ 1.10 (https://github.com/JiaoLaboratory/CRAQ,
# MIT License, Copyright (c) 2023 JiaoLaboratory; see ../LICENSE.CRAQ), used by
# asmqc M9. Changed (marked "asmqc" below): the per-base depth table is
# streamed to its readers instead of written, and the BAM filter, the clip and
# indel scans and the depth stream run per segment of whole sequences in
# parallel. Everything else as in CRAQ.

src=`cd $(dirname $0); pwd -P`
#echo "$src"
pipline=$(basename $0)
LRname="LR"
minclip_num=2
lhe_cutoff_left=0.4
lhe_cutoff_right=0.6
lrbk_cutoff=0.65
LRavg_depth=100
max_depratio=0.05
mapquality=20
minindel=40
next_clip_dis=50000
x="map-hifi";
t=5
report_SNV="F"

for com in perl samtools minimap2 
do
        mg=$(command -v $com)
        if [ "$mg" == "" ]
        then
                echo -e "\n\tError: Command $com is NOT in you PATH. Please check.\n"
                exit 1
        fi
done


Usage="\nUsage:\n\t$pipline -g  Genome.fa -z  Genome.fa.size -1 SMS_sorted.bam -m minclip_num -q mapq -f lhe_cutoff_left -h lhe_cutoff_right -r lrbk_cutoff \nor\t$pipline -g  Genome.fasta -z  Genome.fasta.size -1 SMS.fa.gz -x map-hifi -m minclip_num -q mapq -f lhe_cutoff_left -h lhe_cutoff_right -r lrbk_cutoff -n max_depratio"

while getopts "a:g:x:z:1:d:m:q:f:h:d:r:t:v:" opt
do
    case $opt in
        g)      ref_fa=$OPTARG ;;
        z)      ref_fa_size=$OPTARG;;
        1)      inquery=$OPTARG;;
	m)	minclip_num=$OPTARG;;
	a)	LRavg_depth=$OPTARG;;
	d)	next_clip_dis=$OPTARG;;
	f)	lhe_cutoff_left=$OPTARG;;
	h)	lhe_cutoff_right=$OPTARG;;
	r)	lrbk_cutoff=$OPTARG;;	
	q)	mapquality=$OPTARG;;
	n)	max_depratio=$OPTARG;;
	v)	report_SNV=$OPTARG;;
	t)      t=$OPTARG;;
	x)	x=$OPTARG;;
        ?)
        echo ":| WARNING: Unknown option. Ignoring: Exiting!"
        exit 1;;
    esac
done
  

 if [ ! -e "$ref_fa" ]
then
       echo -e "\n\tgenome.fa is not found,  please check the README.md for the requirements of input files!\n\t$Usage \n"
       exit 1
fi

if [ ! -e "$ref_fa_size" ]
then
       echo -e "\n\tGenome.fasta.size is not found,  please check the README.md for the requirements of input files!\n\t$Usage \n"
       exit 1
fi

if [ `echo "$minclip_num < 1"|bc` -eq 1 ] ; then
        echo -e "\n\t  minclip_num ERROR: $minclip_num  Exit !"
        exit 1
fi

if [ `echo "$lhe_cutoff_left < 0"|bc` -eq 1 ] ; then
        echo -e "\n\t min_bkrate ERROR: $lhe_cutoff_left  Exit !"
        exit 1
fi


if [ `echo "$t <= 0"|bc` -eq 1 ] ; then
        echo -e "\n\tthread ERROR: $t  Exit !"
        exit 1
fi

if [ -d "LRout" ];then
   echo -e "Error::  LRout already exists, Exit !"	
	exit 1
fi
mkdir LRout
#echo "$inquery"
inquery_tmp=$(echo $inquery | tr [A-Z] [a-z])
if [[ "$inquery_tmp" =~ (fa$)|(fq$)|(fasta$)|(fastq$)|(fa.gz$)|(fq.gz$)|(fasta.gz$)|(fastq.gz$)|(bam$) ]]; then
        echo "worker_pipeline::"

     if [[ "$inquery" =~ (bam$) ]];then
		echo "Skipping alignment::"
                if [ ! -e "$inquery" ];then
                echo -e "\n\t $inquery is not found,  please check the README.md for the requirements of input files!\n\t$Usage \n"
                exit 1
                fi
                #if [ ! -e $inquery".bai" ];then
                #echo -e "\n\t $inquery".bai" is not found, cannot read index for $inquery \n\t$Usage \n"
                #exit 1
                #fi

     input_bam=$inquery
     echo -e "[M::worker_pipeline:: Filtering bamfiles]"
     # asmqc: filtered per segment of whole sequences (header order), in
     # parallel, then concatenated: the same records in the same order
     mkdir -p LRout/asmqc_parts
     idx_bam=$(readlink -f $input_bam)  # CRAQ links the BAM without its index
     segs=($(perl $src/asmqc_segments.pl $idx_bam $t LRout/asmqc_parts/seg)) || exit 1
     cmds=()
     for s in "${segs[@]}"; do
          cmds+=("samtools view -h -q $mapquality -F 1796 -M -L $s $idx_bam | perl $src/lrsam_cigar_filter.pl - | samtools view -h -S -b - -o ${s%.bed}.filter.bam")
     done
     perl $src/asmqc_par.pl $t "${cmds[@]}" || exit 1
     samtools cat -o LRout/$LRname"_sort.bam" "${segs[@]/%.bed/.filter.bam}" || exit 1
     rm -f "${segs[@]/%.bed/.filter.bam}"
     samtools index -@ $t LRout/$LRname"_sort.bam"
     fi

     if [[ "$inquery_tmp" =~ (fa$)|(fq$)|(fasta$)|(fastq$)|(fa.gz$)|(fq.gz$)|(fasta.gz$)|(fastq.gz$) ]]; then
	mkdir  -p LRout/tmp_bam
	array=(${inquery//,/ })   
	for query in ${array[@]}
         	do
		if [ ! -e "$query" ];then
        	echo -e "\n\t $query is not found,  please check the README.md for the requirements of input files!\n\t$Usage \n"
        	exit 1
		fi
	       done
	 for query in ${array[@]}
 	        do
	 	minimap2 -ax $x -R "@RG\tID:foo\tSM:bar1\tLB:lib" $ref_fa  $query -t $t | samtools view -h  -S -b -@ $t -  -o LRout/tmp_bam/$query"_unsort.bam"
        	samtools sort LRout/tmp_bam/$query"_unsort.bam"  -@ $t -o LRout/tmp_bam/$query"_sort.bam"
        	samtools index LRout/tmp_bam/$query"_sort.bam"
		rm LRout/tmp_bam/$query"_unsort.bam"
		done

      	 echo -e "[M::worker_pipeline:: Merge bamfiles]"
       	 samtools merge -@ $t -f LRout/tmp_bam/tmp_lr_merge.bam  LRout/tmp_bam/*_sort.bam 
	 rm  LRout/tmp_bam/*_sort.bam*
       	 input_bam="LRout/tmp_bam/tmp_lr_merge.bam"
       	 echo -e "[M::worker_pipeline:: Filtering bamfiles]"
         samtools view -h -q $mapquality -F 1796  $input_bam -@ $t | perl $src/lrsam_cigar_filter.pl - | samtools view -h -S -b -@ $t -  -o LRout/$LRname"_sort.bam"
         samtools index LRout/$LRname"_sort.bam"
	rm -r LRout/tmp_bam/
     fi	 
     
     else
     echo -e "Error, query should be fasta, fq or a sorted.bam  \n!"
     exit 1
fi


# asmqc (see header): no per-base depth table (LR_sort.depth) is written. Clip
# and indel sites come from the BAM alone and are extracted first; a
# "samtools depth -a" stream then feeds every former reader of the table at
# once. File names and contents are as in CRAQ 1.10, plus LR_depth.seqs (the
# sequences in the stream) for the region queries in runAQI.sh.
# Per segment of whole sequences, in parallel ($t jobs): the clip and indel
# scans, then the depth stream. Parts are joined in segment (= header =
# depth table) order; asmqc_merge.pl restores the genome-wide layout where a
# script's output is not in table order.
segs=($(perl $src/asmqc_segments.pl LRout/$LRname"_sort.bam" $t LRout/asmqc_parts/seg)) || exit 1
echo -e "[M::worker_pipeline:: Collect potential CRE|H]"
if [ "$report_SNV" != "T" ] ; then
	# renamed to LR_DI.cov.dimin3.tmp below, after "rm LRout/*tmp" (as in CRAQ)
	di_in=LRout/asmqc_DI.cov.dimin3; di_out=LRout/$LRname"_DI.cov.dimin3.tmp.dep"; di_min=3
else
	di_in=LRout/$LRname"_DI.cov"; di_out=LRout/$LRname"_DI.covRate.filter.all"; di_min=1
fi
echo -e "[M::worker_pipeline:: Extract SMS clipping signal]"
cmds=()
for s in "${segs[@]}"; do
	b=${s%.bed}
	cmds+=("samtools view -M -L $s LRout/${LRname}_sort.bam | perl $src/caculate_breakpoint_depth.pl - >$b.clip"
	       "samtools view -M -L $s LRout/${LRname}_sort.bam | perl $src/caculate_clipDI_cov.pl - $di_min >$b.di")
done
perl $src/asmqc_par.pl $t "${cmds[@]}" || exit 1
cat "${segs[@]/%.bed/.clip}" > LRout/$LRname"_clipped.cov"
perl $src/asmqc_merge.pl dici "${segs[@]/%.bed/.di}" > $di_in || exit 1
	perl -alne  'print if($F[3]>='$minclip_num')' LRout/$LRname"_clipped.cov" >LRout/$LRname"_clipped.cov.tmp"

echo -e "[M::worker_pipeline:: Compute effective SMS coverage]"
cmds=()
for s in "${segs[@]}"; do
	b=${s%.bed}
	cmds+=("samtools view -h -M -L $s LRout/${LRname}_sort.bam | samtools depth -a - | perl $src/asmqc_fanout.pl \\
		'perl $src/LReffect_size.pl /dev/stdin $LRavg_depth $max_depratio $b.nonmap >$b.eff' \\
		'perl $src/synthesize_LRbkdep_and_alldep.pl LRout/${LRname}_clipped.cov.tmp /dev/stdin >$b.bk' \\
		'perl $src/synthesize_clipDIcov_and_alldep.pl $di_in /dev/stdin >$b.dd' \\
		'cut -f1 | uniq >$b.seqs'")
done
perl $src/asmqc_par.pl $t "${cmds[@]}" || { echo -e "Error:: long-read depth stream failed, Exit !" >&2; exit 1; }
perl $src/asmqc_merge.pl effsize "${segs[@]/%.bed/.eff}" > LRout/$LRname"_eff.size"
cat "${segs[@]/%.bed/.nonmap}" > LRout/Nonmap.loc
cat "${segs[@]/%.bed/.bk}" > LRout/$LRname"_clip.coverRate"
cat "${segs[@]/%.bed/.dd}" > $di_out
cat "${segs[@]/%.bed/.seqs}" > LRout/$LRname"_depth.seqs"
rm -r LRout/asmqc_parts

echo -e "[M::worker_pipeline:: Collect potential CSE|H]"
	perl -alne  'print if($F[4]< 2*'$LRavg_depth' && $F[3]>='$minclip_num' && $F[3]/$F[4]>'$lhe_cutoff_left')' LRout/$LRname"_clip.coverRate" >LRout/$LRname"_clip.coverRate.filter"

#get putative.HR
	perl -alne 'print if($F[3]/$F[4]<='$lhe_cutoff_right'  )' LRout/$LRname"_clip.coverRate.filter" >LRout/$LRname"_clip.coverRate.filter.SH"
        perl $src/LER_softclip_filter.pl LRout/$LRname"_clipped.cov"  LRout/$LRname"_clip.coverRate.filter.SH" $next_clip_dis >LRout/$LRname"_putative.SH.tmp"
        perl -alne  'print if($F[5]>0.5)' LRout/$LRname"_putative.SH.tmp" |cut  -f -5 >LRout/$LRname"_putative.SH"

#get putative.ER
        perl -alne 'print if($F[3]/$F[4]>='$lrbk_cutoff'  )' LRout/$LRname"_clip.coverRate.filter" >LRout/$LRname"_clip.coverRate.filter.SE"
        perl $src/LER_softclip_filter.pl LRout/$LRname"_clipped.cov"  LRout/$LRname"_clip.coverRate.filter.SE" $next_clip_dis >LRout/$LRname"_putative.SE.tmp"
        perl -alne  'print if($F[5]<0.1)' LRout/$LRname"_putative.SE.tmp" |cut  -f -5 > LRout/$LRname"_putative.SE"

cat LRout/$LRname"_putative.SE" LRout/$LRname"_putative.SH" >LRout/$LRname"_putative.SE.SH"
rm LRout/*tmp  LRout/LR_clip.coverRate.filter.S*

if [ "$report_SNV" != "T" ] ; then
	mv LRout/asmqc_DI.cov.dimin3 LRout/$LRname"_DI.cov.dimin3.tmp"
	perl -alne  'print if($F[5]>='$minindel')' LRout/$LRname"_DI.cov.dimin3.tmp.dep" > LRout/$LRname"_DI.covRate.filter.R"
fi

if [ "$report_SNV" == "T" ] ; then
	perl -alne 'print if($F[5]>= '$minindel')' LRout/$LRname"_DI.covRate.filter.all" >LRout/$LRname"_DI.covRate.filter.R"

	perl -alne 'print if($F[5] < '$minindel')' LRout/$LRname"_DI.covRate.filter.all" >LRout/$LRname"_DI.covRate.filter.SNV"
        perl -alne  '$a=$F[3]/$F[4]; print if($a >=0.9 && $F[7]>0.7  )' LRout/$LRname"_DI.covRate.filter.SNV" |perl -alne 'print "$F[0]\t$F[1]\t$F[3]\t$F[4]\t$F[5]$F[6]"' >LRout/out_final_indel.err
        perl -alne  '$a=$F[3]/$F[4]; print if($a >='$lhe_cutoff_left' && $a <='$lhe_cutoff_right' &&  $F[7]>0.7  )' LRout/$LRname"_DI.covRate.filter.SNV" |perl -alne 'print "$F[0]\t$F[1]\t$F[3]\t$F[4]\t$F[5]$F[6]"' >LRout/out_final_indel.het
	
	perl -alne  'print if($F[5]>=3)' LRout/$LRname"_DI.covRate.filter.all" >LRout/$LRname"_DI.cov.dimin3.tmp.dep"
	rm LRout/$LRname"_DI.covRate.filter.all"	
fi


#echo -e "\n##########################################################################################\n"
echo -e "LR clipping analysis completed. Check current directory LRout for final results!\n"
