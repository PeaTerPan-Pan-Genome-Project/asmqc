# Patched copy of src/runSR.sh of CRAQ 1.10 (https://github.com/JiaoLaboratory/CRAQ,
# MIT License, Copyright (c) 2023 JiaoLaboratory; see ../LICENSE.CRAQ), used by
# asmqc M9. Changed (marked "asmqc" below): the per-base depth table is
# streamed to its readers instead of written, and the BAM filter, the clip scan
# and the depth stream run per segment of whole sequences in parallel.
# Everything else as in CRAQ.

src=`cd $(dirname $0); pwd -P`
SRname="SR"

she_cutoff_left=0.25
she_cutoff_right=0.6
srbk_cutoff=0.75

SRavg_depth=100
mapquality=20
t=5
minclip_num=2

pipline=$(basename $0)

for com in perl minimap2 samtools 
do
        mg=$(command -v $com)
        if [ "$mg" == "" ]
        then
                echo -e "\n\tError: Command $com is NOT in you PATH. Please check.\n"
                exit 1
        fi
done


Usage="\nUsage:\n\t$pipline -g  Genome.fa  -z Genome.fa.size  -1 NGS_sort.bam \n or \t$pipline -g  Genome.fa  -z Genome.fa.size  -1 fq1.gz -2 fq2.gz -f she_cutoff_left -h she_cutoff_right -r srbk_cutoff \n\t[default: -q 20 -f 0.4 -h 0.6 -r 0.75 -m 2 -a 100 -t 5]"

while getopts "g:z:1:2:q:f:h:r:m:a:t:" opt
do
    case $opt in
        g)	ref_fa=$OPTARG ;;
        z)	ref_fa_size=$OPTARG ;;
        1)	query_1=$OPTARG ;;
        2)	query_2=$OPTARG ;;
	q)	mapquality=$OPTARG ;;
	m)	minclip_num=$OPTARG ;;
	f)	she_cutoff_left=$OPTARG ;;
	h)	she_cutoff_right=$OPTARG ;;
	r)	srbk_cutoff=$OPTARG ;;
	a)	SRavg_depth=$OPTARG ;;
        t)	t=$OPTARG ;;
        ?)
        echo ":| WARNING: Unknown option. Ignoring: Exiting!"
        exit 1;;
    esac
done


if [ ! -e "$ref_fa" ];then
        echo -e "\n\tgenome.fa is not found,  please check the README.md for the requirements of input files!\n\t$Usage \n"
        exit 1
        fi

if [ ! -e "$ref_fa_size" ]
then
       echo -e "\n\tGenome.fasta.size is not found,  please check the README.md for the requirements of input files!\n\t$Usage \n" 
       exit 1
fi
if [ `echo "$minclip_num < 0"|bc` -eq 1 ] ; then
        echo -e "\n\tminclip_num ERROR: $minclip_num  Exit !"
        exit 1
fi


if [ `echo "$t <= 0"|bc` -eq 1 ] ; then
        echo -e "\n\tthread ERROR: $t  Exit !"
        exit 1
fi

if [ `echo "$she_cutoff_left < 0"|bc` -eq 1 ] ; then
        echo -e "\n\tshe_cutoff_left ERROR: $she_cutoff_left  Exit !"
        exit 1
fi

query_1_tmp=$(echo $query_1 | tr [A-Z] [a-z])
query_2_tmp=$(echo $query_2 | tr [A-Z] [a-z])
if [[ "$query_1_tmp" =~ (fa$)|(fq$)|(fasta$)|(fastq$)|(fa.gz$)|(fq.gz$)|(fasta.gz$)|(fastq.gz$)|(bam$) ]]; then
	if [ -d "SRout" ];then
        echo -e "Error::  SRout already exists, Exit !"
        exit 1
        fi
        mkdir SRout
	if [[ "$query_1_tmp" =~ (bam$) ]];then
		if [ ! -e "$query_1" ];then
		echo -e "\n\t $query_1 is not found,  please check the README.md for the requirements of input files!\n\t$Usage \n"
                exit 1
                fi
                #if [ ! -e $query_1".bai" ];then
                #echo -e "\n\t $query_1".bai" is not found, cannot read index for $query_1 \n\t$Usage \n"
                #exit 1
                #fi

	        input_bam=$query_1
		echo -e "Skipping alignment::\n[M::worker_pipeline:: Filtering bamfiles]"
	# asmqc: filtered per segment of whole sequences (header order), in
	# parallel, then concatenated: the same records in the same order
	mkdir -p SRout/asmqc_parts
	idx_bam=$(readlink -f $input_bam)  # CRAQ links the BAM without its index
	segs=($(perl $src/asmqc_segments.pl $idx_bam $t SRout/asmqc_parts/seg)) || exit 1
	cmds=()
	for s in "${segs[@]}"; do
		cmds+=("samtools view -h  -F 1796 -M -L $s $idx_bam  -t $t | perl $src/srsam_cigar_filter.pl $mapquality - | samtools view -h -S -b -  -o ${s%.bed}.filter.bam")
	done
	perl $src/asmqc_par.pl $t "${cmds[@]}" || exit 1
	samtools cat -o SRout/$SRname"_sort.bam" "${segs[@]/%.bed/.filter.bam}" || exit 1
	rm -f "${segs[@]/%.bed/.filter.bam}"
		samtools index -@ $t SRout/$SRname"_sort.bam"
	fi
	
	if [[ "$query_1_tmp" =~ (fa$)|(fq$)|(fasta$)|(fastq$)|(fa.gz$)|(fq.gz$)|(fasta.gz$)|(fastq.gz$) ]]; then
		if [ ! -e "$query_1" ];then
                echo -e "\n\t $query_1 is not found,  please check the README.md for the requirements of input files!\n\t$Usage \n"
                exit 1
                fi
		if [ ! -e "$query_2" ];then
                echo -e "worker_pipeline::\nWARNING: Short pair_end read pair_2 is not found, only $query_1 used"
                fi
		if [[  -e "$query_2"   ]];then
			if [[ "$query_2_tmp" =~  (fa$)|(fq$)|(fasta$)|(fastq$)|(fa.gz$)|(fq.gz$)|(fasta.gz$)|(fastq.gz$) ]]; then
			echo -e "worker_pipeline::"
			else 
               		echo -e "\n\t $query_2 should be fastq/fa suffix,  please check the README.md for the requirements of input files!\n\t$Usage \n"
			exit 1
			fi	
                fi

	        echo -e "worker_pipeline:: NGS reads aligning and filtering"
                minimap2 -ax sr -R "@RG\tID:foo\tSM:bar1\tLB:lib" $ref_fa  $query_1 $query_2  -t $t | perl  $src/srsam_cigar_filter.pl $mapquality -  | samtools view -h -F 1796 -S -b -@ $t -  -o SRout/$SRname"_unsort.bam"

		echo -e "[M::worker_pipeline:: Sort bamfiles]"
      		samtools sort SRout/$SRname"_unsort.bam"   -@ $t -o SRout/$SRname"_sort.bam" 
      		rm SRout/$SRname"_unsort.bam" && samtools index SRout/$SRname"_sort.bam" 
	fi
else
     echo -e "Error, query should be fasta, fq or a sorted.bam  \n!"
     exit 1
fi
#echo -e "\n##########################################################################################\n"
# asmqc (see header): no per-base depth table (SR_sort.depth) is written.
# Clip sites come from the BAM alone and are extracted first; a depth stream
# then feeds every former reader of the table at once. Both run per segment of
# whole sequences in parallel ($t jobs); parts are joined in segment (= header
# = depth table) order, or by asmqc_merge.pl where a script's output is not in
# table order. File names and contents are as in CRAQ 1.10, plus
# SR_depth.seqs (the sequences in the stream) and SR_depth.mapq (its MAPQ
# filter) for the region queries in runAQI.sh.
echo "$mapquality" >SRout/$SRname"_depth.mapq"
mkdir -p SRout/asmqc_parts
segs=($(perl $src/asmqc_segments.pl SRout/$SRname"_sort.bam" $t SRout/asmqc_parts/seg)) || exit 1
if (($minclip_num >= 0)); then
     echo -e "[M::worker_pipeline:: Collect potential CRE|RH]"
#######start extract clipped reads
     cmds=()
     for s in "${segs[@]}"; do
          cmds+=("samtools view -q $mapquality -M -L $s SRout/${SRname}_sort.bam | perl $src/caculate_breakpoint_depth.pl - >${s%.bed}.clip")
     done
     perl $src/asmqc_par.pl $t "${cmds[@]}" || exit 1
     cat "${segs[@]/%.bed/.clip}" > SRout/$SRname"_clipped.cov"
     perl -alne  'print if($F[3]>='$minclip_num')' SRout/$SRname"_clipped.cov" >SRout/$SRname"_clipped.cov.tmp"
fi
echo -e "[M::worker_pipeline:: Compute effective NGS coverage]"
cmds=()
for s in "${segs[@]}"; do
     b=${s%.bed}
     c="samtools view -h -q $mapquality -M -L $s SRout/${SRname}_sort.bam | samtools depth -a - | perl $src/asmqc_fanout.pl 'perl $src/SReffect_size.pl /dev/stdin $b.nonmap >$b.eff' 'cut -f1 | uniq >$b.seqs'"
     if (($minclip_num >= 0)); then
          c="$c 'perl $src/synthesize_SRbkdep_and_alldep.pl SRout/${SRname}_clipped.cov.tmp /dev/stdin >$b.bk' 'perl $src/search_dep0.pl /dev/stdin >$b.dep0'"
     fi
     cmds+=("$c")
done
perl $src/asmqc_par.pl $t "${cmds[@]}" || { echo -e "Error:: short-read depth stream failed, Exit !" >&2; exit 1; }
perl $src/asmqc_merge.pl effsize "${segs[@]/%.bed/.eff}" > SRout/$SRname"_eff.size"
cat "${segs[@]/%.bed/.nonmap}" > SRout/Nonmap.loc
cat "${segs[@]/%.bed/.seqs}" > SRout/$SRname"_depth.seqs"
if (($minclip_num >= 0)); then
     cat "${segs[@]/%.bed/.bk}" > SRout/$SRname"_clip.coverRate"
     perl $src/asmqc_merge.pl dep0 "${segs[@]/%.bed/.dep0}" > SRout/$SRname"_putative.RE.0.tmp"
fi
rm -r SRout/asmqc_parts

if (($minclip_num >= 0)); then
     perl -alne  'print if($F[3]>='$minclip_num' && $F[4]>=1 && $F[4] < 2*'$SRavg_depth' && $F[3]/$F[4] >'$she_cutoff_left' )' SRout/$SRname"_clip.coverRate" >SRout/$SRname"_putative.RE.bk.tmp"
#######start extract nonmapping locus
     perl $src/srbk_merge_srnoncov.pl SRout/$SRname"_putative.RE.bk.tmp" SRout/$SRname"_putative.RE.0.tmp" |sort -k 1,1 -k 2,2n >SRout/$SRname"_putative.RE.RH"
     perl -alne  'print if( $F[3]/$F[4] <= '$she_cutoff_right' )' SRout/$SRname"_putative.RE.RH" >SRout/$SRname"_putative.RH"
     perl -alne  'print if( $F[3]/$F[4] >= '$srbk_cutoff' )' SRout/$SRname"_putative.RE.RH" >SRout/$SRname"_putative.RE"

     rm SRout/$SRname"_putative.RE.0.tmp" SRout/$SRname"_putative.RE.bk.tmp"

fi

#if (($minclip_num >=1)); then
#     echo -e "[M::worker_pipeline:: Collect potential CRE|H]"
#        samtools view -@ $t -q $mapquality SRout/$SRname"_sort.bam" | perl $src/caculate_breakpoint_depth.pl -  > SRout/$SRname"_clipped.cov"
#        perl -alne  'print if($F[3]>='$minclip_num')' SRout/$SRname"_clipped.cov" >SRout/$SRname"_clipped.cov.tmp"
#        perl   $src/synthesize_SRbkdep_and_alldep.pl  SRout/$SRname"_clipped.cov.tmp" SRout/$SRname"_sort.depth" >SRout/$SRname"_clip.coverRate"
#        perl -alne  'print if($F[3]>='$minclip_num' && $F[4]>=1 && $F[4] < 2*'$SRavg_depth' && $F[3]/$F[4] >'$she_cutoff_left')' SRout/$SRname"_clip.coverRate" > SRout/$SRname"_clip.coverRate.tmp"
 #       perl -alne  'print if( $F[3]/$F[4] <= '$she_cutoff_right' )' SRout/$SRname"_clip.coverRate.tmp" >SRout/$SRname"_putative.RH"
  #      perl -alne  'print if( $F[3]/$F[4] >= '$srbk_cutoff' )' SRout/$SRname"_clip.coverRate.tmp" >SRout/$SRname"_putative.RE"
   #     cat SRout/$SRname"_putative.RH" SRout/$SRname"_putative.RE" >SRout/$SRname"_putative.RE.RH"
    #    rm SRout/$SRname"_clipped.cov.tmp" SRout/$SRname"_clip.coverRate.tmp"

#fi









#echo -e "\n##########################################################################################\n"
echo -e "SR clipping analysis completed. Check current directory SRout for final results!"
