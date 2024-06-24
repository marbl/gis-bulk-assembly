#!/bin/sh
#
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --time=4-0
#SBATCH --output=./analysis.err
#SBATCH --job-name=va$samp
#

module load samtools
module load mashmap
module load bedtools
module load yak
module load minimap2

set -e
set -x

trap "rm -f analysis.jid" EXIT

cpus=$SLURM_CPUS_PER_TASK

refn="/data/Phillippy/t2t-share/assemblies/release/v2.0/chm13v2.0.fasta"
refc="/data/korens/devel/sg_sandbox/resources/reference.compressed.fasta"

samp=$1
root="/data/Phillippy2/projects/hprc-assemblies"

compleasm="/data/korens/devel/compleasm_kit/compleasm.py"
compledir=`dirname $compleasm`


#mkdir -p $root/assemblies-v2/$samp/$flav/analysis
#cd       $root/assemblies-v2/$samp/$flav/analysis

export REF_CACHE=$root/hprc-cache/samtools

#
#  This is getOutput.sh (it calls getComplete.sh at the end)
#    /data/Phillippy/projects/verkko/julian/HG02809/asm_trio/getOutput.sh

if [ ! -e ../assembly.fasta.fai ]; then
    echo "Index ../assembly.fasta"
    samtools faidx ../assembly.fasta
fi

if [ ! -e telomere ] ; then
    echo "Find telomere."
    /data/korens/devel/utils/telomere/find ../assembly.fasta > telomere
fi

#  Find telomeres
if [ ! -e telomere.bed ] ; then
    echo "Find telomere.bed."
    java -cp /data/korens/devel/utils FindTelomereWindows telomere 99.9 \
        | awk '{if ($NF > 0.5) print $2"\t"$4"\t"$5"\t"$3"\t"$NF}' \
        | sed s/\>//g \
        | bedtools merge -d -500 -i - -c 4 -o distinct \
        | bedtools sort -i - \
        | bedtools merge -i - > telomere.bed
fi


#  Find gaps
if [ ! -e assembly.gaps ] ; then
    echo "Find gaps."
    java -jar -Xmx4g $root/software/merqury/trio/fastaGetGaps.jar ../assembly.fasta assembly.gaps
fi

#
#  This is getComplete.sh
#    /data/Phillippy/projects/verkko/julian/HG02809/asm_trio/getComplete.sh
#

# get list of nodes with gaps
awk < assembly.gaps '{ print $1 }' | sort | uniq > nodes-with-gaps

# exclude those nodes from the full list of scaffolds in the asm
java -cp /data/korens/devel/utils SubFile nodes-with-gaps ../assembly.fasta.fai 0 -1 true | awk '{print $1}' > tmp2

# get list of things with >1 telomere
awk '{print $1}' < telomere.bed \
    | sort \
    | uniq -c \
    | awk '{if ($1 > 1) print $NF}' > tmp

# get list of nodes w/>1 telomere and no gaps
java -cp /data/korens/devel/utils SubFile tmp tmp2 > tmp3

# keep only those that have two telomeres near the ends
grep -w -f tmp telomere.bed | awk '{print $1}' > tmp4

grep -w -f tmp4 ../assembly.fasta.fai  \
    | awk '{print $1"\t"$2}'  \
    | sort -sk1,1 > tmp5

grep -w -f tmp telomere.bed  \
    | sort -sk1,1 > tmp4

join tmp4 tmp5  \
    | awk -v PREV="" '{if ($1 != PREV) { if (PREV != "" && C >=1 && E >=1 ) print PREV; PREV=$1; C=0; E=0 } if ($2 < 51000) C++; if ($3 + 51000 > $NF && $NF > 50000) { E++; } } END { if (C >=1 && E >= 1) print PREV}' > t2t_scfs

grep -w -f tmp3 telomere.bed | awk '{print $1}' > tmp4

grep -w -f tmp4 ../assembly.fasta.fai  \
    | awk '{print $1"\t"$2}'  \
    | sort -sk1,1 > tmp5

grep -w -f tmp3 telomere.bed| sort -sk1,1 > tmp4

join tmp4 tmp5  \
    | awk -v PREV="" '{if ($1 != PREV) { if (PREV != "" && C >=1 && E >=1 ) print PREV; PREV=$1; C=0; E=0 } if ($2 < 51000) C++; if ($3 + 51000 > $NF && $NF > 50000 ) { E++; } } END { if (C >=1 && E >= 1) print PREV}' > t2t_ctgs

echo ""
if [ -s t2t_ctgs ] ; then
    echo Ungapped two telomere: $( cat t2t_ctgs | wc -l )
    grep -w -f t2t_ctgs telomere.bed
else
    echo Ungapped two telomere: NONE
fi

echo ""
if [ -s t2t_scfs ] ; then
    echo Gapped two telomere: $( cat t2t_scfs | wc -l )
    grep -w -f t2t_scfs telomere.bed
else
    echo Gapped two telomere: NONE
fi

#
#  This is get.sh
#    /data/Phillippy/projects/verkko/julian/HG02809/asm_trio/get.sh
#

if [ ! -e assembly-ref.norm.mashmap ]; then
    mashmap -r $refn -q ../assembly.fasta --pi 95 -s 10000 -t 8 -o assembly-ref.norm.mashmap
fi

if [ ! -e assembly-ref.comp.mashmap ]; then
    mashmap -r $refc -q ../../verkko-base/contigs.fasta --pi 95 -s 10000 -t 8 -o assembly-ref.comp.mashmap
fi

g="."
if grep -q NC assembly-ref.norm.mashmap ; then
    g="NC"
fi

cg="."
if grep -q NC assembly-ref.comp.mashmap ; then
    cg="CNC"
fi


if grep -q mat- ../assembly.fasta ; then
    label1="mat-"
    label2="pat-"
elif grep -q h1tg ../assembly.fasta ; then
    label1="h1tg"
    label2="h2tg"
elif grep -q contig- ../assembly.fasta ; then
    label1="contig-"
    label2="none-ignore"
else
    label1="haplotype1-"
    label2="haplotype2-"
fi
echo "$label1 $label2 compNC: $cg regNC: $g"

echo  > nc-to-chr.sed 
#cho >> nc-to-chr.sed 's/NC_052532.1/chr1/g'    #  This seems to be mapping chicken IDs to chromosome numbers
#cho >> nc-to-chr.sed 's/NC_052533.1/chr2/g'
#cho >> nc-to-chr.sed 's/NC_052534.1/chr3/g'
#cho >> nc-to-chr.sed 's/NC_052535.1/chr4/g'
#cho >> nc-to-chr.sed 's/NC_052536.1/chr5/g'
#cho >> nc-to-chr.sed 's/NC_052537.1/chr6/g'
#cho >> nc-to-chr.sed 's/NC_052538.1/chr7/g'
#cho >> nc-to-chr.sed 's/NC_052539.1/chr8/g'
#cho >> nc-to-chr.sed 's/NC_052540.1/chr9/g'
#cho >> nc-to-chr.sed 's/NC_052541.1/chr10/g'
#cho >> nc-to-chr.sed 's/NC_052542.1/chr11/g'
#cho >> nc-to-chr.sed 's/NC_052543.1/chr12/g'
#cho >> nc-to-chr.sed 's/NC_052544.1/chr13/g'
#cho >> nc-to-chr.sed 's/NC_052545.1/chr14/g'
#cho >> nc-to-chr.sed 's/NC_052546.1/chr15/g'
#cho >> nc-to-chr.sed 's/NC_052547.1/chr16/g'
#cho >> nc-to-chr.sed 's/NC_052548.1/chr17/g'
#cho >> nc-to-chr.sed 's/NC_052549.1/chr18/g'
#cho >> nc-to-chr.sed 's/NC_052550.1/chr19/g'
#cho >> nc-to-chr.sed 's/NC_052551.1/chr20/g'
#cho >> nc-to-chr.sed 's/NC_052552.1/chr21/g'
#cho >> nc-to-chr.sed 's/NC_052553.1/chr22/g'
#cho >> nc-to-chr.sed 's/NC_052554.1/chr23/g'
#cho >> nc-to-chr.sed 's/NC_052555.1/chr24/g'
#cho >> nc-to-chr.sed 's/NC_052556.1/chr25/g'
#cho >> nc-to-chr.sed 's/NC_052557.1/chr26/g'
#cho >> nc-to-chr.sed 's/NC_052558.1/chr27/g'
#cho >> nc-to-chr.sed 's/NC_052559.1/chr28/g'
#cho >> nc-to-chr.sed 's/NC_052560.1/chr29/g'
#cho >> nc-to-chr.sed 's/NC_052561.1/chr30/g'
#cho >> nc-to-chr.sed 's/NC_052562.1/chr31/g'
#cho >> nc-to-chr.sed 's/NC_052563.1/chr32/g'
#cho >> nc-to-chr.sed 's/NC_052564.1/chr33/g'
#cho >> nc-to-chr.sed 's/NC_052565.1/chr34/g'
#cho >> nc-to-chr.sed 's/NC_052566.1/chr35/g'
#cho >> nc-to-chr.sed 's/NC_052567.1/chr36/g'
#cho >> nc-to-chr.sed 's/NC_052568.1/chr37/g'
#cho >> nc-to-chr.sed 's/NC_052569.1/chr38/g'
#cho >> nc-to-chr.sed 's/NC_052570.1/chr39/g'
#cho >> nc-to-chr.sed 's/NC_052571.1/chrW/g'
#cho >> nc-to-chr.sed 's/NC_052572.1/chrZ/g'
#cho >> nc-to-chr.sed 's/NC_053523.1/chrM/g'


if [ ! -e assembly.homopolymer-compressed.chr.csv ] ; then
    echo -e "node\tchr" > assembly.homopolymer-compressed.chr.csv
    for i in $( cat assembly-ref.comp.mashmap |awk '{if ($NF > 99) print $6}'|sort |uniq ); do
        echo "Chr $i"
        cat assembly-ref.comp.mashmap | \
            awk '{if ($NF > 99 && $4-$3 > 5000000) print $1"\t"$6"\t"$2}' | \
            grep -w $i |sort -srnk3,3 | \
            awk '{print $1"\t"$2}' | \
            sort | \
            uniq | \
            grep "$cg" | \
            sed -f nc-to-chr.sed >> assembly.homopolymer-compressed.chr.csv
    done

    cat assembly.homopolymer-compressed.chr.csv | awk '{print $1}' | sort | uniq > tmp4

    /data/korens/devel/sg_sandbox/gfacpp/build/neighborhood ../assembly.homopolymer-compressed.noseq.gfa tmp.gfa -n tmp4 --drop-sequence -r 1000

    cat tmp.gfa | grep "^S" |awk '{print $2}' > tmp4

    #second pass to get missing chr w/shorter matches
    cat assembly-ref.comp.mashmap | \
        grep -w -v -f tmp4 | \
        awk '{if ($NF > 99 && $4-$3 > 500000) print $1"\t"$6}' | \
        sort | \
        uniq | \
        grep "$cg" | \
        sed -f nc-to-chr.sed >> assembly.homopolymer-compressed.chr.csv
fi

if [ ! -e translation_hap1 -o ! -e translation_hap2 ] ; then
    cat assembly-ref.norm.mashmap | grep "$label1" | grep $g | awk '{if ($NF > 99 && $4-$3 > 1000000) print $1"\t"$6"\t"$2"\t"$7}' | sort | uniq | sed -f nc-to-chr.sed > translation_hap1
    cat assembly-ref.norm.mashmap | grep "$label2" | grep $g | awk '{if ($NF > 99 && $4-$3 > 1000000) print $1"\t"$6"\t"$2"\t"$7}' | sort | uniq | sed -f nc-to-chr.sed > translation_hap2
fi

if [ ! -e chr_completeness_max_hap1 -o ! -e chr_completeness_max_hap2 ] ; then
    cat translation_hap1|sort -k2,2 | awk '{if ($3 > 15000000) print $0}' | awk -v LAST="" -v S="" '{if (LAST != $2) { if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG; SUM=0; MAX=0; C=0; } LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1} END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}' | awk '{print $1"\t"$4}' | sort -nk1,1 -s > chr_completeness_max_hap1
    cat translation_hap2|sort -k2,2 | awk '{if ($3 > 15000000) print $0}' | awk -v LAST="" -v S="" '{if (LAST != $2) { if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG; SUM=0; MAX=0; C=0; } LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1} END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}' | awk '{print $1"\t"$4}' | sort -nk1,1 -s > chr_completeness_max_hap2
fi

rm -f tmp tmp? tmp.gfa

#
#  /data/walenzbp/hprc/software/marbl_utils/asm_evaluation/getYakStats.sh sh
#

#  Old hacked version of compleasm.
#
#for asm in assembly.haplotype1 assembly.haplotype2 ; do
#    for odb in `cd $root/hprc-cache/busco ; ls -d *odb10` ; do
#        faa="$root/hprc-cache/busco/$odb/refseq_db.faa.gz"
#
#        if [ ! -e "$faa" ] ; then
#            echo "Failed to find '$faa'."
#            exit 1
#        fi
#
#        if [ ! -e $asm.$odb.full_table.tsv ]; then
#            if [ ! -s $asm.$odb.aln.gff ]; then
#                $root/software/miniprot/miniprot -u --outs=0.95 -t$cpus --gff ../$asm.fasta $faa > $asm.$odb.aln.gff
#            fi
#
#            #  This came from minibusco dbf13d032cd6790c7d992f993abc3b604acc5cea
#            #  with a small bugfix.
#            python3 $root/hprc/analyze-miniprot.py \
#                    -g $asm.$odb.aln.gff \
#                    --full_table_file $asm.$odb.full_table.tsv \
#                    --complete_file $asm.$odb.summary.txt
#        fi
#    done
#done


for asm in assembly.haplotype1 assembly.haplotype2 ; do
    for db in primates_odb10 ; do
        if [ ! -e $asm.$db.full_table.tsv -o ! -e $asm.$db.summary.txt ]; then
            $compleasm run -t$cpus -l $db --library_path $compledir/mb_downloads -a ../$asm.fasta -o $asm \
            && \
            mv $asm/summary.txt        $asm.$db.summary.txt \
            && \
            mv $asm/$db/full_table.tsv $asm.$db.full_table.tsv
        fi
    done
done


for asm in assembly assembly.haplotype1 assembly.haplotype2 ; do
    if [ ! -e $asm.yak.qv -a -e $root/hprc-data/$samp/yakmers/ilmn.yak ] ; then
        yak qv       -t $cpus -l 100000 $root/hprc-data/$samp/yakmers/ilmn.yak ../$asm.fasta > $asm.yak.qv
    else
        touch $asm.yak.qv.NO_ILMN_DATA
    fi
    if [ ! -e $asm.yak.trioeavl -a -e $root/hprc-data/$samp/yakmers/mati.yak -a -e $root/hprc-data/$samp/yakmers/pati.yak ] ; then
        yak trioeval -t $cpus           $root/hprc-data/$samp/yakmers/mati.yak \
                                        $root/hprc-data/$samp/yakmers/pati.yak ../$asm.fasta > $asm.yak.trioeval
    else
        touch $asm.yak.trioeval.NO_MATI_or_PATI_DATA
    fi
done


touch analysis.complete
#  'analysis.jid' removed by trap on EXIT.
