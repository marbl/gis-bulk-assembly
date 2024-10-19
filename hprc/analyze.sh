#!/bin/sh
#
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --time=4:00:00
#SBATCH --partition=norm,quick
#SBATCH --output=./analysis.err
#SBATCH --job-name=va${1}
#

module load samtools
module load mashmap
module load mash
module load bedtools
module load yak
module load minimap2
module load seqtk

set -e
set -x

trap "rm -f analysis.jid" EXIT

cpus=$SLURM_CPUS_PER_TASK

refn="/data/Phillippy/t2t-share/assemblies/release/v2.0/chm13v2.0.fasta"
refc="/data/korens/devel/sg_sandbox/resources/reference.compressed.fasta"
marbl_utils="/data/korens/devel/marbl_utils"

samp=$1
root="/data/Phillippy2/projects/hprc-assemblies"
base=`pwd`
base=`dirname $base |awk -F "/" '{print $NF}' |sed s/hi-c-test/base/g | sed s/hi-c/base/g |sed s/trio/base/g |sed s/thic/base/g`

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

sh $marbl_utils/asm_evaluation/getT2T.sh ../assembly.fasta

# remove rDNA and add telomere nodes
if [ ! -e assembly.homopolymer-compressed.add_telo.noseq.gfa ]; then
   ln -s ../assembly.colors.csv
   repeatUnit="/data/Phillippy/references/hg38/rDNA_compressed.fasta"
   cat ../assembly.homopolymer-compressed.gfa |awk '{if (match($1, "^S")) { print ">"$2; print $3}}' | mash sketch -i - -o sketch.msh
   mash screen sketch.msh $repeatUnit | awk '{if ($1 > 0.9 && $4 < 0.05) print $NF}' > rdna.nodes
   python $root/software-v4/verkko/lib/verkko/scripts/remove_nodes_add_telomere.py -t assembly.telomere.bed -g ../assembly.homopolymer-compressed.noseq.gfa -s ../assembly.scfmap -p ../assembly.paths.tsv -o assembly.homopolymer-compressed.add_telo.noseq.gfa -c assembly.colors.add_telo_add_rdna.csv
   python $root/software-v4/verkko/lib/verkko/scripts/remove_nodes_add_telomere.py -r rdna.nodes -t assembly.telomere.bed -g ../assembly.homopolymer-compressed.noseq.gfa -s ../assembly.scfmap -p ../assembly.paths.tsv -o assembly.homopolymer-compressed.add_telo_remove_rdna.noseq.gfa -c assembly.colors.add_telo_add_rdna.csv
   rm -f ./assembly.colors.csv
fi

#
#  This is get.sh
#    /data/Phillippy/projects/verkko/julian/HG02809/asm_trio/get.sh
#

if [ ! -e assembly-ref.norm.mashmap ]; then
    mashmap -r $refn -q ../assembly.fasta --pi 95 -s 10000 -t 8 -o assembly-ref.norm.mashmap
fi

if [ ! -e assembly-ref.comp.mashmap ]; then
    mashmap -r $refc -q ../../$base/contigs.fasta --pi 95 -s 10000 -t 8 -o assembly-ref.comp.mashmap
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
    label1="pat-"
    label2="mat-"
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

# this is chromosome assignment
isXY=`grep chrY translation_hap* |wc -l |awk '{print $1}'`
if [ $isXY -ne 0 ]; then
   yak sexchr -t 8 $root/hprc-data/chrY-no-par.yak $root/hprc-data/chrX-no-par.yak $root/hprc-data/par.yak ../assembly.haplotype1.fasta ../assembly.haplotype2.fasta > assembly.yak.sexchr
fi

if [ ! -e assembly.refOriented.fasta ]; then
   sh $marbl_utils/verkko_helpers/reorientByRef.sh assembly-ref.norm.mashmap > assembly-ref.reorient.tsv
   cat assembly-ref.reorient.tsv |awk '{print $1}' > tmp
   grep -w -v -f tmp ../assembly.fasta.fai |awk '{print $1"\t0\t"$2}' >>  assembly-ref.reorient.tsv
   rm ./tmp

   java -cp /data/korens/devel/utils:. SubFasta  assembly-ref.reorient.tsv ../assembly.fasta >  assembly.refOriented.fasta
   for i in `seq 1 2`; do
      parent=`echo $i |awk '{if ($1 == 2) print "mat"; else print "pat"}'`
      # we have XY then we use the assignment information to make sure chrX is is haplotype 2 (this is checking $6/$5 which is fraction of sex markers is hight and $7/($7+$8) is more Y chr than X markers while $8/($7+$8) is more X than Y
      if [ $isXY -ne 0 ]; then
         if [ $parent = "mat" ]; then
           cat assembly.yak.sexchr |grep "^S" | awk '{if ($6/$5 > 0.9 && $7/($7+$8) > 0.5) print $2}' > ignore.tmp
           cat assembly.yak.sexchr |grep "^S" | awk '{if ($6/$5 > 0.9 && $8/($7+$8) > 0.5) print $2}' > include.tmp
           grep -w -f include.tmp assembly-ref.reorient.tsv > tmp
         elif [ $parent = "pat" ]; then
           cat assembly.yak.sexchr |grep "^S" | awk '{if ($6/$5 > 0.9 && $8/($7+$8) > 0.5) print $2}' > ignore.tmp
           cat assembly.yak.sexchr |grep "^S" | awk '{if ($6/$5 > 0.9 && $7/($7+$8) > 0.5) print $2}' > include.tmp
           grep -w -f include.tmp assembly-ref.reorient.tsv > tmp
         fi
	  else
	     touch ./ignore.tmp
		 > tmp
      fi
      grep "_haplotype$i" assembly-ref.reorient.tsv >> tmp || true
      grep "_$parent" assembly-ref.reorient.tsv >> tmp || true
      grep "_homozygous" assembly-ref.reorient.tsv >> tmp || true
      grep -v chr  assembly-ref.reorient.tsv |grep "$parent" >> tmp || true
      grep -v chr  assembly-ref.reorient.tsv |grep "haplotype$i" >> tmp || true
      cat tmp | sort |uniq | grep -w -v -f ignore.tmp > tmp2
      java -cp /data/korens/devel/utils:. SubFasta tmp2 ../assembly.fasta > assembly.refOriented.haplotype$i.fasta

      rm -f ./tmp? ./ignore.tmp ./include.tmp
   done
fi

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


for asm in assembly.haplotype1 assembly.haplotype2 assembly.unassigned ; do
    for db in primates_odb10 ; do
        if [ ! -e $asm.$db.full_table.tsv -o ! -e $asm.$db.summary.txt ]; then
            $compleasm run -t$cpus -l $db --library_path $compledir/mb_downloads -a ../$asm.fasta -o $asm \
            && \
            mv $asm/summary.txt        $asm.$db.summary.txt \
            && \
            mv $asm/$db/full_table.tsv $asm.$db.full_table.tsv \
            && \
            rm -rf $asm/$db/hmmer_output
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
