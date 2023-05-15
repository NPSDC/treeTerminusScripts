CMD=$1
DIR=$2
# for i in `seq 1 4`; 
# do
    bash $CMD -l $DIR/sample_01_1.fasta -r $DIR/sample_01_2.fasta;
    bash $CMD -l $DIR/sample_02_1.fasta -r $DIR/sample_02_2.fasta;
    rm $DIR/*.fasta
    pigz -1 -p 4 $DIR/sample_01_1_shuffled.fa;
    pigz -1 -p 4 $DIR/sample_01_2_shuffled.fa;
    pigz -1 -p 4 $DIR/sample_02_1_shuffled.fa;
    pigz -1 -p 4 $DIR/sample_02_2_shuffled.fa;
# done
