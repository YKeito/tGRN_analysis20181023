cd /home/yasue/bigdata/yasue/tGRN_Groping/inflation4/motif/multi-fasta/target
ls > /home/yasue/bigdata/yasue/motif/allfasta/inflation4_allfasta.txt
for i in `cat /home/yasue/bigdata/yasue/motif/allfasta/inflation4_allfasta.txt`
do
ame --o ${i}_results --evalue-report-threshold 20.0 --control /home/yasue/bigdata/yasue/tGRN_Groping/inflation4/motif/multi-fasta/control/control_${i} /home/yasue/bigdata/yasue/tGRN_Groping/inflation4/motif/multi-fasta/target/${i} /home/yasue/bigdata/yasue/MEME/Databases/motif_databases/CIS-BP/Arabidopsis_thaliana.meme
done
find /home/yasue/bigdata/yasue/tGRN_Groping/inflation4/motif/multi-fasta/target/ -name "*_results" | xargs -I% mv % /home/yasue/bigdata/yasue/tGRN_Groping/inflation4/motif/motif_results/
