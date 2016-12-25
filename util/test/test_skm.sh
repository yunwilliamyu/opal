export PATH=../ext/gdl-1.1/GDL/bin:../ext/gdl-1.1/GDL/include:$PATH
export LD_LIBRARY_PATH=../ext/gdl-1.1/GDL/lib:$LD_LIBRARY_PATH

# draw fragments of size 20 covering each sequence 0.8 times
../drawfrag -i input/seq.fasta -t input/seq.taxid -l 8 -c 0.1 -o output/frags.fasta -g output/frags.gi2taxid 

# extract fragments taxids
cut -f 2 output/frags.gi2taxid > output/frags.taxid

# convert fragments to k-mers of size 6 in a format compliant with Vowpal Wabbit 
../fasta2skm -i output/frags.fasta -t output/frags.taxid -k 8 -o output/frags.vw -d output/vw-dico.txt -p input/patterns.txt



