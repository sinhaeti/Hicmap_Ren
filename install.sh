cd htslib/
./configure    # Optional, needed for choosing optional functionality
make
cd ../
gcc -I htslib -L htslib src/pair2mates.c src/utils.c htslib/libhts.a -lz -lpthread -o bin/pair2mates
chmod +x bin/hicmap
chmod +x bin/chimeric.pl
chmod +x bin/chipseq
chmod +x bin/atacseq