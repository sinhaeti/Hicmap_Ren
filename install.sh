gcc -I htslib -L htslib src/pair2mates.c src/utils.c htslib/libhts.a -lz -lpthread -o bin/pair2mates
chmod +x bin/hicmap
chmod +x bin/chimeric.pl
tar -xzf data/mm9.MboI.500bp.tar.gz -C data/
tar -xzf data/hg19.Hind3.500bp.tar.gz -C data/
tar -xzf data/mm9.Hind3.500bp.tar.gz -C data/
