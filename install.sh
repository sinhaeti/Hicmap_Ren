gcc -I htslib -L htslib src/pair2mates.c src/utils.c htslib/libhts.a -lz -lpthread -o bin/pair2mates
chmod u+x bin/hicmap
tar -xzf data/mm9.MboI.500bp.tar.gz -C data/