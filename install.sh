#!/usr/bin/env bash
set -e

cd htslib/
make
cd ../
gcc -I htslib -L htslib src/pair2mates.c src/utils.c htslib/libhts.a -lz -lpthread -o bin/pair2mates
chmod +x bin/hicmap
chmod +x bin/chimeric.pl
chmod +x bin/chipseq
chmod +x bin/atacseq
chmod +x bin/hicmap_cutter_sites_filter
chmod +x bin/hicmap_pair_up_filter

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "export PATH=\$PATH:$DIR/bin" >> ~/.bash_profile
bash ~/.bash_profile
