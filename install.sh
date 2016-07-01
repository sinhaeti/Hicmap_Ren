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

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "export PATH=\$PATH:$DIR/bin" >> ~/.bash_profile
bash ~/.bash_profile
