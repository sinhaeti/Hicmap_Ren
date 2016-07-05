#!/usr/bin/env bash
set -e
chmod +x bin/hicmap
chmod +x bin/chimeric.pl
chmod +x bin/hicmap_cutter_sites_filter
chmod +x bin/hicmap_pair_up_filter
chmod +x bin/bedGraphToBigWig

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "export PATH=\$PATH:$DIR/bin" >> ~/.bash_profile
bash ~/.bash_profile
