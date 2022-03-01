#!/usr/local_rwth/bin/zsh

### beginning of executable commands
cd $HOME/Dokumente/gbp_paper/gbp/sim

for run in {1..3}
do
        ./sim ./input/input_files/input_test.json output/tmp
done
