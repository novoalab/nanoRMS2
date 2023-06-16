# dockerfiles

```bash
cd ~/src
tar cpfz ~/public/src/nanoRMS2.tar.gz --exclude=*/{.git,docs,test,notebooks,**__pycache__,*.pyc} nanoRMS2/
# and make chmod a+r

tool=nanorms2 # has to be lowercase
cd ~/src/$tool/dockerfiles/3.6.1
version=$(../run --version 2> /dev/null)
guppyversion=$(basename $PWD)
name=$tool
echo $name:$version-guppy$guppyversion

docker build --pull -t lpryszcz/$name:$version .
docker tag lpryszcz/$name:$version-guppy$guppyversion lpryszcz/$name:latest
## needed to comment out os.setpgrp() lines in /run

docker push lpryszcz/$name:$version-guppy$guppyversion && docker push lpryszcz/$name:latest


# testing
cd ~/test/nanoRMS2/test
acc=oligopl; find $acc -name "*.bam"|xargs rm
docker run --gpus all -u $UID:$GID -v `pwd`:/data lpryszcz/nanorms2 \
     /opt/app/run -f /data/ref/ECOLI.fa -o /data/nanoRMS2/$acc \
     -i /data/$acc/C600_control_WGA/ /data/$acc/C600_native/ \
     -c dna_r9.4.1_450bps_pcr_fast.cfg --host /usr/bin/guppy_basecall_server \
     -b "NC_000913.3:1-100000"
```