# CompleteWGS
Pipeline for combining stLFR and PCR free data

 
1. On a Linux server, install singularity >= 3.8.1 with root on every nodes.
   
2. Download the singularity images (internet connection required) by the following commands:
 
cat <<EOF > CWGS.def
Bootstrap: docker
From: stlfr/cwgs:1.0.3
%post
    cp /90-environment.sh /.singularity.d/env/
EOF
 
singularity build --fakeroot CWGS.sif CWGS.def
# If the singularity doesn't support --fakeroot, you need sudo permission to run this command:
# sudo singularity build CWGS.sif CWGS.def
 
singularity exec -B`pwd -P` CWGS.sif cp -rL /usr/local/bin/CWGS /usr/local/bin/runit /usr/local/app/CWGS/demo .

3. Download the database (internet connection required) by this command:
 
./CWGS -createdb
 
This command will download around 32G data from internet and build index locally, which will occupy another 30G storage.
If users are using MegaBolt or ZBolt nodes, they should create database by ./CWGS -createdb –megabolt .
 
4. Test demo data:
 
cat << EOF > samplelist.txt
sample  stlfr1                      stlfr2                      pcrfree1                 pcrfree2
demo    demo/stLFR_demo_1M_1.fq.gz  demo/stLFR_demo_1M_2.fq.gz  demo/PF_demo_1M_1.fq.gz demo/PF_demo_1M_2.fq.gz
EOF
 
./CWGS samplelist.txt -local
 
Test demo data on clusters by SGE (Sun Grid Engine):
 
./CWGS samplelist.txt --queue mgi.q --project none
 
Test demo data on clusters by SGE (Sun Grid Engine) with MegaBolt/ZBolt nodes:
 
./CWGS samplelist.txt --queue mgi.q --project none --use_megabolt true --boltq fpga.q
