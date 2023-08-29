# SASA calculation


Calculate SASA with Gromacs for one file

```bash
gmx sasa -f .xtc -s Fab.pdb -n .ndx -oa atomarea.xvg -or resarea.xvg -tv volumn.xvg -o area.xvg
```



Calculate SASA with GROMCAS by loop

```bash
for i in *.xtc; # find all the files that end with .xtc
do f="${i%.*}"; # extract the filename without extensions
	echo 1| gmx sasa -f $i -s Fab.pdb -oa "${f}_atomarea.xvg" -or "${f}_resarea.xvg" -tv "${f}_volumn.xvg" -o "${f}_area.xvg"|| break; #calculate SASA for each XTC file
done
```

Reference: 

Loop through all the files with a certain extension: https://stackoverflow.com/questions/14505047/loop-through-all-the-files-with-a-specific-extension

Extract file extension: https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash

https://stackoverflow.com/questions/2664740/extract-file-basename-without-path-and-extension-in-bash



script version - kathleen

```bash
#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=2:00:0
#$ -l mem=2G
#$ -N 
#$ -pe mpi 40
#$ -cwd
#$ -m beas

module load compilers/intel/2018/update3
module load mpi/intel/2018/update3/intel
module load gromacs/2019.3/intel-2018

for i in *.xtc; # find all the files that end with .xtc
do f="${i%.*}"; # extract the filename without extensions
	echo 1| gmx sasa -f $i -s Fab.pdb -oa "${f}_atomarea.xvg" -or "${f}_resarea.xvg" -tv "${f}_volumn.xvg" -o "${f}_area.xvg"|| break; #calculate SASA for each XTC file
done
```



Make the index file

```bash
gmx make_ndx -f Fab.pdb -o index.ndx
```

or 

```bash
gmx make_ndx -f md_0_1.gro -o index.ndx
```

Depends on whether you want to include water and salts as groups



Calculate non-polar SASA with the manually created index file -12/06/23

```bash
for i in *.xtc;
do f="${i%.*}";
echo 15| gmx sasa -f $i -s Fab.pdb -n index.ndx -oa "${f}_atomarea.xvg" -or "${f}_resarea.xvg" -tv "${f}_volumn.xvg" -o "${f}_area.xvg"|| break;
done
```



script

```bash
#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=4:00:0 # 4 hours is sufficient
#$ -l mem=2G
#$ -N SASAnonpolar_54runs_120623
#$ -pe mpi 80
#$ -cwd
#$ -m beas

module load compilers/intel/2018/update3
module load mpi/intel/2018/update3/intel
module load gromacs/2019.3/intel-2018

for i in *.xtc;
do f="${i%.*}";
echo 15| gmx sasa -f $i -s Fab.pdb -n index.ndx -oa "${f}_atomarea.xvg" -or "${f}_resarea.xvg" 
-tv "${f}_volumn.xvg" -o "${f}_area.xvg"|| break;
done
```

