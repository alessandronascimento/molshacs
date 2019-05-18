<h1> MolShaCS </h1>




<h2> Introduction </h2> 
MolShaCS (Molecular Shape and Charge Similarity) is a computational tool designed to evaluate the pairwise similarity among small molecules.

<h2> Details </h2> 
The software uses a Gaussian description of molecular shape and uses the very same concept for charge distribution description.

<h2> Usage </h2> 
MolShaCS can be compiled and used using a Qt graphical user interface or (my preferred method) using command line. For command line, type:

`MolShaCS Molshacs.inp`

Where <code> Molshacs.inp </code> is a text file with the following syntax:

```
refmol_mol2 mol1.mol2
output_prefix MolShaCS 
minimizer nlopt_mma 
align_molecules yes 
timeout 60 
write_coordinates yes 
mol2_aa no 
box_size 30.0 30.0 30.0 
multimol molecules.list 
step 1.0E-5 
tol 1.0E-4 
delta 1.0E-5
```

A file molecules.list should also exist in the directory where MolShaCS is running. This file must be a text file with the path for the comparing molecules, only. For example:

`$ more molecules.list ../mol2/mol1.mol2 ../mol2/mol2.mol2 ../mol2/mol3.mol2 ../mol2/mol4.mol2 ...`

<h2> Downloads </h2>
  
We currently distribute Windows binaries, including the Qt GUI. Also, the source code is available here in GitHub.


<h2> Sample Run </h2> 

As an example on how to use MolShaCS, we will compare aldosterone to a set of FDA approved molecules and take a look at the molecules in the top of the list. First of all, let’s get the molecules 

1. Let’s go to ZINC and download the Drugbank list of approved drugs with 1761 representative molecules. The molecules are provided as MOL2 files through script to download in Linux or Windows. 

2. Supposing you already have the mol2 files named ‘fda80.1.mol2’, fda80.2.mol2’, fda80.3.mol2’, etc, in a separate folder named ‘mol2’, lets generate a molecules.list file. If you use Linux or Cygwin in Windows, this should be easy (see below). These lines will look for the files in the folder mol2 and put them on the list if the file exists. 

```
$ for i in `seq 1 1761`; 
do 
  if [ -e ../mol2/fda80.$i.mol2 ]; 
  then 
    echo ../mol2/fda80.$i.mol2 >> molecules.list ; 
  fi; 
done 
```

3. Let’s get aldosterone from ZINC again, using this link. 

4. Now, let’s prepare an input file for MolShaCS. This file should have the following instructions (see below). Save the file as ‘input.inp’.

```
refmol_mol2 aldo.mol2 
output_prefix MolShaCS
minimizer nlopt_mma_sog
align_molecules yes
timeout 60 
write_coordinates yes
write_coord_threshold 0.85
mol2_aa no
box_size 30.0
multimol molecules.list
step 1.0E-5 
tol 1.0E-4 
delta 1.0E-5
```

5. Ok, we are ready to start the computation with the command <code> $MolShaCS input.inp </code>. Note that the file vdw.param distributed together with MolShaCS should be in the same folder where the computation takes place or you can use the environment variable MOLSHACS_DIR to point to the folder where the file is located. 

6. After a couple of minutes the calculation is done and two files are written: MolShaCS.log and MolShaCS.cc.dat. The latter has the similarities computed for each of the provided molecules. We can rank the results using a bash command again (below) and we will find the top scored molecules:

`$ more MolShaCS.cc.dat | sort –n –r –k 5 | more`
