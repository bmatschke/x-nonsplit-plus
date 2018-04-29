# Slurm scripts

Contains SLURM scripts that were used to run the main code nonsplit.sage in parallel on the super computer [PlaFRIM](https://www.plafrim.fr/).

### Examples

For prime p=7, run *sbatch nonsplit_for_7*, which uses 1 node with 24 CPUs. A testrun output the examplary file *nonsplit-p7-347904-output.txt* that can be found in this folder.

For larger primes, we used a clusters with 4 or 8 nodes, each of them having 24 CPUs. As there was a time limit of 3 days for each run, for p>=59 we had to cut the workload into further "chunks" for p>=59. E.g. for p=59, we run *sbatch nonsplit_for_59_inChunksAndParts_$i*, where $i denotes the chunk's index, here 0 or 1.

### Authors

Aurelien Bajolet, Yuri Bilu, Benjamin Matschke.

### License

Creative Commons BY-NC 4.0.

### Reference

The code is based on the authors' paper [Computing integral points on X_ns^+(p)](https://arxiv.org/abs/1212.0665).
This repository is published on [github](https://github.com/bmatschke/x-nonsplit-plus/).

