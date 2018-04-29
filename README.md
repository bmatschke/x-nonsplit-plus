# X-nonsplit-plus

We compute the integral points of the moduli curve X<sub>ns</sub><sup>+</sup>(p), where p is a prime. In particular we find that for 11 <= p <= 97 the only integral points are CM points, i.e. they correspond to elliptic curves with complex multiplication.
This project is motivated by Serre’s uniformity problem for Galois representations, which would be solved if one could show that for large p the sets X<sub>ns</sub><sup>+</sup>(p)(Q) and X<sub>ns</sub><sup>+</sup>(p)(Q) consist only of cusps and CM points.

### Contents

 - *nonsplit.sage* - Source code for computing integral points on the moduli curve X<sub>ns</sub><sup>+</sup>(p), written in [sage](https://www.sagemath.org/).

 - *data/* - Data files that contain the computed non-cusp integral points.

 - *slurm-scripts/* - Contains slurm scripts to run nonsplit.sage in parallel; see its README.md.

### Authors

Aurelien Bajolet, Yuri Bilu, Benjamin Matschke.

### License

Creative Commons BY-NC 4.0.

### Reference

The code is based on the authors' paper [Computing integral points on X_ns^+(p)](https://arxiv.org/abs/1212.0665).
This repository is published on [github](https://github.com/bmatschke/x-nonsplit-plus/).

