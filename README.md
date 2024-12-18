# QracBounds
This repository contains MATLAB programs used for generating numerical results from [arXiv:2312.14142](https://arxiv.org/abs/2312.14142). 

The files:
- RAC333.m and RAC333_run.m are programs used to generate the upper bounds from Table I for quantum RAC with d=D=n=3.
- SeeSawQRAC.m is a general code for the seesaw algorithm for quantum RAC. 

The prerequisites:
- [YALMIP](https://yalmip.github.io/) toolbox
- An SDP solver, e.g., [mosek](https://www.mosek.com/)
- [QDimSum](https://denisrosset.github.io/qdimsum/) package
