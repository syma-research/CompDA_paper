This is the github repository for reproducing the analyses performed in the paper "Compositional Differential Abundance 
Testing: Defining and Finding a New Type of Health-Microbiome Associations" (Ma, Huttenhower, and Janson, 2024).

* The analyses involved here requires the `CompDA` R package, which can be found in [this companion software repo](https://github.com/syma-research/CompDA).

* Each analysis script is organized by sequence under `mds/`.

* Scripts `mds/1.0` - `mds/4.0` conduct the simulation analyses. Scripts `mds/1.2`, `mds/2.1`, `mds/3.1`, and `mds/4.1`
  summarize findings from the analyses and generate figures. Script `mds/5.0` conducts the real-world (CRC) analysis.

* For the simulation evaluation scripts (mds/1.0-4.0), these were designed to run large-scale benchmarking using the
  Harvard research computing environment (https://www.rc.fas.harvard.edu). These use the `batchtools` R package to
  interact with the job scheduling system.

    * It is unlikely that a regular desktop or laptop computer will have the sufficient computation power to conduct
      the full set of simulations. However, you may choose to select a subset of the simulation scenarios
      (specified with the `tb_job` object in each script) and evaluate them. To this end, the `one_job` function performs
      the relevant simulation analysis for each scenario specification (i.e., one "row" of `tb_job`).

    * If you have access to a high-performance computing environment, you should be able to similarly interact with
      `r_batch_tools` to perform the same analysis.

* `R/` has utility scripts used in the simulations. It also includes the relevant matlab code, directly copied from
  [Lu, Shi, and Li, Biometrics, 2019](https://pubmed.ncbi.nlm.nih.gov/30039859/), to evaluate the performance of the debiased
  generalized linear lasso method described in that paper.
