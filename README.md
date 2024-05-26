This is the github repository for reproducing the analyses performed in the paper "Compositional Differential Abundance 
Testing: Defining and Finding a New Type of Health-Microbiome Associations" ((Ma, Huttenhower, and Janson, 2024).

* Each analysis script is organized by sequence under `mds/`.

* For the simulation evaluation scripts (mds/1.0-4.0), these were designed to run large-scale benchmarking using the
  Harvard research computing environment (https://www.rc.fas.harvard.edu). These use the `batchtools` R package to
  interact with the job scheduling system.

      * It is unlikely that a regular desktop or laptop computer will have the sufficient computation power to con
