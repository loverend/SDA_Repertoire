#!/bin/bash
## Code to run SDA analysis on isotypr metrics
## Run for x number of iterations...! 

## BCR 
## 185 samples in total 
## 503 features
mkdir /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/TensorDecomposition/
for VARIABLE in {1..10}
do
mkdir /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/TensorDecomposition/${VARIABLE}
export OPENMP_NUM_THREADS=${NSLOTS:-1}
/apps/well/sda/1.1/./sda --data /gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/BCR_MATRIX_FOR_TENSOR_PRODUCTIVE_DATA.txt --N 185 --out /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/TensorDecomposition/${VARIABLE} --num_openmp_threads ${NSLOTS:-1} --num_blocks ${NSLOTS:-1} --eigen_parallel true --num_comps 100 --remove_zero_comps true --max_iter 3000 --impute_missing true 
done 

## 196 
## TCRGD
mkdir /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/TensorDecomposition
for VARIABLE in {1..10}
do
mkdir /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/TensorDecomposition/${VARIABLE}
export OPENMP_NUM_THREADS=${NSLOTS:-1}
/apps/well/sda/1.1/./sda --data /gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Summary/isotyper_metrics_filtered_FINAL_ALL_TENSOR_FORMAT_PRODUCTIVE.txt --N 194 --out /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/TensorDecomposition/${VARIABLE} --num_openmp_threads ${NSLOTS:-1} --num_blocks ${NSLOTS:-1} --eigen_parallel true --num_comps 50 --remove_zero_comps true --max_iter 3000 --impute_missing true 
done 

## TCRAB
mkdir /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/TensorDecomposition/
for VARIABLE in {1..10}
do
mkdir /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/TensorDecomposition/${VARIABLE}
export OPENMP_NUM_THREADS=${NSLOTS:-1}
/apps/well/sda/1.1/./sda --data /gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Summary/isotyper_metrics_filtered_FINAL_ALL_TENSOR_FORMAT_PRODUCTIVE.txt --N 193 --out /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/TensorDecomposition/${VARIABLE} --num_openmp_threads ${NSLOTS:-1} --num_blocks ${NSLOTS:-1} --eigen_parallel true --num_comps 50 --remove_zero_comps true --max_iter 3000 --impute_missing true 
done 

