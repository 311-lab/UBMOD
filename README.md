# UBMOD

## Introduction
This repository serves as the source code repository of the Unsaturated zone soil water Balance Model, known as UBMOD. The model a one-dimensional soil water balance model for simulating both downward and upward soil water movement. The model was designed by Wei Mao, Yan Zhu, Tianxing Zhao and Jinzhong yang. The [[Paper](https://www.sciencedirect.com/science/article/pii/S0022169418301562)] can be find here.

## Requirements
Fortran 90

## concept map


## Table of contents
index | name | content
-|-|-
A | Source code  | src|
B |  Demo file   | Rh1D.in, Rh1D.out|
C |  Document    | aaa|

### A. Source code
totally 8 files  
**Main_Waterbalance.f90**     
> The main program  

**parm.f90**  
> The parameters used in the model  

**ETp.f90**  
> The module of evapotranspiration  

**DateF.f90**
> The module about time

**input.f90**
> The input module

**preprocess.f90**
> The preprocess procedures

**majorprocess.f90**
> The major calculation procedures

**output.f90**
> The output module

### B. Demo file
There are totally 3 demoes. The Rh1D.in folder contains the input files, while the Rh1D.out folder contains the output files. The details of the demoes can be found in the technical report.

### C. Document
Technical report.

## Cited
Mao W, Yang J, Zhu Y, et al. An efficient soil water balance model based on hybrid numerical and statistical methods[J]. Journal of hydrology, 2018, 559: 721-735, doi: 10.1016/j.jhydrol.2018.02.074.     

## Contact us
weimao@whu.edu.cn  
2013301580339@whu.edu.cn  
zyan0701@163.com  
