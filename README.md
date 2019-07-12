# UBMOD

## Introduction
This repository serves as the source code repository of the Unsaturated zone soil water Balance Model, known as UBMOD. The model a one-dimensional soil water balance model for simulating both downward and upward soil water movement. The model was designed by Wei Mao, Yan Zhu, Tianxing Zhao and Jinzhong yang. The [[Paper](https://www.sciencedirect.com/science/article/pii/S0022169418301562)] can be find here.  

UBMOD can simulate both upward and downward flow in the unsaturated zone. The accommodation of the model to the upward flow, complex boudnary conditions, and soil heterogeneity conditions have been demonstrated. What's more, the model keeps mass balance well, and can give satisfactory results with coarse spatial and temporal discretization. Only four physical meaning soil hydraulic parameters are essential to calculate the soil water movement. All of these make the model more applicable in practice.

## Requirements
Fortran 90

## concept map
![Concmap](https://github.com/Weiwei-Mao/UBMOD/blob/master/images/UBMOD.png)

There are four major components to describe the soil water movement in UBMOD model, as shown in the figure. Firstly, the vertical soil column is divided into a cascade of “buckets” and each “bucket” corresponds to a soil layer. The “buckets” will be filled to saturation from the top layer to the bottom layer if there is infiltration, which is referred as the allocation of infiltration water. Secondly, when the soil water content exceeds the field capacity, the soil water will move downward driven by the gravitational potential. Thirdly, the source/sink terms are used to account for soil evaporation and crop transpiration. Lastly, we calculate the diffusive movement driven by the matric potential.

## Table of contents
index | name | content
-|-|-
A | Source code  | src|
B |  Demo file   | Rh1D.in, Rh1D.out|
C |  Document    | tech.pdf|

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
