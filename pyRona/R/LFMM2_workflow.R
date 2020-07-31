#!/usr/bin/Rscript
# Copyright 2018-2020 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of pyRona.
# pyRona is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pyRona is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pyRona.  If not, see <http://www.gnu.org/licenses/>.

## Load required libraries
if (!require("LEA")) {install.packages("BiocManager");
                      BiocManager::install("LEA", update=F)}
library("LEA")

if (!require("nFactors")) install.packages("nFactors")
library("nFactors")

## Set variables:
# Path to vcf file to convert to lfmm eg. "~/qsuber_GEO.vcf"
original_vcf = ""

# Path to the converted lfmm file eg. "~/qsuber_GEO.lfmm"
target_lfmm = ""

# Number of "K" to use. You should run the PCA first and then define this number
# eg. 5
PCA_points = 5

# Path to file with environmental variables. The first column should be the name
# of the "population" each sample belongs to. Eg. "~/env_names.txt"
ENV_FILE = ""

# Path to file where the association table will be written to.
# Each column represents a covariate and each line represents a marker.
# Eg. "~/associations.csv"
ASSOCIATION_TABLE = ""

## Functions

# Converts a VCF file to LFMM format
convert_file = function(original_vcf, target_lfmm) {
    vcf2lfmm(original_vcf, target_lfmm)
}


# Performs a preliminary PCA analysis. Use it to determine the best "K"
preliminary_pca = function(genetic_data, PCA_points) {

    pc <- prcomp(genetic_data)
    if (PCA_points == "estimate") {
        PCA_points = as.numeric(as.character(nFactors::nCng(pc$sdev^2, cor=F, details=F)$nFactors))
    }
    plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
    points(PCA_points, pc$sdev[PCA_points]^2, type="h", lwd=3, col="blue")

    return(PCA_points)
}

data_impute = function(genetic_data, K, repetitions=10, method="mode") {
    project.snmf = snmf(genetic_data, K=K, entropy=T, repetitions=repetitions, project="new")
    best = which.min(cross.entropy(project.snmf, K=K))
    impute(project.snmf, genetic_data, method = method, K=K, run = best)

}

# LFMM estimates
lfmm_estimates = function(env_file,
                          geno_data,
                          PCA_points) {

    env_file = env_file_compat(env_file)

    mod = lfmm2(input=genetic_data,
                env=env_file,
                K=PCA_points)
    
    pv <- lfmm2.test(object=mod,
                     input=genetic_data,
                     env=env_file, linear = TRUE)
    
    pvalues <- t(pv$pvalues)

    return(pvalues)
}


# Checks of the env file is made for lfmm1.5 or lfmm2
# If a file's extension does not end with ".env"
# If it is meant for lfmm1.5 make arrangements for it to work with lfmm2
env_file_compat = function(env_file) {
    if (substring(env_file, nchar(env_file) - 3) != ".env") {
        env_vars = read.csv(env_file, sep="\t", header=FALSE)[,-1]
        new_filename = paste(substr(env_file, 0, nchar(env_file) - 4), ".env", sep="")
        write.table(file=new_filename,
                    x=env_vars,
                    sep="\t",
                    col.names=F,
                    row.names=F)
        env_file = new_filename
    }
    
    return(env_file)
}


# Write down the associations table
write_associations_table = function(assoc_table, pvalues){
    write.table(x=pvalues,
                col.names=FALSE,
                row.names=F,
                sep=",",
                file=assoc_table)
}

## Function invocation
#convert_file(original_vcf, target_lfmm)
genetic_data = read.lfmm(target_lfmm)

PCA_points = preliminary_pca(genetic_data, PCA_points)
data_impute(genetic_data=target_lfmm, K=PCA_points)

# Change the genetic data to the imputed data file
genetic_data = read.lfmm(paste(target_lfmm ,"_imputed.lfmm", sep=""))

pvalues = lfmm_estimates(ENV_FILE,
                         genetic_data,
                         PCA_points)

write_associations_table(ASSOCIATION_TABLE, pvalues)
