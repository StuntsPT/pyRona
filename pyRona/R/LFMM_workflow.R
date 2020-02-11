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
if (!require("lfmm")) install.packages("lfmm")
library("lfmm")

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

# Calibration type to use for LFMM tests. Valid values are "gif" or "median+MAD"
CALIBRATION = "gif"

# Estimate type to use - valid values are "ridge" or "lasso"
ESTIMATE = "ridge"

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


# LFMM estimates
lfmm_estimates = function(env_file,
                          geno_data,
                          PCA_points,
                          calibration,
                          estimate) {

    env_vars = read.csv(env_file, sep="\t", header=FALSE)[,-1]
    env_vars = scale(env_vars)

    # Fit an LFMM, i.e, compute B, U, V estimates
    if (estimate == "ridge") {
        mod.lfmm <- lfmm_ridge(Y = geno_data,
                               X = env_vars,
                               K = PCA_points)
    } else if (estimate == "lasso") {
        mod.lfmm <- lfmm_lasso(Y = geno_data,
                               X = env_vars,
                               K = PCA_points,
                               nozero.prop = 0.01)
    } else {
        stop("This script only supports 'ridge' or 'lasso' estimates.")
    }

    # Performs association testing using the fitted model:
    pv <- lfmm_test(Y = geno_data,
                    X = env_vars,
                    lfmm = mod.lfmm,
                    calibrate = calibration)

    pvalues <- pv$calibrated.pvalue

    return(pvalues)
}


# Draw a QQ-plot
draw_qq_plot = function(pvalues) {
    qqplot(rexp(length(pvalues), rate=log(10)),
           -log10(pvalues), xlab="Expected quantile",
           pch = 19, cex = .4)
    abline(0,1)
}


# Write down the assocaitions table
write_associations_table = function(assoc_table, pvalues){
    write.table(x=pvalues,
                col.names=FALSE,
                row.names=F,
                sep=",",
                file=assoc_table)
}

## Function invocation
convert_file(original_vcf, target_lfmm)
genetic_data = read.lfmm(target_lfmm)

PCA_points = preliminary_pca(genetic_data, PCA_points)

pvalues = lfmm_estimates(ENV_FILE,
                         genetic_data,
                         PCA_points,
                         CALIBRATION,
                         ESTIMATE)

draw_qq_plot(pvalues)

write_associations_table(ASSOCIATION_TABLE, pvalues)
