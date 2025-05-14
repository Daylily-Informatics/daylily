./plot_all_three_var2.R all_giab_proc.tsv ./atv2.pdf Fscore 0.995 1

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

geom_smooth() using formula = 'y ~ x'
Warning message:
Using size aesthetic for lines was deprecated in ggplot2 3.4.0.
ℹ Please use linewidth instead.
null device
          1
geom_smooth() using formula = 'y ~ x'
null device
          1

Regression summaries per SNPClass and Pipeline:

Group: All | ILMN (sent+sentd)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.48044 -0.20161  0.05221  0.23649  0.29993

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.487657   0.062216   7.838 4.00e-09 ***
Coverage    0.027041   0.005029   5.377 5.56e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2661 on 34 degrees of freedom
Multiple R-squared:  0.4596,	Adjusted R-squared:  0.4437
F-statistic: 28.91 on 1 and 34 DF,  p-value: 5.565e-06


Group: All | ONT (ont+sentdont)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
      Min        1Q    Median        3Q       Max
-0.036098 -0.019248  0.005519  0.014702  0.026439

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.9234625  0.0117658  78.487  < 2e-16 ***
Coverage    0.0027887  0.0006477   4.306 0.000854 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.02158 on 13 degrees of freedom
Multiple R-squared:  0.5878,	Adjusted R-squared:  0.5561
F-statistic: 18.54 on 1 and 13 DF,  p-value: 0.0008543


Group: All | UG (ug+sentdug)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.06893 -0.02346  0.01385  0.02474  0.04427

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.9015179  0.0167610  53.787  < 2e-16 ***
Coverage    0.0027352  0.0007202   3.798  0.00158 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.03637 on 16 degrees of freedom
Multiple R-squared:  0.4741,	Adjusted R-squared:  0.4413
F-statistic: 14.43 on 1 and 16 DF,  p-value: 0.001579


Group: DEL_50 | ILMN (sent+sentd)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.49200 -0.18016  0.05456  0.20810  0.26516

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.524172   0.057478   9.119 1.17e-10 ***
Coverage    0.025045   0.004646   5.390 5.34e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2459 on 34 degrees of freedom
Multiple R-squared:  0.4608,	Adjusted R-squared:  0.4449
F-statistic: 29.06 on 1 and 34 DF,  p-value: 5.345e-06


Group: DEL_50 | ONT (ont+sentdont)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
      Min        1Q    Median        3Q       Max
-0.095764 -0.036649  0.005563  0.014796  0.104111

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.742405   0.035640  20.831 2.27e-11 ***
Coverage    0.007445   0.001962   3.795  0.00223 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06536 on 13 degrees of freedom
Multiple R-squared:  0.5256,	Adjusted R-squared:  0.4891
F-statistic:  14.4 on 1 and 13 DF,  p-value: 0.002229


Group: DEL_50 | UG (ug+sentdug)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.11536 -0.03798  0.00466  0.02235  0.10983

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.762811   0.030895  24.691 3.63e-14 ***
Coverage    0.004889   0.001327   3.683  0.00201 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06703 on 16 degrees of freedom
Multiple R-squared:  0.4588,	Adjusted R-squared:  0.4249
F-statistic: 13.56 on 1 and 16 DF,  p-value: 0.002015


Group: INS_50 | ILMN (sent+sentd)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.49878 -0.17883  0.05868  0.21499  0.27341

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.524325   0.058455   8.970 1.75e-10 ***
Coverage    0.025061   0.004725   5.304 6.93e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.25 on 34 degrees of freedom
Multiple R-squared:  0.4528,	Adjusted R-squared:  0.4367
F-statistic: 28.13 on 1 and 34 DF,  p-value: 6.926e-06


Group: INS_50 | ONT (ont+sentdont)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
      Min        1Q    Median        3Q       Max
-0.094804 -0.037541  0.008369  0.015437  0.102337

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.748480   0.035362  21.166 1.86e-11 ***
Coverage    0.007268   0.001947   3.734   0.0025 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06486 on 13 degrees of freedom
Multiple R-squared:  0.5174,	Adjusted R-squared:  0.4803
F-statistic: 13.94 on 1 and 13 DF,  p-value: 0.002505


Group: INS_50 | UG (ug+sentdug)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
      Min        1Q    Median        3Q       Max
-0.114858 -0.033935  0.006273  0.021074  0.100683

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.778699   0.028423  27.396 7.14e-15 ***
Coverage    0.004872   0.001221   3.989  0.00106 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06167 on 16 degrees of freedom
Multiple R-squared:  0.4986,	Adjusted R-squared:  0.4673
F-statistic: 15.91 on 1 and 16 DF,  p-value: 0.001056


Group: Indel_50 | ILMN (sent+sentd)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.44675 -0.26802  0.03611  0.27350  0.34874

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.445889   0.067660   6.590 1.49e-07 ***
Coverage    0.029416   0.005469   5.378 5.54e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2894 on 34 degrees of freedom
Multiple R-squared:  0.4597,	Adjusted R-squared:  0.4438
F-statistic: 28.93 on 1 and 34 DF,  p-value: 5.538e-06


Group: Indel_50 | ONT (ont+sentdont)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.12040 -0.04597 -0.01098  0.04883  0.13322

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.727866   0.044962  16.188 5.37e-10 ***
Coverage    0.007319   0.002475   2.957   0.0111 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.08246 on 13 degrees of freedom
Multiple R-squared:  0.4022,	Adjusted R-squared:  0.3562
F-statistic: 8.745 on 1 and 13 DF,  p-value: 0.01112


Group: Indel_50 | UG (ug+sentdug)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
      Min        1Q    Median        3Q       Max
-0.224335 -0.066259  0.000416  0.070670  0.204840

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.682420   0.058010  11.764 2.74e-09 ***
Coverage    0.006644   0.002492   2.666   0.0169 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.1259 on 16 degrees of freedom
Multiple R-squared:  0.3075,	Adjusted R-squared:  0.2642
F-statistic: 7.105 on 1 and 16 DF,  p-value: 0.01692


Group: SNPts | ILMN (sent+sentd)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.47892 -0.20455  0.05217  0.23975  0.30360

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.484316   0.062744   7.719 5.62e-09 ***
Coverage    0.027234   0.005072   5.370 5.69e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2684 on 34 degrees of freedom
Multiple R-squared:  0.4589,	Adjusted R-squared:  0.443
F-statistic: 28.83 on 1 and 34 DF,  p-value: 5.688e-06


Group: SNPts | ONT (ont+sentdont)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
      Min        1Q    Median        3Q       Max
-0.027822 -0.014223  0.007645  0.015335  0.016532

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.9441567  0.0096363  97.979  < 2e-16 ***
Coverage    0.0022793  0.0005305   4.297 0.000868 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.01767 on 13 degrees of freedom
Multiple R-squared:  0.5868,	Adjusted R-squared:  0.555
F-statistic: 18.46 on 1 and 13 DF,  p-value: 0.0008683


Group: SNPts | UG (ug+sentdug)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.06152 -0.02182  0.01171  0.02698  0.03538

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.919245   0.015618  58.859  < 2e-16 ***
Coverage    0.002447   0.000671   3.647  0.00217 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.03389 on 16 degrees of freedom
Multiple R-squared:  0.4539,	Adjusted R-squared:  0.4198
F-statistic:  13.3 on 1 and 16 DF,  p-value: 0.002174


Group: SNPtv | ILMN (sent+sentd)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.47623 -0.20572  0.05005  0.23875  0.30426

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.480920   0.062752   7.664 6.57e-09 ***
Coverage    0.027381   0.005072   5.398 5.22e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.2684 on 34 degrees of freedom
Multiple R-squared:  0.4615,	Adjusted R-squared:  0.4457
F-statistic: 29.14 on 1 and 34 DF,  p-value: 5.223e-06


Group: SNPtv | ONT (ont+sentdont)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
      Min        1Q    Median        3Q       Max
-0.030600 -0.014990  0.007897  0.015914  0.017669

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.941334   0.010136  92.872  < 2e-16 ***
Coverage    0.002380   0.000558   4.266  0.00092 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.01859 on 13 degrees of freedom
Multiple R-squared:  0.5833,	Adjusted R-squared:  0.5512
F-statistic:  18.2 on 1 and 13 DF,  p-value: 0.0009198


Group: SNPtv | UG (ug+sentdug)

Call:
lm(formula = as.formula(paste(plot_var, "~ Coverage")), data = .)

Residuals:
     Min       1Q   Median       3Q      Max
-0.06690 -0.02286  0.01255  0.02718  0.03852

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.9097103  0.0163249  55.725  < 2e-16 ***
Coverage    0.0026503  0.0007014   3.778  0.00165 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.03542 on 16 degrees of freedom
Multiple R-squared:  0.4715,	Adjusted R-squared:  0.4385
F-statistic: 14.28 on 1 and 16 DF,  p-value: 0.001646

# A tibble: 0 × 0
# Rowwise:

---------

Your output is correct and meaningful. Here’s how to interpret it clearly:

Interpretation Summary:
Coverage Effect (All Pipelines, All SNP Classes):

ILMN (sent+sentd): Highly significant positive effect of Coverage (p=5.56e-06).

ONT (ont+sentdont): Also highly significant positive effect (p=0.000854), though smaller in magnitude.

UG (ug+sentdug): Significant positive effect (p=0.00158), similar magnitude to ONT.

Overall Interpretation:
Increasing coverage consistently improves Fscore significantly across all pipelines, but the magnitude of improvement per coverage increment is largest in ILMN.

Pipeline-specific Interpretations by SNP Class:
DEL_50 and INS_50 SNP Classes:

All pipelines show significant positive effects of increased coverage, strongly suggesting these SNP classes are coverage-sensitive across all technologies.

Indel_50 SNP Class:

ILMN shows very strong significant effects (p=5.54e-06), indicating coverage sensitivity.

ONT (p=0.0111) and UG (p=0.0169) effects are smaller, though still significant, suggesting lower sensitivity compared to ILMN.

SNPts and SNPtv Classes:

All pipelines again show strong, significant coverage effects, reinforcing the general conclusion that higher coverage benefits accuracy significantly for these SNP types.

Practical Recommendations:
Prioritize increasing coverage for maximum accuracy gain, particularly if using ILMN data.

For classes like DEL_50, INS_50, and SNP (ts, tv), accuracy improvements from higher coverage are especially consistent and meaningful.

ONT and UG pipelines benefit from higher coverage as well, though the improvement per coverage increment is smaller compared to ILMN.

Possible Further Analysis:
Consider interactions between Coverage and SNP classes or Pipelines to reveal nuanced improvements.

Run more detailed models or non-linear regressions to identify thresholds beyond which accuracy improvements diminish.

In conclusion, your analyses successfully identified clear, statistically relevant effects of coverage increments on accuracy across different SNP classes and pipelines.










======



# A tibble: 0 × 0
# Rowwise:
   Coverage ILMN_Lower ILMN_Upper  UG_Lower  UG_Upper Equivalent
1         5  0.5878162  0.6632398 0.8183403 0.8735216      FALSE
2         6  0.6160144  0.6887676 0.8234933 0.8764475      FALSE
3         7  0.6438211  0.7146869 0.8286134 0.8794063      FALSE
4         8  0.6712044  0.7410294 0.8336961 0.8824025      FALSE
5         9  0.6981456  0.7678142 0.8387365 0.8854409      FALSE
6        10  0.7246415  0.7950442 0.8437289 0.8885274      FALSE
7        11  0.7507059  0.8227057 0.8486668 0.8916683      FALSE
8        12  0.7763665  0.8507711 0.8535433 0.8948707      FALSE
9        13  0.8016609  0.8792027 0.8583506 0.8981423       TRUE
10       14  0.8266313  0.9079582 0.8630802 0.9014915       TRUE
11       15  0.8513208  0.9369946 0.8677237 0.9049269       TRUE
12       16  0.8757698  0.9662716 0.8722724 0.9084572       TRUE
13       17  0.9000146  0.9957527 0.8767180 0.9120904       TRUE
14       18  0.9240870  1.0254063 0.8810534 0.9158339      FALSE
15       19  0.9480138  1.0552054 0.8852728 0.9196933      FALSE
16       20  0.9718176  1.0851275 0.8893727 0.9236723      FALSE
17       21  0.9955174  1.1151537 0.8933517 0.9277722      FALSE
18       22  1.0191288  1.1452682 0.8972111 0.9319917      FALSE
19       23  1.0426647  1.1754582 0.9009546 0.9363270      FALSE
20       24  1.0661360  1.2057129 0.9045879 0.9407727      FALSE
21       25  1.0895517  1.2360232 0.9081181 0.9453213      FALSE
22       26  1.1129192  1.2663816 0.9115535 0.9499648      FALSE
23       27  1.1362448  1.2967820 0.9149027 0.9546945      FALSE
24       28  1.1595338  1.3272189 0.9181743 0.9595017      FALSE
25       29  1.1827908  1.3576878 0.9213767 0.9643782      FALSE
26       30  1.2060195  1.3881851 0.9245176 0.9693162      FALSE



Corrected Interpretation of the Table:
Coverage < 13x:
UG (Ultima) and ILMN (Illumina) confidence intervals do not overlap, and because UG’s intervals are entirely above ILMN’s intervals, UG is statistically outperforming ILMN at these low coverage points.

Coverage from 13x–17x:
Confidence intervals overlap (Equivalent = TRUE), meaning at these points UG and ILMN are statistically indistinguishable (no significant difference).

Coverage ≥ 18x:
Confidence intervals no longer overlap, and ILMN’s confidence interval is now entirely above UG’s intervals. This means ILMN begins to significantly outperform UG at coverage ≥ 18x.

Simple Summary of Findings:
Coverage Level	Which Technology Performs Best?
< 13x	✅ Ultima (UG)
13x–17x	⚖️ Statistically Equivalent
≥ 18x (to 30x)	✅ Illumina (ILMN)