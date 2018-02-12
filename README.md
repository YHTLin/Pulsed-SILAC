# Pulsed-SILAC
Analysis of pulsed-SILAC proteomics data on the effect of mTOR inhibitors in GTML5 mouse cells.

## Motivation
Signaling pathways are responsible for coordinating cell behavior, and when dysregulated they often lead to
disease. The PIK3CA-AKT-mTOR pathway is commonly activated in human cancers. It plays an important
role in cell cycle regulation, metabolism, and protein synthesis. One of its members, the mechanistic target
of rapamycin (mTOR), is an serine/threonine kinase whose over-activation promotes tumor progression.
To reverse the effect of increased mTOR activity, mTOR inhibitors, such as rapamycin (first generation)
and mTOR kinase inhibitors (second generation), are prescribed as cancer therapy. However, over time,
natural selection of tumors results in drug resistance and clinical relapse. To counteract mTOR mutants,
Rodrik-Outmezguine et al. recently developed a third-generation mTOR inhibitor called Rapalink-1 (M1071).
M1071 has been established by Fan et al. as a potent drug for brain tumors owing to its ability to cross the
blood-brain barrier. While the mechanism of action of M1071 is known by design, its downstream effects on
biological processes, such as protein translation, are not well understood. The goal of this study is to identify
proteins whose biosyntheses are sensitive to M1071 in the context of neural cells. We aim to monitor changes
in protein translation using two techniques, ribosome profiling and pulsed SILAC proteomics.

## File Description
+ *pulsed_SILAC_ozlem.R*: R analysis script
+ *pulsed_SILAC.RData*: Saved workspace after running the script above
+ *deployDirectory/app.R*: Source code behind interactive app.

## Usage
An interactive web tool, built using ShinyApp in R, is available for dynamic data visualization. Visit the site [here](https://tony-lin.shinyapps.io/deploydirectory/).
