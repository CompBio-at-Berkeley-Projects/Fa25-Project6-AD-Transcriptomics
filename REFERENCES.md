# References - Alzheimer's Disease Transcriptomics Project

## Primary Dataset Source
**Used for**: Core single-nucleus RNA-seq data (GSE157827) containing 169,496 nuclei from 21 donors (12 AD, 9 control) from human prefrontal cortex. This dataset forms the foundation of all gene expression analyses in this project.

Lau, S.-F., Cao, H., Fu, A. K. Y., & Ip, N. Y. (2020). Single-nucleus transcriptome analysis reveals dysregulation of angiogenic endothelial cells and neuroprotective glia in Alzheimer's disease. *Proceedings of the National Academy of Sciences of the United States of America*, *117*(41), 25800-25809. https://doi.org/10.1073/pnas.2008762117

---

## Database Resources
**Used for**: Access to curated, analysis-ready single-cell and spatial RNA-seq datasets for Alzheimer's disease through the ssREAD portal. Provided standardized metadata and quality-controlled count matrices.

Wang, C., Xiang, Y., Wang, F., Peng, Y., Zhou, Y., Zhang, X., Cheng, X., Yang, B., Qiu, X., Chen, Y., Li, Z., Yu, H., Zhang, Z., Sun, L., Tang, R., Liu, Y., Li, J., Zhao, D., Cheng, X., ... Gu, H. (2024). A single-cell and spatial RNA-seq database for Alzheimer's disease (ssREAD). *Nature Communications*, *15*, Article 4335. https://doi.org/10.1038/s41467-024-48207-6

---

## GWAS and Genetic Risk Studies
**Used for**: Validation of genetic risk factors and identification of GWAS-significant genes (APOE, APP, GRN, ACE, APH1B) analyzed in this project. Provided genetic evidence linking target genes to AD pathology.

Bellenguez, C., Küçükali, F., Jansen, I. E., Kleineidam, L., Moreno-Grau, S., Amin, N., Naj, A. C., Campos-Martin, R., Grenier-Boley, B., Andrade, V., Holmans, P. A., Boland, A., Damotte, V., van der Lee, S. J., Costa, M. R., Kuulasmaa, T., Yang, Q., de Rojas, I., Bis, J. C., ... Lambert, J.-C. (2022). New insights into the genetic etiology of Alzheimer's disease and related dementias. *Nature Genetics*, *54*(4), 412-436. https://doi.org/10.1038/s41588-022-01024-z

Kunkle, B. W., Grenier-Boley, B., Sims, R., Bis, J. C., Damotte, V., Naj, A. C., Boland, A., Vronskaya, M., van der Lee, S. J., Amlie-Wolf, A., Bellenguez, C., Frizatti, A., Chouraki, V., Martin, E. R., Sleegers, K., Badarinarayan, N., Jakobsdottir, J., Hamilton-Nelson, K. L., Moreno-Grau, S., ... Pericak-Vance, M. A. (2019). Genetic meta-analysis of diagnosed Alzheimer's disease identifies new risk loci and implicates Aβ, tau, immunity and lipid processing. *Nature Genetics*, *51*(3), 414-430. https://doi.org/10.1038/s41588-019-0358-2

---

## Computational Methods & Cell-Type Annotation
**Used for**: SCINA (Semi-supervised Category Identification and Assignment) algorithm for marker-based cell-type annotation. This method assigned biological identities (astrocytes, microglia, neurons, etc.) to transcriptomic clusters.

Zhang, Z., Luo, D., Zhong, X., Choi, J. H., Ma, Y., Wang, S., Mahrt, E., Guo, W., Stawiski, E. W., Modrusan, Z., Seshagiri, S., Kapur, P., Hon, G. C., Brugarolas, J., & Wang, T. (2019). SCINA: A semi-supervised subtyping algorithm of single cells and bulk samples. *Genes*, *10*(7), 531. https://doi.org/10.3390/genes10070531

---

## Single-Cell RNA-seq Analysis Software
**Used for**: Seurat v5 R package - the core analytical pipeline for quality control, normalization (LogNormalize), dimensionality reduction (PCA, UMAP), clustering, differential expression analysis, and visualization of single-nucleus RNA-seq data.

Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, *42*(2), 293-304. https://doi.org/10.1038/s41587-023-01767-y

Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W. M., Zheng, S., Butler, A., Lee, M. J., Wilk, A. J., Darby, C., Zager, M., Hoffman, P., Stoeckius, M., Papalexi, E., Mimitou, E. P., Jain, J., Srivastava, A., Stuart, T., Fleming, L. M., Yeung, B., ... Satija, R. (2021). Integrated analysis of multimodal single-cell data. *Cell*, *184*(13), 3573-3587.e29. https://doi.org/10.1016/j.cell.2021.04.048

Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck, W. M., Hao, Y., Stoeckius, M., Smibert, P., & Satija, R. (2019). Comprehensive integration of single-cell data. *Cell*, *177*(7), 1888-1902.e21. https://doi.org/10.1016/j.cell.2019.05.031

---

## Reference-Based Cell-Type Annotation
**Used for**: Azimuth reference-based annotation tool for validation of cell-type assignments and analysis of ACE expression in GSE175814 dataset.

Hao, Y., et al. (2024). *Nature Biotechnology*. (See above)

---

## APOE Biology & Lipid Metabolism
**Used for**: Understanding APOE's role in cholesterol transport, lipid homeostasis, Aβ clearance, and astrocyte-mediated neuroprotection. APOE-ε4 is the strongest genetic risk factor for late-onset AD.

Liu, C.-C., Kanekiyo, T., Xu, H., & Bu, G. (2013). Apolipoprotein E and Alzheimer disease: Risk, mechanisms and therapy. *Nature Reviews Neurology*, *9*(2), 106-118. https://doi.org/10.1038/nrneurol.2012.263

Yamazaki, Y., Zhao, N., Caulfield, T. R., Liu, C.-C., & Bu, G. (2019). Apolipoprotein E and Alzheimer disease: Pathobiology and targeting strategies. *Nature Reviews Neurology*, *15*(9), 501-518. https://doi.org/10.1038/s41582-019-0228-7

Shi, Y., Yamada, K., Liddelow, S. A., Smith, S. T., Zhao, L., Luo, W., Tsai, R. M., Spina, S., Grinberg, L. T., Rojas, J. C., Gallardo, G., Wang, K., Roh, J., Robinson, G., Finn, M. B., Jiang, H., Sullivan, P. M., Baufeld, C., Wood, M. W., ... Holtzman, D. M. (2017). ApoE4 markedly exacerbates tau-mediated neurodegeneration in a mouse model of tauopathy. *Nature*, *549*(7673), 523-527. https://doi.org/10.1038/nature24016

---

## APP Biology & Amyloid Processing
**Used for**: Understanding APP's role as amyloid precursor protein, its processing by secretases (α, β, γ), and generation of Aβ peptides (Aβ40, Aβ42) central to amyloid cascade hypothesis.

O'Brien, R. J., & Wong, P. C. (2011). Amyloid precursor protein processing and Alzheimer's disease. *Annual Review of Neuroscience*, *34*, 185-204. https://doi.org/10.1146/annurev-neuro-061010-113613

Haass, C., Kaether, C., Thinakaran, G., & Sisodia, S. (2012). Trafficking and proteolytic processing of APP. *Cold Spring Harbor Perspectives in Medicine*, *2*(5), a006270. https://doi.org/10.1101/cshperspect.a006270

Zhang, Y.-W., Thompson, R., Zhang, H., & Xu, H. (2011). APP processing in Alzheimer's disease. *Molecular Brain*, *4*, 3. https://doi.org/10.1186/1756-6606-4-3

---

## GRN Biology & Microglial Function
**Used for**: Understanding GRN (Progranulin) as a key regulator of microglial homeostasis, lysosomal function, and neuroinflammation. GRN haploinsufficiency causes frontotemporal dementia (FTD) and is implicated in AD susceptibility.

Kao, A. W., McKay, A., Singh, P. P., Brunet, A., & Huang, E. J. (2017). Progranulin, lysosomal regulation and neurodegenerative disease. *Nature Reviews Neuroscience*, *18*(6), 325-333. https://doi.org/10.1038/nrn.2017.36

Zhou, X., Sun, L., Bastos de Oliveira, F., Qi, X., Brown, W. J., Smolka, M. B., Sun, Y., & Hu, F. (2015). Prosaposin facilitates sortilin-independent lysosomal trafficking of progranulin. *Journal of Cell Biology*, *210*(6), 991-1002. https://doi.org/10.1083/jcb.201502029

Minami, S. S., Min, S.-W., Krabbe, G., Wang, C., Zhou, Y., Asgarov, R., Li, Y., Martens, L. H., Elia, L. P., Ward, M. E., Mucke, L., Farese, R. V., & Gan, L. (2014). Progranulin protects against amyloid β deposition and toxicity in Alzheimer's disease mouse models. *Nature Medicine*, *20*(10), 1157-1164. https://doi.org/10.1038/nm.3672

---

## ACE Biology & Aβ Degradation
**Used for**: Understanding ACE (Angiotensin-Converting Enzyme) as an Aβ-degrading enzyme with vascular expression. ACE activity influences Aβ clearance and cerebrovascular health in AD.

Hemming, M. L., & Selkoe, D. J. (2005). Amyloid β-protein is degraded by cellular angiotensin-converting enzyme (ACE) and elevated by an ACE inhibitor. *Journal of Biological Chemistry*, *280*(45), 37644-37650. https://doi.org/10.1074/jbc.M508460200

Miners, J. S., Ashby, E., Van Helmond, Z., Chalmers, K. A., Palmer, L. E., Love, S., & Kehoe, P. G. (2008). Angiotensin-converting enzyme (ACE) levels and activity in Alzheimer's disease, and relationship of perivascular ACE-1 to cerebral amyloid angiopathy. *Neuropathology and Applied Neurobiology*, *34*(2), 181-193. https://doi.org/10.1111/j.1365-2990.2007.00885.x

Kehoe, P. G., Miners, S., & Love, S. (2009). Angiotensins in Alzheimer's disease – friend or foe? *Trends in Neurosciences*, *32*(12), 619-628. https://doi.org/10.1016/j.tins.2009.07.006

---

## APH1B Biology & γ-Secretase Complex
**Used for**: Understanding APH1B as a core component of the γ-secretase complex responsible for intramembranous cleavage of APP to generate Aβ peptides. APH1B variants influence Aβ42/Aβ40 ratio.

De Strooper, B., Iwatsubo, T., & Wolfe, M. S. (2012). Presenilins and γ-secretase: Structure, function, and role in Alzheimer disease. *Cold Spring Harbor Perspectives in Medicine*, *2*(1), a006304. https://doi.org/10.1101/cshperspect.a006304

Serneels, L., Van Biervliet, J., Craessaerts, K., Dejaegere, T., Horré, K., Van Houtvin, T., Esselmann, H., Paul, S., Schäfer, M. K., Berezovska, O., Hyman, B. T., Sprangers, B., Sciot, R., Moons, L., Jucker, M., Yang, Z., May, P. C., Karran, E., Wiltfang, J., ... De Strooper, B. (2009). γ-Secretase heterogeneity in the Aph1 subunit: Relevance for Alzheimer's disease. *Science*, *324*(5927), 639-642. https://doi.org/10.1126/science.1171176

Bai, X.-C., Yan, C., Yang, G., Lu, P., Ma, D., Sun, L., Zhou, R., Scheres, S. H. W., & Shi, Y. (2015). An atomic structure of human γ-secretase. *Nature*, *525*(7568), 212-217. https://doi.org/10.1038/nature14892

---

## Astrocyte Dysfunction in AD
**Used for**: Understanding reactive astrocyte states (A1 vs A2), loss of homeostatic functions, and astrocyte contributions to neuroinflammation and synaptic dysfunction in AD.

Liddelow, S. A., Guttenplan, K. A., Clarke, L. E., Bennett, F. C., Bohlen, C. J., Schirmer, L., Bennett, M. L., Münch, A. E., Chung, W.-S., Peterson, T. C., Wilton, D. K., Frouin, A., Napier, B. A., Panicker, N., Kumar, M., Buckwalter, M. S., Rowitch, D. H., Dawson, V. L., Dawson, T. M., ... Barres, B. A. (2017). Neurotoxic reactive astrocytes are induced by activated microglia. *Nature*, *541*(7638), 481-487. https://doi.org/10.1038/nature21029

Escartin, C., Galea, E., Lakatos, A., O'Callaghan, J. P., Petzold, G. C., Serrano-Pozo, A., Steinhäuser, C., Volterra, A., Carmignoto, G., Agarwal, A., Allen, N. J., Araque, A., Barbeito, L., Barzilai, A., Bergles, D. E., Bonvento, G., Butt, A. M., Chen, W.-T., Cohen-Salmon, M., ... Verkhratsky, A. (2021). Reactive astrocyte nomenclature, definitions, and future directions. *Nature Neuroscience*, *24*(3), 312-325. https://doi.org/10.1038/s41593-020-00783-4

Habib, N., McCabe, C., Medina, S., Varshavsky, M., Kitsberg, D., Dvir-Szternfeld, R., Green, G., Dionne, D., Nguyen, L., Marshall, J. L., Chen, F., Zhang, F., Kaplan, T., Regev, A., & Schwartz, M. (2020). Disease-associated astrocytes in Alzheimer's disease and aging. *Nature Neuroscience*, *23*(6), 701-706. https://doi.org/10.1038/s41593-020-0624-8

---

## Microglial Dysfunction in AD
**Used for**: Understanding disease-associated microglia (DAM) phenotypes, loss of homeostatic microglial function, chronic neuroinflammation, and impaired phagocytosis in AD.

Keren-Shaul, H., Spinrad, A., Weiner, A., Matcovitch-Natan, O., Dvir-Szternfeld, R., Ulland, T. K., David, E., Baruch, K., Lara-Astaiso, D., Toth, B., Itzkovitz, S., Colonna, M., Schwartz, M., & Amit, I. (2017). A unique microglia type associated with restricting development of Alzheimer's disease. *Cell*, *169*(7), 1276-1290.e17. https://doi.org/10.1016/j.cell.2017.05.018

Krasemann, S., Madore, C., Cialic, R., Baufeld, C., Calcagno, N., El Fatimy, R., Beckers, L., O'Loughlin, E., Xu, Y., Fanek, Z., Greco, D. J., Smith, S. T., Tweet, G., Humulock, Z., Zrzavy, T., Conde-Sanroman, P., Gacias, M., Weng, Z., Chen, H., ... Butovsky, O. (2017). The TREM2-APOE pathway drives the transcriptional phenotype of dysfunctional microglia in neurodegenerative diseases. *Immunity*, *47*(3), 566-581.e9. https://doi.org/10.1016/j.immuni.2017.08.008

---

## Neuroinflammation in AD
**Used for**: Understanding chronic inflammation, cytokine signaling (IL-1β, TNF-α, IL-6), complement activation, and inflammatory cascades in AD pathogenesis.

Heneka, M. T., Carson, M. J., El Khoury, J., Landreth, G. E., Brosseron, F., Feinstein, D. L., Jacobs, A. H., Wyss-Coray, T., Vitorica, J., Ransohoff, R. M., Herrup, K., Frautschy, S. A., Finsen, B., Brown, G. C., Verkhratsky, A., Yamanaka, K., Koistinaho, J., Latz, E., Halle, A., ... Kummer, M. P. (2015). Neuroinflammation in Alzheimer's disease. *The Lancet Neurology*, *14*(4), 388-405. https://doi.org/10.1016/S1474-4422(15)70016-5

Kinney, J. W., Bemiller, S. M., Murtishaw, A. S., Leisgang, A. M., Salazar, A. M., & Lamb, B. T. (2018). Inflammation as a central mechanism in Alzheimer's disease. *Alzheimer's & Dementia: Translational Research & Clinical Interventions*, *4*, 575-590. https://doi.org/10.1016/j.trci.2018.06.014

---

## Oligodendrocyte Biology in AD
**Used for**: Understanding myelin dysfunction, oligodendrocyte vulnerability, and white matter degeneration in AD. APP is highly expressed in oligodendrocytes and may contribute to oligodendrocyte-mediated pathology.

Quintela-López, T., Ortiz-Sanz, C., Serrano-Regal, M. P., Gaminde-Blasco, A., Valero, J., Baleriola, J., Sánchez-Gómez, M. V., Matute, C., & Alberdi, E. (2019). Aβ oligomers promote oligodendrocyte differentiation and maturation via integrin β1 and Fyn kinase signaling. *Cell Death & Disease*, *10*(6), 445. https://doi.org/10.1038/s41419-019-1636-8

Vanzago, G., Thion, M. S., Castoldi, F., & Maggi, R. (2023). Oligodendrocyte progenitor cells in Alzheimer's disease: From physiology to pathology. *Frontiers in Cellular Neuroscience*, *17*, 1125136. https://doi.org/10.3389/fncel.2023.1125136

---

## Differential Expression Analysis Methods
**Used for**: Statistical methods for identifying differentially expressed genes between AD and control samples in single-cell data. MAST and Wilcoxon rank-sum tests were potential methods for DEG analysis.

Finak, G., McDavid, A., Yajima, M., Deng, J., Gersuk, V., Shalek, A. K., Slichter, C. K., Miller, H. W., McElrath, M. J., Prlic, M., Linsley, P. S., & Gottardo, R. (2015). MAST: A flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. *Genome Biology*, *16*, 278. https://doi.org/10.1186/s13059-015-0844-5

Soneson, C., & Robinson, M. D. (2018). Bias, robustness and scalability in single-cell differential expression analysis. *Nature Methods*, *15*(4), 255-261. https://doi.org/10.1038/nmeth.4612

---

## Data Integration & Batch Correction
**Used for**: Methods for integrating multiple single-cell datasets and removing batch effects. Harmony and Seurat CCA integration are standard approaches for combining AD and control samples from different sequencing batches.

Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., Baglaenko, Y., Brenner, M., Loh, P.-R., & Raychaudhuri, S. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*, *16*(12), 1289-1296. https://doi.org/10.1038/s41592-019-0619-0

---

## UMAP Dimensionality Reduction
**Used for**: UMAP (Uniform Manifold Approximation and Projection) algorithm for 2D visualization of high-dimensional single-cell transcriptomic data. UMAP projections were used to visualize cell clusters and gene expression patterns.

McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform manifold approximation and projection for dimension reduction. *arXiv preprint arXiv:1802.03426*. https://arxiv.org/abs/1802.03426

Becht, E., McInnes, L., Healy, J., Dutertre, C.-A., Kwok, I. W. H., Ng, L. G., Ginhoux, F., & Newell, E. W. (2019). Dimensionality reduction for visualizing single-cell data using UMAP. *Nature Biotechnology*, *37*, 38-44. https://doi.org/10.1038/nbt.4314

---

## Quality Control in Single-Cell RNA-seq
**Used for**: Standards and best practices for quality control filtering (mitochondrial percentage, gene counts, UMI counts) to remove low-quality cells and doublets from single-nucleus RNA-seq data.

Luecken, M. D., & Theis, F. J. (2019). Current best practices in single-cell RNA-seq analysis: A tutorial. *Molecular Systems Biology*, *15*(6), e8746. https://doi.org/10.15252/msb.20188746

Ilicic, T., Kim, J. K., Kolodziejczyk, A. A., Bagger, F. O., McCarthy, D. J., Marioni, J. C., & Teichmann, S. A. (2016). Classification of low quality cells from single-cell RNA-seq data. *Genome Biology*, *17*, 29. https://doi.org/10.1186/s13059-016-0888-1

---

## Normalization Methods
**Used for**: LogNormalize method in Seurat for normalizing raw UMI counts to account for sequencing depth differences across cells.

Hafemeister, C., & Satija, R. (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology*, *20*, 296. https://doi.org/10.1186/s13059-019-1874-1

---

## Gene Ontology & Pathway Enrichment
**Used for**: Methods for interpreting gene lists through pathway enrichment analysis (GSEA, Enrichr) to identify biological processes dysregulated in AD.

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S., & Mesirov, J. P. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences of the United States of America*, *102*(43), 15545-15550. https://doi.org/10.1073/pnas.0506580102

Chen, E. Y., Tan, C. M., Kou, Y., Duan, Q., Wang, Z., Meirelles, G. V., Clark, N. R., & Ma'ayan, A. (2013). Enrichr: Interactive and collaborative HTML5 gene list enrichment analysis tool. *BMC Bioinformatics*, *14*, 128. https://doi.org/10.1186/1471-2105-14-128

---

## Spatial Transcriptomics Technology
**Used for**: Understanding spatial transcriptomics technologies (Visium, MERFISH) that preserve spatial context of gene expression in tissue sections, enabling mapping of gene expression to neuroanatomical regions.

Rao, A., Barkley, D., França, G. S., & Yanai, I. (2021). Exploring tissue architecture using spatial transcriptomics. *Nature*, *596*(7871), 211-220. https://doi.org/10.1038/s41586-021-03634-9

Chen, W.-T., Lu, A., Craessaerts, K., Pavie, B., Sala Frigerio, C., Corthout, N., Qian, X., Laláková, J., Kühnemund, M., Voytyuk, I., Wolfs, L., Mancuso, R., Salta, E., Balusu, S., Snellinx, A., Munck, S., Jurek, A., Fernandez Navarro, J., Saido, T. C., ... De Strooper, B. (2020). Spatial transcriptomics and in situ sequencing to study Alzheimer's disease. *Cell*, *182*(4), 976-991.e19. https://doi.org/10.1016/j.cell.2020.06.038

---

## Single-Nucleus vs Single-Cell RNA-seq
**Used for**: Methodological justification for using single-nucleus RNA-seq (snRNA-seq) instead of single-cell RNA-seq for frozen postmortem human brain tissue. SnRNA-seq better preserves cell-type diversity and reduces dissociation artifacts.

Bakken, T. E., Hodge, R. D., Miller, J. A., Yao, Z., Nguyen, T. N., Aevermann, B., Barkan, E., Bertagnolli, D., Casper, T., Dee, N., Garren, E., Goldy, J., Graybuck, L. T., Kroll, M., Lasken, R. S., Lathia, K., Parry, S., Rimorin, C., Scheuermann, R. H., ... Lein, E. S. (2018). Single-nucleus and single-cell transcriptomes compared in matched cortical cell types. *PLOS ONE*, *13*(12), e0209648. https://doi.org/10.1371/journal.pone.0209648

Thrupp, N., Sala Frigerio, C., Wolfs, L., Skene, N. G., Fattorelli, N., Poovathingal, S., Fourne, Y., Matthews, P. M., Theys, T., Mancuso, R., de Strooper, B., & Fiers, M. (2020). Single-nucleus RNA-seq is not suitable for detection of microglial activation genes in humans. *Cell Reports*, *32*(13), 108189. https://doi.org/10.1016/j.celrep.2020.108189

---

## 10x Genomics Technology
**Used for**: Description of 10x Genomics Chromium platform for high-throughput droplet-based single-nucleus RNA-seq used to generate the GSE157827 dataset.

Zheng, G. X. Y., Terry, J. M., Belgrader, P., Ryvkin, P., Bent, Z. W., Wilson, R., Ziraldo, S. B., Wheeler, T. D., McDermott, G. P., Zhu, J., Gregory, M. T., Shuga, J., Montesclaros, L., Underwood, J. G., Masquelier, D. A., Nishimura, S. Y., Schnall-Levin, M., Wyatt, P. W., Hindson, C. M., ... Bielas, J. H. (2017). Massively parallel digital transcriptional profiling of single cells. *Nature Communications*, *8*, 14049. https://doi.org/10.1038/ncomms14049

---

## Cerebrovascular Dysfunction in AD
**Used for**: Understanding blood-brain barrier (BBB) breakdown, endothelial cell dysfunction, pericyte loss, and cerebral amyloid angiopathy (CAA) in AD. ACE expression in vascular cells is relevant to cerebrovascular pathology.

Sweeney, M. D., Sagare, A. P., & Zlokovic, B. V. (2018). Blood-brain barrier breakdown in Alzheimer disease and other neurodegenerative disorders. *Nature Reviews Neurology*, *14*(3), 133-150. https://doi.org/10.1038/nrneurol.2017.188

Nelson, A. R., Sweeney, M. D., Sagare, A. P., & Zlokovic, B. V. (2016). Neurovascular dysfunction and neurodegeneration in dementia and Alzheimer's disease. *Biochimica et Biophysica Acta*, *1862*(5), 887-900. https://doi.org/10.1016/j.bbadis.2015.12.016

---

## Amyloid Cascade Hypothesis
**Used for**: Theoretical framework explaining AD pathogenesis through sequential steps: APP processing → Aβ production → oligomer formation → plaque deposition → tau pathology → neurodegeneration.

Hardy, J., & Selkoe, D. J. (2002). The amyloid hypothesis of Alzheimer's disease: Progress and problems on the road to therapeutics. *Science*, *297*(5580), 353-356. https://doi.org/10.1126/science.1072994

Karran, E., Mercken, M., & De Strooper, B. (2011). The amyloid cascade hypothesis for Alzheimer's disease: An appraisal for the development of therapeutics. *Nature Reviews Drug Discovery*, *10*(9), 698-712. https://doi.org/10.1038/nrd3505

---

## Tau Pathology
**Used for**: Understanding tau hyperphosphorylation, neurofibrillary tangle formation, and tau's relationship to APP/Aβ pathology in AD.

Iqbal, K., Liu, F., Gong, C.-X., & Grundke-Iqbal, I. (2010). Tau in Alzheimer disease and related tauopathies. *Current Alzheimer Research*, *7*(8), 656-664. https://doi.org/10.2174/156720510793611592

Guo, T., Noble, W., & Hanger, D. P. (2017). Roles of tau protein in health and disease. *Acta Neuropathologica*, *133*(5), 665-704. https://doi.org/10.1007/s00401-017-1707-9

---

## Synaptic Dysfunction in AD
**Used for**: Understanding early synaptic loss, synaptic dysfunction preceding neuronal death, and correlation of cognitive decline with synaptic markers rather than plaque burden.

Selkoe, D. J. (2002). Alzheimer's disease is a synaptic failure. *Science*, *298*(5594), 789-791. https://doi.org/10.1126/science.1074069

Hong, S., Beja-Glasser, V. F., Nfonoyim, B. M., Frouin, A., Li, S., Ramakrishnan, S., Merry, K. M., Shi, Q., Rosenthal, A., Barres, B. A., Lemere, C. A., Selkoe, D. J., & Stevens, B. (2016). Complement and microglia mediate early synapse loss in Alzheimer mouse models. *Science*, *352*(6286), 712-716. https://doi.org/10.1126/science.aad8373

---

## Clustering Algorithms
**Used for**: Graph-based clustering algorithms (Louvain, Leiden) used in Seurat to identify transcriptionally distinct cell populations.

Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics: Theory and Experiment*, *2008*(10), P10008. https://doi.org/10.1088/1742-5468/2008/10/P10008

Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: Guaranteeing well-connected communities. *Scientific Reports*, *9*, 5233. https://doi.org/10.1038/s41598-019-41695-z

---

## Human Brain Cell Atlas
**Used for**: Reference atlases of human brain cell types for validation of cell-type annotations and understanding normal brain cellular diversity.

Hodge, R. D., Bakken, T. E., Miller, J. A., Smith, K. A., Barkan, E. R., Graybuck, L. T., Close, J. L., Long, B., Johansen, N., Penn, O., Yao, Z., Eggermont, J., Höllt, T., Levi, B. P., Shehata, S. I., Aevermann, B., Beller, A., Bertagnolli, D., Brouner, K., ... Lein, E. S. (2019). Conserved cell types with divergent features in human versus mouse cortex. *Nature*, *573*(7772), 61-68. https://doi.org/10.1038/s41586-019-1506-7

---

## AD Neuropathology Standards
**Used for**: Understanding Braak staging, CERAD criteria, and neuropathological diagnostic standards for AD used to classify postmortem samples.

Braak, H., & Braak, E. (1991). Neuropathological stageing of Alzheimer-related changes. *Acta Neuropathologica*, *82*(4), 239-259. https://doi.org/10.1007/BF00308809

Hyman, B. T., Phelps, C. H., Beach, T. G., Bigio, E. H., Cairns, N. J., Carrillo, M. C., Dickson, D. W., Duyckaerts, C., Frosch, M. P., Masliah, E., Mirra, S. S., Nelson, P. T., Schneider, J. A., Thal, D. R., Thies, B., Trojanowski, J. Q., Vinters, H. V., & Montine, T. J. (2012). National Institute on Aging–Alzheimer's Association guidelines for the neuropathologic assessment of Alzheimer's disease. *Alzheimer's & Dementia*, *8*(1), 1-13. https://doi.org/10.1016/j.jalz.2011.10.007

---

## Marker Gene Databases
**Used for**: CellMarker, PanglaoDB, and other curated databases of cell-type-specific marker genes used to define SCINA signatures for astrocytes (GFAP, AQP4), microglia (CX3CR1, TMEM119), neurons (SYP, RBFOX3), etc.

Zhang, X., Lan, Y., Xu, J., Quan, F., Zhao, E., Deng, C., Luo, T., Xu, L., Liao, G., Yan, M., Ping, Y., Li, F., Shi, A., Bai, J., Zhao, T., Li, X., & Xiao, Y. (2019). CellMarker: A manually curated resource of cell markers in human and mouse. *Nucleic Acids Research*, *47*(D1), D721-D728. https://doi.org/10.1093/nar/gky900

Franzén, O., Gan, L.-M., & Björkegren, J. L. M. (2019). PanglaoDB: A web server for exploration of mouse and human single-cell RNA sequencing data. *Database*, *2019*, baz046. https://doi.org/10.1093/database/baz046

---

## Reproducibility & Open Science
**Used for**: Best practices for computational reproducibility, version control (Git/GitHub), and transparent reporting of bioinformatics workflows.

Sandve, G. K., Nekrutenko, A., Taylor, J., & Hovig, E. (2013). Ten simple rules for reproducible computational research. *PLOS Computational Biology*, *9*(10), e1003285. https://doi.org/10.1371/journal.pcbi.1003285

Wilson, G., Bryan, J., Cranston, K., Kitzes, J., Nederbragt, L., & Teal, T. K. (2017). Good enough practices in scientific computing. *PLOS Computational Biology*, *13*(6), e1005510. https://doi.org/10.1371/journal.pcbi.1005510

---

## R Programming & Bioconductor
**Used for**: R statistical programming language and Bioconductor project providing core infrastructure for genomic data analysis.

R Core Team. (2023). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing. https://www.R-project.org/

Huber, W., Carey, V. J., Gentleman, R., Anders, S., Carlson, M., Carvalho, B. S., Bravo, H. C., Davis, S., Gatto, L., Girke, T., Gottardo, R., Hahne, F., Hansen, K. D., Irizarry, R. A., Lawrence, M., Love, M. I., MacDonald, J., Obenchain, V., Oleś, A. K., ... Morgan, M. (2015). Orchestrating high-throughput genomic analysis with Bioconductor. *Nature Methods*, *12*(2), 115-121. https://doi.org/10.1038/nmeth.3252

---

## GitHub Repository
**Used for**: Version control and code/documentation sharing for reproducible research.

GitHub Repository: https://github.com/bhavnam/Fa25-Project6-AD-Transcriptomics

---

*Last updated: December 8, 2025*
*Total references: 63*
