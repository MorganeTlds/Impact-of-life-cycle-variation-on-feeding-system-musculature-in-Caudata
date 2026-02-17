# Impact of Life Cycle Variation on Feeding System Musculature in *Caudata*
**2026‑02‑17 Morgane Taillades**

## General Information

**Project:** Impact of Life Cycle Variation on Feeding System Musculature in Caudata

**Description**  
This project investigates how differences in life cycles, ecological transitions, morphological transformations, and adult habitat influence jaw and hyoid musculature involved in feeding across salamanders.  
A total of 25 species representing different developmental strategies were dissected, from fully metamorphosed taxa to facultative and fully paedomorphic species.

For each species, we performed detailed dissections of the feeding system and quantified for each muscle of the system:
- Muscle volume  
- Fiber length
- PCSA (physiological cross sectional area)

Our results reveal differences in muscle organization associated with ecological strategies, with additional strong effects of head size. Consistent patterns emerge depending on ecological transitions, morphological transformations, and adult habitat types.


## Data Description

[Raw_Data.xlsx](Raw_Data.xlsx)

Excel file for the raw data for each dissected specimen:
- Specimen ID (if there was one)
- Family,
- Species  
- Life cycle   
- Adult habitat (aquatic, terrestrial, semi-aquatic)  
- Ecological transition (yes/no)  
- Morphological transformation (yes/no)  
- Morphometric measurements (mm): SVL, HL, HW, HH, LJ  
- Individual muscle weights (mg)  
- Muscle fiber lengths (mm)  
- PCSA values (calculated following the manuscript formula,cm²)

[Mean_by_species.xlsx](Mean_by_species.xlsx)

Excel file for species-level averages (when multiple individuals shared the same developmental and ecological strategies):
- Family
- Species  
- Life cycle strategy  
- Adult habitat  
- Ecological transition (T_eco)  
- Morphological transformation (T_morpho)  
- Morphometric measurements (mm)  
- Muscle weights (mg)  
- Muscle PCSA (cm²)

[BigTree.tre](BigTree.tre)

Time-calibrated salamander phylogeny from Stewart & Wiens (2025).

[PhySalamanders.nex](PhySalamanders.nex)

Nexus file, subset of species used in this study, extracted from *BigTree.tre*, with a polytomy added for *Ambystoma mexicanum metamorph* and *Ambystoma andersoni metamorph*.

## R Scripts

[Phylogeny.R](Phylogeny.R)

R script which permits to generate **PhySalamanders.nex** from **BigTree.tre**.

[Analysis.R](Analysis.R)

R script used for all statistical analyses (requires **Mean_by_species.xlsx**).

**Main workflow:**
1. Creation of dataframe *df*; each individual muscle was grouped by his functional group 
2. Log10 transformation of variables  
3. Loading of phylogeny (PhySalamanders.nex)  

**Three ecological variables (Ecological transition, Morphological transformation and Adult habitat) can be tested. 
Select the appropriate one at lines 69–71.**

**Analytical structure (for both muscle volume and PCSA):**
- Evolutionary model testing  
  - OU (Ornstein–Uhlenbeck) currently performs best  
- Phylogenetic PCA  
- Phylogenetic MANCOVA  
- Phylogenetic ANCOVA  

### **PCA_plot.R**

[PCA_plot.R](PCA_plot.R)

R script which generates visualizations of the PCA phylomorphospace.

## Abbreviations

### Morphometric measurements
| Abbreviation | Meaning |
|-------------|---------|
| **SVL** | Snout-vent length |
| **HL** | Head length |
| **HW** | Head width |
| **HH** | Head height |
| **LJ** | Lower jaw length |

### Muscle parameters
| Abbreviation | Meaning |
|-------------|---------|
| **pds** | Weight |
| **pcsa** | Physiological cross-sectional area |
| **lg** | Fiber length |
|**Vol**|Muscle volume|
|**pcsa**|Muscle PCSA|

### Individual muscles
| Code | Muscle Name |
|------|-------------|
| MQP | quadratopectoralis |
| MIM | intermandibularis |
| MIH | interhyoideus |
| MIM_MIH | fused intermandibularis + interhyoideus |
| MDM | depressor mandibulae |
| MDMa | anterior depressor mandibulae |
| MDMp | posterior depressor mandibulae |
| MGG | genioglossus |
| MGH | geniohyoideus |
| MSbH | subhyoideus |
| MSbR | subarcualis rectus |
| MRC | rectus cervicis |
| MRCsup | superficial rectus cervicis |
| MRCpro | profundus rectus cervicis |
| MAME | adductor mandibulae externus |
| MAMEsup | superficial adductor mandibulae externus |
| MAMEpro | profundus adductor mandibulae externus |
| MAMIproI | adductor mandibulae internus profundus I |
| MAMIproII | adductor mandibulae internus profundus II |
| MAMIsup | adductor mandibulae internus superficialis |
| MAMP | adductor mandibulae posterior |
| MPt | pterygoideus |
| MBh | branchiohyoideus |
| Mbhint | branchiohyoideus internus |

### Functional muscle group
|   Code   | Muscle group Name |
|----------|-------------|
| HyoLev   | Hyoid levators |
| HyoPro   | Hyoid protractors |
| HyoRet   | Hyoid retractors |
| MouthOp  | Mouth opening muscles |
| MouthClo | Mouth closing muscles |

### Evolutionary models
| Code | muscle name Name |
|------|-------------|
| BM      | Brownian Motion |
| OU      | Ornstein–Uhlenbeck |
| EB      | Early Burst |
| LAMBDA  | Pagel’s λ model |


## Reference

Stewart, A. A., & Wiens, J. J. (2025). *A time‑calibrated salamander phylogeny including 765 species and 503 genes.*  
**Molecular Phylogenetics and Evolution**, 204, 108272.  
https://doi.org/10.1016/j.ympev.2024.108272
