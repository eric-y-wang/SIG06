# SIG06

### <ins>Goal<ins>
Screen single ligands and identifiy candidates for combinatorial studies in activated CD4 T cells

### <ins>Experimental Design<ins>
**Cells**: TCRb+CD4+CD44+CD62L- naive CD4 T cells were sorted from 9W old B6 mice straight from JAX and kept in autoclaved cages on DND. Mice were harvested and split into 3 biological replicates, 3F and 2M mice per group. Mice were from 5 different cages. Each group consisted of mice from each cage (to minimize cage effects influencing samples). Cells were isolated from skin-draining LNs, mesenteric LNs, and spleen. Cells were activated for 24h with 2.5 ug/mL plate-bound anti-CD3, 3 ug/mL soluble anti-CD28, and 30 U/mL hIL2. After 24h, cells were lifted from plate, washed once, and rested in cRMPI for 12h.

**Treatment**: Cells were treated in 96-well plates at 100,000 cells per well in 100 ul cRPMI. Cells were spinfected with virus or recombinant proteins 2000g x 30m and then placed into TC incubator for 6h (kept in viral supernatent). Edge-wells were filled with PBS during incubation. Cells were then lysed and RNA was extracted. 

**Sequencing**: RNA-seq libraries were prepared using BRB-seq with no normalization of RNA quantity. RNA-seq libraries were processed in 4 separate plates.

### <ins>Analysis<ins>
#### <ins>Identify robust clusters of ligand-driven transcriptional phenotypes.<ins>
Ultimately, the goal of this experiment is to identify single ligands that would be appropriate to use in combinatorial studies. It would be valuable to identify groups of ligands that have shared transcriptional programs and to examine overlap between these clusters. It would also be interesting to examine how these clusters correspond to receptor identity and known biological functions.
##### Files:

#### <ins>Decompose transcriptional phenotypes<ins>
Ideally, we can decompose the transcriptional signatures in some manner to identify co-regulated gene modules that are biologically interpretable. I guess this would fall into the class of "factor analysis" although a lot of those tools are developed for scRNA-seq and I'm not sure how well they work on the more limited sample size of bulk RNA-seq.
* **PCA**: Perhaps we can just use the principal components since by definition they represent orthogonal groups of genes. Need to focus on "coherent" PCs as some of these will capture variation due to technical factors. Will probably need to filter on HVGs too.

#### <ins>Identify covariates driving variation between groups.<ins>
There seems to consistently be some type of batch effect driving divergence of linker-only samples and variation in ligand-treated samples with weaker perturbation effects. For ligands with strong perturbations, I think the strength of transcriptonal changes masks these group to group variation. I can imagine a number of potential reasons for this.
* **sequencing depth (can be fixed)**
  * Because of the nature of BRB-seq, because I didn't do any RNA normalization prior to input into the library given the equal cell numbers in each well, and because the different library pools weren't sequenced equally there's some level of variability in the sequencing depth across given plates and samples. This can be seen in **/processing/count_matrix_generation.Rmd**.
* **true biological variability**
* **differences in RNA-seq plate to plate library preparation**

##### Files:  

#### <ins>Examine believability of "type 1" signature across diverse ligands.<ins>

##### Files:  

#### <ins>Compare concordance between recombinant protein and viral ligand stimulation.<ins>

##### Files:

