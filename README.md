# ferroviz3
R package to perform ferroptosis analyses on limma output cross with mitochondria and metabolism databases


## code to install rcd package
```r
library(devtools)
install_github("cdesterke/ferroviz3")
```


## load data and volcanoplot on limma output
```r
library(ferroviz3)
data(resulths)
data(ferrdb)
volcano_ferroptosis(
  resulths,
  ferrdb,
  fc_thresh = 0.5,
  p_thresh = 0.05,
  plot_title = "DEGs and ferroptosis",
  base_size = 14,
  label_size = 4
)

```

![volcano](https://github.com/cdesterke/ferroviz3/blob/main/volcano.png)


## barplot of ferroptosis related features
```r
barploths(resulths, ferrdb, fc = 0.25, size = 16)

```
![bar](https://github.com/cdesterke/ferroviz3/blob/main/barplot.png)


## table of significant ferroptosis related features
```r
table<-export_ferroptosis_table(resulths, ferrdb, fc = 0.25)
head(table)
write.csv(table,file="ferroptosis_features.csv",row.names=F)
```
![table](https://github.com/cdesterke/ferroviz3/blob/main/table.png)


## heatmap significant ferroptosis versus mitochondria
```r
data(mito)
heatmap_ferro_mito(resulths, mito, fc = 0.1)
```
![mito](https://github.com/cdesterke/ferroviz3/blob/main/mito.png)

## metabolism KEGG abbreviation table
```r
kegg_table <- as.data.frame(table_kegg_abbrev(kegg2))
print(kegg_table)
```
![abb](https://github.com/cdesterke/ferroviz3/blob/main/abb.png)

## heatmap significant ferroptosis versus metabolism KEGG
```r
load(kegg2)
heatmap_ferro_kegg(resulths, kegg2, fc = 0.1)
```
![kegg](https://github.com/cdesterke/ferroviz3/blob/main/kegg.png)


# REFERENCES

> Zhou N, Peng L, Luo Q, Yin T, Sun H, Zhang Y, Shi X, Peng X, Bao J, Ning Y, Yuan X. FerrDb V3: expanding the manually curated resource for regulators and disease associations from ferroptosis to regulated cell death. Nucleic Acids Res. 2026 Jan 6;54(D1):D572-D582. doi: 10.1093/nar/gkaf1119. PMID: 41171133; PMCID: PMC12807728.

> Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Res. 2000 Jan 1;28(1):27-30. doi: 10.1093/nar/28.1.27. PMID: 10592173; PMCID: PMC102409.

> Calvo SE, Clauser KR, Mootha VK. MitoCarta2.0: an updated inventory of mammalian mitochondrial proteins. Nucleic Acids Res. 2016 Jan 4;44(D1):D1251-7. doi: 10.1093/nar/gkv1003. Epub 2015 Oct 7. PMID: 26450961; PMCID: PMC4702768.
