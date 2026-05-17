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
