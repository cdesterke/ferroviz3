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
load(resulths)
load(ferrdb)
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
barploths(resulths, ferrdb, fc = 0.25, size = 16)

```
