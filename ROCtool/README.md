# How to Make Use of ROCtool

This is a short file on how to best make use of ROCtool, a metric for the accuracy of DNase1-seq and ATAC-seq footprinting data.

## Input

Pass the command

```
sh ROCtool.sh -fp your_footprint_file_here.bed > your_AUC_file_here.txt
```

## Output

The footprint file will be processed, editing the contents of the ChIP-validation folder. Then, ROCtool.py will be called on the ChIP-validation folder and the areas under the sensitivity-specificity curve (AUC) for each TF in the folder will be returned as a .txt file.

The format will be

```
ENSG1 AUC: 0.8
ENSG2 AUC: 0.9
...
Mean AUC: 0.77829412
```

So, if you just want the mean value, you can write in python:

```
openfile = open('your_AUC_file_here.txt','r')
lastline = openfile.read().split('\n')[-1]
openfile.close()
meanAUC = float(lastline.split(' AUC: ')[1])
print(meanAUC)
```