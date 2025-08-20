#### Source Dir Tree

After you download the .tar.gz package from the URL, unpack it and reorganize the extracted files into the following directory structure:

```shell
.
├── Atlas.Model
│   ├── model.pkl
│   └── top10_markers.tsv
├── genome
│   └── CR-GRCh38
├── genomeIndex
│   ├── hg38.exclusive.idx
│   └── hg38_te.exclusive.idx
├── softwares
│   ├── bcsubset
│   └── scTE
└── whitelist
    ├── 3M-february-2018.txt
    ├── 737K-april-2014_rc.txt
    ├── 737K-august-2016.txt
    └── whitelist.md
```
Note: you can install softwares at other path, you just need update the config of "tool.tsv" to regist bcsubset / scTE
