# sc-flow

## broad order of steps:

In general, a typical scRNA-seq workflow will go through these steps (in order):
1. align
2. preprocess
3. cluster
4. plot
5. diff\_exp

Each of these corresponds to an argument for the ```--workflow``` command in the sc-flow pipeline:
1. align: ```--workflow align```
