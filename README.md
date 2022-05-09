# Alpha/Beta Cell Model
This python package implements computational model, implementing mechanisms of glucose-dependent glucagon and insulin secretion from pancreatic alpha and beta cells, respectively. Metabolism is modeled as a joint class (with different parameters for each cell) and depends on the models by Grubelnik et al. (2020)[^grubelnik]. The glucose-dependent intracellular cAMP concentration is modeled according to Zmazek et al. (2021)[^zmazek2021]. The glucagon secretion depends on the previously developed electrophysiological model by Montefusco & Pedersen (2015)[^montefusco]. Modifications to the cAMP signaling and secretion components of the model additionally implement the amino acid (AA) concentration-dependent glucagon secretion[^zmazek2022].
  [^montefusco]: Montefusco, F., & Pedersen, M. G. (2015). Mathematical modelling of local calcium and regulated exocytosis during inhibition and stimulation of glucagon secretion from pancreatic alpha‐cells. The Journal of physiology, 593(20), 4519-4530.
  [^grubelnik]: Grubelnik, V., Zmazek, J., Markovič, R., Gosak, M., & Marhl, M. (2020). Modelling of energy-driven switch for glucagon and insulin secretion. Journal of Theoretical Biology, 493, 110213.
  [^zmazek2021]: Zmazek, J., Grubelnik, V., Markovič, R., & Marhl, M. (2021). Role of cAMP in double switch of glucagon secretion. Cells, 10(4), 896.
  [^zmazek2022]: Zmazek, J., Grubelnik, V., Markovič, R., & Marhl, M. (2022). Modeling the Amino Acid Effect on Glucagon Secretion from Pancreatic Alpha Cells. Metabolites, 12(4), 348..

## Requirements
This package requires Python 3 installed with additional packages numpy, scipy, matplotlib, cython, multiprocessing, and pathlib.

## Installation
Clone package from GitHub (git must be installed):
```bash
git clone https://github.com/janzmazek/Alpha-Beta-Model/
```
Compile cython module:
```bash
cd Alpha-Beta-Model/CellModel 
python setup_montefusco.py build_ext --inplace
cd ../..
```
Run package:
```bash
python Alpha-Beta-Model
```
