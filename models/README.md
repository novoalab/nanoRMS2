# nanoRMS2 models

This directory lists models trained by nanoRMS2.
Those models can be used with nanoRMS2 via
```python
src/encode_mods.py -m model.lzma
```

Available models:
- human_METTL3_WT_vs_KO.lzma - models (for 328 7-mers) trained with Mettl3 WT and KO rep1 (from  [Pratanwanich et al. Nat Biotechnol (2021)](https://doi.org/10.1038/s41587-021-00949-w))


If you wish to explore individual models, first install sklearn
```bash
pip install -U joblib scikit-learn==1.0.1
```
and then you can load all models
```python
import joblib
fn = "models/human_METTL3_WT_vs_KO.lzma"
mer2clf = joblib.load(fn)
len(mer2clf)
print(mer2clf.keys())
```
