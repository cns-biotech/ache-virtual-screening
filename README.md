# Acetylcholinesterase inhibitors in ChEMBL: does similarity to donepezil predict activity?

I collected every IC50 measurement for human acetylcholinesterase (AChE) available in ChEMBL, cleaned the data, calculated physicochemical descriptors, and checked whether molecules that resemble donepezil are also more potent. I then trained a Random Forest model on pIC50. On a scaffold split, which tests the model on chemistry it has not seen, it reaches R² = 0.49 (0.74 on a random split).

## Question

AChE is a therapeutic target in Alzheimer's disease. Patients lose cholinergic neurons, and inhibiting AChE raises acetylcholine levels, which partly makes up for that loss. Donepezil is an approved AChE inhibitor.

I wanted to know three things:

1. How potent are the molecules reported against AChE in public data?
2. How well does structural similarity to donepezil tell us whether a molecule is potent?
3. Can a model predict potency for chemistry it has not seen?

## Data

- **Source:** ChEMBL, target `CHEMBL220` (human acetylcholinesterase), activity type IC50
- **Accessed:** 7 November 2026 via chembl_webresource_client, database release ChEMBL_37 (released 1 May 2026).
- **Starting point:** 10,079 IC50 records

### Curation

| Step | Records left |
|---|---|
| Initial IC50 records for CHEMBL220 | 10,079 |
| Units standardized to nM | 8,723 |
| Exact values only (relation "=") | 7,827 |
| No missing values; values > 0 | 7,825 |
| One value per molecule (median) | 6,287 molecules |
| Salts and mixtures removed, values re-aggregated | 6,174 molecules |
| pIC50 between 3 and 11 | 6,044 molecules |

Notes on each step:

- **Units.** 1,265 records had no standardized unit (no value reported). 83 records were given in µg/mL and 8 had malformed unit strings. I excluded the µg/mL and malformed records: converting µg/mL to nM needs the molecular weight of every compound, and these 91 records are under 1% of the data. 
- **Censored values.** I excluded 859 records with >, 12 with <, 2 with ~, 2 with >> and 21 with no relation recorded: 896 in total, all censored values with no measured IC50.
- **Salts and mixtures.** 1,035 records contained several fragments. I kept the largest one. After converting to canonical SMILES, 113 of these molecules merged with molecules that were already in the dataset.
- **Suspicious values.** 17 different molecules shared the same IC50 of 0.000005 nM (5 fM). That is far below the picomolar potencies reported even for the most potent AChE inhibitors, so I treated it as a curation artifact. These records fall outside the pIC50 range below and were removed at that step.
- **pIC50 range.** I chose the range 3 to 11 by hand and report it here for transparency. It removed 130 molecules (6,174 to 6,044).

## Methods

`pIC50 = 9 - log10(IC50 in nM)`

Descriptors calculated with RDKit:

| Descriptor | What it is |
|---|---|
| MW | Molecular weight in daltons (Da) |
| LogP | Calculated octanol-water partition coefficient (RDKit's Crippen method), a measure of lipophilicity |
| HBD | Number of hydrogen bond donors (N-H and O-H groups) |
| HBA | Number of hydrogen bond acceptors (N and O atoms) |
| TPSA | Topological polar surface area in Å², used to estimate membrane permeability |

**Lipinski's rule of five** is a rough filter for drug-likeness. Poor oral absorption is more likely when a molecule breaks more than one of these limits: MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10. I report both the strict count (no violations) and Lipinski's own criterion (at most one violation).

**Morgan fingerprints** (radius 2, 2048 bits, also called ECFP4) describe each molecule as a 2048-bit vector that encodes the neighborhood of every atom up to 2 bonds away.

**Tanimoto similarity to donepezil** (CHEMBL502) is the number of fingerprint bits two molecules share, divided by the number of bits set in either one. It runs from 0 (nothing in common) to 1 (identical fingerprints).

## Statistics

**Spearman correlation between Tanimoto similarity to donepezil and pIC50:** rho = 0.13, p = 1e-23. The correlation is statistically significant but weak: similarity explains under 2 percent of the variance in rank.

**Mann-Whitney U test,** molecules with similarity 0.4 or above (n = 144) against the rest (n = 5,900): p = 2.8e-16, difference in medians 0.96 log units, about ninefold in potency.

With a sample of this size almost any difference reaches significance, so the effect size is reported alongside the p-value.

**Tests were chosen for the data:** Spearman rather than Pearson because only a monotonic relationship is expected and the data contain outliers; Mann-Whitney rather than a t-test because group sizes differ greatly and normality is not assumed.

## Results

Median pIC50 is 5.85 (IC50 about 1.4 µM), so most compounds reported against AChE are moderately potent at best. Defining an active compound as pIC50 >= 7 (IC50 <= 100 nM), 1,397 of 6,044 molecules (23.1 percent) are active. The threshold is my own and is stated here so the number can be interpreted.

Median TPSA is 56.8 Å², below the value of 90 often used as a rough ceiling for blood-brain barrier penetration. This matters for a CNS target.

Molecules with Tanimoto similarity 0.4 or above to donepezil (n = 144) have a median pIC50 of 6.77 against 5.81 for the rest.

Donepezil itself is present in the dataset with a similarity of 1.0 and pIC50 about 7.5, which serves as a sanity check.

![Distributions](fig_distributions.png)

**Lipinski violations by criterion:** LogP above 5 in 2,076 molecules (34.3 percent), MW above 500 in 1,278 (21.1 percent), HBD above 5 in 34, HBA above 10 in 36. Lipophilicity is the dominant reason for failure, hydrogen bonding is almost never the limiting factor.

**Distribution of violations:** 0 in 3,713 molecules, 1 in 1,277, 2 in 1,025, 3 in 19, 4 in 10.

Under the strict reading (no violations at all) 3,713 of 6,044 molecules pass (61.4 percent). Under Lipinski’s original formulation, which allows a single violation, 4,990 molecules pass (82.6 percent). Both numbers are reported; the second follows the rule as published.

Similarity to donepezil alone is a weak predictor of potency across the dataset (rho = 0.13), but compounds with high similarity are enriched in potent molecules.

## QSAR model

I trained a Random Forest regressor on pIC50 using the 2,048-bit Morgan fingerprint plus the five descriptors above (2,053 features; 300 trees, `min_samples_leaf=2`, fixed seed). Tanimoto similarity to donepezil was deliberately left out of the features: it is a distance to one particular molecule, not a property of the compound, and it would not exist in a real screening campaign.

### Scaffold split

Molecules were grouped by their Bemis-Murcko scaffold (2,539 scaffolds for 6,044 molecules). Whole scaffold groups were assigned to one side of the split, largest group first, so that no scaffold occurs on both sides: 4,835 training and 1,209 test molecules, zero shared scaffolds.

| Split | R² on test |
|---|---|
| Random 80/20 | 0.742 |
| Scaffold 80/20 | 0.489 |

Same model, same features, same seed. The random split overstates performance by 0.25 in R². ChEMBL contains large series of analogs — 144 molecules here reduce to a bare benzene scaffold, 77 to the tacrine core — and under a random split part of a series trains the model while the rest tests it, so the model is rewarded for recognizing a series it has already seen. The scaffold number describes performance on chemistry the model has not seen, and that is the number I report. The two groups of 144 are unrelated: only N molecules belong to both.

On the scaffold split, RMSE is 1.00 log units against 1.40 for a baseline that always predicts the mean training pIC50: a 29 percent reduction in error.

The scaffold split is also less strict than it sounds. Murcko scaffolds are compared as SMILES strings, so two close analogs whose rings differ by a single atom receive different scaffolds and can still end up on opposite sides of the split. The true generalization gap is therefore at least as large as the table suggests.

Predicted against true values shows the flattening typical of an averaging regressor: low values are overpredicted, high values underpredicted. Compounds with a true pIC50 above 9 are predicted at 6 to 8. The model is usable for ranking candidates, not for quoting an absolute potency. Spearman correlation between predicted and true pIC50 on the test set is rho = 0.696, so the ranking is largely preserved even though absolute values are compressed.

### What the model uses

One fingerprint bit, bit 1048, carries 10.9 percent of the total feature importance, two and a half times the next feature. Decoding it with RDKit's bit-info map shows a 4-aminoquinoline motif — the core shared by tacrine and its analogs, the first drug class approved for Alzheimer's disease. 675 of the 6,044 molecules contain it, with a median pIC50 of 7.42 against 5.66 for the rest. The most important feature therefore corresponds to a known pharmacophore, which is consistent with the model learning real structure-activity relationships, although it may also reflect the size of this series. This last comparison is descriptive and uses the whole dataset, training molecules included, so it is not an independent validation.

The remaining top features are a mixture of further bits and the descriptors LogP, MW and TPSA.

![QSAR](fig_qsar.png)

A model that does not see similarity to donepezil reaches R² = 0.49 on unseen scaffolds, so structure carries much more information about potency than similarity to one drug.

## Limitations

- IC50 values come from different assays and laboratories, so they are not strictly comparable.
- Similarity does not guarantee activity. Small structural changes can cause large drops in potency (activity cliffs).
- Excluding censored values (">") removes mostly weak or inactive compounds, so the data lean toward measurable activity.
- ChEMBL contains many series of analogs built around known drugs. Similarity to donepezil may partly reflect which compounds chemists chose to make and test.
- I used a single similarity cut-off (0.4) chosen by hand; a different threshold would change the size of the similar group.
- Keeping only the largest fragment is wrong for genuine two-component mixtures.
- Impurity-based feature importances in a Random Forest favor continuous features over binary bits, so the position of the descriptors in the ranking should be read with care.
- No hyperparameter search was performed and there is no validation set; the test set was used once, for the numbers reported here.
- Every number comes from a single run with one seed, so I cannot give a confidence interval.

## Reproducing

All results use random_state = 0. Run with scikit-learn 1.6.1 and rdkit 2026.09.1

Open ache_virtual_screening.ipynb in Google Colab and run all cells. ChEMBL is updated regularly, so the numbers may differ slightly from the ones above.
