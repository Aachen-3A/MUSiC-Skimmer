Information on datasets & run ranges from:
   * https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters?rev=60#EB_or_EE_Xtals_with_large_laser
   * https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions?rev=403#2012_A_B_C_D_datasets_re_reco_pr
   * and DAS queries.

Rereco with CMSSW_5_3_3_patch1:

Run A:
   /SingleElectron/Run2012A-13Jul2012-v1/AOD
   runs: 190456 -- 193621

Run B:
   /SingleElectron/Run2012B-13Jul2012-v1/AOD
   runs: 193834 -- 196531

Run C:
   /SingleElectron/Run2012C-24Aug2012-v1/AOD
   runs: 198022 -- 198523


Additional PromptRecos:

Run C v1:
DON'T SKIM THIS DATASET!
The addtional runs in this dataset (compared to the Run C ReReco)
are not included in the golden JSON. They are only high pile-up runs
We don't need them, since we have the Aug 24 ReReco.
   /SingleElectron/Run2012C-PromptReco-v1/AOD
   runs: 197770 -- 198913

Run C v2:
This dataset has some additional runs worth skimming.
   /SingleElectron/Run2012C-PromptReco-v2/AOD
   runs: 198934 -- 203755

Run D:
   /SingleElectron/Run2012D-PromptReco-v1/AOD
   runs: 203773 -- 209465

