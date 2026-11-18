# IBDverse annotation vocabulary vs our states (name-level only)

Source: `IBDverse_multitissue_annotation.csv` (EBI BioStudies S-BSST2922, CC0): 102 clusters, categories T (22), Myeloid (18), B (10), Plasma (4), plus epithelial and mesenchymal. It lists cluster -> category -> cell-type name only; no expression data.

`ibdverse_vocabulary_check.csv` is my reading of the names, not a computed comparison. It shows which of our high/medium-confidence states have a same-name, related-name or no counterpart in IBDverse's list. 3 of our 34 high/medium states have no counterpart by name: MAIT/Th17-like, mature DC (LAMP3+) and neutrophil. The others have a same or related name.

Limits: names only (no cell-level or marker-level check); the table is the multi-tissue atlas table (terminal ileum, rectum, blood) and may differ from the 57 cell types of the terminal-ileum paper; IBDverse is Crohn's vs healthy, not inflamed vs non-inflamed. Our states are unchanged.
