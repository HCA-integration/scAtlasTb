import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd



def trajectory_evaluation(
        adata,
        output_type='embed',
        batch_key=None,
        label_key=None,
        celltypes_list={},
        GT_traj_key='CCF score',
        **kwargs
    ):

    from sctram.api._lower_level import TrajectoryEvaluationAPI
    from sctram.input import InputTrajectories


    emb ='X_emb' if output_type == 'embed' else 'X_pca'

    print(f"Evaluating trajectory preservation for cell types: {list(celltypes_list.values())} ...", flush=True)

    # Check CCF score (or ground truth trajectory) exists
    if GT_traj_key not in adata.obs.columns:
        raise ValueError(f"{GT_traj_key} column is missing from adata.obs. Please ensure it is present for trajectory evaluation.")
    
    
    # --------------------------------------------------
    # Precompute GT locations for all cell types
    # --------------------------------------------------

    celltype_key = list(celltypes_list.keys())[0]
    celltypes=list(celltypes_list.values())[0]

    
    locations = {}
    for ct in celltypes:
        locations[ct]=[]
        adata_sub=adata[adata.obs[celltype_key]==ct]
        locations[ct].extend(list(adata_sub.obs[GT_traj_key].unique()))
    
    
    
    results = []     
    for ct in celltypes:
        print(ct)
        
        adata_ct = adata[adata.obs[celltype_key]==ct]
    
        if len(adata_ct)>10:
    
            adata_emb = ad.AnnData(X=adata_ct.obsm[emb].copy(), obs=adata_ct.obs.copy())
    
            sc.pp.neighbors(adata_emb, use_rep='X', n_neighbors=15)
            
    
            locations_ct_ordered=sorted(locations[ct])
            gt_traj = {}
            gt_traj[ct] = [
                    (str(locations_ct_ordered[i]), str(locations_ct_ordered[i+1]), {"distance": float(locations_ct_ordered[i+1]) - float(locations_ct_ordered[i])})
                    for i in range(len(locations_ct_ordered) - 1)
                ]
            
            api = TrajectoryEvaluationAPI(
                adata=adata_emb,
                input_trajectory=InputTrajectories(gt_traj).get_trajectory(ct, include_additional_nodes=False),
                labels_obs=GT_traj_key,
                root_label=sorted(adata_emb.obs[GT_traj_key].unique())[0],
                logger_level="DEBUG"
            )
    
            api.evaluate_embedding(metrics = [
                "neighborhood_preservation_score"
            ])
            api.evaluate_pseudotime(metrics=[
                            "spearman_correlation",
                            "kendall_correlation"
                                                ])
            api.evaluate_adjacency(metrics=[
                            "persistence_diagram_distance",
                            "average_shortest_path_difference"
                                                    ])
                    
        
            df = api.get_all_results()
            df['cell type'] = ct
    
            results.append(df)
    
    results_all = pd.concat(results, ignore_index=True)
    
    # ------------------------------------------------------
    # Normalize metrics independently
    # ------------------------------------------------------
    normalized_scores = []
    
    for metric in results_all["metric"].unique():
    
        df_metric = results_all[results_all["metric"] == metric].copy()
    
        vals = df_metric["score"].values.astype(float)
    
        # Min-max normalization
        vmin = np.nanmin(vals)
        vmax = np.nanmax(vals)
    
        if vmax == vmin:
            norm = np.ones_like(vals)
        else:
            norm = (vals - vmin) / (vmax - vmin)
    
        # Flip if lower is better
        lower_better_metrics = [
            "persistence_diagram_distance"
        ]
        if metric in lower_better_metrics:
            norm = 1 - norm
    
        df_metric["normalized_score"] = norm
    
        normalized_scores.append(df_metric)
    
    df_norm = pd.concat(normalized_scores)
    
    # ------------------------------------------------------
    # Aggregate composite score
    # ------------------------------------------------------
    df_agg = (
        df_norm
        .groupby(["cell type"])["normalized_score"]
        .mean()
        .reset_index()
    )
    
        
    
    metrics_names = [f"Traj pres score:{df_agg['cell type'].iloc[i]}" for i in range(len(df_agg))]
    scores = df_agg['normalized_score'].tolist()
    
    
    
    return scores, metrics_names






        
        