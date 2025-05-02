#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["apply_t1_labels_transfo":"$params.apply_t1_labels_transfo",
                "output_dir":"$params.output_dir",
                "run_commit":"$params.run_commit",
                "use_commit2":"$params.use_commit2",
                "b_thr":"$params.b_thr",
                "nbr_dir":"$params.nbr_dir",
                "ball_stick":"$params.ball_stick",
                "para_diff":"$params.para_diff",
                "perp_diff":"$params.perp_diff",
                "iso_diff":"$params.iso_diff",
                "no_pruning":"$params.no_pruning",
                "no_remove_loops":"$params.no_remove_loops",
                "no_remove_outliers":"$params.no_remove_outliers",
                "min_length":"$params.min_length",
                "max_length":"$params.max_length",
                "loop_max_angle":"$params.loop_max_angle",
                "outlier_threshold":"$params.outlier_threshold",
                "length_weighting":"$params.length_weighting",
                "sh_basis":"$params.sh_basis",
                "run_afd_rd":"$params.run_afd_rd",
                "min_lesion_vol":"$params.min_lesion_vol",
                "use_similarity_metric":"$params.use_similarity_metric",
                "nbr_subjects_for_avg_connections":"$params.nbr_subjects_for_avg_connections",
                "processes_register":"$params.processes_register",
                "processes_commit":"$params.processes_commit",
                "processes_afd_rd":"$params.processes_afd_rd",
                "processes_avg_similarity":"$params.processes_avg_similarity",
                "processes_connectivity":"$params.processes_connectivity",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "Run Connectivity Construction"
log.info "============================="
log.info ""

log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "[Inputs]"
log.info "Root: $params.input"
log.info "Template: $params.template"
log.info "Labels list: $params.labels_list"
log.info "Labels image prefix: $params.labels_img_prefix"
log.info "Output directory: $params.output_dir"
log.info ""

log.info "Options"
log.info "======="
log.info "Apply transformation: $params.apply_t1_labels_transfo"
log.info "Run COMMIT: $params.run_commit"
log.info "Use COMMIT2: $params.use_commit2"
log.info "bval tolerance: $params.b_thr"
log.info "Nbr directions: $params.nbr_dir"
log.info "Ball & Stick: $params.ball_stick"
log.info "Parallel diffusion: $params.para_diff"
log.info "Perpendicular diffusion: $params.perp_diff"
log.info "Isotropic diffusion: $params.iso_diff"
log.info "NO pruning: $params.no_pruning"
log.info "NO loops removal: $params.no_remove_loops"
log.info "NO remove outliers: $params.no_remove_outliers"
log.info "Minimal length: $params.min_length"
log.info "Maximal length: $params.max_length"
log.info "Angle loops removal: $params.loop_max_angle"
log.info "Outliers removal threshold: $params.outlier_threshold"
log.info "Run AFD & RD: $params.run_afd_rd"
log.info "Length weighting: $params.length_weighting"
log.info "SH basis: $params.sh_basis"
log.info "Minimum lesion volume: $params.min_lesion_vol"
log.info "Use similarity metric: $params.use_similarity_metric"
log.info "Average from N subjects: $params.nbr_subjects_for_avg_connections"
log.info ""

log.info "Number of processes per tasks"
log.info "============================="
log.info "Template registration: $params.processes_register"
log.info "COMMIT: $params.processes_commit"
log.info "AFD & RD computation: $params.processes_afd_rd"
log.info "Average / Similarity: $params.processes_avg_similarity"
log.info "Compute Connectivity: $params.processes_connectivity"
log.info ""

root = file(params.input)
/* Watch out, files are ordered alphabetically in channel */
Channel
    .fromPath("$root/**/*t1.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{in_t1;subjects_for_count}

Channel
    .fromPath("$root/**/*$params.labels_img_prefix*labels.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .set{in_labels}

in_transfo = Channel
    .fromFilePairs("$root/**/{*0GenericAffine.mat,*1Warp.nii.gz}",
                    size: 2,
                    maxDepth:1,
                    flat: true) {it.parent.name}

in_tracking = Channel
    .fromFilePairs("$root/**/{*tracking*.*,}",
                    size: -1,
                    maxDepth:1) {it.parent.name}

Channel
    .fromPath("$root/**/*fodf.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{fodf_for_afd_rd;fodf_for_count}

Channel
    .fromFilePairs("$params.input/**/*lesion_mask.nii.gz",
        size: -1) { it.parent.name }
    .set{lesion_for_lesion_load}

Channel.fromPath(file(params.template))
    .into{template_for_registration;template_for_transformation_data;template_for_transformation_metrics;template_for_transformation_lesions}

Channel.fromPath(file(params.labels_list))
    .into{labels_list_for_compute;labels_list_for_visualize}

in_opt_metrics = Channel
    .fromFilePairs("$root/**/metrics/*.nii.gz",
                    size: -1,
                    maxDepth:2) {it.parent.parent.name}

in_dwi_data = Channel
    .fromFilePairs("$root/**/{*dwi.bval,*dwi.bvec,*dwi.nii.gz,*peaks.nii.gz}",
                    size: 4,
                    maxDepth:1,
                    flat: true) {it.parent.name}

(dwi_for_count, data_for_commit) = in_dwi_data
    .map{sid, bval, bvec, dwi, peaks ->
        [tuple(sid, dwi),
        tuple(sid, bval, bvec, dwi, peaks)]}
    .separate(2)

subjects_for_count.count().into{ number_subj_for_null_check; number_subj_for_compare_dwi; number_subj_for_compare_fodf; number_subj_for_compare_similarity}
dwi_for_count.count().into{ dwi_for_null_check; dwi_for_compare }
fodf_for_count.count().into{ fodf_for_null_check; fodf_for_compare }

number_subj_for_null_check
.subscribe{a -> if (a == 0)
    error "Error ~ No subjects found. Please check the naming convention, your --input path."}

number_subj_for_compare_similarity
.subscribe{a -> if (a < params.nbr_subjects_for_avg_connections && params.use_similarity_metric)
    error "Error --nbr_subjects_for_avg_connections is higher than the number of subjects. Please modify it or use --he number of subjects. Please modify it or use --use_similarity_metric"}


run_commit = params.run_commit
dwi_for_null_check
.subscribe{a -> if (a == 0 && params.run_commit)
    run_commit = false}
    log.warn "Warning ~ No DWI or peaks found. COMMIT will not be run."

run_afd_rd = params.run_afd_rd
fodf_for_null_check
.subscribe{a -> if (a == 0 && params.run_afd_rd)
    run_afd_rd = false}
    log.warn "Warning ~ No FODF found. AFD & RD will not be run."

number_subj_for_compare_dwi
    .concat(dwi_for_compare)
    .toList()
    .subscribe{a, b -> if (a != b && b > 0)
    error "Error ~ Mismatch between the number of subjects and DWI"}

number_subj_for_compare_fodf
    .concat(fodf_for_compare)
    .toList()
    .subscribe{a, b -> if (a != b && b > 0)
    error "Error ~ Mismatch between the number of subjects and FODF"}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

if (!params.apply_t1_labels_transfo) {
    in_t1
        .set{ori_anat}
    in_labels
        .set{ori_labels}
    anat_for_transformation = Channel.empty()
}
else {
    in_t1
        .join(in_labels)
        .join(in_transfo)
        .set{anat_for_transformation}
    ori_anat = Channel.empty()
    ori_labels = Channel.empty()
}

process Transform_T1_Labels {
    cpus 1
    memory '2 GB'

    input:
    set sid, file(anat), file(labels), file(mat), file(warp) from anat_for_transformation

    output:
    set sid, "${sid}__labels_warped_int16.nii.gz" into transformed_labels
    set sid, "${sid}__t1_warped.nii.gz" into transformed_anat

    script:
    """
    antsApplyTransforms -d 3 -i $anat -r $warp -o "${sid}__t1_warped.nii.gz" -t $warp $mat -n Linear
    antsApplyTransforms -d 3 -i $labels -r $warp -o labels_warped.nii.gz -t $warp $mat -n NearestNeighbor
    scil_volume_math.py convert labels_warped.nii.gz "${sid}__labels_warped_int16.nii.gz" --data_type int16
    """
}

ori_anat
    .concat(transformed_anat)
    .into{anat_for_registration;anat_for_concatenate;anat_for_metrics;anat_for_lesions}
ori_labels
    .concat(transformed_labels)
    .into{labels_for_transformation;labels_for_decompose}

in_tracking.set{tracking_for_decompose}

tracking_for_decompose
    .join(labels_for_decompose)
    .set{tracking_labels_for_decompose}

process Decompose_Connectivity {
    cpus 1
    memory { 7.B * trackings.size() }

    input:
    set sid, file(trackings), file(labels) from tracking_labels_for_decompose

    output:
    set sid, "${sid}__decompose.h5" into h5_for_commit, h5_for_skip_commit

    script:
    no_pruning_arg = ""
    if (params.no_pruning) {
        no_pruning_arg = "--no_pruning"
    }
    no_remove_loops_arg = ""
    if (params.no_remove_loops) {
        no_remove_loops_arg = "--no_remove_loops"
    }
    no_remove_outliers_arg = ""
    if (params.no_pruning) {
        no_remove_outliers_arg = "--no_pruning"
    }
    no_remove_outliers_arg = ""
    if (params.no_remove_outliers) {
        no_remove_outliers_arg = "--no_remove_outliers"
    }
    """
    scil_tractogram_segment_connections_from_labels.py $trackings $labels "${sid}__decompose.h5" --no_remove_curv_dev \
        $no_pruning_arg $no_remove_loops_arg $no_remove_outliers_arg --min_length $params.min_length \
        --max_length $params.max_length --loop_max_angle $params.loop_max_angle \
        --outlier_threshold $params.outlier_threshold
    """
}

data_for_commit
    .join(h5_for_commit)
    .set{data_tracking_for_commit}

process Run_COMMIT {
    cpus params.processes_commit
    memory params.commit_memory_limit

    input:
    set sid, file(bval), file(bvec), file(dwi), file(peaks), file(h5) from data_tracking_for_commit

    output:
    set sid, "${sid}__results_bzs/"
    set sid, "${sid}__decompose_commit.h5" into h5_for_afd_rd, h5_for_skip_afd_rd

    when:
    run_commit

    script:
    ball_stick_arg=""
    perp_diff_arg=""
    if (params.ball_stick) {
        ball_stick_arg="--ball_stick"
    }
    else {
        perp_diff_arg="--perp_diff $params.perp_diff"
    }
    if (params.use_commit2) {
    """
    scil_tractogram_commit.py $h5 $dwi $bval $bvec "${sid}__results_bzs/" --ball_stick --commit2 --in_peaks $peaks \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir \
        --para_diff $params.para_diff $perp_diff_arg --iso_diff $params.iso_diff
    mv "${sid}__results_bzs/commit_2/decompose_commit.h5" ./"${sid}__decompose_commit.h5"
    """
    }
    else {
    """
    scil_tractogram_commit.py $h5 $dwi $bval $bvec "${sid}__results_bzs/" --in_peaks $peaks \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir $ball_stick_arg \
        --para_diff $params.para_diff $perp_diff_arg --iso_diff $params.iso_diff
    mv "${sid}__results_bzs/commit_1/decompose_commit.h5" ./"${sid}__decompose_commit.h5"
    """
    }
}

if (!run_commit) {
    h5_for_skip_commit
        .into{h5_for_afd_rd;h5_for_skip_afd_rd}
}

h5_for_afd_rd
    .join(fodf_for_afd_rd)
    .set{h5_fodf_for_afd_rd}

process Compute_AFD_RD {
    cpus params.processes_afd_rd
    memory '2 GB'

    input:
    set sid, file(h5), file(fodf) from h5_fodf_for_afd_rd

    output:
    set sid, "${sid}__decompose_afd_rd.h5" into h5_for_transformation

    when:
    run_afd_rd

    script:
    length_weighting_arg = ""
    if (params.length_weighting) {
        length_weighting_arg = "--length_weighting"
    }
    """
    scil_bundle_mean_fixel_afd.py $h5 $fodf "${sid}__decompose_afd_rd.h5" $length_weighting_arg \
        --sh_basis $params.sh_basis --processes $params.processes_afd_rd
    """
}

anat_for_registration
    .combine(template_for_registration)
    .set{anats_for_registration}
process Register_Anat {
    cpus params.processes_register
    memory '2 GB'

    input:
    set sid, file(native_anat), file(template) from anats_for_registration

    output:
    set sid, "${sid}__output0GenericAffine.mat", "${sid}__output1Warp.nii.gz", "${sid}__output1InverseWarp.nii.gz"  into transformation_for_data, transformation_for_metrics, transformation_for_lesions
    file "${sid}__outputWarped.nii.gz"
    script:
    """
    antsRegistrationSyNQuick.sh -d 3 -m ${native_anat} -f ${template} -n ${params.processes_register} -o "${sid}__output" -t s
    """
}

in_opt_metrics
    .transpose()
    .concat(anat_for_metrics)
    .groupTuple()
    .flatMap{ sid, metrics -> metrics.collect{ [sid, it] } }
    .combine(transformation_for_metrics, by: 0)
    .combine(template_for_transformation_metrics)
    .set{metrics_transformation_for_metrics}
process Transform_Metrics {
    cpus 1
    memory '2 GB'

    input:
    set sid, file(metric), file(transfo), file(warp), file(inverse_warp), file(template) from metrics_transformation_for_metrics

    output:
    set sid, "*_mni.nii.gz" into metrics_for_compute

    script:
    """
    antsApplyTransforms -d 3 -i $metric -r $template -t $warp $transfo -o ${metric.getSimpleName()}_mni.nii.gz
    """
}

lesion_for_lesion_load
    .join(transformation_for_lesions)
    .combine(template_for_transformation_lesions)
    .set{lesions_transformation_for_data}

process Transform_Lesions {
    cpus 1

    input:
    set sid, file(lesion), file(transfo), file(warp), file(inverse_warp), file(template) from lesions_transformation_for_data

    output:
    set sid, "lesion_mask_mni.nii.gz" into lesions_for_compute

    script:
    """
    antsApplyTransforms -d 3 -i $lesion -r $template -t $warp $transfo -o ${lesion.getSimpleName()}_mni.nii.gz -n NearestNeighbor
    scil_volume_math.py convert lesion_mask_mni.nii.gz lesion_mask_mni.nii.gz --data_type uint8 -f
    """
}

if (!run_afd_rd) {
    h5_for_skip_afd_rd
        .set{h5_for_transformation}
}
h5_for_transformation
    .join(labels_for_transformation)
    .join(transformation_for_data)
    .combine(template_for_transformation_data)
    .set{labels_tracking_transformation_for_data}
process Transform_Data {
    cpus 1
    memory '2 GB'

    input:
    set sid, file(h5), file(labels), file(transfo), file(warp), file(inverse_warp), file(template) from labels_tracking_transformation_for_data

    output:
    set sid, "${sid}__decompose_warped_mni.h5", "${sid}__labels_warped_mni_int16.nii.gz" into h5_labels_for_compute
    file "${sid}__decompose_warped_mni.h5" into h5_for_similarity

    script:
    """
    scil_tractogram_apply_transform_to_hdf5.py $h5 $template ${transfo} "${sid}__decompose_warped_mni.h5" --inverse --in_deformation $inverse_warp
    antsApplyTransforms -d 3 -i $labels -r $template -t $warp $transfo -n NearestNeighbor -o labels_mni.nii.gz
    scil_volume_math.py convert labels_mni.nii.gz "${sid}__labels_warped_mni_int16.nii.gz" --data_type int16
    """
}


h5_for_similarity
    .take(params.nbr_subjects_for_avg_connections)
    .collect()
    .set{all_h5_for_similarity}

process Average_Connections {
    cpus params.processes_avg_similarity
    memory '2 GB'
    publishDir = "$params.avg_conn_output_dir"

    input:
    file(all_h5) from all_h5_for_similarity

    output:
    file "avg_per_edges/" into edges_for_similarity

    when:
    params.use_similarity_metric

    script:
    """
    scil_connectivity_hdf5_average_density_map.py $all_h5 avg_per_edges/ --binary --processes $params.processes_avg_similarity
    """
}

metrics_for_compute
    .concat(lesions_for_compute)
    .groupTuple()
    .set{all_metrics_for_compute}

if (params.use_similarity_metric) {
    h5_labels_for_compute
        .join(all_metrics_for_compute)
        .combine(edges_for_similarity)
        .combine(labels_list_for_compute)
        .set{h5_labels_similarity_list_for_compute}
    h5_labels_list_for_compute = Channel.empty()
}
else {
    h5_labels_for_compute
        .join(all_metrics_for_compute)
        .combine(labels_list_for_compute)
        .set{h5_labels_list_for_compute}
    h5_labels_similarity_list_for_compute = Channel.empty()
}

process Compute_Connectivity_with_similiarity {
    cpus params.processes_connectivity
    memory '2 GB'
    publishDir = {"${params.output_dir}/$sid/Compute_Connectivity"}

    input:
    set sid, file(h5), file(labels), file(metrics), file(avg_edges), file(labels_list) from h5_labels_similarity_list_for_compute

    output:
    set sid, "*.npy" into matrices_for_visualize_with_similarity, matrices_w_similarity_for_merge

    script:
    String metrics_list = metrics.join(", ").replace(',', '')
    """
    metrics_args=""
    lesion_args=""
    for metric in $metrics_list; do
        base_name=\$(basename \${metric/_mni/})
        base_name=\${base_name/_warped/}
        base_name=\${base_name/"${sid}__"/}
        if [[ \$metric == lesion_mask_mni.nii.gz ]]; then
            lesion_args="--lesion_load \$metric ./"
        else
            metrics_args="\${metrics_args} --metrics \${metric} \$(basename \$base_name .nii.gz).npy"
        fi
    done

    scil_connectivity_compute_matrices.py $h5 $labels --force_labels_list $labels_list \
        --volume vol.npy --streamline_count sc.npy \
        --length len.npy --similarity $avg_edges sim.npy \$metrics_args \
        --density_weighting --no_self_connection \
        --include_dps dps_matrices \$lesion_args --min_lesion_vol $params.min_lesion_vol \
        --processes $params.processes_connectivity

    rm rd_fixel.npy -f
    scil_connectivity_normalize.py sc.npy sc_edge_normalized.npy \
        --parcel_volume $labels $labels_list
    scil_connectivity_normalize.py vol.npy sc_vol_normalized.npy \
        --parcel_volume $labels $labels_list
    """
}

process Compute_Connectivity_without_similiarity {
    cpus params.processes_connectivity
    memory '2 GB'
    publishDir = {"${params.output_dir}/$sid/Compute_Connectivity"}

    input:
    set sid, file(h5), file(labels), file(metrics), file(labels_list) from h5_labels_list_for_compute

    output:
    set sid, "*.npy" into matrices_for_visualize_without_similarity, matrices_wo_similarity_for_merge

    script:
    String metrics_list = metrics.join(", ").replace(',', '')
    """
    metrics_args=""
    lesion_args=""
    for metric in $metrics_list; do
        base_name=\$(basename \${metric/_mni/})
        base_name=\${base_name/_warped/}
        base_name=\${base_name/"${sid}__"/}
        if [[ \$metric == lesion_mask_mni.nii.gz ]]; then
            lesion_args="--lesion_load \$metric ./"
        else
            metrics_args="\${metrics_args} --metrics \${metric} \$(basename \$base_name .nii.gz).npy"
        fi
    done

    scil_compute_connectivity.py $h5 $labels --force_labels_list $labels_list \
        --volume vol.npy --streamline_count sc.npy \
        --length len.npy \$metrics_args --density_weighting \
        --no_self_connection --include_dps dps_matrices \$lesion_args \
        --processes $params.processes_connectivity

    rm rd_fixel.npy -f
<<<<<<< HEAD
    scil_connectivity_normalize.py sc.npy sc_edge_normalized.npy \
        --parcel_volume $labels $labels_list
    scil_connectivity_normalize.py vol.npy sc_vol_normalized.npy \
=======
    scil_normalize_connectivity.py sc.npy sc_parcel_vol_normalized.npy \
>>>>>>> master
        --parcel_volume $labels $labels_list
    scil_normalize_connectivity.py sc.npy sc_bundle_vol_normalized.npy \
        --bundle_volume vol.npy
    """
}
 
matrices_w_similarity_for_merge
  .mix(matrices_wo_similarity_for_merge)
  .set{matrices_for_connectivity_in_csv}

process Connectivity_in_csv {
    cpus 1
    memory '2 GB'
    publishDir = {"${params.output_dir}/$sid/Compute_Connectivity"}

    input:
    set sid, file(matrices) from matrices_for_connectivity_in_csv

    output:
    file "*csv"

    script:
    String matrices_list = matrices.join("\",\"")
    """
    #!/usr/bin/env python3
    import numpy as np
    import os, sys

    for data in ["$matrices_list"]:
      fmt='%1.8f'
      if data == 'sc.npy':
        fmt='%i'

      curr_data = np.load(data)
      np.savetxt(data.replace(".npy", ".csv"), curr_data, delimiter=",", fmt=fmt)
    """
}

matrices_for_visualize_with_similarity
    .concat(matrices_for_visualize_without_similarity)
    .combine(labels_list_for_visualize)
    .set{matrices_labels_list_for_visualize}

process Visualize_Connectivity {
    cpus 1
    memory '2 GB'

    input:
    set sid, file(matrices), file(labels_list) from matrices_labels_list_for_visualize

    output:
    set sid, "*.png"

    script:
    String matrices_list = matrices.join(", ").replace(',', '')
    """
    for matrix in $matrices_list; do
        scil_viz_connectivity.py \$matrix \${matrix/.npy/_matrix.png} --labels_list $labels_list --name_axis \
            --display_legend --histogram \${matrix/.npy/_hist.png} --nb_bins 50 --exclude_zeros --axis_text_size 5 5
    done
    """
}
