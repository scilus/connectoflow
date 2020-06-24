#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["tracking_ext":"$params.tracking_ext",
                "length_weighting":"$params.length_weighting",
                "sh_basis":"$params.sh_basis",
                "run_commit":"$params.run_commit",
                "b_thr":"$params.b_thr",
                "nbr_dir":"$params.nbr_dir",
                "ball_stick":"$params.ball_stick",
                "para_diff":"$params.para_diff",
                "perp_diff":"$params.perp_diff",
                "iso_diff":"$params.iso_diff",
                "register_processes":"$params.register_processes",
                "processes_commit":"$params.processes_commit",
                "processes_afd_rd":"$params.processes_afd_rd",
                "processes_avg_similarity":"$params.processes_avg_similarity",
                "compute_connectivity":"$params.compute_connectivity"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "Run Connectivity Construction"
log.info "========================="
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
log.info "Root: $params.root"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info ""

root = file(params.root)
/* Watch out, files are ordered alphabetically in channel */
in_data = Channel
    .fromFilePairs("$root/**/{*fa.nii.gz,*labels.nii.gz,}",
                    size: 2,
                    maxDepth:1,
                    flat: true) {it.parent.name}

in_tracking = Channel
    .fromFilePairs("$root/**/{*tracking*.*,}",
                    size: -1,
                    maxDepth:1) {it.parent.name}

fodf_for_afd_rd = Channel
    .fromFilePairs("$root/**/*fodf.nii.gz",
                    size: 1,
                    maxDepth:1,
                    flat: true) {it.parent.name}

Channel.fromPath(file(params.template))
    .into{template_for_transformation;template_for_transformation_data;template_for_transformation_metrics}

Channel.fromPath(file(params.labels_list))
    .into{labels_list_for_compute;labels_list_for_visualize}

in_opt_metrics = Channel
    .fromFilePairs("$root/**/metrics/*.nii.gz",
                    size: -1,
                    maxDepth:2) {it.parent.parent.name}

in_opt_data = Channel
    .fromFilePairs("$root/**/{*dwi.bval,*dwi.bvec,*dwi.nii.gz,*peaks.nii.gz}",
                    size: 4,
                    maxDepth:1,
                    flat: true) {it.parent.name}


(anat_for_transformation, anat_for_concatenate, anat_for_ic, labels_for_transformation, labels_for_decompose) = in_data
    .map{sid, anat, labels ->
        [tuple(sid, anat),
        tuple(sid, anat),
        tuple(sid, anat),
        tuple(sid, labels),
        tuple(sid, labels)]}
    .separate(5)

(data_for_commit) = in_opt_data
    .map{sid, bval, bvec, dwi, peaks -> 
        [tuple(sid, bval, bvec, dwi, peaks)]}
    .separate(1)

in_tracking
    .join(anat_for_concatenate)
    .set{trackings_anat_for_concatenate}

process Concatenate_Tracking {
    cpus 1

    input:
    set sid, file(trackings), file(ref) from trackings_anat_for_concatenate

    output:
    set sid, "${sid}__tracking_union.trk" into tracking_for_ic

    script:
    """
    scil_streamlines_math.py concatenate $trackings ${sid}__tracking_union.trk --ignore_invalid --reference $ref
    """
}

tracking_for_ic
    .join(anat_for_ic)
    .set{tracking_anat_for_ic}

process Remove_IC {
    cpus 1

    input:
    set sid, file(tracking), file(ref) from tracking_anat_for_ic

    output:
    set sid, "${sid}__tracking_ic.trk" into tracking_for_commit, tracking_for_skip

    script:
    """
    scil_remove_invalid_streamlines.py $tracking "${sid}__tracking_ic.trk" --cut --remove_single --remove_overlapping
        --reference $ref
    """
}

data_for_commit
    .join(tracking_for_commit)
    .set{data_tracking_for_commit}

process Run_COMMIT {
    cpus params.processes_commit

    input:
    set sid, file(bval), file(bvec), file(dwi), file(peaks), file(tracking) from data_tracking_for_commit

    output:
    set sid, "${sid}__results_bzs/"
    set sid, "${sid}__essential_tractogram.trk" into tracking_for_decompose

    when:
    params.run_commit

    script:
    """
    ball_stick_arg=""
    perp_diff=""
    if $params.ball_stick; then
        ball_stick_arg="--ball_stick"
    else
        perp_diff="--perp_diff $params.perp_diff"
    fi
    scil_compute_streamlines_density_map.py $tracking tracking_mask.nii.gz --binary
    scil_run_commit.py $tracking $dwi $bval $bvec ${sid}__results_bzs/ --in_peaks $peaks --in_tracking_mask tracking_mask.nii.gz
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir \$ball_stick_arg
        --para_diff $params.para_diff \$perp_diff --iso_diff $params.iso_diff
    mv ${sid}__results_bzs/essential_tractogram.trk ./"${sid}__essential_tractogram.trk"
    """
}

if (!params.run_commit) {
    tracking_for_skip
        .set{tracking_for_decompose}
}

tracking_for_decompose
    .join(labels_for_decompose)
    .set{tracking_labels_for_decompose}

process Decompose_Connectivity {
    cpus 1

    input:
    set sid, file(tracking), file(labels) from tracking_labels_for_decompose

    output:
    set sid, "${sid}__decompose.h5" into h5_for_afd_rd, h5_for_skip

    script:
    """
    scil_decompose_connectivity.py $tracking $labels ${sid}__decompose.h5
    """
}

h5_for_afd_rd
    .join(fodf_for_afd_rd)
    .set{h5_fodf_for_afd_rd}

process Compute_AFD_RD {
    cpus params.processes_afd_rd

    input:
    set sid, file(h5), file(fodf) from h5_fodf_for_afd_rd

    output:
    set sid, "${sid}__decompose_afd_rd.h5" into h5_for_transformation

    when:
    params.run_afd_rd

    script:
    """
    length_weighting_arg=""
    if $params.length_weighting; then
        length_weighting_arg="--length_weighting"
    fi
    scil_compute_mean_fixel_afd_from_hdf5.py $h5 $fodf "${sid}__decompose_afd_rd.h5" \$length_weighting_arg
        --sh_basis $params.sh_basis --processes $params.processes_afd_rd
    """
}

anat_for_transformation
    .combine(template_for_transformation)
    .set{anats_for_transformation}
process Register_Anat {
    cpus params.register_processes
    input:
    set sid, file(native_anat), file(template) from anats_for_transformation

    output:
    set sid, "${sid}__output0GenericAffine.mat", "${sid}__output1Warp.nii.gz", "${sid}__output1InverseWarp.nii.gz"  into transformation_for_data, transformation_for_metrics
    file "${sid}__outputWarped.nii.gz"
    script:
    """
    antsRegistrationSyNQuick.sh -d 3 -m ${native_anat} -f ${template} -n ${params.register_processes} -o ${sid}__output -t s
    """ 
}

in_opt_metrics
    .flatMap{ sid, metrics -> metrics.collect{ [sid, it] } }
    .combine(transformation_for_metrics, by: 0)
    .combine(template_for_transformation_metrics)
    .set{metrics_transformation_for_metrics}
process Transform_Metrics {
    cpus 1

    input:
    set sid, file(metric), file(transfo), file(warp), file(inverse_warp), file(template) from metrics_transformation_for_metrics

    output:
    set sid, "*_warped.nii.gz" into metrics_for_compute

    script:
    """
    antsApplyTransforms -d 3 -i $metric -r $template -t $warp $transfo -o ${metric.getSimpleName()}_warped.nii.gz
    """
}

if (!params.run_afd_rd) {
    h5_for_skip
        .set{h5_for_transformation}
}

h5_for_transformation
    .join(labels_for_transformation)
    .join(transformation_for_data)
    .combine(template_for_transformation_data)
    .set{labels_tracking_transformation_for_data}
process Transform_Data {
    input:
    set sid, file(h5), file(labels), file(transfo), file(warp), file(inverse_warp), file(template) from labels_tracking_transformation_for_data

    output:
    set sid, "${sid}__decompose_warped.h5", "${sid}__labels_warped_int16.nii.gz" into h5_labels_for_compute,lol
    file "${sid}__decompose_warped.h5" into h5_for_similarity

    script:
    """
    scil_apply_transform_to_hdf5.py $h5 $template ${transfo} "${sid}__decompose_warped.h5" --inverse --in_deformation $inverse_warp --cut_invalid
    antsApplyTransforms -d 3 -i $labels -r $template -t $warp $transfo -n NearestNeighbor -o ${sid}__labels_warped.nii.gz
    scil_image_math.py convert ${sid}__labels_warped.nii.gz ${sid}__labels_warped_int16.nii.gz --data_type int16
    """
}


h5_for_similarity
    .collect()
    .set{all_h5_for_similarity}

process Average_Connections {
    cpus params.processes_avg_similarity
    publishDir = {"./results_conn/$task.process"}

    input:
    file(all_h5) from all_h5_for_similarity

    output:
    file "avg_per_edges/" into edges_for_similarity

    when:
    params.use_similarity_metric

    script:
    """
    scil_compute_hdf5_average_density_map.py $all_h5 avg_per_edges/ --binary --processes $params.processes_avg_similarity
    """
}

metrics_for_compute
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
    cpus params.compute_connectivity
    publishDir = {"./results_conn/$sid/Compute_Connectivity"}

    input:
    set sid, file(h5), file(labels), file(metrics), file(avg_edges), file(labels_list) from h5_labels_similarity_list_for_compute

    output:
    set sid, "*.npy" into matrices_for_visualize_with_similarity

    script:
    String metrics_list = metrics.join(", ").replace(',', '')
    """
    metrics_args=""
    for metric in $metrics_list; do
        metrics_args="\${metrics_args} --metrics \${metric} \$(basename \${metric/_warped/} .nii.gz).npy" 
    done

    scil_compute_connectivity.py $h5 $labels --force_labels_list $labels_list --volume vol.npy --streamline_count sc.npy
        --length len.npy --similarity $avg_edges sim.npy \$metrics_args --density_weighting --no_self_connection
        --include_dps ./ --processes $params.compute_connectivity
    scil_normalize_connectivity.py sc.npy sc_edge_normalized.npy --parcel_volume $labels $labels_list
    scil_normalize_connectivity.py vol.npy sc_vol_normalized.npy --parcel_volume $labels $labels_list
    """
}

process Compute_Connectivity_without_similiarity {
    cpus params.compute_connectivity
    publishDir = {"./results_conn/$sid/Compute_Connectivity"}

    input:
    set sid, file(h5), file(labels), file(metrics), file(labels_list) from h5_labels_list_for_compute

    output:
    set sid, "*.npy" into matrices_for_visualize_without_similarity

    script:
    String metrics_list = metrics.join(", ").replace(',', '')
    """
    metrics_args=""
    for metric in $metrics_list; do
        metrics_args="\${metrics_args} --metrics \${metric} \$(basename \${metric/_warped/} .nii.gz).npy" 
    done

    scil_compute_connectivity.py $h5 $labels --force_labels_list $labels_list --volume vol.npy --streamline_count sc.npy
        --length len.npy \$metrics_args --density_weighting --no_self_connection --include_dps ./
        --processes $params.compute_connectivity
    scil_normalize_connectivity.py sc.npy sc_edge_normalized.npy --parcel_volume $labels $labels_list
    scil_normalize_connectivity.py vol.npy sc_vol_normalized.npy --parcel_volume $labels $labels_list
    """
}

matrices_for_visualize_with_similarity
    .concat(matrices_for_visualize_without_similarity)
    .combine(labels_list_for_visualize)
    .set{matrices_labels_list_for_visualize}

process Visualize_Connectivity {
    input:
    set sid, file(matrices), file(labels_list) from matrices_labels_list_for_visualize

    output:
    set sid, "*.png"

    script:
    String matrices_list = matrices.join(", ").replace(',', '')
    """
    for matrix in $matrices_list; do
        scil_visualize_connectivity.py \$matrix \${matrix/.npy/_matrix.png} --labels_list $labels_list --name_axis
            --display_legend --histogram \${matrix/.npy/_hist.png} --nb_bins 50 --exclude_zeros --axis_text_size 5 5
    done
    """
}
