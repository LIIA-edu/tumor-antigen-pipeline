
nextflow.enable.dsl=2

include { ANNDATA_OVERVIEW } from './modules/anndata_overview'
include { DATA_CATEGORIZATION } from './modules/data_categorization'
include { SPECIFICITY } from './modules/specificity'
include { PLOT_SPECIFICITY } from './modules/plot_specificity'
include { SET_ANTIGENS_TYPE } from './modules/set_antigens_type'
include { PLOT_ANTIGENS_TYPE } from './modules/plot_antigens_type'
include { NON_MALIGNANTS_HEALTHY_TISSUES } from './modules/non_malignants_healthy_tissues'
include { CANDIDATES_PER_TUMORS } from './modules/candidates_per_tumors'
include { PLOT_CANDIDATES_PER_TUMORS } from './modules/plot_candidates_per_tumors'
include { HETEROGENEITY_CANDIDATES } from './modules/heterogeneity_candidates'
include { PLOT_HETEROGENEITY_CANDIDATES } from './modules/plot_heterogeneity_candidates'
include { HPA_TISSUE_SPECIFICITY } from './modules/hpa_tissue_specificity'
include { ENRICHMENT_MODULES } from './modules/enrichment_modules'
include { PRIORITIZE_CANDIDATES } from './modules/prioritize_candidates'
include { PLOT_EPITOPES } from './modules/plot_epitopes'


workflow {

    // ========================================================
    // INPUT
    // ========================================================
    
    input_ch = Channel.fromPath(
        params.input,
        checkIfExists: true
    )
    // ========================================================
    // DATA CATEGORIZATION
    // ========================================================

    DATA_CATEGORIZATION(
        input_ch
    )

    // ========================================================
    // ANNDATA OVERVIEW
    // ========================================================

    ANNDATA_OVERVIEW(input_ch)

    // ========================================================
    // SPECIFICITY
    // ========================================================

    if (params.specificity == 'hpa') {

        SPECIFICITY(input_ch)

        PLOT_SPECIFICITY(
            SPECIFICITY.out.summary
        )

    } else {
        error "Unsupported specificity method: ${params.specificity}"
    }

    // ========================================================
    // HPA TISSUE SPECIFICITY
    // ========================================================

    HPA_TISSUE_SPECIFICITY(
        SPECIFICITY.out.enriched
    )

    // ========================================================
    // ANTIGEN TYPE
    // ========================================================

    SET_ANTIGENS_TYPE(
        SPECIFICITY.out.classified
    )

    PLOT_ANTIGENS_TYPE(
        SET_ANTIGENS_TYPE.out.classified_labeled
    )

    // ========================================================
    // NON-MALIGNANT CELLS / HEALTHY TISSUES
    // ========================================================

    NON_MALIGNANTS_HEALTHY_TISSUES(
        SET_ANTIGENS_TYPE.out.classified_labeled
    )

    // ========================================================
    // PAN-CANCER ANALYSIS
    // ========================================================

    if (params.pan_cancer) {

        // ----------------------------------------------------
        // VALIDATE PARAMETERS
        // ----------------------------------------------------

        if (!params.epitopes) {
            error "The --epitopes parameter is required when --pan_cancer is enabled."
        }

        if (!params.gene_name_mapping) {
            error "The --gene_name_mapping parameter is required when --pan_cancer is enabled."
        }

        // ----------------------------------------------------
        // EPITOPE INPUT
        // ----------------------------------------------------

        epitopes_ch = Channel.fromPath(
            params.epitopes,
            checkIfExists: true
        )

        // ----------------------------------------------------
        // GENE NAME MAPPING
        // ----------------------------------------------------

        gene_name_mapping_ch = Channel.value(
            file(params.gene_name_mapping, checkIfExists: true)
        )

        // ----------------------------------------------------
        // CANDIDATES PER TUMOR
        // ----------------------------------------------------

        CANDIDATES_PER_TUMORS(
            input_ch,
            gene_name_mapping_ch,
            SPECIFICITY.out.classified
        )

        // ----------------------------------------------------
        // PLOT CANDIDATES PER TUMOR
        // ----------------------------------------------------

        PLOT_CANDIDATES_PER_TUMORS(
            CANDIDATES_PER_TUMORS.out.binary_matrix,
            epitopes_ch
        )

        // ----------------------------------------------------
        // PLOT EPITOPES
        // ----------------------------------------------------

        PLOT_EPITOPES(
            epitopes_ch,
            SPECIFICITY.out.enriched
        )

        // ----------------------------------------------------
        // HETEROGENEITY CANDIDATES
        // ----------------------------------------------------

        heterogeneity_input_ch = Channel.fromPath(
            params.heterogeneity_input ?: params.input,
            checkIfExists: true
        )

        HETEROGENEITY_CANDIDATES(
            heterogeneity_input_ch,
            SPECIFICITY.out.classified,
            gene_name_mapping_ch
        )

        // ----------------------------------------------------
        // PLOT HETEROGENEITY CANDIDATES
        // ----------------------------------------------------
        
        PLOT_HETEROGENEITY_CANDIDATES(
            HETEROGENEITY_CANDIDATES.out.modules_by_disease,
            HETEROGENEITY_CANDIDATES.out.mean_expression_by_disease,
            HETEROGENEITY_CANDIDATES.out.pct_expression_by_disease,
            epitopes_ch
        )
        // ----------------------------------------------------
        // FUNCTIONAL ENRICHMENT
        // ----------------------------------------------------
        
        ENRICHMENT_MODULES(
            HETEROGENEITY_CANDIDATES.out.modules_by_disease
        )

        // ----------------------------------------------------
        // PRIORITIZE CANDIDATES
        // ----------------------------------------------------

        PRIORITIZE_CANDIDATES(
            SPECIFICITY.out.enriched,
            HPA_TISSUE_SPECIFICITY.out.tissue_specificity,
            epitopes_ch,
            CANDIDATES_PER_TUMORS.out.binary_matrix,
            SET_ANTIGENS_TYPE.out.merged_summed_genes
        )

    }
}
