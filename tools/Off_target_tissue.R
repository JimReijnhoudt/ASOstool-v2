library(jsonlite)

parse_PAtlas <- function(ENSMBL_id, target_tissue){
  tissue_types <- list(
    "Brain" = c("brain_RNA_amygdala",
                "brain_RNA_basal_ganglia", 
                "brain_RNA_cerebellum",
                "brain_RNA_cerebral_cortex", 
                "brain_RNA_choroid_plexus",
                "brain_RNA_hippocampal_formation",
                "brain_RNA_hypothalamus",
                "brain_RNA_medulla_oblongata", 
                "brain_RNA_midbrain",
                "brain_RNA_pons", 
                "brain_RNA_spinal_cord",
                "brain_RNA_thalamus"),
    "Eye" = c("t_RNA_retina"), 
    "Endrocrine tissue" = c("t_RNA_thyroid_gland",
                            "t_RNA_parathyroid_gland", 
                            "t_RNA_adrenal_gland",
                            "t_RNA_pituitary_gland"), 
    "Respiratory system" = c("t_RNA_lung"),
    "Proximal digestive tract" = c("t_RNA_salivary_gland",
                                   "t_RNA_esophagus",
                                   "t_RNA_tongue"), 
    "Gastrointestinal tract" = c("t_RNA_stomach_1",
                                 "t_RNA_duodenum", 
                                 "t_RNA_small_intestine",
                                 "t_RNA_colon",
                                 "t_RNA_rectum"),
    "Liver & galbladder" = c("t_RNA_gallbladder"), 
    "Pancreas" = c("t_RNA_pancreas"), 
    "Kidney & urinary bladder" = c("t_RNA_kidney",
                                   "t_RNA_urinary_bladder"),
    "Male tissues" = c("t_RNA_testis",
                       "t_RNA_epididymis",
                       "t_RNA_seminal_vesicle",
                       "t_RNA_prostate"), 
    "Female tissues" = c("t_RNA_vagina",
                         "t_RNA_ovary",
                         "t_RNA_fallopian_tube",
                         "t_RNA_endometrium_1",
                         "t_RNA_cervix",
                         "t_RNA_breast"), 
    "Muscle tissues" = c("t_RNA_heart_muscle",
                         "t_RNA_smooth_muscle",
                         "t_RNA_skeletal_muscle"), 
    "Connective & soft tissues" = c("t_RNA_adipose_tissue"), 
    "Skin" = c("t_RNA_skin_1"), 
    "Bone marrow & lymphoid tissues" = c("t_RNA_appendix",
                                         "t_RNA_spleen",
                                         "t_RNA_lymph_node",
                                         "t_RNA_tonsil",
                                         "t_RNA_bone_marrow",
                                         "t_RNA_thymus")
                       )

  api_columns <- "g,eg,t_RNA__tau,t_RNA_liver,"
  api_options <- "compress=no&format=json"
  base_url <- "https://www.proteinatlas.org/api/search_download.php?search="
  tissue_columns <- paste(tissue_types[[target_tissue]], collapse=",")

  full_url <- paste0(base_url, ENSMBL_id, "&columns=", api_columns, tissue_columns, "&", api_options)
  res <- fromJSON(full_url)
  res <- res %>%
    mutate(across(4:ncol(.), ~ as.numeric(as.character(.)))) %>%
    rowwise() %>%
    mutate(total_expression = sum(c_across(4:ncol(.)), na.rm = TRUE)) %>%
    ungroup()
  return(res)
}
