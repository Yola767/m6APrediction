#' @import randomForest
NULL

#' Convert DNA 5-mer strings into a one-hot encoded data frame
#'
#' This internal helper takes a character vector of DNA 5-mers (all of equal length)
#' and returns a data frame where each nucleotide position is a factor column
#' with levels `"A"`, `"T"`, `"C"`, `"G"`.  The resulting columns are named
#' `nt_pos1` … `nt_pos5`.
#'
#' @param dna_strings Character vector containing DNA 5-mer sequences.
#'
#' @return A data frame with 5 factor columns (`nt_pos1` to `nt_pos5`).
#'   Each column has levels `c("A", "T", "C", "G")`.
#'
#' @keywords internal
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Predict m6A modification for multiple sequences
#'
#' Given a fitted randomForest model and a data frame of features,
#' this function returns the original features together with two new columns:
#' `predicted_m6A_prob` (probability of being Positive) and
#' `predicted_m6A_status` (`"Positive"` or `"Negative"`).
#'
#' @param ml_fit A fitted `randomForest` model that was trained with a
#'   factor response having levels `"Negative"` and `"Positive"`.
#' @param feature_df Data frame containing the required feature columns
#'   (see **Details**).  The column `DNA_5mer` must hold character strings
#'   of length 5.
#' @param positive_threshold Numeric scalar (default `0.5`).  Probability
#'   above this threshold is classified as `"Positive"`.
#'
#' @details The input `feature_df` **must** contain the following columns
#'   (exact names):
#'   \itemize{
#'     \item `gc_content`
#'     \item `RNA_type`      (factor with levels `mRNA`, `lincRNA`,
#'                           `lncRNA`, `pseudogene`)
#'     \item `RNA_region`    (factor with levels `CDS`, `intron`,
#'                           `3'UTR`, `5'UTR`)
#'     \item `exon_length`
#'     \item `distance_to_junction`
#'     \item `evolutionary_conservation`
#'     \item `DNA_5mer`      (character, length 5)
#'   }
#'
#' @return A data frame with the original columns (except `DNA_5mer`,
#'   which is replaced by one-hot columns `nt_pos1` … `nt_pos5`) plus
#'   two new columns:
#'   \describe{
#'     \item{predicted_m6A_prob}{Rounded probability (4 decimals) of class
#'                               `"Positive"`}
#'     \item{predicted_m6A_status}{Factor with levels `"Negative"`,
#'                                 `"Positive"`}
#'   }
#'   The column order is fixed as shown in the source code.
#'   #'
#' @examples
#' # 加载包内示例数据
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' input_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#'
#' # 批量预测
#' predictions <- prediction_multiple(rf_model, input_df)
#' head(predictions)
#'
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region",
                  "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  #new
  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  dna_onehot <- dna_encoding(feature_df$DNA_5mer)
  feature_df <- cbind(feature_df[ , !colnames(feature_df) %in% "DNA_5mer"], dna_onehot)
  prob_pos <- predict(ml_fit, newdata = feature_df, type = "prob")[,"Positive"]
  status <- ifelse(prob_pos > positive_threshold, "Positive", "Negative")
  status <- factor(status, levels = c("Negative", "Positive"))
  feature_df$predicted_m6A_prob <- round(prob_pos, 4)
  feature_df$predicted_m6A_status <- status
  feature_df <- feature_df[ , c("gc_content", "RNA_type", "RNA_region",
                                "exon_length", "distance_to_junction",
                                "evolutionary_conservation",
                                "nt_pos1","nt_pos2","nt_pos3","nt_pos4","nt_pos5",
                                "predicted_m6A_prob", "predicted_m6A_status")]
  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}


#' Predict m6A modification for a single sequence
#'
#' Convenience wrapper around [prediction_multiple()] for a single set
#' of feature values supplied as separate arguments.
#'
#' @param ml_fit Fitted `randomForest` model (same as in
#'   [prediction_multiple()]).
#' @param gc_content Numeric. GC content of the region.
#' @param RNA_type Character (will be converted to factor). One of
#'   `"mRNA"`, `"lincRNA"`, `"lncRNA"`, `"pseudogene"`.
#' @param RNA_region Character (will be converted to factor). One of
#'   `"CDS"`, `"intron"`, `"3'UTR"`, `"5'UTR"`.
#' @param exon_length Numeric. Length of the exon.
#' @param distance_to_junction Numeric. Distance to nearest splice junction.
#' @param evolutionary_conservation Numeric. Conservation score.
#' @param DNA_5mer Character. 5-mer DNA sequence centered on the site.
#' @param positive_threshold Numeric (default `0.5`). Classification threshold.
#'
#' @return A **named vector** with two elements:
#'   \item{predicted_m6A_prob}{Rounded probability of `"Positive"` (4 decimals)}
#'   \item{predicted_m6A_status}{Character `"Positive"` or `"Negative"`}
#'
#' @examples
#' # 加载模型
#' rf_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'
#' # 单个位点预测
#' prediction_single(
#'   ml_fit = rf_model,
#'   gc_content = 0.48,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 120,
#'   distance_to_junction = 30,
#'   evolutionary_conservation = 0.92,
#'   DNA_5mer = "GACTA"
#' )
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5){
  one_row <- data.frame(
    gc_content = gc_content, RNA_type = RNA_type, RNA_region = RNA_region,
    exon_length = exon_length, distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer, stringsAsFactors = FALSE)
  one_row$RNA_type <- factor(one_row$RNA_type,
                             levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  one_row$RNA_region <- factor(one_row$RNA_region,
                               levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  res_df <- prediction_multiple(ml_fit, one_row, positive_threshold)
  returned_vector <- c(
    predicted_m6A_prob = res_df$predicted_m6A_prob[1],
    predicted_m6A_status = as.character(res_df$predicted_m6A_status[1])
  )
  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}
