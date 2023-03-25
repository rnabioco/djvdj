
# # LINEAGE TREES
# # * how to get the correct alignments for sequence/germline inputs?
# # * with unaligned sequences, get warning from formatClones() and error from
# #   getTress()
# # * need to mask D segment in germline sequence
#
# # Parse bam file
# dat_dir  <- "inst/extdata/splen/BL6_BCR/"
# bam_file <- here::here(dat_dir, "concat_ref.bam")
# fa_file  <- here::here(dat_dir, "concat_ref.fasta")
#
# bam_info <- Rsamtools::scanBam(bam_file)[[1]]
#
# bam_info <- tibble::tibble(
#   contig_id = bam_info$qname,
#   ref       = bam_info$rname,
#   cigar     = bam_info$cigar,
#   seq       = as.character(bam_info$seq)
# )
#
# bam_info <- bam_info %>%
#   dplyr::filter(grepl("_consensus_[0-9]+$", contig_id)) %>%
#   mutate(
#     clonotype_id = str_extract(ref, "^clonotype[0-9]+")
#   )
#
# # Parse germline ref fasta
# fa <- Biostrings::readDNAStringSet(fa_file)
#
# fa_info <- as.character(fa)
#
# fa_info <- tibble::tibble(
#   ref = names(fa_info),
#   germline = unname(fa_info)
# )
#
# # Parse other info
# vdj <- import_vdj(
#   vdj_dir = dat_dir,
#   data_cols = "raw_consensus_id",
#   filter_paired = TRUE
# ) %>%
#   fetch_vdj() %>%
#   filter(chains == "IGH", n_chains == 2) %>%
#   distinct(
#     clonotype_id, exact_subclonotype_id, raw_consensus_id,
#     v_gene, j_gene, cdr3_nt_length, chains
#   )
#
# # Merge
# vdj <- bam_info %>%
#   left_join(fa_info, by = "ref") %>%
#   inner_join(vdj, by = c("clonotype_id", contig_id = "raw_consensus_id"))
#
# # Format airrClone object for lineage construction
# airr <- vdj %>%
#   dplyr::rename(
#     sequence_alignment        = "seq",
#     germline_alignment_d_mask = "germline",
#     v_call                    = "v_gene",
#     j_call                    = "j_gene",
#     junction_length           = "cdr3_nt_length",
#     clone_id                  = "clonotype_id",
#     sequence_id               = "contig_id",
#     minseq                    = 1
#   ) %>%
#   dowser::formatClones(
#     subclone = "exact_subclonotype_id",
#     minseq = 1
#   )
#
# airr <- dowser::getTrees(airr, build = "pml")
