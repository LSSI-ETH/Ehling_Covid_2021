###
# Helper functions for data loading, clonotyping and extraction of most expanded sequences
###

.load_datasets <- function(file, sep, dataset_id, short=TRUE, filter=TRUE){

  dataset <- read.delim(file, sep = sep, stringsAsFactors = FALSE)

  # Remove allele information

  dataset$v_call <- gsub(as.character(dataset$v_call),
                         pattern = '\\*\\d\\d', replacement = '')
  dataset$v_call <- unlist(lapply(strsplit(dataset$v_call, ','), function(x) x[1]))

  dataset$j_call <- gsub(as.character(dataset$j_call),
                         pattern = '\\*\\d\\d', replacement = '')
  dataset$j_call <- unlist(lapply(strsplit(dataset$j_call, ','), function(x) x[1]))

  # Only keep productive sequences

  if(filter){
    dataset <- dataset[dataset$productive == 'True', ]
    dataset <- dataset[dataset$is_kept, ]
  }

  # Only keep CDR3 and V/J genes

  if(short){
    keep_cols <- c('junction_aa', 'j_call', 'v_call')
    dataset <- dataset[, names(dataset) %in% keep_cols]
  }

  # Compute CDR3 length and drop short sequences (below 3)

  dataset$cdr3_len <- nchar(as.character(dataset$junction_aa))
  dataset <- dataset[dataset$cdr3_len > 3, ]

  # Duplicate rows for multiple V/J gene assignments
  dataset <- distinct(dataset)
  dataset$dataset_id <- dataset_id

  return(dataset)
}

# Function to compute clonotypes
.clonotyping <- function(data, linkage, threshold){

  data <- separate_rows(data, 'v_call', sep = ',')
  data <- separate_rows(data, 'j_call', sep = ',')

  idx_list <- split(seq_len(nrow(data)),
                    list(data$v_call,data$j_call, data$cdr3_len))

  idx_list <- idx_list[lapply(idx_list, length) > 0]

  message('Start building clonotype groups...')

  cltyps <- 0

  msg <- paste('0 /',length(idx_list),'combinations scanned -', cltyps,
               'clonotypes found', sep = ' ')

  message(msg, "\r", appendLF=FALSE)
  flush.console()

  clust_list <- vector(mode='list', length = length(idx_list))
  dist_list <- vector(mode='list', length = length(idx_list))

  for(i in 1:length(idx_list)){

    if(length(data[idx_list[[i]], ]$junction_aa) > 1){
      dist_mat <- stringdistmatrix(a = as.character(data[idx_list[[i]], ]$junction_aa),
                                   b = as.character(data[idx_list[[i]], ]$junction_aa),
                                   method = 'hamming')
      dists <- apply(dist_mat, 1, function(x) min(x[x > 0]))
      clust <- hclust(as.dist(dist_mat), method = linkage)
      threshold_group <- data[idx_list[[i]], ]$cdr3_len[1]*threshold
      clust <- cutree(clust, h = threshold_group)
    } else {
      clust <- 1
      dists <- data[idx_list[[i]], ]$cdr3_len[1]
    }

    dist_list[[i]] <- dists
    clust_list[[i]] <- clust + cltyps
    cltyps <- max(clust_list[[i]])
    msg <- paste(i, '/',length(idx_list),'combinations scanned -', cltyps,
                 'clonotypes found', sep = ' ')
    message(msg, "\r", appendLF=FALSE)
    flush.console()

  }

  data$clonotype <- 0
  data$nn <- 0

  data[unlist(idx_list), ]$clonotype <- unlist(clust_list)
  data[unlist(idx_list), ]$nn <- unlist(dist_list)

  return(data)
}
#nn how far next neighbor: inf: clonotype is empty (nur eine sequenz) --> no neighbor --> set inf.
#clonotype: clusternummer vom clonotpying


.combine_chains <- function(x){
  ids <- unique(x$cell_id)

  x$combined_chain <- NA_character_

  for(idx in c(1:length(ids))){

    meta_x <- x[x$cell_id%in%ids[idx],]
    meta_x_ighv <- unique(meta_x$sequence_vdj[meta_x$chain == 'IGH'])
    meta_x_ighl <- unique(meta_x$sequence_vdj[meta_x$chain%in%c('IGL','IGK')])

    if((length(meta_x_ighl) == 1) & (length(meta_x_ighv) == 1)){

      x[x$cell_id%in%ids[idx], ]$combined_chain <- paste(meta_x_ighv, meta_x_ighl, sep='---')

    }


  }

  return(x)
}


#based on 10x cell umis count nb of cells
.count_cells <- function(x){

  seq_chain <- unique(x$combined_chain)
  x$cells_with_combined_chain <- 0

  for(idx in c(1:length(seq_chain))){
    x$cells_with_combined_chain[x$combined_chain==seq_chain[idx]] <- length(unique(x$cell_id[x$combined_chain==seq_chain[idx]]))

  }
  return(x)
}
#cellrangeer check consensus count.. different compared to bulk therefore focus on nb of cells with combined chain..


.get_most_expanded_IgGclone <- function(x){

  meta_x <- x[grepl(x$c_call, pattern = 'IGHG'), ]
  max_seq <- meta_x$combined_chain[which.max(meta_x$cells_with_combined_chain)]
  idx_igh <- which((x$combined_chain == max_seq) & (x$chain == 'IGH'))[1]
  idx_igl <- which((x$combined_chain == max_seq) & (x$chain%in%c('IGL', 'IGK')))[1]

  #return(x[c(idx_igh, idx_igl), ])
  return(x[c(idx_igh), ])

}

.match_clonotypes <- function(x, clonotyped_data){

  x$clonotype_large <- 0

  for(idx in c(1:length(x$cell_id))){

    #for current cell_id in single dataset find mclonotype of same cell_id in global clonotype data
    #in case cell id was kicked out (due to not paired, ie only light chain was present) this value is set to NA
    x$clonotype_large[idx] <- clonotyped_data$clonotype[which(clonotyped_data$cell_id == x$cell_id[idx])][1]

  }
  return(x)
}


.select_candidates <- function(PATH="igblast_results/"){


  dataset_paths <- list.files(PATH, full.names = TRUE, recursive = TRUE, pattern = 'combined')
  names(dataset_paths) <- unlist(lapply(strsplit(dataset_paths, split = '[/-]'),
                                        function(x) grep(pattern = 'combined', x=x,
                                                         value = TRUE)))

  datasets <- list()
  for(i in c(1:length(dataset_paths))){

    datasets[[i]] <- .load_datasets(dataset_paths[i], sep = '\t', dataset = names(dataset_paths)[i], short = FALSE)
    names(datasets)[i] <- names(dataset_paths)[i]

  }

  datasets <- lapply(datasets, .combine_chains)


  datasets_paired <- lapply(datasets, function(x) x[!is.na(x = x$combined_chain), ])
  datasets_paired <- lapply(datasets_paired, .count_cells)
  datasets_paired <- lapply(datasets_paired, .clonotyping, linkage='single', threshold = 0.2)


  #row bind into df that contains all patients and cleanup and colmerge H/L
  datasets_paired <- bind_rows(datasets_paired)
  datasets_paired <- datasets_paired[!is.na(datasets_paired$sequence_id), ]
  mask <- datasets_paired$chain=='IGH'
  datasets_paired <- merge(datasets_paired[mask, ], datasets_paired[!mask, ],
                           by = c('dataset_id', 'cell_id', 'combined_chain'))

  names(datasets_paired) <- gsub(names(datasets_paired), pattern = 'clonotype.x', replacement = 'clonotype.heavy')
  names(datasets_paired) <- gsub(names(datasets_paired), pattern = '\\.x', replacement = '')
  names(datasets_paired) <- gsub(names(datasets_paired), pattern = '\\.y', replacement = '.light')

  #clonotyping across all datasets.
  datasets_paired <- .clonotyping(datasets_paired, 'single', 0.2)

  datasets_out<-list()
  for(i in 1:length(unique(datasets_paired$dataset_id))){
    curr_id<-unique(datasets_paired$dataset_id)[i]
    curr_df<-datasets_paired[datasets_paired$dataset_id==curr_id,]
    curr_df_igg<-curr_df[curr_df$c_call %in% c("IGHG","IGHG1","IGHG2","IGHG3","IGHG4"),]
    cat(max(curr_df_igg$cells_with_combined_chain),"\n")
    datasets_out[[i]]<-.get_most_expanded_IgGclone(curr_df)
  }

  datasets_out <- bind_rows(datasets_out)


  TopClone_IgG_PerPatient <-datasets_out






  # Clonotype

  #bind rows before count cells (after combination but before non-HL-matched excluded)
  clonotyped_data <- bind_rows(datasets)
  clonotyped_data <- .clonotyping(clonotyped_data, 'single', 0.2)

  #calculate percent similarity (for ones without neighbor set to 100%
  clonotyped_data$nn[is.infinite(clonotyped_data$nn)] <- nchar(as.character(clonotyped_data$junction_aa[is.infinite(clonotyped_data$nn)]))
  clonotyped_data$nn_norm <- clonotyped_data$nn/nchar(as.character(clonotyped_data$junction_aa))
  #exclude light chain clonotypes.
  clonotyped_data <- clonotyped_data[!grepl(clonotyped_data$v_call, pattern = 'IGKV|IGHL'), ]

  # go back. do clonotyping per dataset (begfore bind row,)
  datasets <- lapply(datasets, .clonotyping, linkage='single', threshold = 0.2)
  #match within dataset clonotype to across datasets clonotype
  datasets <- lapply(datasets, .match_clonotypes, clonotyped_data=clonotyped_data)


  #combine and cleanup for only paired IGHG..
  datasets_combined <- bind_rows(datasets)

  datasets_combined$clonotype_large[!datasets_combined$chain=='IGH'] <- NA
  datasets_combined$clonotype[!datasets_combined$chain=='IGH'] <- NA
  datasets_combined$igg_clonotype <- datasets_combined$clonotype
  datasets_combined$igg_clonotype[!grepl(datasets_combined$c_call, pattern = 'IGHG')] <- NA
  datasets_combined$igg_clonotype[is.na(datasets_combined$combined_chain)] <- NA

  #data summary after cleanup
  dataset_summary <- datasets_combined %>% group_by(dataset_id) %>% summarise('# Cells' = n_distinct(cell_id), '# Unique Clones (At least one chain)' = n_distinct(raw_clonotype_id), '# Unique paired clones' = n_distinct(combined_chain), '#Unique Heavy chain clonotypes' = n_distinct(clonotype), '# Unique paired IgG clonotypes' = n_distinct(igg_clonotype))
  dataset_summary <- melt(as.data.table(dataset_summary), id.vars = 'dataset_id', value.name = 'Count', variable.name = 'Subset')
  dataset_summary$dataset_id <- gsub(dataset_summary$dataset_id, pattern = '_combined_db', replacement = '')
  labels_summary <- c(paste(rep('G'), c(1:8)), paste(rep('H'), c(1:8)))



  ## Chain plots


  #group by dataset_id and clonotype
  #count distinct cell_id
  dataset_summary2 <- datasets_combined %>% distinct(junction_aa, cell_id, .keep_all = TRUE) %>% group_by(dataset_id, clonotype) %>%
    summarise('junction_aa' = paste(unique(junction_aa),collapse=","),'Expansion' = n_distinct(cell_id), 'IgM' = sum(grepl(c_call, pattern='IGHM')),
              'IgG' = sum(grepl(c_call, pattern='IGHG')), 'IgA' = sum(grepl(c_call, pattern='IGHA')),
              'IgD' = sum(grepl(c_call, pattern='IGHD')), 'IgE' = sum(grepl(c_call, pattern='IGHE')),
              'Missing' = sum(grepl(c_call, pattern = 'None')))  %>%
    arrange(desc(Expansion)) %>%
    arrange(dataset_id) %>%
    drop_na() %>%
    mutate(Rank = seq(n()))



  # Get TopClonotypes
  dataset_summary_top <- dataset_summary2 %>% group_by(dataset_id)

  #reformat in preparation of subsetting
  dataset_summary_top <- melt(as.data.table(dataset_summary_top),
                              id.vars = c('dataset_id', 'clonotype', 'Expansion', 'Rank'), measure.vars = c('IgA', 'IgD', 'IgG', 'IgM', 'Missing', 'IgE'),
                              value.name = 'Count', variable.name = 'Isotype')

  dataset_summary_top$dataset_id <- as.factor(dataset_summary_top$dataset_id)
  #levels(dataset_summary_top$dataset_id) <-  c(paste(rep('G'), c(1:8)), paste(rep('H'), c(1:8)))

  #subset to 4 datasets: with highest quality H1-H4.
  # take top 10 IgG with count>1
  #note ranking is based on global count not IgG count
  #not only ten per since some have same count.
  dataset_summaryIgG_top <- dataset_summary_top %>%
    filter(Isotype=='IgG', Count > 1, dataset_id%in%c('137988_combined_db', '137989_combined_db',
                                                      '137990_combined_db', '137991_combined_db')) %>%
    group_by(dataset_id) %>% top_n(10)


  #resplit dataset paired
  dataset_paired_list<-list()
  for(i in 1:length(unique(datasets_paired$dataset_id))){
    curr_id<-unique(datasets_paired$dataset_id)[i]
    curr_df<-datasets_paired[datasets_paired$dataset_id==curr_id,]
    dataset_paired_list[[curr_id]]<-curr_df
  }

  #top IgG clonotypes --. for each top clonotype extract most expanded IgG clone. (if more than one equally expanded pick first in line)
  datasets_out2 <- list()
  for(irow in c(1:length(dataset_summaryIgG_top$dataset_id))){

    df1 <- datasets[[dataset_summaryIgG_top$dataset_id[irow]]]
    df2 <- dataset_paired_list[[dataset_summaryIgG_top$dataset_id[irow]]]

    df1_cdr3s <- df1$junction_aa[df1$clonotype == dataset_summaryIgG_top$clonotype[irow]]
    df2_cellIds <- df2$cell_id[df2$junction_aa%in%df1_cdr3s]
    df2 <- df2[df2$cell_id%in%df2_cellIds, ]

    datasets_out2[[irow]] <- .get_most_expanded_IgGclone(df2)

  }

  #datasets_out = top 16 clones
  #datsets_out2 = top clones for each top clonotype.
  datasets_out2 <- bind_rows(datasets_out2)
  #kick out clones that are already covered in top 16 clones
  datasets_out2 <- datasets_out2[!datasets_out2$junction_aa%in%datasets_out$junction_aa, ]
  ##kick out clones with missing sequence id (cellranger.. sanity check that data is complete)
  datasets_out2 <- datasets_out2[!is.na(datasets_out2$sequence_id), ]



  #50 top clonotypes --> 50 clones based on clonotypes need to be paired (some get lost)
  #add 16 top clones that are not in top clonotypes (some overlap there)
  #--> 36
  Top36Set <- paste(datasets_out2$junction_aa, datasets_out2$junction_aa.light, sep = '-')





  #Second set
  # Get all duplicated cell ids with shared sequence and unique V-Gene assigment and print their overlap
  flagged_cell_ids <- datasets_paired %>% group_by(sequence, cell_id) %>% summarise(Count = n(),
                                                                                    VH = length(unique(v_gene)),
                                                                                    VL = length(unique(v_gene.light)))
  flagged_cell_ids <- flagged_cell_ids$cell_id[flagged_cell_ids$Count>1 & flagged_cell_ids$VH==1 & flagged_cell_ids$VL==1]


  # Remove flagged ids
  datasets_paired <- datasets_paired[!datasets_paired$cell_id%in%flagged_cell_ids, ]
  datasets_paired_df<-as.data.frame(datasets_paired)
  library(data.table)
  list_datasets_covid<-split(datasets_paired,list(datasets_paired$dataset_id))
  names(list_datasets_covid)

  idx_list <- split(seq_len(nrow(data)),
                    list(data$v_call,data$j_call, data$cdr3_len))


  # Largest clonotypes overall

  paired_summary <- datasets_paired %>% group_by(clonotype)  %>%
    summarise('Expansion' = n_distinct(cell_id), 'IgM' = sum(grepl(c_call, pattern='IGHM')),
              'IgG' = sum(grepl(c_call, pattern='IGHG')), 'IgA' = sum(grepl(c_call, pattern='IGHA')),
              'IgD' = sum(grepl(c_call, pattern='IGHD')), 'IgE' = sum(grepl(c_call, pattern='IGHE')),
              'Missing' = sum(grepl(c_call, pattern = 'None')))  %>%
    arrange(desc(Expansion)) %>%
    drop_na() %>%
    mutate(Rank = seq(n()))

  #extract clonotypes that have IgG and at least 2 cells
  #(majority clonotype does not need to be IgG)
  clonotypes_to_express <- paired_summary$clonotype[paired_summary$IgG > 0 & paired_summary$Expansion > 1]


  #load top 16 clones and top 36 clonotype clones
  first_batch <- TopClone_IgG_PerPatient
  second_batch <- Top36Set


  #kick out clonotypes_to_express that are in either of the two
  clonotypes_to_express <- clonotypes_to_express[!clonotypes_to_express%in%datasets_paired$clonotype[datasets_paired$cell_id%in%first_batch$cell_id]]
  clonotypes_to_express <- clonotypes_to_express[!clonotypes_to_express%in%datasets_paired$clonotype[datasets_paired$cell_id%in%second_batch$cell_id]]
  clonotypes_to_express <- clonotypes_to_express[c(1:96)]

  Top96Clonontype_IgG <- list()

  for(idx in c(1:length(clonotypes_to_express))){

    top_mAb <- datasets_paired[datasets_paired$clonotype==clonotypes_to_express[idx], ] %>%
      group_by(combined_chain) %>%
      summarise(Expansion = length(unique(cell_id))) %>%
      arrange(desc(Expansion))

    Top96Clonontype_IgG[[idx]] <- datasets_paired[datasets_paired$combined_chain%in%top_mAb$combined_chain[1], ][1,]

  }

  Top96Clonontype_IgG <- bind_rows(Top96Clonontype_IgG)


  cands<-list()
  cands[["Top36Clonontype_IgG"]]<-Top36Set
  cands[["Top96Clonontype_IgG"]]<-Top96Clonontype_IgG

  return(cands)


}
