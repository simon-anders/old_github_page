library( readr )
library( dplyr )
library( tidyr )
library( reshape2 )
library( pheatmap )

tbl <- 
   read_csv( "~/work/AML_RNA/voom_cor_analysis_protein_coding.csv" ) %>%
   rename( ensg = X )

tbl <- tbl %>%
   mutate( adj.P.Val = p.adjust( P.Value, "BH" ) )

load( "~/work/AML_RNA/dg.rda" )

dgi <- 7  # index to dg; 6 is immunomodulators

sigGenes <- 
   tbl %>%
   inner_join( dg[[dgi]] ) %>%
   filter( adj.P.Val < .05, abs(logFC) > .2 ) %>%
   distinct( ensg )
sigGenes

lfcM <- 
   tbl %>%
   inner_join( dg[[dgi]] ) %>%
   inner_join( sigGenes ) %>%
   acast( Gene_name ~ Drug.name, value.var = "logFC" )

pheatmap( lfcM )
   
lfcMc <- lfcM[ rev(hclust(as.dist(1-cor(t(lfcM))))$order), hclust(dist(t(lfcM)))$order ]

dsrt <- read.csv( "~/work/AML_RNA/AML_DSS_52S.csv", row.names=2 )[,-1][colnames(lfcMc),]

cnts <- read.csv( "~/work/AML_RNA/AML_read_count_75S_15062015.csv", row.names=1 )
lexpr <- log2( t( t(cnts) / DESeq2::estimateSizeFactorsForMatrix(cnts) ) + .1 )

genes <- inner_join( 
   tibble( Gene_name = rownames(lfcMc) ), 
   tbl %>% distinct( ensg, Gene_name ) )

samples <- intersect( colnames(cnts), colnames(dsrt) )

smplTbl <- read.csv( "~/work/AML_RNA/target_80S.csv" )

f = file( "AML/input.js", "w" )
writeLines( "genes = ", f )
writeLines( toJSON( rownames(lfcMc) ), f )
writeLines( "\ndrugs = ", f )
writeLines( toJSON( colnames(lfcMc) ), f )
writeLines( "\nlfc = ", f )
writeLines( toJSON( lfcMc ), f )
writeLines( "\nexprs = ", f )
writeLines( toJSON( lexpr[ genes$ensg, samples ] ), f )
writeLines( "\ndss = ", f )
writeLines( toJSON( as.matrix(dsrt[ , samples ]) ), f )
writeLines( "\nsmplInfo = ", f )
writeLines( toJSON( inner_join( smplTbl, tibble( Sample_ID = samples ) ) ), f )
close(f)
