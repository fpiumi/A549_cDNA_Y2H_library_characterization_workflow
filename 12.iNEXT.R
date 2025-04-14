library(iNEXT)

# saturation curve calculation:

iNEXT_genome_splice_A549 <- iNEXT(genome_splice_counts, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)


# plot:

iNEXT_genome_splice_plot <- ggiNEXT(iNEXT_genome_splice_A549, type=1, se=TRUE, grey=FALSE) + theme_minimal()