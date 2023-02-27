Input_Paths = "/AWS/directories/we/need/to/pull/from"
`%notin%` <- Negate(`%in%`)

# Function for calculating PRS elements
calculatePRSElement = function(dosages, effectSizes) {
    prsElements = effectSizes * dosages
    return(prsElements)
}

# Load dosages.Dosage is a weighted metric for whatever allele is designated the alternate allele
dosages = read.table(
    paste0(Input_Paths, "/prs_dosages_outer_test.tsv"),
    header = T
)

# Load effect sizes. Effect size should be provided by the paper/ papers in question.
effects = read.table(
    paste0(Input_Paths, "/prs_effect_sizes_outer_test.tsv"),
    header = T
)

# Merge
merged = merge(
    effects,
    dosages,
    by = "ID",
    sort = F
)

# Release mem
rm("dosages")
rm("effects")

# Remove variants with NAs
merged = na.omit(merged)

# Remove variants with mismatched alleles
remove = merged[(merged$ALT.x != merged$REF.y) && (merged$ALT.x != merged$ALT.y),]$ID
merged = merged[!(merged$ID %in% remove),]

# Calculate PRS elements
merged[merged['ALT_EFFECT_SIZE'] < 0, c(7:ncol(merged))] = 2 - merged[merged['ALT_EFFECT_SIZE'] < 0, c(7:ncol(merged))]

merged['ALT_EFFECT_SIZE'] = abs(merged['ALT_EFFECT_SIZE'])

prsElements = apply(merged[7:ncol(merged)], 2, calculatePRSElement, effectSizes = merged$ALT_EFFECT_SIZE)

rownames(prsElements) = merged$ID

prsElements = as.data.frame(t(prsElements))


#Calculates the final PRS.That is done by summing the (effect size x dosage) for all SNPs included in the PRS for each individual and then dividing by the total number of SNPs (i.e., you find the average effect size x dosage across all SNPs included in the PRS for each individual). This average is the PRS for the individual.
prs = data.frame(
  ID = rownames(prsElements),
  PRS = rep(0, nrow(prsElements))
)
prs["PRS"] = apply(prsElements[colnames(prsElements) %in% variants$ID], 1, sum) / nrow(variants)




