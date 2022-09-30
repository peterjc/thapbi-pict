# Export Figure 3 data
#
# Assume's you are in RStudio and have run the Figure 3 section in
# https://github.com/alexpiper/HemipteraMetabarcodingMS/blob/master/hemiptera_metabarcoding.Rmd
#
# We're just aggregating the data already filtered and compiled, and exporting as TSV.
# To simplify linking the tables, aggregating at genus level to pool the two Acizzia spp. etc.

require(comprehenr)

export_ALL <- df_exp[c("pool_comp", "Genus", "Abundance")]
export_ALL <- export_ALL[order(export_ALL$Genus, export_ALL$pool_comp),]
names(export_ALL)[3] <- "Expected"

export_COI <- df_coi[c("pool_comp", "Genus", "Abundance")]
export_COI <- export_COI[order(export_COI$Genus, export_COI$pool_comp),]
export_COI <- aggregate(. ~  pool_comp + Genus, data = export_COI, sum)
assertthat::are_equal(export_COI$pool_comp, export_ALL$pool_comp)
assertthat::are_equal(export_COI$Genus, export_ALL$Genus)
export_ALL["COI"] <- export_COI$Abundance

export_18S <- df_18s[c("pool_comp", "Genus", "Abundance")]
export_18S <- export_18S[order(export_18S$Genus, export_18S$pool_comp),]
export_18S <- aggregate(. ~  pool_comp + Genus, data = export_18S, sum)
assertthat::are_equal(export_18S$pool_comp, export_ALL$pool_comp)
assertthat::are_equal(export_18S$Genus, export_ALL$Genus)
export_ALL["18S"] <- export_18S$Abundance

export_12S <- df_12s[c("pool_comp", "Genus", "Abundance")]
export_12S <- export_12S[order(export_12S$Genus, export_12S$pool_comp),]
export_12S <- aggregate(. ~  pool_comp + Genus, data = export_12S, sum)
assertthat::are_equal(export_12S$pool_comp, export_ALL$pool_comp)
assertthat::are_equal(export_12S$Genus, export_ALL$Genus)
export_ALL["12S"] <- export_12S$Abundance

export_ALL <- export_ALL[order(export_ALL$pool_comp, export_ALL$Genus),]
# Relabel to match desired captions by converting pool_comp which is e.g. 4.5
# for individuals class 4 (i.e. from 100, 250, 500 and 1000) and pool 5:
names(export_ALL)[1] <- "Caption"
export_ALL[1] <- to_vec(for(px in export_ALL[1]) paste(c("100", "250", "500", "1000")[round(10*(px-floor(px)))], "Pool", floor(px)))

write.table(export_ALL, file="figure3.tsv", sep="\t", quote=FALSE, row.names=FALSE)
