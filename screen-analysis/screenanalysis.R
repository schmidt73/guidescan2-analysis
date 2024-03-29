library(data.table)
library(plyr)
library(qvalue)
library(ggplot2)


# load data

design <- fread("datafiles/a549_screen.csv")

foldername <- "datafiles/"
filenames <- c("GS1.sgrna_summary.txt", "GS2.sgrna_summary.txt",
               "GS3.sgrna_summary.txt", "GS4.sgrna_summary.txt",
               "GS5.sgrna_summary.txt", "GS6.sgrna_summary.txt")

res <- sapply(paste0(foldername, filenames), fread, simplify = FALSE)
names(res) <- c("GuideScan2", "Root2016", "Moffat2015", "Elling2020",
                "Bassik2017", "Sabatini2015")

design.guidescan.control <- fread(paste0(foldername, "GS1_library_control.csv"),
                                  header = FALSE)[, .(V1, V2)]
colnames(design.guidescan.control) <- c("Identifier", "sgRNA")
design.updated <- merge(design, design.guidescan.control, by = "sgRNA",
                        all.x = TRUE, suffixes = c("", ".y"))
design.updated[!is.na(Identifier.y), Identifier := Identifier.y]

# merge all datasets
res <- as.data.table(ldply(res, .id = "lib"))
# remove gRNAs with very low read count
res <- res[control_count > 100]
# reorder library names for consistenty with other plots in the analysis
res$lib <- factor(res$lib, levels = c("GuideScan2", "Moffat2015", "Bassik2017",
                                      "Root2016", "Elling2020", "Sabatini2015"))

# this table is a simple merge
res.simple <- merge(res, design.updated, by.x = "sgrna", by.y = "Identifier")

# but for the gRNAs that appear in multiple libraries, we want their association
# with each of these libraries preserved, i.e. screen results duplicated for
# each occurrence of gRNA in any of the libraries
res.withgrna <- merge(res, design.updated[, .(Identifier, sgRNA)],
                      by.x = "sgrna", by.y = "Identifier")
res <- merge(res.withgrna, design, by = "sgRNA", all.y = TRUE,
             suffixes = c("", ".y"))
res <- res[!is.na(LFC)]
res[, gene := Gene]
res[Library == "Guidescan", Library := "GuideScan2"]
res$Library <- factor(res$Library, levels = c("GuideScan2", "Moffat2015", "Bassik2017",
                                      "Root2016", "Elling2020", "Sabatini2015"))
saveRDS(res, "datafiles/results-crispr-screen.rds")
# library: library name assigned to a gRNA sequence, unique for the gRNA and
# multiple copies of the same gRNA across libraries
# library_design: library name for a gRNA at the time of design, will be
# different for multiple compies of the same gRNA included in different libraries
fwrite(res[, .(sgRNA, Identifier, library = lib, library_design = Library,
               gene = Gene, type = Type,
               control_count, treatment_count, LFC)],
       "datafiles/results-crispr-screen.csv")
# refer to the variable "Library" in subsequent analysis for comparing gRNA
# across libraries in different arms of the experiment


# analysis of essential genes

res.ess <- res[Type == "essential_gene_targeting"]
res.ess[, .N, by = c("gene", "Library")]
# test signficance for each library
res.ess.testsignif <- res.ess[, .(p = wilcox.test(LFC, alternative = "less")[[3]]),
                              by = c("gene", "Library")]
res.ess.testsignif$Library <- factor(res.ess.testsignif$Library,
                                     levels = c("GuideScan2", "Moffat2015", "Bassik2017",
                                                "Root2016", "Elling2020", "Sabatini2015"))
res.ess.testsignif[, .N, by = Library]
res.ess.testsignif[p < 0.1][, .N, by = Library]
# Library  N
# 1:     Root2016 54
# 2: Sabatini2015 77
# 3:   Elling2020 80
# 4:   Moffat2015 60
# 5:   GuideScan2 83
# 6:   Bassik2017 76
ggplot(res.ess.testsignif, aes(x = Library, fill = p < 0.1)) +
    geom_bar() +
    scale_fill_manual(values = c("gray", "orange")) +
    labs(x = "", y = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/fraction-sig-essential.pdf", width = 3, height = 4)

# how many genes for each library?
res.ess[, .N, by = c("Library", "gene")][, mean(N) , by = Library]
# all libraries have 6 gRNAs or more per gene except Root2016, which has 4 gRNAs
# subsample GuideScan to 4 gRNAs per gene and do the test
set.seed(0)
values <- sapply(1:100,
                 function(i)
                 res.ess[Library == "GuideScan2"][,
        .(p = wilcox.test(sample(LFC, min(4, .N)), alternative = "less")[[3]]),
        by = gene][p < 0.1][, .N])
mean(values); sd(values)
# [1] 61.88
# [1] 2.345337

# this selects 85 out of 100 genes that are consistently strongly affected by all libraries:
essential.genes <- res.ess[, .(mean.lfc = mean(LFC)), by = gene][mean.lfc < -1][, gene]
res.ess.selected <- res.ess[gene %in% essential.genes]


# boxplots for essential and non-essential genes

table(res$Type)
res$Type <- factor(res$Type,
    levels = c("essential_gene_targeting", "non_essential_gene_targeting"))

ggplot(res, aes(x = Library, y = LFC, fill = Type)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.2) +
    scale_fill_manual(values = c("orange", "lightblue"),
                      labels = c("essential", "nonessential")) +
    guides(fill=guide_legend(title="")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(y = "gRNA LFC", x = "") +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/results-boxplot-notches-grna.pdf", width = 4, height = 4)

ggplot(res[Type == "essential_gene_targeting"], aes(x = Library, y = LFC, fill = Type)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.2) +
    scale_fill_manual(values = c("orange"), labels = c("essential")) +
    guides(fill=guide_legend(title="")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(y = "gRNA LFC", x = "") +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/results-boxplot-notches-grna-essential.pdf", width = 3.5, height = 4)

ggplot(res[Type == "non_essential_gene_targeting"], aes(x = Library, y = LFC, fill = Type)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.2) +
    scale_fill_manual(values = c("lightblue"), labels = c("nonessential")) +
    guides(fill=guide_legend(title="")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(y = "gRNA LFC", x = "") +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/results-boxplot-notches-grna-nonessential.pdf", width = 3.5, height = 4)

ggplot(res[, .(mean.lfc = mean(LFC)), by = c("Gene", "Library", "Type")],
       aes(x = Library, y = mean.lfc, fill = Type)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.2) +
    scale_fill_manual(values = c("orange", "lightblue"),
                      labels = c("essential", "nonessential")) +
    guides(fill=guide_legend(title="")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(y = "gene mean LFC", x = "") +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/results-boxplot-notches.pdf", width = 4, height = 4)

ggplot(res[, .(mean.lfc = mean(LFC)), by = c("Gene", "Library", "Type")][Type == "essential_gene_targeting"],
       aes(x = Library, y = mean.lfc, fill = Type)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.2) +
    scale_fill_manual(values = c("orange"), labels = c("essential")) +
    guides(fill=guide_legend(title="")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(y = "gene mean LFC", x = "") +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/results-boxplot-notches-essential.pdf", width = 3, height = 4)

ggplot(res[, .(mean.lfc = mean(LFC)), by = c("Gene", "Library", "Type")][Type == "non_essential_gene_targeting"],
       aes(x = Library, y = mean.lfc, fill = Type)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.2) +
    scale_fill_manual(values = c("lightblue"), labels = c("nonessential")) +
    guides(fill=guide_legend(title="")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(y = "gene mean LFC", x = "") +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/results-boxplot-notches-nonessential.pdf", width = 3, height = 4)


# analysis of non-essential genes

res.noness <- res[Type == "non_essential_gene_targeting"]
res.noness[, .N, by = c("gene", "Library")]

lfc.threshold <- -2
for (libname in c("Root2016", "Elling2020", "Bassik2017", "Moffat2015",
                  "Sabatini2015")) {
    res.noness.lib <-
        res.noness[Library %in% c(libname, "GuideScan2")]
    profiled.genes <- intersect(res.noness.lib[Library == libname, gene],
                                res.noness.lib[Library == "GuideScan2", gene])
    res.noness.lib <- res.noness.lib[gene %in% profiled.genes]
    res.noness.lib$Library <- droplevels(factor(res.noness.lib$Library))
    res.noness.lib.gene <- res.noness.lib[, .(min.lfc = min(LFC)),
                                          by = c("gene", "Library")]
    contable <- res.noness.lib.gene[,
                                    addmargins(table(min.lfc < lfc.threshold, Library))]
    pval <- phyper(contable[2, 2], contable[2, 3], contable[1, 3], contable[3, 2],
                   lower.tail = FALSE)
    # pval <- fisher.test(contable[1:2, 1:2])$p.value
    ggplot(res.noness.lib.gene, aes(x = Library, fill = min.lfc < lfc.threshold)) +
        geom_bar() +
        scale_fill_manual(values = c("gray", "red")) +
        labs(x = "", y = "") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90)) +
        annotate("text", x = 1.5, y = (contable[3, 1] + contable[3, 2]) / 2,
                 label = ifelse(pval < 0.05,
                                sprintf("p < %.1e", pval),
                                sprintf("p = %.1e", pval))) +
        guides(fill=guide_legend(title=sprintf("nonessential\ngenes\nwith gRNA\nLFC < %s",
                                               lfc.threshold)))
    ggsave(sprintf("plots/fraction-negative-GuideScan2-%s.pdf", libname),
           width = 2.65, height = 4)
}


# show dependence on specificity:
# turns out every library has gRNAs with specificity < 0.1, and such gRNAs
# are more likely to be false positives i.e. have low negative logFC
ggplot(res.noness, aes(x = Library, y = LFC, fill = Specificity < 0.05)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.2) +
    scale_fill_manual(values = c("gray", "red"), name = "specificity < 0.05") +
    labs(x = "", y = "gRNA LFC") + # , title = "non-essential genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/results-nonessential-boxplot-specificity.pdf",
       width = 4, height = 4)

ggplot(res.noness, aes(x = Library, fill = Specificity < 0.05)) +
    geom_bar() +
    scale_fill_manual(values = c("gray", "red")) +
    labs(x = "", y = "gRNA LFC") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/count-nonessential-specificity-all.pdf",
       width = 4, height = 4)
ggplot(res.noness[Specificity < 0.05], aes(x = Library)) +
    geom_bar(fill = "red") +
    labs(x = "", y = "gRNA count, specificity < 0.05") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/count-nonessential-specificity.pdf",
       width = 2, height = 4)



# analyze GuideScan2 non-targeting controls

guidescan.updated <- "../fromMinsi/essential_genome/GS1_with_control/GS1_all.sgrna_summary.txt"
guidescan.updated <- fread(guidescan.updated)
guidescan.updated <- guidescan.updated[control_count > 100]
guidescan <- res[Library == "GuideScan2"][, .(sgrna, gene = Gene, control_count,
                                     treatment_count, LFC)]
guidescan.updated <- guidescan.updated[, .(sgrna, gene = Gene, control_count,
                                           treatment_count, LFC)]
guidescan.merged <- merge(guidescan, guidescan.updated, by = "sgrna", all = TRUE)
guidescan.merged[!is.na(LFC.x)][!is.na(LFC.y)][, cor(LFC.x, LFC.y, method = "spearman")]
# shows correlation of 1

guidescan.updated <- merge(guidescan.updated, design.updated,
                           by.x = "sgrna", by.y = "Identifier")
table(guidescan.updated$Type)
guidescan.updated$lib <- "GuideScan2"

guidescan.updated$Type <-
    factor(guidescan.updated$Type,
           levels = c("essential_gene_targeting", "non_essential_gene_targeting",
                              "non_targeting_control", "safe_targeting_control"))
ggplot(guidescan.updated,
       aes(x = lib, y = LFC, fill = Type)) +
    geom_boxplot(notch = TRUE, outlier.size = 0.2) +
    scale_fill_manual(values = c("orange", "lightblue", "gray60", "gray80"),
                      labels = c("essential", "nonessential",
                                 "nontargeting", "safetargeting")) +
    guides(fill=guide_legend(title="")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(y = "gRNA LFC", x = "") +
    theme(axis.text.x = element_text(angle = 90))
ggsave("plots/results-boxplot-notches-grna-guidescan.pdf", width = 3.5, height = 4)

saveRDS(guidescan.updated,
        "datafiles/results-crispr-screen-guidescan-withcontrols.rds")
fwrite(guidescan.updated[, .(sgRNA, Identifier = sgrna, library = lib, gene = Gene,
                             type = Type,
               control_count, treatment_count, LFC)],
       "datafiles/results-crispr-screen-guidescan-withcontrols.csv")


sessionInfo()

