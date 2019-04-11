
# library imports
library(tidyverse)
library(scales)
library(limma)
library(edgeR)
library(psych)

# get the default plot width and height
width <- options()$repr.plot.width
height <- options()$repr.plot.height

# load the IRS-normalized data and check the table
data_import <- read_tsv("labeled_grouped_protein_summary_TMT_9_131N_IRS_normalized.txt", guess_max = 10326)

# the "Filter" column flags contams and decoys
# the "Missing" column flags proteins without reporter ion intensities (full sets missing)
# the prepped table from pandas is sorted so these are the upper rows
data_all <- filter(data_import, is.na(Filter), is.na(Missing))

# save gene names for edgeR so we can double check that results line up
accessions <- data_all$Accession

# see how many rows in the table
nrow(data_all)

# we want to get the SL normed columns, and subsetted by condition
sl_all <- data_all %>%
  select(starts_with("SLNorm"))
sl_HeH <- sl_all %>% select(contains("_HeH_"))
sl_ETV6 <- sl_all %>% select(contains("_ETV6-RUNX1_"))

# and the IRS normed columns by condition
irs_all <- data_all %>%
  select(starts_with("IRSNorm"))
irs_HeH <- irs_all %>% select(contains("_HeH_"))
irs_ETV6 <- irs_all %>% select(contains("_ETV6-RUNX1_"))

# and collect the pooled channels before and after IRS
sl_pool <- sl_all %>% select(contains("pool"))
irs_pool <- irs_all %>% select(contains("pool"))

# multi-panel scatter plot grids from the psych package
pairs.panels(log2(sl_pool), lm = TRUE, main = "Pooled Std before IRS")
pairs.panels(log2(irs_pool), lm = TRUE, main = "Pooled Std after IRS")

# multi-panel scatter plot grids
heh_sample <- sample(1:18, 5)
pairs.panels(log2(sl_HeH[heh_sample]), lm = TRUE, main = "HeH before IRS (random 5)")
pairs.panels(log2(irs_HeH[heh_sample]), lm = TRUE, main = "HeH after IRS (same 5)")

# multi-panel scatter plot grids
etv6_sample  <- sample(1:9, 5)
pairs.panels(log2(sl_ETV6[etv6_sample]), lm = TRUE, main = "ETV6-RUNX1 before IRS (random 5)")
pairs.panels(log2(irs_ETV6[etv6_sample]), lm = TRUE, main = "ETV6-RUNX1 after IRS (same 5)")

# get the biological sample data into a DGEList object
group = c(rep('HeH', 18), rep('ETV6', 9))
y_sl <- DGEList(counts = cbind(sl_HeH, sl_ETV6), group = group, genes = accessions)
y_irs <- DGEList(counts = cbind(irs_HeH, irs_ETV6), group = group, genes = accessions)

# run TMM normalization (also includes a library size factor)
y_sl <- calcNormFactors(y_sl)
y_irs <- calcNormFactors(y_irs)

# set some colors by condition
colors = c(rep('red', 18), rep('blue', 9))

# check the clustering
plotMDS(y_sl, col = colors, main = "SL: all samples")
plotMDS(y_irs, col = colors, main = "IRS: all samples")

# we do not want the technical replicates in the mix for dispersion estimates
irs <- cbind(irs_HeH, irs_ETV6)

# load a new DGEList object (need to update the groups)
y <- DGEList(counts = irs, group = group, genes = accessions) # group was set above
y <- calcNormFactors(y)

# see what the normalization factors look like
y$samples

# Compute the normalized intensities (start with the IRS data)
# sample loading adjusts each channel to the same average total
lib_facs <- mean(colSums(irs)) / colSums(irs)

# print("Sample loading normalization factors")
print("Library size factors")
round(lib_facs, 4)

# the TMM factors are library adjustment factors (so divide by them)
norm_facs <- lib_facs / y$samples$norm.factors

# print these final correction factors
print("Combined (lib size and TMM) normalization factors")
round(norm_facs, 4)

# compute the normalized data as a new data frame
irs_tmm <- sweep(irs, 2, norm_facs, FUN = "*")
colnames(irs_tmm) <- str_c(colnames(irs), "_TMMnorm") # add suffix to col names
# head(results) # check that the column headers are okay

long_results <- gather(irs_tmm, key = "sample", value = "intensity") %>%
  mutate(log_int = log10(intensity)) %>%
  extract(sample, into = 'group', ".*?_(.*?)_", remove = FALSE)
head(long_results)

ggplot(long_results, aes(x = sample, y = log_int, fill = group)) +
  geom_boxplot(notch = TRUE) +
  coord_flip() + 
  ggtitle("edgeR normalized data")

# look at normalized intensity distributions for each sample
boxplot(log10(irs_tmm), col = colors,
        xlab = 'TMT samples', ylab = 'log10 Intensity', 
        main = 'edgeR normalized data', notch = TRUE)

ggplot(long_results, aes(x = log_int, color = sample)) +
  geom_density() +
  guides(color = FALSE) +
  ggtitle("edgeR normalized data (with legend is too busy)")

# we can compare CVs before and after IRS
sl <- cbind(sl_HeH, sl_ETV6)

# save column indexes for different conditions (indexes to data_raw frame)
# these make things easier (and reduce the chance for errors)
HeH <- 1:18
ETV6 <- (1:9) + 18

# create a CV computing function
CV <- function(df) {
    ave <- rowMeans(df)
    sd <- apply(df, 1, sd)
    cv <- 100 * sd / ave
}

# put CVs in data frames to simplify plots and summaries
cv_frame <- data.frame(HeH_sl = CV(sl[HeH]), HeH_final = CV(irs_tmm[HeH]), 
                       ETV6_sl = CV(sl[ETV6]), ETV6_final = CV(irs_tmm[ETV6]))


# see what the median CV values are
medians <- apply(cv_frame, 2, FUN = median)
print("Median CVs by condition, before/after IRS (%)")
round(medians, 1)

# see what the CV distibutions look like
# need long form for ggplot
long_cv <- gather(cv_frame, key = "condition", value = "cv") %>%
  extract(condition, into = 'group', "(.*?)_+", remove = FALSE)

# traditional boxplots
cv_plot <- ggplot(long_cv, aes(x = condition, y = cv, fill = group)) +
  geom_boxplot(notch = TRUE) +
  ggtitle("CV distributions")

# vertical orientation
cv_plot

# horizontal orientation
cv_plot + coord_flip()

# density plots
ggplot(long_cv, aes(x = cv, color = condition)) +
  geom_density() +
  coord_cartesian(xlim = c(0, 150)) +
  ggtitle("CV distributions")

# compute dispersions and plot BCV
y <- estimateDisp(y)
plotBCV(y, main = "BCV plot of IRS normed, TMM normed, all 27")

# the exact test object has columns like fold-change, CPM, and p-values
et <- exactTest(y, pair = c("HeH", "ETV6"))

# this counts up, down, and unchanged genes (proteins) at 10% FDR
summary(decideTestsDGE(et, p.value = 0.10))

# the topTags function adds the BH FDR values to an exactTest data frame 
# make sure we do not change the row order (the sort.by parameter)!
topTags(et, n = 25)
tt <- topTags(et, n = Inf, sort.by = "none")
tt <- tt$table    # tt is a list. We just need the "table" data frame

# make an MD plot (like MA plot)
plotMD(et, p.value = 0.10)
abline(h = c(-1, 1), col = "black")

# check the p-value distribution
ggplot(tt, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(et$table$PValue, breaks = 100, 
                                    plot = FALSE)$counts[26:100])) +
  ggtitle("HeH vs ETV6 PValue distribution")

# get the averages within each condition 
# results already has the normalized data in its left columns
tt$ave_HeH <- rowMeans(irs_tmm[HeH])
tt$ave_ETV6 <- rowMeans(irs_tmm[ETV6])

# add the cadidate status column
tt <- tt %>%
  mutate(candidate = cut(FDR, breaks = c(-Inf, 0.01, 0.05, 0.10, 1.0),
  labels = c("high", "med", "low", "no")))

tt %>% count(candidate)  # count candidates

ggplot(tt, aes(x = logFC, fill = candidate)) +
  geom_histogram(binwidth=0.1, color = "black") +
  facet_wrap(~candidate) +
  coord_cartesian(xlim = c(-4, 4)) +
  ggtitle("HeH vs ETV6-RUNX1 logFC distributions by candidate")

# ================= reformat edgeR test results ================================

collect_results <- function(df, tt, x, xlab, y, ylab) {
    # Computes new columns and extracts some columns to make results frame
        # df - data in data.frame
        # tt - top tags table from edgeR test
        # x - columns for first condition
        # xlab - label for x
        # y - columns for second condition
        # ylab - label for y
        # returns a new dataframe
    
    # condition average vectors
    ave_x <- rowMeans(df[x])
    ave_y <- rowMeans(df[y])
    
    # FC, direction, candidates
    fc <- ifelse(ave_y > ave_x, (ave_y / ave_x), (-1 * ave_x / ave_y))
    direction <- ifelse(ave_y > ave_x, "up", "down")
    candidate <- cut(tt$FDR, breaks = c(-Inf, 0.01, 0.05, 0.10, 1.0), 
                     labels = c("high", "med", "low", "no"))
    
    # make data frame
    temp <- cbind(df[c(x, y)], data.frame(logFC = tt$logFC, FC = fc, 
                                          PValue = tt$PValue, FDR = tt$FDR, 
                                          ave_x = ave_x, ave_y = ave_y, 
                                          direction = direction, candidate = candidate, 
                                          Acc = tt$genes)) 
    
    # fix column headers for averages
    names(temp)[names(temp) %in% c("ave_x", "ave_y")]  <- str_c("ave_", c(xlab, ylab))    
    
    temp # return the data frame
}

# get the results
results <- collect_results(irs, tt, HeH, "HeH", ETV6, "ETV6")

transform <- function(results, x, y) {
    # Make data frame with some transformed columns
        # results - results data frame
        # x - columns for x condition
        # y - columns for y condition
        # return new data frame
    df <- data.frame(log10((results[x] + results[y])/2), 
                     log2(results[y] / results[x]), 
                     results$candidate,
                     -log10(results$FDR))
    colnames(df) <- c("A", "M", "candidate", "P")
    
    df # return the data frame
}

MA_plots <- function(results, x, y, title) {
    # makes MA-plot DE candidate ggplots
        # results - data frame with edgeR results and some condition average columns
        # x - string for x-axis column
        # y - string for y-axis column
        # title - title string to use in plots
        # returns a list of plots 
    
    # uses transformed data
    temp <- transform(results, x, y)
    
    # 2-fold change lines
    ma_lines <- list(geom_hline(yintercept = 0.0, color = "black"),
                     geom_hline(yintercept = 1.0, color = "black", linetype = "dotted"),
                     geom_hline(yintercept = -1.0, color = "black", linetype = "dotted"))

    # make main MA plot
    ma <- ggplot(temp, aes(x = A, y = M)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scale_y_continuous(paste0("logFC (", y, "/", x, ")")) +
        scale_x_continuous("Ave_intensity") +
        ggtitle(title) + 
        ma_lines
    
    # make separate MA plots
    ma_facet <- ggplot(temp, aes(x = A, y = M)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scale_y_continuous(paste0("log2 FC (", y, "/", x, ")")) +
        scale_x_continuous("log10 Ave_intensity") +
        ma_lines +
        facet_wrap(~ candidate) +
        ggtitle(str_c(title, " (separated)"))

    # make the plots visible
    print(ma)
    print(ma_facet)
}    

scatter_plots <- function(results, x, y, title) {
    # makes scatter-plot DE candidate ggplots
        # results - data frame with edgeR results and some condition average columns
        # x - string for x-axis column
        # y - string for y-axis column
        # title - title string to use in plots
        # returns a list of plots
    
    # 2-fold change lines
    scatter_lines <- list(geom_abline(intercept = 0.0, slope = 1.0, color = "black"),
                          geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted"),
                          geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted"),
                          scale_y_log10(),
                          scale_x_log10())

    # make main scatter plot
    scatter <- ggplot(results, aes_string(x, y)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        ggtitle(title) + 
        scatter_lines

    # make separate scatter plots
    scatter_facet <- ggplot(results, aes_string(x, y)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scatter_lines +
        facet_wrap(~ candidate) +
        ggtitle(str_c(title, " (separated)")) 

    # make the plots visible
    print(scatter)
    print(scatter_facet)
}

volcano_plot <- function(results, x, y, title) {
    # makes a volcano plot
        # results - a data frame with edgeR results
        # x - string for the x-axis column
        # y - string for y-axis column
        # title - plot title string
    
    # uses transformed data
    temp <- transform(results, x, y)
    
    # build the plot
    ggplot(temp, aes(x = M, y = P)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        xlab("log2 FC") +
        ylab("-log10 FDR") +
        ggtitle(str_c(title, " Volcano Plot"))
}

# make the DE plots
MA_plots(results, "ave_HeH", "ave_ETV6", "HeH vs ETV6/RUNX1")
scatter_plots(results, "ave_HeH", "ave_ETV6", "HeH vs ETV6/RUNX1")
volcano_plot(results, "ave_HeH", "ave_ETV6", "HeH vs ETV6/RUNX1")

# ============== individual protein expression plots ===========================

# function to extract the identifier part of the accesssion
get_identifier <- function(accession) {
    identifier <- str_split(accession, "\\|", simplify = TRUE)
    identifier[,3]
}

set_plot_dimensions <- function(width_choice, height_choice) {
    options(repr.plot.width=width_choice, repr.plot.height=height_choice)
}

plot_top_tags <- function(results, nleft, nright, top_tags) {
    # results should have data first, then test results (two condition summary table)
    # nleft, nright are number of data points in each condition
    # top_tags is number of up and number of down top DE candidates to plot
    # get top ipregulated
    up <- results %>% 
        filter(logFC >= 0) %>%
        arrange(FDR)
    up <- up[1:top_tags, ]
    
    # get top down regulated
    down <- results %>% 
        filter(logFC < 0) %>%
        arrange(FDR)
    down <- down[1:top_tags, ]
    
    # pack them into one data frame
    proteins <- rbind(up, down)
        
    color = c(rep("red", nleft), rep("blue", nright))
    for (row_num in 1:nrow(proteins)) {
        row <- proteins[row_num, ]
        vec <- as.vector(unlist(row[1:(nleft + nright)]))
        names(vec) <- colnames(row[1:(nleft + nright)])
        title <- str_c(get_identifier(row$Acc), ", int: ", scientific(mean(vec), 2), 
                       ", p-val: ", scientific(row$FDR, digits = 3), 
                       ", FC: ", round(row$FC, digits = 1))
        barplot(vec, col = color, main = title)
    }    
}


# set plot size, make plots, reset plot size
set_plot_dimensions(6, 4)                      
plot_top_tags(results, length(HeH), length(ETV6), 25)
set_plot_dimensions(width, height)

write.table(results, "IRS_R_pools_results.txt", sep = "\t",
           row.names = FALSE, na =  " ")

sessionInfo()


