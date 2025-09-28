library(readr)
library(dplyr)
library(stringr)

# ---- Step 1: Load accession + length ----
acc_len <- read_tsv("acc_len.tsv", col_names = c("accession", "length"), col_types = "ci") %>%
  mutate(accession_short = str_replace(accession, "\\.\\d+$", ""))

# ---- Step 2: Load accession + taxid ----
acc_taxid <- read_tsv("acc_taxid.tsv", col_names = c("accession", "taxid"), col_types = "cc")

# ---- Step 3: Merge length + taxid via short accession ----
acc_info <- acc_len %>%
  inner_join(acc_taxid, by = c("accession_short" = "accession")) %>%
  select(accession, length, taxid)

# Load lineage.tsv with 3 columns (taxid, names, taxids)
lineage <- read_tsv("lineage.tsv", col_names = c("taxid", "lineage_names", "lineage_taxids"), col_types = "ccc")



lineage_parsed <- lineage %>%
  mutate(
    lineage_split = str_split(lineage_names, ";\\s*"),

    # Find last binomial name (Genus species) as species
    species = vapply(lineage_split, function(x) {
      hits <- x[str_detect(x, "^[A-Z][a-z]+ [a-z]+$")]
      if (length(hits) > 0) tail(hits, 1) else NA_character_
    }, character(1)),

    # Genus is first word of species (e.g., Arabidopsis thaliana → Arabidopsis)
    genus = word(species, 1),

    # Extract first matching family/order
    family = vapply(lineage_split, function(x) {
      fx <- x[str_detect(x, "aceae$")]
      if (length(fx) > 0) fx[1] else NA_character_
    }, character(1)),

    order = vapply(lineage_split, function(x) {
      ox <- x[str_detect(x, "ales$")]
      if (length(ox) > 0) ox[1] else NA_character_
    }, character(1))
  ) %>%
  select(taxid, species, genus, family, order)

# Save to TSV
write_tsv(lineage_parsed, "lineage_parsed.tsv")

# ---- Step 5: Join all tables ----
final <- acc_len %>%
  mutate(acc_short = str_replace(accession, "\\.\\d+$", "")) %>%
  inner_join(acc_taxid, by = c("acc_short" = "accession")) %>%
  inner_join(lineage_parsed, by = "taxid") %>%
  select(accession, length, taxid, species, genus, family, order)


# ---- Step 6: Deduplicate by genus, keeping the longest sequence ----
final_dedup <- final %>%
  group_by(genus) %>%
  arrange(desc(length)) %>%
  slice(1) %>%
  ungroup()

# ---- Step 7: Write to file ----
write_tsv(final_dedup, "taxonomy_deduplicated.tsv")

