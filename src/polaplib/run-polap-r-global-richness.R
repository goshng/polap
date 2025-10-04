#!/usr/bin/env Rscript

################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("taxize"))
suppressPackageStartupMessages(library("rgbif"))
suppressPackageStartupMessages(library("ggplot2"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

# https://chatgpt.com/share/6739a37d-0c7c-800e-9473-f748e3c5acee

# Parse the main subcommand
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  subcommand <- "sample"
  subcommand_args <- c()
  # stop("No subcommand provided. Use 'analyze' or 'report'.")
} else {
  subcommand <- args[1]
  subcommand_args <- args[-1] # Remaining arguments for the subcommand
}


# Define option lists for each subcommand
option_list_analyze <- list(
  make_option(c("-t", "--taxon"),
    action = "store",
    help = "Organelle annotation table",
    metavar = "<TAXON>"
  ),
  make_option(c("-o", "--output"),
    action = "store",
    help = "Output contig seeds filename"
  )
)

option_list_report <- list(
  make_option(c("-r", "--report"),
    type = "character",
    help = "Report file"
  ),
  make_option(c("-s", "--summary"),
    action = "store",
    help = "Include summary"
  ),
  make_option(c("-o", "--output"),
    action = "store",
    help = "Output contig seeds filename"
  )
)

option_list_sample <- list(
  make_option(c("-r", "--richness"),
    type = "character",
    help = "Species richness file"
  ),
  make_option(c("-s", "--seed"),
    type = "integer",
    action = "store",
    help = "Random number seed"
  ),
  make_option(c("-t", "--taxonomy"),
    action = "store",
    help = "Include summary"
  ),
  make_option(c("-o", "--output"),
    action = "store",
    help = "Output contig seeds filename"
  )
)

# Use `switch` to handle subcommands
parser <- switch(subcommand,
  "analyze" = OptionParser(option_list = option_list_analyze),
  "report" = OptionParser(option_list = option_list_report),
  "sample" = OptionParser(option_list = option_list_sample),
  stop(paste("Unknown subcommand:", subcommand))
)

# Parse options
options <- parse_args(parser, args = subcommand_args)

# Execute the corresponding logic
switch(subcommand,
  "analyze" = {
    # cat("Running 'analyze' with the following options:\n")
    # print(options)
    if (is_null(options$output)) {
      input_dir0 <- file.path(".")
      input1 <- file.path("Magnoliopsida")
      output1 <- file.path(input_dir0, "out.txt")
      options <- parse_args(parser, args = c("--taxon", input1, "-o", output1))
    }
    print(options)

    # Add analyze-specific logic here
    # Search for the taxon
    taxon <- name_backbone(name = args1$taxon)

    # Retrieve children taxa (e.g., families)
    children <- name_usage(key = taxon$usageKey, data = "children")

    # Initialize a data frame for results
    species_richness <- data.frame(
      Family = character(),
      TotalSpecies = integer()
    )

    # Loop through families and count occurrences
    for (child in children$data$canonicalName) {
      # Get the taxon key for the child
      child_info <- name_backbone(name = child)

      # Retrieve the occurrence count
      if (!is.null(child_info$usageKey)) {
        occurrences <- occ_search(taxonKey = child_info$usageKey, limit = 0)
        occurrence_count <- ifelse(!is.null(occurrences$meta$count),
          occurrences$meta$count, 0
        )

        # Only add to the data frame if occurrences are valid
        species_richness <- rbind(
          species_richness,
          data.frame(
            Family = child,
            TotalOccurrences = occurrence_count
          )
        )

        if (occurrence_count > 0) {
        } else {
          cat("No occurrences found for:", child, "\n")
        }
      } else {
        cat("No usageKey for:", child, "\n")
      }
    }

    # Save the results to a CSV file
    write_csv(species_richness, "species_richness.csv")
  },
  "sample" = {
    # cat("Running 'analyze' with the following options:\n")
    # print(options)
    if (is_null(options$output)) {
      input_dir0 <- file.path(".")
      input1 <- file.path(input_dir0, "taxonomy_ids_from_fasta.txt")
      input2 <- file.path(input_dir0, "species_richness.csv")
      output1 <- file.path(input_dir0, "sampled_accessions.csv")
      options <- parse_args(parser, args = c(
        "-t", input1,
        "-r", input2,
        "-o", output1,
        "-s", 123
      ))
    }
    print(options)

    # Set random seed
    print(typeof(options$seed))
    if (!is_null(options$seed)) {
      set.seed(options$seed)
    }

    # Load species richness data
    species_richness <- read_csv(options$richness)

    # Load taxonomy IDs data
    taxonomy_ids <- read_tsv(options$taxonomy,
      col_names = c(
        "Accession", "TaxonomyID", "Species",
        "orderColumn", "orderID", "order",
        "familyColumn", "familyID", "family"
      )
    )

    # TaxonomyID and class, order, family

    # Join taxonomy data with order information
    joined_data <- taxonomy_ids %>%
      distinct(TaxonomyID, .keep_all = TRUE)
    # %>%
    #   inner_join(species_richness, by = c("order" = "Family"))

    # Randomly sample two species per order and one accession per species
    sampled_data <- joined_data %>%
      group_by(order) %>%
      sample_n(2, replace = TRUE) %>%
      group_by(TaxonomyID) %>%
      sample_n(1) %>%
      ungroup()

    # Save sampled data
    write_csv(sampled_data, options$output)
  },
  "report" = {
    if (is_null(options$output)) {
      input_dir0 <- file.path(".")
      input1 <- file.path("Magnoliopsida")
      output1 <- file.path(input_dir0, "out.txt")
      options <- parse_args(parser, args = c("--summary", input1, "-o", output1))
    }
    print(options)
    # Add report-specific logic here
    # Example list of NCBI Taxonomy IDs for your mitochondrial genomes
    # taxonomy_ids <- c(3702, 39947, 4577) # Replace with your actual Taxonomy IDs
    taxonomy_ids <- scan("taxa.txt")
    taxonomy_ids <- unique(taxonomy_ids)

    # Get taxonomy classification for each Taxonomy ID
    tax_data <- classification(taxonomy_ids, db = "ncbi")

    # Extract the family-level taxonomy for each ID
    tax_summary <- data.frame(
      TaxID = taxonomy_ids,
      Family = sapply(tax_data, function(x) {
        family_row <- x[x$rank == "order", ]
        if (nrow(family_row) > 0) {
          return(family_row$name)
        } else {
          return(NA)
        }
      })
    )

    # Example global species richness data by family (replace with actual data)
    # global_richness <- data.frame(
    #   Family = c("Brassicaceae", "Poaceae"),
    #   TotalSpecies = c(3500, 12000)
    # )
    # To restore the data
    global_richness <- read_csv("species_richness.csv")

    # Count the number of species with mitochondrial genomes for each family
    tax_summary$SpeciesCount <- 1 # Each row represents one species
    summary_by_family <- aggregate(SpeciesCount ~ Family, data = tax_summary, sum)

    # Merge with global species richness
    comparison <- merge(global_richness, summary_by_family, by = "Family", all.x = TRUE)
    comparison$SpeciesCount[is.na(comparison$SpeciesCount)] <- 0 # Fill missing counts with 0

    # Calculate proportions
    comparison$Proportion <- comparison$SpeciesCount / comparison$TotalSpecies


    # Bar plot showing proportions of explored taxa
    p1 <- ggplot(comparison, aes(x = reorder(Family, Proportion), y = Proportion)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(
        title = "Underexplored Taxa for Plant Mitochondrial Genomes",
        x = "Family",
        y = "Proportion of Species Explored"
      )
    ggsave("p1.pdf", p1)
  }
)
