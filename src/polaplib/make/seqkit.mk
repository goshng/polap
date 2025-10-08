# polaplib/make/seqkit.mk
FASTQ_LIST    ?= $(shell git ls-files '*/*/summary-data/l.fq*' 2>/dev/null)
DATASET_TSV   ?= md/dataset-summary.tsv

.PHONY: dataset.summary
dataset.summary: $(DATASET_TSV)

$(DATASET_TSV): $(FASTQ_LIST)
	@command -v seqkit >/dev/null 2>&1 || { echo "[ERR] seqkit not found"; exit 2; }
	@mkdir -p $(dir $@)
	@echo -e "species\ttotal_bases\tread_count\tmean_length\tN50\tavg_qual\tgc_content" > "$@"
	@for fq in $(FASTQ_LIST); do \
	  sp=$$(echo "$$fq" | awk -F/ '{print $$1}'); \
	  line=$$(seqkit stats -T -a "$$fq" | tail -n +2); \
	  IFS=$$'\t' read -r f fmt typ num sum min avg max n50 q20 q30 gc avgQ <<<"$$line"; \
	  echo -e "$$sp\t$$sum\t$$num\t$$avg\t$$n50\t$$avgQ\t$$gc" >> "$@"; \
	done
