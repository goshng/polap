# 1) One mapping only (reads vs reads) — make once.
# HiFi: -x ava-pb ; ONT: -x ava-ont
minimap2 -x ava-pb -t 32 --secondary=yes -N 50 --mask-level 0.50 \
	reads.fq.gz reads.fq.gz >allvsall.paf

# 2) Provide your seed → eads grouping (seed_groups.tsv):
# seedA   read123
# seedA   read987
# seedB   read555
# ...

# 3) Run the sampler
./polap-py-seed-bridge-sampler.py allvsall.paf seed_groups.tsv \
	-o bridge_out \
	--min-alen 800 --min-ident 0.88 --min-weight 0.12 \
	--walks-per-pair 200 --max-hops 30 \
	--temperature 0.25 --min-edge-w-walk 0.12 --simple --restart 0.02 \
	--seed 2025 --min-count 10
