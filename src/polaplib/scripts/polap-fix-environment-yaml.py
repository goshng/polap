#!/usr/bin/env python3
# Version: v0.2.0
import sys, yaml, os

if len(sys.argv) < 3:
    print("Usage: fix-yaml <environment.yml> <env_name> [channels...]", file=sys.stderr)
    sys.exit(2)
yml, env_name = sys.argv[1], sys.argv[2]
channels = sys.argv[3].split() if len(sys.argv) > 3 else ["conda-forge", "bioconda"]
with open(yml) as f:
    data = yaml.safe_load(f)
if not isinstance(data, dict):
    data = {}
data["name"] = env_name
if not data.get("channels"):
    data["channels"] = channels
# remove empty dependencies
deps = data.get("dependencies", [])
data["dependencies"] = [d for d in deps if d]
tmp = yml + ".tmp"
with open(tmp, "w") as f:
    yaml.safe_dump(data, f, sort_keys=False)
os.replace(tmp, yml)
