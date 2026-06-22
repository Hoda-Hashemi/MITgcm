#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
repo_dir=$(cd "$script_dir/../../.." && pwd -P)

exec "$repo_dir/jobs/clean_run.sh" "$@"
