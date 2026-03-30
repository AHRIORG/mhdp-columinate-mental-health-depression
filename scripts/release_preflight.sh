#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(git rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -z "${ROOT_DIR}" ]]; then
  echo "ERROR: run this inside the git repository."
  exit 2
fi

cd "${ROOT_DIR}"

usage() {
  cat <<'EOF'
Usage:
  scripts/release_preflight.sh --staged
  scripts/release_preflight.sh --all
  scripts/release_preflight.sh <path> [<path> ...]

Checks for release blockers:
  - blocked data/object extensions
  - .DS_Store and _quarantine paths
  - high-confidence secret patterns
  - private machine path markers in non-doc files
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

tmpfile="$(mktemp)"
trap 'rm -f "${tmpfile}"' EXIT

collect_paths() {
  if [[ "${#@}" -eq 0 || "${1:-}" == "--staged" ]]; then
    git diff --cached --name-only
    return
  fi

  if [[ "${1:-}" == "--all" ]]; then
    git ls-files
    git ls-files --others --exclude-standard
    return
  fi

  for path in "$@"; do
    if [[ -e "${path}" ]]; then
      if [[ -d "${path}" ]]; then
        find "${path}" -type f
      else
        printf '%s\n' "${path}"
      fi
    else
      printf '%s\n' "${path}"
    fi
  done
}

collect_paths "$@" | sed '/^$/d' | sort -u > "${tmpfile}"

if [[ ! -s "${tmpfile}" ]]; then
  echo "No files selected for preflight scan."
  echo "Tip: stage files first, or run with --all / explicit paths."
  exit 1
fi

echo "Preflight files:"
cat "${tmpfile}"
echo

fail=0
ALLOWED_RDS_PATH="scripts/obj00/data_management/data_examples/1_staging_snippets/derived_models/irt_joint_models.rds"

echo "Check 1/5: blocked object extensions..."
blocked_object_files=()
while IFS= read -r f; do
  [[ -z "${f}" ]] && continue
  blocked_object_files+=("${f}")
done < <(rg '\.(RData|rda|rds|qs|fst|feather|parquet|pkl|joblib|sqlite)$' "${tmpfile}" || true)
if [[ "${#blocked_object_files[@]}" -gt 0 ]]; then
  disallowed_object_files=()
  for f in "${blocked_object_files[@]}"; do
    if [[ "${f}" == "${ALLOWED_RDS_PATH}" ]]; then
      continue
    fi
    disallowed_object_files+=("${f}")
  done

  if [[ "${#disallowed_object_files[@]}" -gt 0 ]]; then
    printf '%s\n' "${disallowed_object_files[@]}"
    echo "FAIL: blocked data/object extension found."
    fail=1
  fi
fi
echo

echo "Check 2/5: blocked housekeeping and quarantine paths..."
if rg -n '(^|/)\.DS_Store$' "${tmpfile}"; then
  echo "FAIL: .DS_Store file detected."
  fail=1
fi
if rg -n '(^|/)_quarantine(/|$)' "${tmpfile}"; then
  echo "FAIL: _quarantine path detected."
  fail=1
fi
echo

existing_files=()
while IFS= read -r f; do
  if [[ -f "${f}" ]]; then
    existing_files+=("${f}")
  fi
done < "${tmpfile}"
if [[ "${#existing_files[@]}" -eq 0 ]]; then
  echo "No existing files to content-scan (delete-only or missing paths)."
  if [[ "${fail}" -ne 0 ]]; then
    exit 1
  fi
  echo "PASS: path checks complete."
  exit 0
fi

echo "Check 3/5: high-confidence secret patterns..."
if rg -n -I \
  -e 'ghp_[A-Za-z0-9]{36,}' \
  -e 'github_pat_[A-Za-z0-9_]{20,}' \
  -e 'AKIA[0-9A-Z]{16}' \
  -e '-----BEGIN (RSA|EC|OPENSSH|DSA|PRIVATE) KEY-----' \
  "${existing_files[@]}"; then
  echo "FAIL: potential secret detected."
  fail=1
fi
echo

echo "Check 4/5: private path markers in non-doc files..."
non_doc_files=()
for f in "${existing_files[@]}"; do
  if [[ "${f}" == docs/* || "${f}" == "README.md" ]]; then
    continue
  fi
  if [[ "${f}" == ".github/workflows/release-guard.yml" || "${f}" == "scripts/release_preflight.sh" ]]; then
    continue
  fi
  non_doc_files+=("${f}")
done
if [[ "${#non_doc_files[@]}" -gt 0 ]]; then
  if rg -n -I \
    -e '/Users/' \
    -e 'OneDrive' \
    -e '_private_use' \
    -e 'file://' \
    -e '[A-Za-z]:\\' \
    "${non_doc_files[@]}"; then
    echo "FAIL: machine-specific/private path marker detected in non-doc files."
    fail=1
  fi
else
  echo "No non-doc files selected."
fi
echo

echo "Check 5/5: website-render spot check for local path leaks..."
if [[ -d website ]]; then
  if rg -n -I -g '*.html' -e '/Users/|OneDrive|_private_use|file://' website; then
    echo "FAIL: local path marker found in rendered website HTML."
    fail=1
  fi
else
  echo "No website directory present; skipped."
fi
echo

if [[ "${fail}" -ne 0 ]]; then
  echo "Release preflight FAILED."
  exit 1
fi

echo "Release preflight PASSED."
