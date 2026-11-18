"""Step 1 - Download GSE266616 per-sample 10x matrices straight from NCBI GEO and verify them.

For the 33 Crohn's disease inflamed / non-inflamed samples (sample_metadata_inf_vs_non.csv) this script
  1. downloads <GSM>_<name>.tar.gz from the GEO FTP (over https), with retries;
  2. checks the file size against GEO's own listing (filelist.csv from the series /suppl/ folder);
  3. reads every archive end-to-end (gzip CRC + tar structure) and unpacks the 10x triplet;
  4. checks matrix.mtx.gz header dimensions against the number of lines in features.tsv.gz / barcodes.tsv.gz;
  5. optionally (--audit-dir) compares SHA-256 with a pre-existing copy of the same files, to test whether
     that copy is byte-identical to what GEO serves today.
Output: <outdir>/raw/*.tar.gz, <outdir>/mtx/<sample_folder>/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz},
        <outdir>/download_audit.csv
Stdlib only, so it runs on any Python 3.
"""
import argparse, csv, gzip, hashlib, os, sys, tarfile, time, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

ap = argparse.ArgumentParser()
ap.add_argument("--meta", default="sample_metadata_inf_vs_non.csv")
ap.add_argument("--filelist", default="filelist.csv")          # GEO /suppl/filelist.txt converted to csv (Name, Size)
ap.add_argument("--outdir", default="data")
ap.add_argument("--audit-dir", default=None, help="folder with pre-existing copies of the tarballs to compare by SHA-256")
ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()

def rd(path):
    with open(path, newline="") as f: return list(csv.DictReader(f))

meta = rd(a.meta)
geo_size = {r["Name"]: int(r["Size"]) for r in rd(a.filelist)}
raw_dir, mtx_dir = os.path.join(a.outdir, "raw"), os.path.join(a.outdir, "mtx")
os.makedirs(raw_dir, exist_ok=True); os.makedirs(mtx_dir, exist_ok=True)

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(8 << 20), b""): h.update(ch)
    return h.hexdigest()

def download(url, dst, expect, tries=5):
    if os.path.exists(dst) and os.path.getsize(dst) == expect: return "cached"
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
            with urllib.request.urlopen(req, timeout=180) as r, open(dst + ".part", "wb") as f:
                while True:
                    ch = r.read(1 << 20)
                    if not ch: break
                    f.write(ch)
            os.replace(dst + ".part", dst); return "downloaded"
        except Exception as e:
            err = repr(e); time.sleep(3 * (i + 1))
    return f"FAILED: {err}"

def count_lines_gz(p):
    n = 0
    with gzip.open(p, "rt") as f:
        for _ in f: n += 1
    return n

def process(r):
    name = r["tarball"]; gsm = r["gsm"]; folder = r["sample_folder"]
    url = r["suppl"].replace("ftp://", "https://"); dst = os.path.join(raw_dir, name)
    out = dict(gsm=gsm, sample=folder, geo_listed_size=geo_size.get(name), status="")
    out["status"] = download(url, dst, geo_size.get(name, -1))
    if out["status"].startswith("FAILED"): return out
    out["fresh_size"] = os.path.getsize(dst); out["size_matches_geo"] = out["fresh_size"] == out["geo_listed_size"]
    # --- archive integrity + unpack
    d = os.path.join(mtx_dir, folder); os.makedirs(d, exist_ok=True)
    try:
        with tarfile.open(dst, "r:gz") as t:
            for m in t.getmembers():
                if m.isfile():
                    m.name = os.path.basename(m.name); t.extract(m, d)
        out["archive_ok"] = True
    except Exception as e:
        out["archive_ok"] = False; out["archive_err"] = repr(e); return out
    # --- matrix dimensions vs barcode / feature files
    try:
        with gzip.open(os.path.join(d, "matrix.mtx.gz"), "rt") as f:
            line = f.readline()
            while line.startswith("%"): line = f.readline()
            nr, nc, nnz = map(int, line.split())
        nfeat = count_lines_gz(os.path.join(d, "features.tsv.gz")); nbc = count_lines_gz(os.path.join(d, "barcodes.tsv.gz"))
        out.update(mtx_rows=nr, mtx_cols=nc, nnz=nnz, n_features_file=nfeat, n_barcodes_file=nbc,
                   dims_consistent=(nr == nfeat and nc == nbc))
    except Exception as e:
        out["dims_consistent"] = False; out["dims_err"] = repr(e)
    # --- byte-for-byte comparison with an existing copy
    if a.audit_dir:
        old = os.path.join(a.audit_dir, name)
        out["sha256_fresh"] = sha256(dst)
        if os.path.exists(old):
            out["sha256_existing"] = sha256(old); out["identical_to_existing"] = out["sha256_fresh"] == out["sha256_existing"]
        else:
            out["identical_to_existing"] = None
    return out

rows = []
with ThreadPoolExecutor(a.threads) as ex:
    futs = {ex.submit(process, r): r["gsm"] for r in meta}
    for f in as_completed(futs):
        o = f.result(); rows.append(o); print({k: o.get(k) for k in ("gsm", "sample", "status", "size_matches_geo", "archive_ok", "dims_consistent", "identical_to_existing")}, flush=True)

cols = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("gsm", "sample"), k))
with open(os.path.join(a.outdir, "download_audit.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); [w.writerow(r) for r in sorted(rows, key=lambda r: r["gsm"])]
n = len(rows)
print(f"SUMMARY n={n} downloaded_ok={sum(not r['status'].startswith('FAILED') for r in rows)} "
      f"size_ok={sum(bool(r.get('size_matches_geo')) for r in rows)} archive_ok={sum(bool(r.get('archive_ok')) for r in rows)} "
      f"dims_ok={sum(bool(r.get('dims_consistent')) for r in rows)} identical_to_existing={sum(r.get('identical_to_existing') is True for r in rows)}")
