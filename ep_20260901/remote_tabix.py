"""Minimal pure-Python remote tabix client for EBI GWAS Catalog harmonised files.
Reads .tbi index, fetches BGZF byte-ranges over HTTP, decompresses with zlib.
No pysam / htslib required.
"""
import gzip
import io
import json
import struct
import time
import urllib.request

CHUNK = 65536  # max BGZF compressed block size


def _http_get(url, start=None, end=None, retries=4, timeout=90):
    """Return bytes for url via curl; optional inclusive byte range [start, end]."""
    import subprocess

    last = None
    for attempt in range(retries):
        cmd = ["curl", "-sS", "--max-time", str(timeout)]
        if start is not None:
            cmd += ["-r", f"{start}-{end if end is not None else ''}"]
        cmd.append(url)
        try:
            r = subprocess.run(cmd, capture_output=True, timeout=timeout + 15)
            if r.returncode in (0, 28) and r.stdout:  # 28 = our max-time cut
                return r.stdout
            if r.returncode == 0 and not r.stdout and start is None:
                return b""
            last = f"curl rc={r.returncode} err={r.stderr.decode()[:200]}"
        except Exception as e:  # noqa: BLE001
            last = e
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"GET failed {url} range=({start},{end}): {last}")


def bgzf_decompress(buf):
    """Decompress a byte string of concatenated BGZF/gzip members."""
    out = bytearray()
    pos = 0
    import zlib
    n = len(buf)
    while pos < n:
        d = zlib.decompressobj(31)
        try:
            out += d.decompress(buf[pos:])
        except zlib.error:
            break
        consumed = n - pos - len(d.unused_data)
        if consumed <= 0:
            break
        pos += consumed
    return bytes(out)


class RemoteTabix:
    def __init__(self, tsv_url, tbi_path):
        self.tsv_url = tsv_url
        raw = bgzf_decompress(open(tbi_path, "rb").read())
        self._parse_index(raw)

    def _parse_index(self, raw):
        assert raw[:4] == b"TBI\x01", "not a tabix index"
        (n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm) = struct.unpack(
            "<8i", raw[4:36]
        )
        self.format, self.col_seq, self.col_beg = fmt, col_seq, col_beg
        self.meta, self.skip = chr(meta), skip
        names = raw[36 : 36 + l_nm].decode().split("\x00")
        self.ref_names = [x for x in names if x]
        off = 36 + l_nm
        self.refs = {}
        for i in range(n_ref):
            (n_bin,) = struct.unpack("<i", raw[off : off + 4])
            off += 4
            bins = {}
            for _ in range(n_bin):
                bin_id, n_chunk = struct.unpack("<Ii", raw[off : off + 8])
                off += 8
                chunks = []
                for _ in range(n_chunk):
                    cb, ce = struct.unpack("<QQ", raw[off : off + 16])
                    off += 16
                    chunks.append((cb, ce))
                bins[bin_id] = chunks
            (n_intv,) = struct.unpack("<i", raw[off : off + 4])
            off += 4 + 8 * n_intv  # skip linear index (bin-based query only)
            self.refs[self.ref_names[i]] = bins

    @staticmethod
    def reg2bins(beg, end):
        """Candidate bins for zero-based half-open [beg, end)."""
        end -= 1
        bins = [0]
        bins += range(1 + (beg >> 26), 1 + (end >> 26) + 1)
        bins += range(9 + (beg >> 23), 9 + (end >> 23) + 1)
        bins += range(73 + (beg >> 20), 73 + (end >> 20) + 1)
        bins += range(585 + (beg >> 17), 585 + (end >> 17) + 1)
        bins += range(4681 + (beg >> 14), 4681 + (end >> 14) + 1)
        return bins

    def _fetch_text(self, chunks):
        """Fetch compressed ranges, decompress, return (start_voff, text)."""
        # merge compressed byte ranges
        spans = sorted((cb >> 16, (ce >> 16) + CHUNK) for cb, ce in chunks)
        merged = []
        for s, e in spans:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        min_voff = min(cb for cb, _ in chunks)
        pieces = []
        for s, e in merged:
            pieces.append(bgzf_decompress(_http_get(self.tsv_url, s, e - 1)))
        return min_voff, b"".join(pieces)

    def query(self, chrom, beg, end):
        """Yield parsed field lists for lines in [beg, end) (0-based half-open)."""
        if chrom not in self.refs:
            return []
        cand = set(self.reg2bins(beg, end))
        chunks = []
        for b in cand:
            chunks.extend(self.refs[chrom].get(b, []))
        if not chunks:
            return []
        min_voff, text = self._fetch_text(chunks)
        start_inner = min_voff & 0xFFFF
        data = text[start_inner:]
        rows = []
        for line in data.decode(errors="replace").split("\n"):
            if not line or line[0] == self.meta:
                continue
            f = line.split("\t")
            try:
                pos0 = int(f[self.col_beg - 1]) - 1  # generic tabix: 1-based col
            except (IndexError, ValueError):
                continue
            if beg <= pos0 < end and f[self.col_seq - 1] == chrom:
                rows.append(f)
        return rows


def get_header(tsv_url, meta="#"):
    """Fetch first BGZF block and return header column names."""
    text = bgzf_decompress(_http_get(tsv_url, 0, CHUNK - 1)).decode(errors="replace")
    for line in text.split("\n"):
        if line.startswith(meta):
            return line.lstrip(meta).strip().split("\t")
        return line.strip().split("\t")  # no comment char: first line is header
    return []


if __name__ == "__main__":
    import sys

    print(json.dumps(get_header(sys.argv[1])))
