HEADER = "#OTU ID"
FORBIDDEN_IN_SAMPLE = set(";= \t")


def read_fasta(path):
    label, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if label is not None:
                    yield label, "".join(chunks)
                label, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if label is not None:
        yield label, "".join(chunks)


def write_fasta(records, path):
    with open(path, "w") as f:
        for label, seq in records:
            f.write(f">{label}\n{seq}\n")


def seq_id(label):
    """'Uniq1;size=10;' -> 'Uniq1'."""
    fields = label.split()
    return fields[0].split(";")[0] if fields else ""


def label_size(label, default=1):
    """The `size=` annotation of a label, or `default` if it has none."""
    fields = label.split()
    for field in fields[0].split(";")[1:] if fields else []:
        if field.startswith("size="):
            return int(field[5:])
    return default


def members(fasta_path, sample, sizes=None):
    """Records of `fasta_path` labelled "<id>;size=<n>;sample=<sample>;"."""
    if not sample or FORBIDDEN_IN_SAMPLE & set(sample):
        raise ValueError(
            f"Sample name {sample!r} can't be used in sequence labels;"
            " it must not be empty or contain ';', '=' or whitespace."
        )
    for label, seq in read_fasta(fasta_path):
        sid = seq_id(label)
        n = sizes.get(sid, 0) if sizes is not None else label_size(label)
        if n > 0:
            yield f"{sid};size={n};sample={sample};", seq


def read_table(path):
    """Return (samples, {seq_id: [count per sample]})."""
    with open(path) as f:
        header = f.readline().rstrip("\r\n").split("\t")
        if header[0] != HEADER:
            raise ValueError(f"Unexpected OTU table header in {path}")
        samples = header[1:]
        rows = {}
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if fields[0]:
                rows[fields[0]] = [int(x) for x in fields[1:]]
    return samples, rows
