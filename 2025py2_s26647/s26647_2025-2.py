#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key, min_len, max_len):
        Entrez.email, Entrez.api_key, Entrez.tool = email, api_key, 'BioScriptEx10'
        self.min_len, self.max_len = min_len, max_len

    def search(self, taxid):
        term = f"txid{taxid}[Organism]"
        handle = Entrez.esearch(db="nucleotide", term=term, usehistory="y", retmax=0)
        res = Entrez.read(handle)
        self.env, self.key, self.count = res["WebEnv"], res["QueryKey"], int(res["Count"])
        return self.count

    def fetch_and_filter(self, max_records=200):
        records = []
        step = 100
        for start in range(0, min(self.count, max_records), step):
            handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                                   retstart=start, retmax=step,
                                   webenv=self.env, query_key=self.key)
            for r in SeqIO.parse(handle, "gb"):
                l = len(r.seq)
                if self.min_len <= l <= self.max_len:
                    records.append((r.id, l, r.description))
        return records

    def save_csv(self, data, filename):
        df = pd.DataFrame(data, columns=["Accession", "Length", "Description"])
        df.sort_values("Length", ascending=False, inplace=True)
        df.to_csv(filename, index=False)
        return df

    def plot_lengths(self, df, output):
        plt.figure(figsize=(10, 6))
        plt.plot(df["Accession"], df["Length"], marker='o')
        plt.xticks(rotation=90, fontsize=8)
        plt.ylabel("Sequence Length")
        plt.title("GenBank Sequence Lengths")
        plt.tight_layout()
        plt.savefig(output)

def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--email", required=True)
    p.add_argument("--api_key", required=True)
    p.add_argument("--taxid", required=True)
    p.add_argument("--min_len", type=int, default=0)
    p.add_argument("--max_len", type=int, default=1e6)
    p.add_argument("--limit", type=int, default=200)
    args = p.parse_args()

    r = NCBIRetriever(args.email, args.api_key, args.min_len, args.max_len)
    print(f"Searching TaxID {args.taxid}...")
    if not r.search(args.taxid):
        print("No records found.")
        return
    print("Fetching & filtering...")
    data = r.fetch_and_filter(args.limit)
    if not data:
        print("No records within length range.")
        return
    csv_file = f"taxid_{args.taxid}_filtered.csv"
    png_file = f"taxid_{args.taxid}_plot.png"
    df = r.save_csv(data, csv_file)
    r.plot_lengths(df, png_file)
    print(f"Saved CSV: {csv_file}")
    print(f"Saved plot: {png_file}")

if __name__ == "__main__":
    main()
