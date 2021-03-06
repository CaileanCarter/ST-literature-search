import re
from sys import argv

import numpy as np
import pandas as pd

from litsearch import PubMed, ask_email


def find_terms(text):
    p = re.compile(r"Salmonella enterica [Ss]erovar (\w+)")
    try:
        n = p.search(text).group(1)
        return n
    except AttributeError:
        return np.NaN


def search(df, search: str, top=10) -> pd.DataFrame:
    """Search for a given serovar"""
    first_letter, remaining = search[0], search[1:]
    match = f"[{first_letter.upper()}{first_letter.lower()}]{remaining}"
    b = df["Serovar"].str.findall(f"{match}")
    b = b[b.astype(str) != '[]']
    result = df[["Title", "Authors", "Journal/Book", "Publication Year", "Serovar", "Times cited"]]
    result = result.loc[b.index]
    if isinstance(top, int):
        return result.head(top)
    else:
        return result


def Serovars(df, top=10):
    print(df["Serovar"].value_counts()[:top])


def main(csv_file):

    lit = pd.read_csv(csv_file, index_col=0)                        # make DataFrame
    lit = lit.drop(["PMCID", "NIHMS ID", "First Author"], axis=1)
    lit["Serovar"] = lit["Title"].apply(find_terms)                      # Identify sequence type from title
    lit = lit[lit["Serovar"].notna()]                                    # Filter out rows without mention of ST
    lit = PubMed.add_times_cited(lit)

    return lit


if __name__ == "__main__":

    if argv[1] == "new":
        fp = argv[2]
        ask_email()
        df = main(fp)
        df.to_csv(fp)


    elif argv[1] == "find":
        fp = argv[2]
        df = pd.read_csv(fp, index_col=0)
        term = argv[3]
        print(search(df, term))


    elif argv[1] == "summary":
        ask_email()
        pmid = argv[2]
        PubMed.long_summary(pmid)

    
    elif argv[1] == "serovars":
        fp = argv[2]
        df = pd.read_csv(fp, index_col=0)

    else:
        raise TypeError("No recognised argument given for first positional argument.")