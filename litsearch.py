"""
These are PubMed search results for 'escherichia coli sequence type' from 2000 to 2021, 
with the first 10,000 of 14,383 results saved to file 'escherichi-set.csv'.

"""

#TODO: make webpage as UI using flask?

import re
from collections import Counter
from sys import argv

import numpy as np
import pandas as pd
from Bio import Entrez


class PubMed:

    @staticmethod
    def long_summary(pmid: str):
        print(Entrez.efetch(db="pubmed", id=pmid, retmode="text", rettype="gb").read())


    @staticmethod
    def times_cited(pmid: [str, list]) -> list:
        result = Entrez.read(Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pubmed_citedin"))
        
        if isinstance(pmid, str):       
            papers = result[0]['LinkSetDb'][0]['Link']
            print("Cited by ", len(papers), " papers.")
            return papers

        elif isinstance(pmid, list):
            citation_count = []
            for index, _ in enumerate(pmid):
                try:
                    citation_count.append(len(result[index]['LinkSetDb'][0]['Link']))
                except IndexError:
                    citation_count.append(0)
            return citation_count
        
        else:
            raise TypeError("Expected string or list.")


    @staticmethod
    def add_times_cited(df) -> pd.DataFrame:
        citation_count = PubMed.times_cited(list(df.index))
        timescited_df = pd.DataFrame(citation_count, index=df.index, columns=["Times cited"])
        df = df.merge(timescited_df, left_index=True, right_index=True)
        df = df.sort_values(by=["Times cited"], ascending=False)
        return df


def uniformST(st: str) -> str:
    """Transform sequence type string into ST###"""
    return re.sub("ST |[Ss]equence [Tt]ype ", "ST", st)


def find_terms(text: str) -> list:
    result = re.findall(r"ST\d+|ST \d+|[Ss]equence [Tt]ype \d+", text)
    return result


def countST(df: pd.DataFrame) -> Counter:
    """Find occurences of sequence types as a Counter dict."""
    STcount = Counter()
    for STs in df["ST"]:
        for st in STs.split(", "):
            STcount[st] += 1
    return STcount


def frequent_journals(df, top=15):
    print(df["Journal/Book"].value_counts().head(top))


def frequent_pub_years(df, top=15):
    print(df["Publication Year"].value_counts().head(top))


def frequent_ST(df, top=15):
    result = countST(df)
    return result.most_common(top)
    

def plotSTasPie(df: pd.DataFrame):
    result = countST(df)
    STdf = pd.DataFrame.from_dict(result, orient="index")
    STdf.loc["Other"] = [len(STdf[STdf == 1].dropna())]
    STdf.loc[STdf[0] != 1].plot.pie(y=0)


def search(df, search: str, top=10) -> pd.DataFrame:
    """Search for a given sequence type (as 'ST###')"""
    b = df["ST"].str.findall(f"{search},|{search}$")
    b = b[b.astype(str) != '[]']
    result = df.loc[b.index]
    if isinstance(top, int):
        return result.head(top)
    else:
        return result


def ask_email():
    try:
        Entrez.email = argv[3]
    except IndexError:
        Entrez.email = input("Please provide email address for PubMed access: ")


def main(csv_file):
    lit = pd.read_csv(csv_file, index_col=0)                        # make DataFrame
    lit["ST"] = lit["Title"].apply(find_terms)                      # Identify sequence type from title
    lit["ST"] = lit["ST"].apply(lambda x : uniformST(", ".join(x))) # Shorten "sequence type" to ST
    lit["ST"] = lit["ST"].replace("", np.NaN)                       # Identify empty values in ST col
    lit = lit[lit["ST"].notna()]                                    # Filter out rows without mention of ST
    lit = PubMed.add_times_cited(lit)                               # Fetch times each article has been cited
    return lit


if __name__ == "__main__":

    if argv[2].endswith(".csv"):
        fp = argv[2]
        skip_arg_two = False
    else:
        fp = "escherichi-set.csv"
        skip_arg_two = True


    if argv[1] == "new":
        ask_email()
        df = main(fp)
        df.to_csv(fp)

    elif argv[1] == "find":
        df = pd.read_csv(fp, index_col=0)
        if not skip_arg_two:
            term = argv[3]
        else:
            term = argv[2]
        print(search(df, term))


    elif argv[1] == "summary":
        ask_email()
        if not skip_arg_two:
            pmid = argv[3]
        else:
            pmid = argv[2]
        PubMed.long_summary(pmid)

    

    
