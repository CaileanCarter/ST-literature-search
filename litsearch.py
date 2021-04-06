"""
These are PubMed search results for 'escherichia coli sequence type' from 2000 to 2021, 
with the first 10,000 of 14,383 results saved to file 'escherichi-set.csv'.

"""

import re
from flashtext import KeywordProcessor
import pandas as pd 
# import csv
from collections import Counter
from Bio import Entrez
import numpy as np 


class PubMed:

    def __init__(self, email=None, API_key=None):
        Entrez.email = email
        Entrez.api_key = API_key


    def long_summary(self, pmid: str):
        print(Entrez.efetch(db="pubmed", id=pmid, retmode="text", rettype="gb").read())


    def citedby(self, pmid: (str|list)) -> list:
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

    
    def add_times_cited(self, df) -> pd.DataFrame:
        citation_count = self.citedby(list(df.index))
        bar = pd.DataFrame(citation_count, index=df.index, columns=["Times cited"])
        df = df.merge(bar, left_index=True, right_index=True)
        df = df.sort_values(by=["Times cited"], ascending=False)
        return df



def uniformST(st: str) -> str:
    return re.sub("ST |[Ss]equence [Tt]ype ", "ST", st)


def countST(df) -> Counter:
    STcount = Counter()
    for STs in df["ST"]:
        for st in STs:
            st = uniformST(st)
            STcount[st] += 1
    # STcount.most_common(20)
    return STcount


def frequent_journals(df, top=15):
    print(df["Journal/Book"].value_counts().head(top))


def frequent_pub_years(df, top=15):
    print(df["Publication Year"].value_counts().head(top))


def find_terms(text: str) -> list:
    result = re.findall(r"ST\d+|ST \d+|[Ss]equence [Tt]ype \d+", text)
    return result


def plotSTaspie(df: Counter):
    STdf = pd.DataFrame.from_dict(df, orient="index")
    STdf.loc["Other"] = [28]
    STdf.loc[STdf[0] != 1].plot.pie(y=0)


def find_ST(df, search: str):
    b = df.str.findall(f"{search},|{search}$")
    return b[b.astype(str) != '[]']


def main(csv_file="escherichi-set.csv"):

    lit = pd.read_csv(csv_file, index_col=0)
    lit["ST"] = lit["Title"].apply(find_terms)
    lit["ST"] = lit["ST"].apply(lambda x : uniformST(", ".join(x)))


    hasST = lit["ST"].replace("", np.NaN).dropna() # dataframe with only ST in title
