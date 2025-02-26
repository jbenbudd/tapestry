# backend/pubmed.py

from Bio import Entrez

Entrez.email = "benbud7@gmail.com"

def search_pubmed(query, max_results=10, start_date=None, end_date=None):
    try:
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "datetype": "pdat",  # Use publication date
        }
        if start_date:
            params["mindate"] = start_date
        if end_date:
            params["maxdate"] = end_date

        handle = Entrez.esearch(**params)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"Error during PubMed search: {e}")
        return []

def fetch_pubmed_details(id_list):
    if not id_list:
        return []

    ids = ",".join(id_list)
    try:
        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Error during PubMed fetch: {e}")
        return []

    articles = []
    for article in records.get("PubmedArticle", []):
        try:
            citation = article["MedlineCitation"]
            article_info = citation["Article"]
            title = article_info.get("ArticleTitle", "No title available")
            # Extract abstract
            abstract_list = article_info.get("Abstract", {}).get("AbstractText", [])
            abstract = abstract_list[0] if abstract_list else "No abstract available"
            pmid = citation.get("PMID", "N/A")

            # Extract authors (join LastName and Initials for each author)
            authors = []
            for author in article_info.get("AuthorList", []):
                last = author.get("LastName", "")
                initials = author.get("Initials", "")
                full = f"{last} {initials}".strip()
                if full:
                    authors.append(full)
            authors_str = ", ".join(authors) if authors else "No authors available"

            # Extract journal name
            journal = article_info.get("Journal", {}).get("Title", "No journal available")

            # Extract publication date from JournalIssue; fall back to MedlineDate if needed
            pub_date = "No publication date available"
            journal_issue = article_info.get("Journal", {}).get("JournalIssue", {})
            if "PubDate" in journal_issue:
                pub_date_data = journal_issue["PubDate"]
                if "Year" in pub_date_data:
                    pub_date = pub_date_data["Year"]
                elif "MedlineDate" in pub_date_data:
                    pub_date = pub_date_data["MedlineDate"]

            # Construct a link to PubMed for this article
            link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        except Exception as e:
            print(f"Error parsing article: {e}")
            continue

        articles.append({
            "id": pmid,
            "title": title,
            "abstract": abstract,
            "authors": authors_str,
            "journal": journal,
            "pub_date": pub_date,
            "link": link
        })

    return articles

def search_and_fetch(query, max_results=10, start_date=None, end_date=None):
    id_list = search_pubmed(query, max_results, start_date, end_date)
    articles = fetch_pubmed_details(id_list)
    return articles