from Bio import Entrez

# Set your email address (NCBI requires this)
Entrez.email = "benbud7@gmail.com"

def search_pubmed(query, max_results=10):
    """
    Search PubMed for the given query and return a list of PubMed IDs.
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"Error during PubMed search: {e}")
        return []

def fetch_pubmed_details(id_list):
    """
    Fetch details for a list of PubMed IDs.
    """
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

    # Parse and extract desired fields (title, abstract, PMID)
    articles = []
    for article in records.get("PubmedArticle", []):
        try:
            citation = article["MedlineCitation"]
            article_info = citation["Article"]
            title = article_info.get("ArticleTitle", "No title available")
            # The abstract may not always be present
            abstract_list = article_info.get("Abstract", {}).get("AbstractText", [])
            abstract = abstract_list[0] if abstract_list else "No abstract available"
            pmid = citation.get("PMID", "N/A")
        except Exception as e:
            print(f"Error parsing article: {e}")
            continue

        articles.append({
            "id": pmid,
            "title": title,
            "abstract": abstract
        })

    return articles

def search_and_fetch(query, max_results=10):
    """
    Convenience function to perform a search and fetch details.
    """
    id_list = search_pubmed(query, max_results)
    articles = fetch_pubmed_details(id_list)
    return articles